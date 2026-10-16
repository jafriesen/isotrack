#define IsoTrack_cxx
#include "IsoTrack.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TStopwatch.h>
#include <TDatime.h>
#include <TMath.h>
#include <TSystem.h>
#include <cstdlib>

#include <TProfile2D.h>
#include <TProfile3D.h>
#include <TH2D.h>
#include <TH3D.h>

#include <set>
#include <map>
#include <iostream>

// Introduce _gainCorrectionRetriever global pointer
//#include "gainCorrections.C"

// Sunanda's CalibCorr.C for gain and phi symmetry corrections
// as well as thresholds and puFactor
#include "CalibCorr.C"

double puFactorCopy(int ieta, double pmom, double eHcal, double ediff, bool debug = false);

double getIsoTrackCorr(int run, int ieta, int depth);

/*
// Also defined in CalibCorr.C
void unpackDetId(unsigned int detId, int& subdet, int& zside, int& ieta, int& iphi, int& depth) {
  // The maskings are defined in DataFormats/DetId/interface/DetId.h
  //                      and in DataFormats/HcalDetId/interface/HcalDetId.h
  // The macro does not invoke the classes there and use them
  subdet = ((detId >> 25) & (0x7));
  if ((detId & 0x1000000) == 0) {
    ieta = ((detId >> 7) & 0x3F);
    zside = (detId & 0x2000) ? (1) : (-1);
    depth = ((detId >> 14) & 0x1F);
    iphi = (detId & 0x3F);
  } else {
    ieta = ((detId >> 10) & 0x1FF);
    zside = (detId & 0x80000) ? (1) : (-1);
    depth = ((detId >> 20) & 0xF);
    iphi = (detId & 0x3FF);
  }
}
*/

// From Long Wang, Re: ROOT file with HCalRespCorrs for comparison, 2024/2/1
double etaVal(int ieta) {
  //if (ieta>=0) ++ieta; // patch iEta!=0 (Mikko's addition)
  double etavl(0.);
  if (ieta <= -24)
    etavl = .1695*ieta + 1.9931;
  else if (ieta <= -1)
    etavl = .0875*ieta + .0489;
  else if  (ieta < 24)
    etavl = .0875*ieta - .0489;
  else
    etavl = .1695*ieta - 1.9931;
  return etavl;
} // etaVal


int delta_ieta(int ieta1, int ieta2) {
  if (ieta1>0 && ieta2>0) return +(ieta1-ieta2);
  if (ieta1<0 && ieta2<0) return -(ieta1-ieta2);
  if (ieta1<0 && ieta2>0) return +(ieta1-ieta2)+1;
  if (ieta1>0 && ieta2<0) return -(ieta1-ieta2)+1;

  return 0;
}

// Estimate effective area based on half granularity on iphi after |ieta|>=21
double areaScale3x5(int ieta, int iphi) {

  double areaScale1  = 5./(9.-5.); // = 1.25, |ieta|<=19
  double areaScale2  = 3./(5.-3.); // = 1.50, |ieta|>=22
  double areaScale3o = 3./(5.-3.); // = 1.50, |ieta|==20, odd iphi
  double areaScale3e = 2./(4.-2.); // = 1.00, |ieta|==20, even iphi
  if (abs(ieta)<=19) return areaScale1;
  if (abs(ieta)==20 && iphi%2==1) return (2./3.*areaScale1+1./3.*areaScale3o);
  if (abs(ieta)==20 && iphi%2==0) return (2./3.*areaScale1+1./3.*areaScale3e);
  if (abs(ieta)==21) return (1./3.*areaScale1+2./3.*areaScale2);
  if (abs(ieta)>=22) return areaScale2;

  return 1.5;
}

// Override new developments
bool useClassic = false;//true;
bool useSunanda = false;//true;

// Apply gain corrections from Yildiray
bool correctGains = false;//true; // v31 off

// Apply phi asymmetry corrections
bool correctPhis = false;//true; // v30 off

// Apply IsoTrack HCAL correction (closure)
bool correctHCAL = false;//true;

// Apply correct thresholds to RecHits after updating scale
bool correctCuts = true; // v32 off

// Test t_detIds3/getWeight bug
bool testWeightBug = false;//true;

// Window to calculate arithmetic mean
double wmin = 0.15;//0.30;//0.15;//0.30;//0.60;//0.30;//0.15;
double wmax = 1.85;//1.70;//1.85;//3.00;//1.20;//3.00;//1.85

// IsoTrack per depth and/or single depth
bool doPerDepth = true;
bool doSingleDepth = true;

// Limited size regions to include containment
bool enforce5x5 = false;
bool enforce3x3 = false;
bool enforce3x5 = true;//false;//true;

// Propagate updated single depths (gains and/or size) to classic single depth
bool updateSingleDepth = true;

// Check if (classic, non-updated) esum and branches are consistent
bool checkConsistency = true;

// Fill cluster shapes without limited-size mask
bool freePassP3 = false;

// Merge depths 1+2
//bool mergeDepths1and2 = false;//true;

TH2D *_h2prevcorr(0);
void IsoTrack::Loop()
{
//   In a ROOT session, you can do:
//      root> .L IsoTrack.C
//      root> IsoTrack t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   cout << "Calling IsoTrack::Loop()..." << endl << flush;
   
   if (correctGains || correctPhis || correctCuts || correctHCAL)    {
     if (correctGains) cout << "Correcting gain on the fly" << endl;
     if (correctPhis)  cout << "Correcting phi asymmetry on the fly" << endl;
     if (correctHCAL)  cout << "Correcting HCAL on the fly" << endl;
     if (correctCuts)  cout << "Correcting hit thresholds on the fly,"
			    << " on the _newly corrected_ energy" << endl;
     if (doPerDepth && updateSingleDepth && doSingleDepth)
       cout << "  with single depths updated" << endl;
   }

   if (doPerDepth)      cout << "Filling profiles per depth" << endl;
   if (doSingleDepth)   cout << "Filling profiles without depth" << endl;

   cout << "Analysis mode: ";
   if (useSunanda)      cout << "useSunanda";
   else if (useClassic) cout << "useClassic";
   else {
     if (enforce5x5)    cout << "enforce5x5";
     if (enforce3x3)    cout << "enforce3x3";
     if (enforce3x5)    cout << "enforce3x5";
   }
   cout << endl;

   // Set up timers to see progress. Full 2024 takes about 25 minutes
   TStopwatch fulltime, laptime;
   fulltime.Start();
   TDatime time_now, estimated_completion_time;
   int nlap(0);
   
   fChain->SetBranchStatus("*",0);

   fChain->SetBranchStatus("t_p",1);
   fChain->SetBranchStatus("t_eHcal",1);
   fChain->SetBranchStatus("t_ieta",1);
   // MIP correction
   fChain->SetBranchStatus("t_eMipDR",1); // also track ID
   // PU
   fChain->SetBranchStatus("t_eHcal10",1);
   fChain->SetBranchStatus("t_eHcal30",1);
   // track ID
   fChain->SetBranchStatus("t_selectTk",1);
   fChain->SetBranchStatus("t_qltyMissFlag",1);
   fChain->SetBranchStatus("t_hmaxNearP",1);
   fChain->SetBranchStatus("t_mindR1",1);

   // Process depths (heavy!)
   if (doPerDepth) {
     // Sanity checks
     fChain->SetBranchStatus("t_iphi",1);
     // Core energy per depth
     fChain->SetBranchStatus("t_DetIds",1);
     fChain->SetBranchStatus("t_HitEnergies",1);
     // PU per depth (just first cone now)
     fChain->SetBranchStatus("t_DetIds1",1);
     fChain->SetBranchStatus("t_HitEnergies1",1);
     fChain->SetBranchStatus("t_DetIds3",1);
     fChain->SetBranchStatus("t_HitEnergies3",1);
   }
   //fChain->SetBranchStatus("t_DetIds",1);

   // Do gain corrections
   if (correctGains || correctPhis) {
     fChain->SetBranchStatus("t_Run",1);
   }

   // Study AlCaRaw trigger bias: L1 pT and duplicate events
   if (true) {
     fChain->SetBranchStatus("t_nVtx");
     fChain->SetBranchStatus("t_nTrk");
     fChain->SetBranchStatus("t_goodPV");
     fChain->SetBranchStatus("t_rhoh");
     
     fChain->SetBranchStatus("t_l1pt");
     fChain->SetBranchStatus("t_phi");
     fChain->SetBranchStatus("t_pt");
     fChain->SetBranchStatus("t_Run",1);
     fChain->SetBranchStatus("t_Event",1);
   }

   // Load gain corrections and activate _gainCorrectionRetriever global pointer
   //gainCorrections();

   // Initialize gain corrector from Sunanda
   CalibDuplicate *cDuplicate_(0);
   if (correctGains) {
     string gainfilename = "textfiles/GainFactors2024New.txt";
     cout << "Read gain corrections from file '"<<gainfilename<<"'\n";
     cDuplicate_ = new CalibDuplicate(gainfilename.c_str(), 3, false);
   }

   // Initialize phi asymm correctors from Sunanda. File lists:
   // run_begin run_end run_textfilename
   CalibCorr *cFactor_(0);
   if (correctPhis) {
     string filelistname = "textfiles/PhiSym2024.txt";
     cout << "Read phi asymmetry corrections from file list '"<<filelistname<<"'...\n";
     cFactor_ = new CalibCorr(filelistname.c_str(), 5, false);
     /*
     ifstream fphi(phyasymfilelistname.c_str());
     int run_begin, run_end;
     string filename;
     while (fphi >> run_begin >> run_end >> filename) {
       cout << " runs ["<<run_begin<<","<<run_end<<"] file '"<<filename<<"'\n";
       cFactor_ = new CalibCorr(run_textfilename, 5, false); // phiasymm
     }
     */
   }

   CalibThreshold *cCut_(0);
   if (correctCuts) {
     // PF thresholds must match the data-taking year: textfiles/PFCuts<year>.txt,
     // year taken from the run number of the first entry (or ISOTRACK_PFCUTS).
     string pfcuts = (getenv("ISOTRACK_PFCUTS") ? getenv("ISOTRACK_PFCUTS") : "");
     if (pfcuts=="") {
       fChain->GetEntry(0);
       int year = (t_Run<355100 ? 0 : t_Run<366442 ? 2022 : t_Run<378981 ? 2023 :
		   t_Run<390000 ? 2024 : t_Run<400000 ? 2025 : 2026);
       pfcuts = Form("textfiles/PFCuts%d.txt",year);
     }
     cout << "Read HCalPFcuts from file " << pfcuts << " (first run " << t_Run << ")" << endl;
     if (gSystem->AccessPathName(pfcuts.c_str())) {
       cout << "ERROR: PF-cuts file " << pfcuts << " not found" << endl << flush;
       exit(1);
     }
     cCut_ = new CalibThreshold(pfcuts.c_str());
   }
   
   Long64_t nentries = fChain->GetEntries();//fChain->GetEntriesFast();

   cout << "Opening output file IsoTrack.root" << endl << flush;
   cout << "Starting loop over "<<nentries<<" events\n";
   TFile *fout = new TFile("IsoTrack.root","RECREATE");

   // AlCaRaw studies
   TH2D *h2nvtx, *h2ntrk, *h2goodpv, *h2rhoh, *h2phi;
   h2nvtx = new TH2D("h2nvtx",";ieta;t_nVtx", 61,-30.5,30.5, 200,0,200);
   h2ntrk = new TH2D("h2ntrk",";ieta;t_nTrk", 61,-30.5,30.5, 400,0,400);
   h2goodpv = new TH2D("h2goodpv",";ieta;t_goodPV", 61,-30.5,30.5, 200,0,200);
   h2rhoh = new TH2D("h2rhoh",";ieta;t_rhoh", 61,-30.5,30.5, 400,0,100);
   h2phi = new TH2D("h2phi",";ieta;t_phi", 61,-30.5,30.5, 72, 0, TMath::TwoPi());
   
   TH2D *h2p, *h2pt, *h2l1pt, *h2ehcal, *h2ehcal10, *h2ehcal30;
   h2p = new TH2D("h2p",";ieta;t_p", 61,-30.5,30.5, 400,0,200);
   h2pt = new TH2D("h2pt",";ieta;t_pt", 61,-30.5,30.5, 400,0,200);
   h2l1pt = new TH2D("h2l1pt",";ieta;t_l1pt", 61,-30.5,30.5, 400,0,200);
   h2ehcal = new TH2D("h2ehcal",";ieta;eHcal", 61,-30.5,30.5, 400,0,200);
   h2ehcal10 = new TH2D("h2ehcal10",";ieta;eHcal10", 61,-30.5,30.5, 400,0,200);
   h2ehcal30 = new TH2D("h2ehcal30",";ieta;eHcal30", 61,-30.5,30.5, 400,0,200);

   map<Int_t, map<Int_t, Bool_t> > mdup;
   TH1D *h1dup, *h1all;
   h1dup = new TH1D("h1dup",";ieta;", 61,-30.5,30.5);
   h1all = new TH1D("h1all",";ieta;", 61,-30.5,30.5);
   
   // Depth-independent results
   TH2D *h2raw, *h2pu1, *h2pu3, *h2mip, *h2c;
   h2raw = new TH2D("h2raw",";ieta;t_eHcal/(t_p-eMIP)", 61,-30.5,30.5, 400,0,4);
   h2pu1 = new TH2D("h2pu1",";ieta;(t_eHcal10-t_eHcal)/(t_p-eMIP)",
		    61,-30.5,30.5, 400,0,4);
   h2pu3 = new TH2D("h2pu3",";ieta;(t_eHcal30-t_eHcal10)/(t_p-eMIP)",
		    61,-30.5,30.5, 400,0,4);
   h2mip = new TH2D("h2mip",";ieta;t_eMipDR/t_p", 61,-30.5,30.5, 400,0,4);
   h2c = new TH2D("h2c",";ieta;(t_eHcal-ePU)/(t_p-eMIP)",
		  61,-30.5,30.5, 400,0,4);

   // Means
   string sw = Form("(t_p-eMIP) [%1.2f,%1.2f]",wmin,wmax);
   const char *cw = sw.c_str();
   TProfile *praw, *ppu1, *ppu3, *pmip, *pc, *pc_even, *pc_odd;
   praw = new TProfile("praw",Form(";ieta;t_eHcal/%s",cw),
		       61,-30.5,30.5);
   ppu1 = new TProfile("ppu1",Form(";ieta;(t_eHcal10-t_eHcal)/%s",cw),
		       61,-30.5,30.5);
   ppu3 = new TProfile("ppu3",Form(";ieta;(t_eHcal30-t_eHcal10)/%s",cw),
		       61,-30.5,30.5);
   pmip = new TProfile("pmip",Form(";ieta;t_eMipDR/%s",cw),
		       61,-30.5,30.5);   
   pc = new TProfile("pc",Form(";ieta;(t_eHcal-ePU)/%s",cw),
		     61,-30.5,30.5);
   pc_even = new TProfile("pc_even",Form(";ieta;(t_eHcal-ePU)/%s",cw),
			  61,-30.5,30.5);
   pc_odd = new TProfile("pc_odd",Form(";ieta;(t_eHcal-ePU)/%s",cw),
			 61,-30.5,30.5);

   // Same exercise now per depth
   // Depth=0 for ECAL (MIP or early shower)
   TH3D *h3raw, *h3pu1, *h3c;
   h3raw = new TH3D("h3raw",";ieta;depth; t_eHcal/(t_p-eMIP)",
		    61,-30.5,30.5, 10,-0.5,9.5, 400,0,4);
   h3pu1 = new TH3D("h3pu1",";ieta;depth;(t_eHcal10-t_eHcal)/(t_p-eMIP)",
		    61,-30.5,30.5, 10,-0.5,9.5, 400,0,4);
   h3c = new TH3D("h3c",";ieta;depth;(t_eHcal-ePU)/(t_p-eMIP)",
		  61,-30.5,30.5, 10,-0.5,9.5, 400,0,4);

   // Previous correction
   _h2prevcorr = new TH2D("h2prevcorr",";ieta;depth;Correction Factor",
			  61,-30.5,30.5, 10,-0.5,9.5);
   
   // Means
   TProfile2D *p2raw, *p2pu1, *p2c, *p2c_even, *p2c_odd;
   p2raw = new TProfile2D("p2raw",Form(";ieta;depth;t_eHcal/%s",cw),
			  61,-30.5,30.5, 10,-0.5,9.5);
   p2pu1 = new TProfile2D("p2pu1",Form(";ieta;depth;(t_eHcal10-t_eHcal)/%s",cw),
			  61,-30.5,30.5, 10,-0.5,9.5);
   p2c = new TProfile2D("p2c",Form(";ieta;depth;(t_eHcal-ePU)/%s",cw),
			61,-30.5,30.5, 10,-0.5,9.5);
   p2c_even = new TProfile2D("p2c_even",Form(";ieta;depth;(t_eHcal-ePU)/%s",cw),
			     61,-30.5,30.5, 10,-0.5,9.5);
   p2c_odd = new TProfile2D("p2c_odd",Form(";ieta;depth;(t_eHcal-ePU)/%s",cw),
			    61,-30.5,30.5, 10,-0.5,9.5);

   // Depth Covariances
   TProfile3D *p3c, *p3c_even, *p3c_odd;
   p3c = new TProfile3D("p3c",Form(";ieta;depth;depth;(t_eHcal-ePU)/%s",cw),
			61,-30.5,30.5, 10,-0.5,9.5, 10,-0.5,9.5);
   p3c_even = new TProfile3D("p3c_even",
			     Form(";ieta;depth;depth;(t_eHcal-ePU)/%s",cw),
			     61,-30.5,30.5, 10,-0.5,9.5, 10,-0.5,9.5);
   p3c_odd = new TProfile3D("p3c_odd",
			    Form(";ieta;depth;depth;(t_eHcal-ePU)/%s",cw),
			    61,-30.5,30.5, 10,-0.5,9.5, 10,-0.5,9.5);

   // Checking shape and energy distribution of isolation cone
   TH2D *h2draw_hb, *h2dpu1_hb;
   TH2D *h2draw_he1, *h2dpu1_he1, *h2draw_he2, *h2dpu1_he2;
   h2draw_hb = new TH2D("h2draw_hb",";#Delta|i#eta|;#Delta|i#phi|;",
			21,-10.5,10.5, 21,-10.5,10.5);
   h2dpu1_hb = new TH2D("h2dpu1_hb",";#Delta|i#eta|;#Delta|i#phi|;",
			21,-10.5,10.5, 21,-10.5,10.5);
   h2draw_he1 = new TH2D("h2draw_he1",";#Delta|i#eta|;#Delta|i#phi|;",
			 21,-10.5,10.5, 21,-10.5,10.5);
   h2dpu1_he1 = new TH2D("h2dpu1_he1",";#Delta|i#eta|;#Delta|i#phi|;",
			 21,-10.5,10.5, 21,-10.5,10.5);
   h2draw_he2 = new TH2D("h2draw_he2",";#Delta|i#eta|;#Delta|i#phi|;",
			 21,-10.5,10.5, 21,-10.5,10.5);
   h2dpu1_he2 = new TH2D("h2dpu1_he2",";#Delta|i#eta|;#Delta|i#phi|;",
			 21,-10.5,10.5, 21,-10.5,10.5);
   
   TProfile2D *p2draw_hb, *p2dpu1_hb;
   TProfile2D *p2draw_he1, *p2dpu1_he1, *p2draw_he2, *p2dpu1_he2;
   p2draw_hb= new TProfile2D("p2draw_hb",";#Deltai#eta;#Deltai#phi;"
			     "E_{calo}",21,-10.5,10.5, 21,-10.5,10.5);
   p2dpu1_hb = new TProfile2D("p2dpu1_hb",";#Deltai#eta;#Deltai#phi;"
			      "E_{calo}",21,-10.5,10.5, 21,-10.5,10.5);
   p2draw_he1 = new TProfile2D("p2draw_he1",";#Deltai#eta;#Deltai#phi;"
			       "E_{calo}",21,-10.5,10.5, 21,-10.5,10.5);
   p2dpu1_he1 = new TProfile2D("p2dpu1_he1",";#Deltai#eta;#Deltai#phi;"
			       "E_{calo}",21,-10.5,10.5, 21,-10.5,10.5);
   p2draw_he2 = new TProfile2D("p2draw_he2",";#Deltai#eta;#Deltai#phi;"
			       "E_{calo}",21,-10.5,10.5, 21,-10.5,10.5);
   p2dpu1_he2 = new TProfile2D("p2dpu1_he2",";#Deltai#eta;#Deltai#phi;"
			       "E_{calo}",21,-10.5,10.5, 21,-10.5,10.5);

   // Systematic study of pileup
   fout->mkdir("pileup");
   fout->cd("pileup");
   map<int, TProfile3D*> mp3_0;
   map<int, TProfile3D*> mp3_1;
   map<int, TProfile3D*> mp3_3;
   for (int ieta = -29; ieta != 29+1; ++ieta) {
     TProfile3D *p3_0 = new TProfile3D(Form("p3_0_ieta%d",ieta),";#Delta|i#eta|;#Delta|i#phi|;depth", 21,-10.5,10.5, 25,-12.5,12.5, 7,0.5,7.5);
     mp3_0[ieta] = p3_0;
     TProfile3D *p3_1 = new TProfile3D(Form("p3_1_ieta%d",ieta),";#Delta|i#eta|;#Delta|i#phi|;depth", 21,-10.5,10.5, 25,-12.5,12.5, 7,0.5,7.5);
     mp3_1[ieta] = p3_1;
     TProfile3D *p3_3 = new TProfile3D(Form("p3_3_ieta%d",ieta),";#Delta|i#eta|;#Delta|i#phi|;depth", 21,-10.5,10.5, 25,-12.5,12.5, 7,0.5,7.5);
     mp3_3[ieta] = p3_3;
   } // for ietax

   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;

      // Progress indicator on the screen
      if (jentry%100000==0) cout << "." << flush;

      // Timers to estimate finish time
      if (jentry==0) { fulltime.Start(); laptime.Start(); nlap = jentry;}
      else if (jentry==100000 || jentry==1000000 || jentry%5000000==0 ||
	       jentry==nentries-1) {

	nlap = jentry-nlap;
	double events_per_second = (double)jentry / fulltime.RealTime();
	double laps_per_second = (double)nlap / laptime.RealTime();
	double entries_to_go = nentries - jentry;
        cout << Form("\nEstimated remaining runtime:  %1.0f sec. "
                     " (%1.0f sec. from last %d)\n",
		     entries_to_go / events_per_second,
		     entries_to_go / laps_per_second, nlap)
             << flush;

	time_now.Set();
        estimated_completion_time.Set(time_now.Convert() + entries_to_go / events_per_second);
        cout << "Estimated completion time (from all events): " << estimated_completion_time.AsSQLString() << endl;
	estimated_completion_time.Set(time_now.Convert() + entries_to_go / laps_per_second);
        cout << "Estimated completion time (from last lap): " << estimated_completion_time.AsSQLString() << endl;
        //
        laptime.Reset(); nlap = jentry;
	fulltime.Continue();
	laptime.Continue();
      } // timers

      // Speedup for PU MC with lots of track below t_p<40
      b_t_p->GetEntry(ientry); //read only this branch
      if (t_p<40 || t_p>60) continue;
      
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      
      bool ok = ((t_selectTk) && (t_qltyMissFlag) && (t_hmaxNearP < 10.0) && (t_eMipDR < 1.0) && (t_mindR1 > 1.0));
      if (!ok) continue;
      
      //double eta = etaVal(t_ieta);
      double esum(0), esum1(0), esum3(0);
      if (!doPerDepth) {
	esum = t_eHcal; esum1 = t_eHcal10; esum3 = t_eHcal30;
      }
      if (doPerDepth) {
	double e[10] =  {0,0,0,0,0, 0,0,0,0,0};
	double e1[10] = {0,0,0,0,0, 0,0,0,0,0};
	double e3[10] = {0,0,0,0,0, 0,0,0,0,0};
	e[0] = t_eMipDR;
	int subdet, zside, ieta, iphi, depth;

	// Keep track of detIds seen so far in this event
	set<unsigned int> ids;
	
	// Core energies
	for (unsigned int idet = 0; idet != t_DetIds->size(); ++idet) {
	  unpackDetId((*t_DetIds)[idet], subdet, zside, ieta, iphi, depth);
	  ieta *= zside;
	  unsigned int id = (*t_DetIds)[idet];
	  double edet = (*t_HitEnergies)[idet];
	  //if (correctGains) edet *= _gainCorrectionRetriever->getCorrection(t_Run, ieta, depth);
	  if (correctGains) edet *= cDuplicate_->getCorr(t_Run, ieta, depth);
	  if (correctPhis)  edet *= cFactor_->getCorr(t_Run, id);
	  if (correctHCAL)  edet *= getIsoTrackCorr(t_Run, ieta, depth);
	  // Apply cut on _corrected_ energy
	  //if (correctCuts)  edet *= (edet>threshold(id, 4) ? 1 : 0);
	  if (correctCuts)  edet *= (edet>cCut_->threshold(id) ? 1 : 0);
 
	  /*
	  if (fabs(ieta-t_ieta)>3) { // looking into 5x5? sometimes off-center 
	    cout << "jentry = " << jentry
		 << ", ieta = " << ieta << ", iphi = " << iphi
		 << ", t_ieta = " << t_ieta << ", edet = " << edet << endl;
	    continue;
	  }
	  if (depth<0 || depth>9) {
	    cout << "jentry = " << jentry
		 << ", ieta = " << ieta << ", iphi = " << iphi
		 << ", depth = " << depth << ", edet = " << edet << endl;
	    continue;
	  }
	  */
	  //int dieta = (ieta-t_ieta)*TMath::Sign(1,t_ieta);
	  //if (ieta*t_ieta<0) dieta += TMath::Sign(1,t_ieta);
	  int dieta = delta_ieta(ieta,t_ieta);
	  //double diphi = (abs(iphi)-abs(t_iphi));
	  int diphi = iphi-t_iphi;
	  if (abs(diphi)>36) diphi += (diphi>0 ? -72 : +72); // iphi is cyclic

	  bool is5x5 = (fabs(dieta)<3 && fabs(diphi)<3);
	  bool is3x3 = (fabs(dieta)<2 && fabs(diphi)<2);
	  bool is3x5 = (fabs(dieta)<2 && fabs(diphi)<3);
	  bool noForce = ((!enforce5x5 && !enforce3x3 && !enforce3x5) ||
			  useClassic || useSunanda);
	  
	  bool pass = ((is5x5 && enforce5x5) || (is3x3 && enforce3x3) ||
		       (is3x5 && enforce3x5) || noForce);
	
	  if (pass) {

	    esum += edet;
	    e[depth] += edet;

	    if (abs(ieta)<=16) {
	      h2draw_hb->Fill(dieta, diphi);
	      p2draw_hb->Fill(dieta, diphi, edet / (t_p - t_eMipDR));
	    }
	    else if (abs(ieta)<=23) {
	      h2draw_he1->Fill(dieta, diphi);
	      p2draw_he1->Fill(dieta, diphi, edet / (t_p - t_eMipDR));
	    }
	    else {
	      h2draw_he2->Fill(dieta, diphi);
	      p2draw_he2->Fill(dieta, diphi, edet / (t_p - t_eMipDR));
	    }
	  } // is5x5

	  TProfile3D *p3_0 = mp3_0[t_ieta];
	  if (p3_0 && (pass || freePassP3)) {
	    p3_0->Fill(dieta, diphi, depth, edet / (t_p - t_eMipDR));
	  }
	  ids.insert((*t_DetIds)[idet]);
	} // for idet

	// Cone energies for PU estimates (small cone)
	for (unsigned int idet = 0; idet != t_DetIds1->size(); ++idet) {
	  unpackDetId((*t_DetIds1)[idet], subdet, zside, ieta, iphi, depth);
	  ieta *= zside;
	  unsigned int id = (*t_DetIds1)[idet];
	  double edet = (*t_HitEnergies1)[idet];
	  //if (correctGains) edet *= _gainCorrectionRetriever->getCorrection(t_Run, ieta, depth);
	  if (correctGains) edet *= cDuplicate_->getCorr(t_Run, ieta, depth);
	  if (correctPhis)  edet *= cFactor_->getCorr(t_Run, id);
	  if (correctHCAL)  edet *= getIsoTrackCorr(t_Run, ieta, depth);
	  // Apply cut on _corrected_ energy
	  //if (correctCuts)  edet *= (edet>threshold(id, 4) ? 1 : 0);
	  if (correctCuts)  edet *= (edet>cCut_->threshold(id) ? 1 : 0);
	  
	  //int dieta = (ieta-t_ieta)*TMath::Sign(1,t_ieta);
	  //if (ieta*t_ieta<0) dieta += TMath::Sign(1,t_ieta);
	  int dieta = delta_ieta(ieta,t_ieta);
	  int diphi = iphi-t_iphi;
	  if (abs(diphi)>36) diphi += (diphi>0 ? -72 : +72); // iphi is cyclic

	  bool is5x5 = (fabs(dieta)<3 && fabs(diphi)<3);
	  bool is3x3 = (fabs(dieta)<2 && fabs(diphi)<2);
	  bool is3x5 = (fabs(dieta)<2 && fabs(diphi)<3);
	  bool noForce = ((!enforce5x5 && !enforce3x3 && !enforce3x5) ||
			  useClassic || useSunanda);
	  
	  bool pass = ((is5x5 && enforce5x5) || (is3x3 && enforce3x3) ||
		       (is3x5 && enforce3x5) || noForce);
	  
	  //if (ids.find((*t_DetIds1)[idet])==ids.end()) {
	  TProfile3D *p3_1 = mp3_1[t_ieta];
	  if (p3_1 && (pass || freePassP3)) {
	    p3_1->Fill(dieta, diphi, depth, edet / (t_p - t_eMipDR));
	  }
	  //ids.insert((*t_DetIds1)[idet]);
	  //}

	  if (pass) {
	    esum1 += edet;
	    e1[depth] += edet;
	  }
	} // for idet1

        // Cone energies for PU estimates (large cone)
	for (unsigned int idet = 0; idet != t_DetIds3->size(); ++idet) {
	  unpackDetId((*t_DetIds3)[idet], subdet, zside, ieta, iphi, depth);
	  ieta *= zside;
	  unsigned int id = (*t_DetIds3)[idet];
	  double edet = (*t_HitEnergies3)[idet];
	  //if (correctGains) edet *= _gainCorrectionRetriever->getCorrection(t_Run, ieta, depth);
	  if (correctGains) edet *= cDuplicate_->getCorr(t_Run, ieta, depth);
	  if (correctPhis)  edet *= cFactor_->getCorr(t_Run, id);
	  if (correctHCAL)  edet *= getIsoTrackCorr(t_Run, ieta, depth);
	  if (testWeightBug) edet *= cDuplicate_->getWeight(id);
	  // Apply cut on _corrected_ energy
	  //if (correctCuts)  edet *= (edet>threshold(id, 4) ? 1 : 0);
	  if (correctCuts)  edet *= (edet>cCut_->threshold(id) ? 1 : 0);
	  
	  /*
	  if (fabs(ieta-t_ieta)>4) { // looking into 7x7?
	    cout << "jentry = " << jentry
		 << ", ieta1 = " << ieta << ", iphi = " << iphi
		 << ", t_ieta = " << t_ieta << ", edet = " << edet << endl;
	    continue;
	  }
	  if (depth<0 || depth>9) {
	    cout << "jentry = " << jentry
		 << ", ieta1 = " << ieta << ", iphi = " << iphi
		 << ", depth = " << depth << ", edet = " << edet << endl;
	    continue;
	  }
	  */

	  //int dieta = (ieta-t_ieta)*TMath::Sign(1,t_ieta);
	  //if (ieta*t_ieta<0) dieta += TMath::Sign(1,t_ieta);
	  int dieta = delta_ieta(ieta,t_ieta);
	  int diphi = iphi-t_iphi;
	  if (abs(diphi)>36) diphi += (diphi>0 ? -72 : +72); // iphi is cyclic
	  
	  //bool is5x9 = (fabs(dieta)<3 && fabs(diphi)<5);
	  //bool is3x9 = (fabs(dieta)<2 && fabs(diphi)<5);

	  bool is5x5 = (fabs(dieta)<3 && fabs(diphi)<3);
	  bool is3x3 = (fabs(dieta)<2 && fabs(diphi)<2);
	  bool is3x5 = (fabs(dieta)<2 && fabs(diphi)<3);
	  
	  bool doForce = ((enforce5x5 || enforce3x3 || enforce3x5) &&
			  !useClassic && !useSunanda);
	  bool is5x9 = (fabs(dieta)<3 && fabs(diphi)<7 &&
			fabs(diphi)!=3 && fabs(diphi)!=4);
	  //bool is3x9 = (fabs(dieta)<2 && fabs(diphi)<7 &&
	  //		fabs(diphi)!=2 && fabs(diphi)!=3);
	  bool is3x9 = (fabs(dieta)<2 && fabs(diphi)<5 &&
	  		fabs(diphi)!=2);
	  bool is3x7 = (fabs(dieta)<2 && fabs(diphi)<5);
	  
	  bool pass = ((is5x9 && enforce5x5) || (is3x9 && enforce3x3) ||
		       (is3x7 && enforce3x5) || !doForce);

	  //bool veto = (((is5x5 && enforce5x5) || (is3x3 && enforce3x3) ||
	  //		(is3x5 && enforce3x5)) && doForce);
	  //pass = (pass && !veto);
	  
	  if (pass) {
	    
	    esum3 += edet;
	    e3[depth] += edet;

	    // Control plots of pileup subtraction
	    if (abs(ieta)<=16) {
	      h2dpu1_hb->Fill(dieta, diphi);
	      p2dpu1_hb->Fill(dieta, diphi, edet / (t_p - t_eMipDR));
	    }
	    else if (abs(ieta)<=23) {
	      h2dpu1_he1->Fill(dieta, diphi);
	      p2dpu1_he1->Fill(dieta, diphi, edet / (t_p - t_eMipDR));
	    }
	    else {
	      h2dpu1_he2->Fill(dieta, diphi);
	      p2dpu1_he2->Fill(dieta, diphi, edet / (t_p - t_eMipDR));
	    }
	  } // is5x9

	  //if (ids.find((*t_DetIds3)[idet])==ids.end()) {
	  TProfile3D *p3_3 = mp3_3[t_ieta];
	  if (p3_3 && (pass || freePassP3)) {
	    p3_3->Fill(dieta, diphi, depth, edet / (t_p - t_eMipDR));
	  }
	  //ids.insert((*t_DetIds3)[idet]);
	  //}
	} // for idet

	bool doForce = ((enforce5x5 || enforce3x3 || enforce3x5) &&
			!useClassic && !useSunanda);
	//if (doForce) { // with veto enabled
	//esum3 += esum;
	//for (int i = 0; i != 10; ++i) e3[i] += e[i];
	//}
	
	if (updateSingleDepth) {
	  t_eHcal = esum;
	  t_eHcal10 = esum1;
	  t_eHcal30 = esum3;
	}
	
	// Verify sums
	if (checkConsistency && !doForce && fabs(esum-t_eHcal)>1) {
	  cout << "jentry = " << jentry
	       << ", t_ieta = " << t_ieta << ", t_iphi = " << t_iphi
	       << ", esum = " << esum << ", t_eHcal = " << t_eHcal << endl;
	  continue;
	}
	if (checkConsistency && !doForce && fabs(esum1-t_eHcal10)>1) {
	  cout << "jentry = " << jentry
	       << ", t_ieta = " << t_ieta << ", t_iphi = " << t_iphi
	       << ", esum1 = " << esum1 << ", t_eHcal1 = " << t_eHcal10 << endl;
	  continue;
	}
	if (checkConsistency && !doForce && fabs(esum3-t_eHcal30)>1) {
	  cout << "jentry = " << jentry
	       << ", t_ieta = " << t_ieta << ", t_iphi = " << t_iphi
	       << ", esum3 = " << esum3 << ", t_eHcal3 = " << t_eHcal30 << endl;
	  continue;
	}

	// Merge depths 1+2 by replacing them both with half-average
	// This will preserve the array sizes for later
	// => matrix was singular, try setting e[1] to eps and e[2] to sum
	/*
	if (mergeDepths1and2) {
	  //double e12_avg = 0.5*(e[1]+e[2]);
	  //double pu12_avg = 0.5*(e1[1]+e1[2]);
	  //e[1] = e[2] = e12_avg;
	  //e1[1] = e1[2] = pu12_avg;
	  double eps = 0.01; // small non-zero energy offset
	  double e12_sum = (e[1]+e[2]);
	  double pu12_sum = (e1[1]+e1[2]);
	  e[1] = eps; e[2] = e12_sum;
	  e1[1] = 0; e1[2] = pu12_sum;
	}
	*/
	
	// Fill histograms and profiles per depth
	const int ndepth = 10;
	double vrc[ndepth];
	//double rcsum = rc; // precalculated (or sum explicitly)
	double areaScale(1.0);
	if (enforce5x5) areaScale = 0.5;//(5.*5.)/(5.*9.-5.*5.);
	if (enforce3x3) areaScale = 0.5;//(3.*3.)/(3.*9.-3.*3.);
	if (enforce3x5) areaScale = areaScale3x5(t_ieta, t_iphi);
	double epusum = (esum3 - esum) * areaScale;
	double rcsum = (esum - epusum) / (t_p - t_eMipDR);
	for (int i = 0; i != ndepth; ++i) {

	  double depth = i;
	  double epu1 = (i==0 ? 0 : e1[i]-e[i]);
	  //double epu3 = (t_eHcal30-t_eHcal10)*0.5;
	  double epu = (i==0 ? 0 : e3[i]-e[i]) * areaScale;
	  double rc = (e[i]-epu) / (t_p - t_eMipDR);

	  if (useClassic) {
	    rc = (e[i]-epu1) / (t_p - t_eMipDR);
	    rcsum = (t_eHcal-(t_eHcal10-t_eHcal)) / (t_p - t_eMipDR);
	    epu = epu1;
	  }
	  
	  if (useSunanda) {
	    // esum=etot in Sunanda's code, t_p=pmom
	    double ediff = esum3 - esum1;
	    double ehcal = esum * puFactor(+8, t_ieta, t_p, esum, ediff);
	    //double ehcal = esum * puFactorCopy(t_ieta, t_p, esum, ediff);
	    // ehcal = esum - k * ediff
	    double k = (ediff>0 ? (esum - ehcal) / ediff : 1);
	    double epu3 = k * (e3[i] - e1[i]);
	    rc = (e[i] - epu3) / (t_p - t_eMipDR);
	    epu = k * ediff;
	    rcsum = (esum - epu) / (t_p - t_eMipDR);
	  }

	  h3raw->Fill(t_ieta, depth, e[i] / (t_p - t_eMipDR));
	  h3pu1->Fill(t_ieta, depth, epu / (t_p - t_eMipDR));
	  h3c->Fill(t_ieta, depth, rc);
	  
	  if (rcsum>wmin && rcsum<wmax) {
	    p2raw->Fill(t_ieta, depth, e[i] / (t_p - t_eMipDR));
	    p2pu1->Fill(t_ieta, depth, epu / (t_p - t_eMipDR));
	    p2c->Fill(t_ieta, depth, rc);
	    if (jentry%2==0) p2c_even->Fill(t_ieta, depth, rc);
	    if (jentry%2==1) p2c_odd->Fill(t_ieta, depth, rc);
	    vrc[i] = rc;
	  }
	  else {
	    vrc[i] = 0.;
	  }
	} // for i
	
	// Calculate covariance
	for (int i = 0; i != ndepth; ++i) {
	  for (int j = 0; j != ndepth; ++j) {
	    if (rcsum>wmin && rcsum<wmax) {
	      p3c->Fill(t_ieta, i, j, vrc[i]*vrc[j]);
	      if (jentry%2==0) p3c_even->Fill(t_ieta, i, j, vrc[i]*vrc[j]);
	      if (jentry%2==1) p3c_odd->Fill(t_ieta, i, j, vrc[i]*vrc[j]);
	    }
	  } // for i
	} // for j

      } // doPerDepth

      
      if (doSingleDepth) {

	double areaScale(1.0);
	if (enforce5x5) areaScale = 0.5;//(5.*5.)/(5.*9.-5.*5.);
	if (enforce3x3) areaScale = 0.5;//(3.*3.)/(3.*9.-3.*3.);
	if (enforce3x5) areaScale = areaScale3x5(t_ieta, t_iphi);
	double epu1 = (t_eHcal10-t_eHcal);
	double epu3 = (t_eHcal30-t_eHcal10)*0.5;
	double epu = (esum3 - esum) * areaScale;
	double rc = (esum-epu) / (t_p - t_eMipDR);

	if (useClassic) {
	  epu = epu1;
	  rc = (t_eHcal-epu1) / (t_p - t_eMipDR);
	}

	if (useSunanda) {
	  // esum=etot in Sunanda's code, t_p=pmom
	  double ediff = esum3 - esum1;
	  double ehcal = esum * puFactor(+8, t_ieta, t_p, esum, ediff);
	  //double ehcal = esum * puFactorCopy(t_ieta, t_p, esum, ediff);
	  // ehcal = esum - k * ediff
	  double k = (ediff>0 ? (esum - ehcal) / ediff : 1);
	  epu = k * ediff;
	  rc = (esum - epu) / (t_p - t_eMipDR);
	}

	bool isDuplicate(false);
	if (mdup[t_Run][t_Event]==1) isDuplicate = true;
	mdup[t_Run][t_Event] = 1;

	h2nvtx->Fill(t_ieta, t_nVtx);
	h2ntrk->Fill(t_ieta, t_nTrk);
	h2goodpv->Fill(t_ieta, t_goodPV);
	h2rhoh->Fill(t_ieta, t_rhoh);
	h2phi->Fill(t_ieta, t_phi);
      
	h2p->Fill(t_ieta, t_p);
	h2pt->Fill(t_ieta, t_pt);
	h2l1pt->Fill(t_ieta, t_l1pt);
	if (!isDuplicate) {
	  h2ehcal->Fill(t_ieta, t_eHcal*t_pt/t_p);
	  h2ehcal10->Fill(t_ieta, t_eHcal10*t_pt/t_p);
	  h2ehcal30->Fill(t_ieta, t_eHcal30*t_pt/t_p);
	}
	if (isDuplicate) h1dup->Fill(t_ieta);
	h1all->Fill(t_ieta);

	h2raw->Fill(t_ieta, t_eHcal / (t_p - t_eMipDR));
	h2pu1->Fill(t_ieta, epu / (t_p - t_eMipDR));
	h2pu3->Fill(t_ieta, epu3 / (t_p - t_eMipDR));
	h2mip->Fill(t_ieta, t_eMipDR / (t_p - t_eMipDR));
	h2c->Fill(t_ieta, rc);
	
	if (rc>wmin && rc<wmax) {
	  praw->Fill(t_ieta, t_eHcal / (t_p - t_eMipDR));
	  ppu1->Fill(t_ieta, epu / (t_p - t_eMipDR));
	  ppu3->Fill(t_ieta, epu3 / (t_p - t_eMipDR));
	  pmip->Fill(t_ieta, t_eMipDR / (t_p - t_eMipDR));
	  pc->Fill(t_ieta, rc);
	  if (jentry%2==0) pc_even->Fill(t_ieta, rc);
	  if (jentry%2==1) pc_odd->Fill(t_ieta, rc);
	}
      }

   } // for jentry

   fout->Write();
   fout->Close();

   cout << "\nClosed output file IsoTrack.root" << endl << flush;
   cout << "IsoTrack::Loop() finished.\n" << endl << flush;
} // Loop


// Copied from CalibCorr.C from Sunanda Banerjee
// type == 8: Run3 MAHI. Deleted all the rest of the code
//double puFactor(int type, int ieta, double pmom, double eHcal, double ediff, bool debug = false) {
double puFactorCopy(int ieta, double pmom, double eHcal, double ediff, bool debug) {
  
  double fac(1.0);
  if (debug)
    std::cout << "Input Type " << 8 << " ieta " << ieta << " pmon " << pmom << " E " << eHcal << ":" << ediff;
  
  int jeta = std::abs(ieta);
  double d2p = (ediff / pmom);
  const double DELTA_CUT = 0.03;
  
  // type == 8: Run3 MAHI. Deleted all the rest
  //} else {  // Mahi 22pu (Jan, 2022)
  const double CONST_COR_COEF[6] = {0.995902, 0.991240, 0.981019, 0.788052, 0.597956, 0.538731};
  const double LINEAR_COR_COEF[6] = {-0.0540563, -0.104361, -0.215936, -0.147801, -0.160845, -0.154359};
  const double SQUARE_COR_COEF[6] = {0, 0, 0.0365911, 0.0161266, 0.0180053, 0.0184295};
  const int PU_IETA_1 = 7;
  const int PU_IETA_2 = 16;
  const int PU_IETA_3 = 25;
  const int PU_IETA_4 = 26;
  const int PU_IETA_5 = 27;
  unsigned icor = (unsigned(jeta >= PU_IETA_1) +
		   unsigned(jeta >= PU_IETA_2) +
		   unsigned(jeta >= PU_IETA_3) +
		   unsigned(jeta >= PU_IETA_4) +
		   unsigned(jeta >= PU_IETA_5));
  double deltaCut = (icor > 2) ? 1.0 : DELTA_CUT;
  if (d2p > deltaCut)
    fac = (CONST_COR_COEF[icor] + LINEAR_COR_COEF[icor] * d2p + SQUARE_COR_COEF[icor] * d2p * d2p);
  if (debug)
    std::cout << " d2p " << d2p << ":" << DELTA_CUT << " coeff " << icor << ":" << CONST_COR_COEF[icor] << ":"
	      << LINEAR_COR_COEF[icor] << ":" << SQUARE_COR_COEF[icor] << " Fac " << fac;
//}
//}
  if (fac < 0 || fac > 1)
    fac = 0;
  if (debug)
    std::cout << " Final factor " << fac << std::endl;
  return fac;
} // puFactorCopy

/*
void test() {
  // PU correction only for loose isolation cut
  double ehcal = (((rcorForm_ == 3) && (cFactor_ != nullptr))
                      ? (etot * cFactor_->getCorr(entry))
                      : ((puCorr_ == 0) ? etot
                                        : ((puCorr_ < 0) ? (etot * puFactor(-puCorr_, t_ieta, pmom, etot, ediff))
                                                         : puFactorRho(puCorr_, t_ieta, t_rhoh, etot))));
}
*/

double _isotrack[59][7]; bool _isotrack_isset(false);
double getIsoTrackCorr(int run, int ieta, int depth) {

  // Load corrections from file. Make smarter and using run range in the future
  if (_isotrack_isset==false) {
    //TFile *f = new TFile("rootfiles/CorrFact_hybrid_lxplus_v24_24CDEFGHI.root","READ"); // gains+phi+cuts+puFactor, no prevcorr
    //TFile *f = new TFile("rootfiles/CorrFact_hybrid_lxplus_v26_24CDEFGHI.root","READ"); // gains+phi+cuts+puFactor+prevcorr_v24
    TFile *f = new TFile("rootfiles/CorrFact_hybrid_lxplus_v27_24CDEFGHI.root","READ"); // gains+phi+cuts+puFactor+prevcorr_v26
    assert(f && !f->IsZombie());
    for (int depth = 1; depth != 8; ++depth) {
      TH1D *h = (TH1D*)f->Get(Form("hf_dd_%d",depth)); assert(h);
      for (int ieta = -29; ieta != 30; ++ieta) {
	int i = (ieta+29);
	int j = (depth-1);
	int k = h->FindBin(ieta);
	double corr = h->GetBinContent(k);
	_isotrack[i][j] = (corr>0.1 ? corr : 1.0);
      } // for ieta
    } // for depth
    f->Close();
    _isotrack_isset = true;

    // Print corrections
    cout << endl << "  double isotrack[59][7] = {\n";
    for (int i = 0; i != 59; ++i) {
      cout << "    {";
      for (int j = 0; j != 7; ++j) {
	double corr = _isotrack[i][j];
	cout << Form("%s%1.3f", j==0 ? "" : ", ", corr);
	if (_h2prevcorr) _h2prevcorr->SetBinContent(i+2, j+2, corr);
      } // for j
      cout << "}" << (i==29 ? "\n" : ",\n");
    } // for i
    cout << "  };\n" << endl;
  } //_isotrack_isset==false

  if (ieta<-29 || ieta>29 || ieta==0) return 1.0;
  int i = (ieta+29); i = max(0,min(58,i));
  int j = (depth-1); j = max(0,min(6,j));

  return _isotrack[i][j];
} // getIsoTrackCorr
