#define IsoTrackV2_cxx
#include "IsoTrackV2.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TStopwatch.h>
#include <TDatime.h>

#include <TProfile2D.h>
#include <TProfile3D.h>
#include <TH2D.h>
#include <TH3D.h>

#include <set>
#include <iostream>
#include <cmath>

// Introduce _gainCorrectionRetriever global pointer
//#include "gainCorrections.C"

// Sunanda's CalibCorr.C for gain and phi symmetry corrections
// as well as thresholds and puFactor
#include "CalibCorr.C"

double getIsoTrackV2Corr(int run, int ieta, int depth);

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
double etaValV2(int ieta) {
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


int delta_ietaV2(int ieta1, int ieta2) {
  /*
  // ieta1-ieta2, mirrored on negative side
  if (ieta1>0 && ieta2>0) return +(ieta1-ieta2);
  if (ieta1<0 && ieta2<0) return -(ieta1-ieta2);
  // ieta=0 does not exist so skip it. Use ieta2 as reference side
  if (ieta1<0 && ieta2>0) return +(ieta1-ieta2)+1;
  if (ieta1>0 && ieta2<0) return -(ieta1-ieta2)+1;
  */
  // Only remove ieta==0
  if ((ieta1*ieta2)>0) return (ieta1-ieta2);
  if (ieta1>0 && ieta2<0) return (ieta1-ieta2)-1;
  if (ieta1<0 && ieta2>0) return (ieta1-ieta2)+1;
  
  return 0;
}

int delta_iphiV2(int iphi1, int iphi2) {

  int diphi = iphi1-iphi2;
  if (abs(diphi)>36) diphi += (diphi>0 ? -72 : +72); // iphi is cyclic
  return diphi;
}

// Estimate effective area based on half granularity on iphi after |ieta|>=21
double areaScale3x5V2(int ieta, int iphi) {

  //double areaScale1  = 5./(9.-5.); // = 1.25, |ieta|<=19
  //double areaScale2  = 3./(5.-3.); // = 1.50, |ieta|>=22
  //double areaScale3o = 3./(5.-3.); // = 1.50, |ieta|==20, odd iphi
  //double areaScale3e = 2./(4.-2.); // = 1.00, |ieta|==20, even iphi
  //if (abs(ieta)<=19) return areaScale1;
  //if (abs(ieta)==20 && iphi%2==1) return (2./3.*areaScale1+1./3.*areaScale3o);
  //if (abs(ieta)==20 && iphi%2==0) return (2./3.*areaScale1+1./3.*areaScale3e);
  //if (abs(ieta)==21) return (1./3.*areaScale1+2./3.*areaScale2);
  //if (abs(ieta)>=22) return areaScale2;

  //if (abs(ieta)<=20) return areaScale1;
  //if (abs(ieta)>=21) return areaScale2;

  // areaScale = 2*sideband / core (#RecHits in |delta_ieta|<=2, 3-4
  double areaScale1  = 5./4.; // =1.25, |ieta|<=20
  double areaScale2o = 3./2.; // =1.50, |ieta|>=21, odd iphi
  double areaScale2e = 2./2.; // =1.00, |ieta|>=21, even iphi
  if (abs(ieta)<=20) return areaScale1;
  if (abs(ieta)>=21 && iphi%2==1) return areaScale2o;
  if (abs(ieta)>=21 && iphi%2==0) return areaScale2e;
  return 1.5;
}

// Apply gain corrections from Yildiray
bool correctGainsV2 = true;

// Apply phi asymmetry corrections
bool correctPhisV2 = true;

// Apply correct thresholds to RecHits
bool correctCutsV2 = true;

// Apply IsoTrackV2 HCAL correction (closure)
bool correctHCALV2 = false;//true;

// Window to calculate arithmetic mean
double wminV2 = 0.15;
double wmaxV2 = 1.85;

TH2D *_h2prevcorrV2(0);
void IsoTrackV2::Loop()
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

   cout << "Calling IsoTrackV2::Loop()..." << endl << flush;
   
   if (correctGainsV2) cout << "Correcting gain on the fly\n";
   else                cout << "NOT correcting gain on the fly\n";
   if (correctPhisV2)  cout << "Correcting phi asymmetry on the fly\n";
   else                cout << "NOT correcting phi asymmetry on the fly\n";
   if (correctCutsV2)  cout << "Correcting hit thresholds on the fly\n";
   else                cout << "NOT correcting hit thresholds on the fly\n";
   if (correctHCALV2)  cout << "Correcting HCAL on the fly\n";
   else                cout << "NOT correcting HCAL on the fly\n";

   // Set up timers to see progress. Full 2024 takes about 25 minutes
   TStopwatch fulltime, laptime;
   fulltime.Start();
   TDatime time_now, estimated_completion_time;
   int nlap(0);
   
   fChain->SetBranchStatus("*",0);

   // Basic coordinates and properties
   fChain->SetBranchStatus("t_Run",1);
   fChain->SetBranchStatus("t_p",1);
   fChain->SetBranchStatus("t_eHcal",1);
   fChain->SetBranchStatus("t_ieta",1);
   fChain->SetBranchStatus("t_iphi",1);

   // MIP correction
   fChain->SetBranchStatus("t_eMipDR",1); // also track ID

   // Track ID
   fChain->SetBranchStatus("t_selectTk",1);
   fChain->SetBranchStatus("t_qltyMissFlag",1);
   fChain->SetBranchStatus("t_hmaxNearP",1);
   fChain->SetBranchStatus("t_mindR1",1);

   // RecHits to reconstructed 3x5 core and 3x2 sideband
   fChain->SetBranchStatus("t_DetIds3",1);
   fChain->SetBranchStatus("t_HitEnergies3",1);

   // Load gain corrections and activate _gainCorrectionRetriever global pointer
   //gainCorrections();

   // Initialize gain corrector from Sunanda
   CalibDuplicate *cDuplicate_(0);
   if (correctGainsV2) {
     string gainfilename = "textfiles/GainFactors2024New.txt";
     cout << "Read gain corrections from file '"<<gainfilename<<"'\n";
     cDuplicate_ = new CalibDuplicate(gainfilename.c_str(), 3, false);
   }

   // Initialize phi asymm correctors from Sunanda. File lists:
   // run_begin run_end run_textfilename
   CalibCorr *cFactor_(0);
   if (correctPhisV2) {
     string filelistname = "textfiles/PhiSym2024.txt";
     cout << "Read phi asymmetry corrections from file list '"<<filelistname<<"'...\n";
     cFactor_ = new CalibCorr(filelistname.c_str(), 5, false);
   }
   

   Long64_t nentries = fChain->GetEntries();//fChain->GetEntriesFast();

   cout << "Opening output file IsoTrackV2.root" << endl << flush;
   cout << "Starting loop over "<<nentries<<" events\n";
   TFile *fout = new TFile("IsoTrackV2.root","RECREATE");

   // Number of "channels" for each ieta to keep track of:
   const int nband = 2; // (core,side)
   const int nwidth = 3; // (delta_ieta=-1,0,1)
   const int ndepth = 7; // depths 1-7
   const int nch = nband*nwidth*ndepth;
   
   // Depth means
   string sw = Form("(t_p-eMIP) [%1.2f,%1.2f]",wminV2,wmaxV2);
   const char *cw = sw.c_str();
   TProfile2D *p2, *p2_odd, *p2_evn;
   p2 = new TProfile2D("p2",Form(";ieta;channel;t_eHcal/%s",cw),
		       58,-29,29, nch,0,nch);
   p2_odd = (TProfile2D*)p2->Clone("p2_odd");
   p2_evn = (TProfile2D*)p2->Clone("p2_evn");
   
   // Depth covariances
   TProfile3D *p3, *p3_odd, *p3_evn;
   p3 = new TProfile3D("p3",Form(";ieta;channel;channel;t_eHcal/%s",cw),
		       58,-29,29, nch,0,nch, nch,0,nch);
   p3_odd = (TProfile3D*)p3->Clone("p3_odd");
   p3_evn = (TProfile3D*)p3->Clone("p3_evn");
   
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

      // Setup regions for (core,side)x(delta_ieta=-1,0,+1)x(depth) measurements
      double esum[nband][nwidth][ndepth];
      for (int i = 0; i != nband; ++i) { // (core,side)
	for (int j = 0; j != nwidth; ++j) {
	  for (int k = 0; k != ndepth; ++k) {
	    esum[i][j][k] = 0;
	  }
	} // for j
      } // for i

      // Sum up RecHits in core and sideband for means and correlations
      int subdet, zside, ieta, iphi, depth;
      for (unsigned int idet = 0; idet != t_DetIds3->size(); ++idet) {

	unpackDetId((*t_DetIds3)[idet], subdet, zside, ieta, iphi, depth);
	ieta *= zside;
	unsigned int id = (*t_DetIds3)[idet];
	double edet = (*t_HitEnergies3)[idet];
	
	//if (correctGainsV2) edet *= _gainCorrectionRetriever->getCorrection(t_Run, ieta, depth);
	if (correctCutsV2)  edet *= (edet>threshold(id, 3) ? 1 : 0);
	if (correctGainsV2) edet *= cDuplicate_->getCorr(t_Run, ieta, depth);
	if (correctPhisV2)  edet *= cFactor_->getCorr(t_Run, id);
	if (correctHCALV2)  edet *= getIsoTrackV2Corr(t_Run, ieta, depth);
	
	int dieta = delta_ietaV2(ieta,t_ieta);
	int diphi = delta_iphiV2(iphi,t_iphi);
	  
	bool is3x5core = (abs(dieta)<2 && abs(diphi)<3);
	bool is3x2side = (abs(dieta)<2 && abs(diphi)>=3 && abs(diphi)<5);

	if (is3x5core) {
	  esum[0][dieta+1][depth-1] += edet;
	}
	if (is3x2side) {
	  esum[1][dieta+1][depth-1] += edet * areaScale3x5V2(ieta, iphi);
	}
      } // for idet

      // Calculate flattened response for windowing later
      double rcsum(0);
      for (int j = 0; j != nwidth; ++j) {
	for (int k = 0; k != ndepth; ++k) {
	  //double epusum = (esum3 - esum) * areaScale;
	  //double rcsum = (esum - epusum) / (t_p - t_eMipDR);
	  double epu = esum[1][j][k];
	  double etot = esum[0][j][k];
	  double rc = (etot - epu) / (t_p - t_eMipDR);
	  rcsum += rc;
	} // for k
      } // for j
	  
      // Fill only results for compatible calo-track pairs
      if (rcsum>wminV2 && rcsum<wmaxV2) {

	// Fast memory access by using fact that arrays are consecutive blocks
	// Just need pointer to first element
	double *p = &(esum[0][0][0]);
	// Then can access any element with pointer incrementation:
	// int ch = (band,i)*nwidth*ndepth + (width,j)*ndepth + (depth,k);
	// double etot = *(p+ch); // = esum[i][j][k]

	int ieta = (t_ieta>0 ? t_ieta-1 : t_ieta);
	
	// Fill mean profile
	for (int ch = 0; ch != nch; ++ch) {
	  double rc = (*(p+ch)) / (t_p - t_eMipDR);
	  p2->Fill(ieta, ch, rc);
	  if (jentry%2==0) p2_odd->Fill(ieta, ch, rc);
	  if (jentry%2==1) p2_evn->Fill(ieta, ch, rc);
	} // for ch
      
	// Fill covariance profile
	for (int ch1 = 0; ch1 != nch; ++ch1) {
	  for (int ch2 = 0; ch2 != nch; ++ch2) {
	    double cov = (*(p+ch1)) * (*(p+ch2)) / pow(t_p - t_eMipDR,2);
	    p3->Fill(ieta, ch1, ch2, cov);
	    if (jentry%2==0) p3_odd->Fill(ieta, ch1, ch2, cov);
	    if (jentry%2==1) p3_evn->Fill(ieta, ch1, ch2, cov);
	  } // for ch2
	} // for ch1
      } // response window
	
   } // for jentry

   fout->Write();
   fout->Close();

   cout << "\nClosed output file IsoTrackV2.root" << endl << flush;
   cout << "IsoTrackV2::Loop() finished.\n" << endl << flush;
} // Loop


// Copied from CalibCorr.C from Sunanda Banerjee
// type == 8: Run3 MAHI. Deleted all the rest of the code
//double puFactor(int type, int ieta, double pmom, double eHcal, double ediff, bool debug = false) {
double puFactorCopyV2(int ieta, double pmom, double eHcal, double ediff, bool debug) {
  
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

double _isotrackV2[59][7]; bool _isotrackV2_isset(false);
double getIsoTrackV2Corr(int run, int ieta, int depth) {

  // Load corrections from file. Make smarter and using run range in the future
  if (_isotrackV2_isset==false) {
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
	_isotrackV2[i][j] = (corr>0.1 ? corr : 1.0);
      } // for ieta
    } // for depth
    f->Close();
    _isotrackV2_isset = true;

    // Print corrections
    cout << endl << "  double isotrackV2[59][7] = {\n";
    for (int i = 0; i != 59; ++i) {
      cout << "    {";
      for (int j = 0; j != 7; ++j) {
	double corr = _isotrackV2[i][j];
	cout << Form("%s%1.3f", j==0 ? "" : ", ", corr);
	if (_h2prevcorrV2) _h2prevcorrV2->SetBinContent(i+2, j+2, corr);
      } // for j
      cout << "}" << (i==29 ? "\n" : ",\n");
    } // for i
    cout << "  };\n" << endl;
  } //_isotrackV2_isset==false

  if (ieta<-29 || ieta>29 || ieta==0) return 1.0;
  int i = (ieta+29); i = max(0,min(58,i));
  int j = (depth-1); j = max(0,min(6,j));

  return _isotrackV2[i][j];
} // getIsoTrackCorrV2
