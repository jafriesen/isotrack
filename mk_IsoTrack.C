#include "IsoTrack.h"
#include "rebinProfiles.C"
#include "drawIsoTrack.C"
#include "drawCorrFact.C"
#include "drawCovariance.C"
#include "hybridCorrFact.C"
#include "drawPileup.C"
// Forward declarations (dont' work here?)
//void rebinProfiles();
//void drawIsoTrack();
//void drawCorrFact();
//void drawCovariance();
//void hybridCorrFact();
//void drawPileup();

//R__LOAD_LIBRARY(IsoTrack.C+g)
R__LOAD_LIBRARY(IsoTrack_C.so)
R__LOAD_LIBRARY(drawIsoTrack_C.so)
R__LOAD_LIBRARY(compareCorrFact_C.so)
R__LOAD_LIBRARY(drawPileup_C.so)
R__LOAD_LIBRARY(rebinProfiles_C.so)
R__LOAD_LIBRARY(drawCovariance_C.so)
R__LOAD_LIBRARY(hybridCorrFact_C.so)

void mk_IsoTrack(string era = "24F", string version = "local", string filelist = "") {

  string path = gSystem->pwd();
  cout << "Current path: " << path << endl << flush;
  //bool runLXPLUS = (path=="/afs/cern.ch/user/v/voutila/scratch0/IsoTrack");
  bool runLXPLUS = TString(path.c_str()).Contains("/afs/cern.ch/");

  cout << "Era "<<era<<", version "<<version<<endl<<flush;

  TString t(era.c_str());
  TChain *c(0);
  if (t.Contains("MC"))
    c = new TChain("hcalIsoTrkAnalyzer/CalibTree");
  else
    c = new TChain("hcalIsoTrackAnalyzer/CalibTree");
  
  const char *ce = era.c_str();
  const char *cv = version.c_str();
  
  if (filelist!="") {
    // Explicit file list (used by the automator; works outside /afs)
    cout << "Reading input files from " << filelist << "..." << endl << flush;
    ifstream fin(filelist.c_str());
    if (!fin.is_open()) {
      cout << "ERROR: cannot open file list " << filelist << endl << flush;
      exit(1);
    }
    string s;
    while (fin >> s) c->AddFile(s.c_str());
  }
  else if (runLXPLUS) {
    // Copy files to lxplus and back:
    // rsync -rutP *.C *.h 24FEB.txt voutila@lxplus.cern.ch:~/scratch0/IsoTrack/
    // rsync -rutP voutila@lxplus.cern.ch:~/scratch0/IsoTrack/IsoTrack.root ./IsoTrack_lxplus.root
    cout << "Running on lxplus..." << endl << flush;
    ifstream fin(Form("%sE%s.txt",ce,era=="24F"?"B":"A"));
    string s;
    while (fin >> s) c->AddFile(s.c_str());
  }
  else {
    cout << "Running on a local test file..." << endl << flush;
    c->AddFile("../data/IsoTrack/EGamma0_40to60_v4.root");
  }

  // Create IsoTrack.root output file for interactive analysis (lxplus only)

  // Uncomment on lxplus

  IsoTrack it(c);
  //if (runLXPLUS) // don't redo histograms if running locally
  it.Loop(); // uncomment on lxplus
  //exit(0);

  // Copy this file over to era+version and backup copy
  // uncomment on lxplus
  gROOT->ProcessLine(Form(".! cp -p IsoTrack.root IsoTrack_%s_%s.root",cv,ce));

  gROOT->ProcessLine(Form(".! cp -p IsoTrack_%s_%s.root rootfiles/"
  			  "IsoTrack_%s_%s_orig.root",cv,ce,cv,ce));

  // Rebin ieta for more stable depths
  //gROOT->ProcessLine(".L rebinProfiles.C+g");
  rebinProfiles("abs",era,version);
  //rebinProfiles("wide",era,version);
  
  // Solve corrections, IsoTrack.root -> CorrFact*.root
  gROOT->ProcessLine(".! mkdir pdf");
  gROOT->ProcessLine(".! mkdir rootfiles");
  //gROOT->ProcessLine(".L drawIsoTrack.C+g");
  //drawIsoTracks("", era, version);

  //gROOT->ProcessLine(".L drawPileup.C+g");
  drawPileup();

  
  // Even-odd
  //drawIsoTracks("_even", era, version);
  //drawIsoTracks("_odd", era, version);
  //drawIsoTracks("", era, version);
  
  //exit(0);
  // draw more variants ((wide,abs) x (odd,even)
  drawIsoTrack(era,version); // crashing?

  // Even/odd tests of uncertainties for CorrFact*.root
  //gROOT->ProcessLine(".L drawCorrFact.C+g");
  drawCorrFact("",era,version);
  //drawCorrFact("_wide",era,version);

  // Control plots of depth covariance in IsoTrack.root
  gROOT->ProcessLine(".! mkdir pdf/drawCovariance");
  gROOT->ProcessLine(Form(".! mkdir pdf/drawCovariance/%s_%s",cv,ce));
  //gROOT->ProcessLine(".L drawCovariance.C+g");
  drawCovariance(0,era,version);

  // Hybridize correction: use per-depth corrections vs |eta| for x2 statistics,
  // and depth-independent for asymmetry (and time-dependence)
  if (era=="24CDEFGHI") {
    //gROOT->ProcessLine(".L hybridCorrFact.C);
    hybridCorrFact("rootfiles/CorrFact_hybrid_lxplus_v32_24CDEFGHI.root",
		   "rootfiles/CorrFact_lxplus_v32_24CDEFGHI.root",
		   "rootfiles/CorrFact_abs_lxplus_v32_24CDEFGHI.root");
    hybridCorrFact("rootfiles/CorrFact_even_hybrid_lxplus_v32_24CDEFGHI.root",
		   "rootfiles/CorrFact_even_lxplus_v32_24CDEFGHI.root",
		   "rootfiles/CorrFact_even_abs_lxplus_v32_24CDEFGHI.root");
    hybridCorrFact("rootfiles/CorrFact_odd_hybrid_lxplus_v32_24CDEFGHI.root",
		   "rootfiles/CorrFact_odd_lxplus_v32_24CDEFGHI.root",
		   "rootfiles/CorrFact_odd_abs_lxplus_v32_24CDEFGHI.root");
    
    string version2 = "hybrid_"+version;
    drawCorrFact("",era,version2);
  }
} // mk_IsoTrack
