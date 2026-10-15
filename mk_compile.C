// Run this script to compile CondFormats libraries. After this can easily run 
// root -l -b -q mk_GamHistosFill.C
// using R__LOAD_LIBRARY to load *.so
{

  // For IsoTrack code (v6.30/04)
  gROOT->ProcessLine(".L IsoTrack.C+g");
  gROOT->ProcessLine(".L drawIsoTrack.C+g");
  gROOT->ProcessLine(".L compareCorrFact.C+g");
  gROOT->ProcessLine(".L drawPileup.C+g");
  gROOT->ProcessLine(".L rebinProfiles.C+g");
  gROOT->ProcessLine(".L drawCovariance.C+g");
  gROOT->ProcessLine(".L hybridCorrFact.C+g");

  gROOT->ProcessLine(".L IsoTrackV2.C+g");

}
