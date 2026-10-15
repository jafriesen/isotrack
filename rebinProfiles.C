// Purpose: quick rebinning of ieta and depth to new histogram,
//          also allowing to reduce profile range
#ifndef __rebinProfiles_C__
#define __rebinProfiles_C__
#include <TFile.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TProfile3D.h>

#include <iostream>
#include <cmath>

void rebinProfileCustom(TProfile* p, TProfile* p_new);
void rebinProfile2DCustom(TProfile2D* p2, TProfile2D* p2_new);
void rebinProfile3DCustom(TProfile3D* p3, TProfile3D* p3_new);

bool doAbsEta(false);
void rebinProfiles(string mode = "wide",
		   string era = "24CDEF", string version="") {

  cout << "Calling rebinProfiles(\""<<mode<<"\")..." << endl << flush;

  const char *ce = era.c_str();
  const char *cv = version.c_str();
  TFile *f = new TFile(Form("IsoTrack_%s_%s.root",cv,ce),"READ");
  assert(f && !f->IsZombie());

  TFile *fout = new TFile(Form("IsoTrack_%s_%s.root",cv,ce),"UPDATE");
  //TFile *fout = new TFile("IsoTrack_wide.root","UPDATE");
  assert(fout && !fout->IsZombie());

  cout << "Input file " << f->GetName() << endl << flush;
  cout << "Output file " << fout->GetName() << " (UPDATE)" << endl << flush;

  // Define new ieta binning
  double vieta_wide[] =
  //{-30.5, -28.5, -24.5, -19.5, -17.5, -15.5, -12.5, -6.5, -0.5,
  //0.5, 6.5, 12.5, 15.5, 17.5, 19.5, 24.5, 28.5, 30.5}; // v1
    {-30.5, -28.5, -24.5, -20.5, -17.5, -15.5, -12.5, -9.5, -6.5, -3.5, -0.5,
     0.5, 3.5, 6.5, 9.5, 12.5, 15.5, 17.5, 19.5, 24.5, 28.5, 30.5}; // v2,4
  const int nieta_wide = sizeof(vieta_wide)/sizeof(vieta_wide[0])-1;
  double vieta_abs[]
    {-30.5, -29.5, -28.5, -27.5, -26.5, -25.5, -24.5, -23.5, -22.5, -21.5,
     -20.5, -19.5, -18.5, -17.5, -16.5, -15.5, -14.5, -13.5, -12.5, -11.5,
     -10.5, -9.5, -8.5, -7.5, -6.5, -5.5, -4.5, -3.5, -2.5, -1.5, -0.5,
     0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5,
     11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5,
     21.5, 22.5, 23.5, 24.5, 25.5, 26.5, 27.5, 28.5, 29.5, 30.5}; // v3
  const int nieta_abs = sizeof(vieta_abs)/sizeof(vieta_abs[0])-1;
  doAbsEta = true; // v3,4 (list -ieta also, but they will be empty

  const char *cm = mode.c_str();
  double *vieta(0);
  int nieta(0);
  if (mode=="wide") { vieta = &(vieta_wide[0]); nieta = nieta_wide; }
  if (mode=="abs")  { vieta = &(vieta_abs[0]);  nieta = nieta_abs; }
  
  //const int nieta = sizeof(vieta)/sizeof(vieta[0])-1;

  // Define new depth binning
  double vdepth[] =
    {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5};
  const int ndepth = sizeof(vdepth)/sizeof(vdepth[0])-1;

  // Process 1D profiles
  TProfile *p, *p_even, *p_odd, *p_new, *p_even_new, *p_odd_new;
  p = (TProfile*)f->Get("pc"); assert(p);
  p_even = (TProfile*)f->Get("pc_even"); assert(p_even);
  p_odd = (TProfile*)f->Get("pc_odd"); assert(p_odd);
  p_new = new TProfile(Form("pc_%s",cm),
		       ";ieta;(t_eHcal-ePU)/(t_p-eMIP) [0.15,1.85]",
  		       nieta, vieta);//61,-30.5,30.5);
  p_even_new = (TProfile*)p_new->Clone(Form("pc_even_%s",cm));
  p_odd_new = (TProfile*)p_new->Clone(Form("pc_odd_%s",cm));
  rebinProfileCustom(p, p_new);
  rebinProfileCustom(p_even, p_even_new);
  rebinProfileCustom(p_odd, p_odd_new);

  // Process 2D profiles
  TProfile2D *p2, *p2_even, *p2_odd, *p2_new, *p2_even_new, *p2_odd_new;
  p2 = (TProfile2D*)f->Get("p2c"); assert(p2);
  p2_even = (TProfile2D*)f->Get("p2c_even"); assert(p2_even);
  p2_odd = (TProfile2D*)f->Get("p2c_odd"); assert(p2_odd);
  p2_new = new TProfile2D(Form("p2c_%s",cm),
			  ";ieta;depth;(t_eHcal-ePU)/(t_p-eMIP)"
			  " [0.15,1.85]", //61,-30.5,30.5, 10,-0.5,9.5);
			  nieta, vieta, ndepth, vdepth);
  p2_even_new = (TProfile2D*)p2_new->Clone(Form("p2c_even_%s",cm));
  p2_odd_new = (TProfile2D*)p2_new->Clone(Form("p2c_odd_%s",cm));
  rebinProfile2DCustom(p2, p2_new);
  rebinProfile2DCustom(p2_even, p2_even_new);
  rebinProfile2DCustom(p2_odd, p2_odd_new);

  // Profess 3D profiles
  TProfile3D *p3, *p3_even, *p3_odd, *p3_new, *p3_even_new, *p3_odd_new;
  p3 = (TProfile3D*)f->Get("p3c"); assert(p3);
  p3_even = (TProfile3D*)f->Get("p3c_even"); assert(p3_even);
  p3_odd = (TProfile3D*)f->Get("p3c_odd"); assert(p3_odd);
  p3_new = new TProfile3D(Form("p3c_%s",cm),";ieta;depth;depth;"
			  "(t_eHcal-ePU)/(t_p-eMIP) [0.15,1.85]",
			  //61,-30.5,30.5, 10,-0.5,9.5, 10,-0.5,9.5);
			  nieta, vieta, ndepth, vdepth, ndepth, vdepth);
  p3_even_new = (TProfile3D*)p3_new->Clone(Form("p3c_even_%s",cm));
  p3_odd_new = (TProfile3D*)p3_new->Clone(Form("p3c_odd_%s",cm));
  rebinProfile3DCustom(p3, p3_new);
  rebinProfile3DCustom(p3_even, p3_even_new);
  rebinProfile3DCustom(p3_odd, p3_odd_new);

  fout->cd();
  p_new->Write(Form("pc_%s",cm), TObject::kOverwrite);
  p_even_new->Write(Form("pc_even_%s",cm), TObject::kOverwrite);
  p_odd_new->Write(Form("pc_odd_%s",cm), TObject::kOverwrite);
  p2_new->Write(Form("p2c_%s",cm), TObject::kOverwrite);
  p2_even_new->Write(Form("p2c_even_%s",cm), TObject::kOverwrite);
  p2_odd_new->Write(Form("p2c_odd_%s",cm), TObject::kOverwrite);
  p3_new->Write(Form("p3c_%s",cm), TObject::kOverwrite);
  p3_even_new->Write(Form("p3c_even_%s",cm), TObject::kOverwrite);
  p3_odd_new->Write(Form("p3c_odd_%s",cm), TObject::kOverwrite);
  fout->Close();

  cout << "Finished rebinProfiles." << endl << flush;
}

void rebinProfileCustom(TProfile* p, TProfile* p_new) {
  // Iterate through the bins of the original TProfile3D
  for (int ix = 1; ix <= p->GetNbinsX(); ix++) {
                
    // Original bin index
    int bin_orig = ix; // p->GetBin(ix, iy);
      
    // Get bin centers from original profile
    double x_center = p->GetXaxis()->GetBinCenter(ix);

    if (doAbsEta) {
      x_center = fabs(x_center);
    }
    
    // Find corresponding bin in the new profile
    int new_ix = p_new->GetXaxis()->FindBin(x_center);
      
    // Get the new bin index
    int bin_new = new_ix; // p_new->GetBin(new_ix, new_iy);
      
    // profile keeps track of sumy, sumwy, sumwy2, sumw2
    // sumw=fArray, sumwy=fBinEntries.fArray, 
    // sumwy2 = fBinSumw2.fArray, sumw2 = fSum2.fArray
    // GetBinContent = sumwy/sumw
    // https://root-forum.cern.ch/t/copy-entries-of-tprofile/11828
      
    // Add content and entries from the original bin to the new one
    (*p_new)[bin_new] = (*p)[bin_orig] + (*p_new)[bin_new]; // sumwy
    (*p_new->GetSumw2())[bin_new] = (*p->GetSumw2())[bin_orig] +
      (*p_new->GetSumw2())[bin_new]; // sumwy2
    p_new->SetBinEntries(bin_new, p->GetBinEntries(bin_orig) + 
			  p_new->GetBinEntries(bin_new)); // sumw
      
    // Copy (if needed) bin sum of weight square
    if (p->GetBinSumw2()->fN > bin_orig) {
      p_new->Sumw2();
      (*p_new->GetBinSumw2())[bin_new] = (*p->GetBinSumw2())[bin_orig] +
	(*p_new->GetBinSumw2())[bin_new]; // sum2
    }
      
    // Accumulate overall profile entries
    p_new->SetEntries(p_new->GetEntries() + p->GetEntries());
  } // for ix
} // rebinProfileCustom

void rebinProfile2DCustom(TProfile2D* p2, TProfile2D* p2_new) {
  // Iterate through the bins of the original TProfile3D
  for (int ix = 1; ix <= p2->GetNbinsX(); ix++) {
    for (int iy = 1; iy <= p2->GetNbinsY(); iy++) {
                
      // Original bin index
      int bin_orig = p2->GetBin(ix, iy);
      
      // Get bin centers from original profile
      double x_center = p2->GetXaxis()->GetBinCenter(ix);
      double y_center = p2->GetYaxis()->GetBinCenter(iy);
      
      if (doAbsEta) {
	x_center = fabs(x_center);
	y_center = fabs(y_center);
      }

      // Find corresponding bin in the new profile
      int new_ix = p2_new->GetXaxis()->FindBin(x_center);
      int new_iy = p2_new->GetYaxis()->FindBin(y_center);
      
      // Get the new bin index
      int bin_new = p2_new->GetBin(new_ix, new_iy);
      
      // profile keeps track of sumw, sumwz, sumwz2, sumw2
      // sumw=fArray, sumwz=fBinEntries.fArray, 
      // sumwz2 = fBinSumw2.fArray, sumw2 = fSum2.fArray
      // GetBinContent = sumwz/sumw
      // https://root-forum.cern.ch/t/copy-entries-of-tprofile/11828
      
      // Add content and entries from the original bin to the new one
      (*p2_new)[bin_new] = (*p2)[bin_orig] + (*p2_new)[bin_new]; // sumwz
      (*p2_new->GetSumw2())[bin_new] = (*p2->GetSumw2())[bin_orig] +
	(*p2_new->GetSumw2())[bin_new]; // sumwz2
      p2_new->SetBinEntries(bin_new, p2->GetBinEntries(bin_orig) + 
			    p2_new->GetBinEntries(bin_new)); // sumw
      
      // Copy (if needed) bin sum of weight square
      if (p2->GetBinSumw2()->fN > bin_orig) {
	p2_new->Sumw2();
	(*p2_new->GetBinSumw2())[bin_new] = (*p2->GetBinSumw2())[bin_orig] +
	  (*p2_new->GetBinSumw2())[bin_new]; // sum2
      }
      
      // Accumulate overall profile entries
      p2_new->SetEntries(p2_new->GetEntries() + p2->GetEntries());
    } // for iy
  } // for ix
} // rebinProfile2DCustom

void rebinProfile3DCustom(TProfile3D* p3, TProfile3D* p3_new) {
  // Iterate through the bins of the original TProfile3D
  for (int ix = 1; ix <= p3->GetNbinsX(); ix++) {
    for (int iy = 1; iy <= p3->GetNbinsY(); iy++) {
      for (int iz = 1; iz <= p3->GetNbinsZ(); iz++) {
                
	// Original bin index
	int bin_orig = p3->GetBin(ix, iy, iz);
                
	// Get bin centers from original profile
	double x_center = p3->GetXaxis()->GetBinCenter(ix);
	double y_center = p3->GetYaxis()->GetBinCenter(iy);
	double z_center = p3->GetZaxis()->GetBinCenter(iz);

	if (doAbsEta) {
	  x_center = fabs(x_center);
	  y_center = fabs(y_center);
	  z_center = fabs(z_center);
	}
	      
	// Find corresponding bin in the new profile
	int new_ix = p3_new->GetXaxis()->FindBin(x_center);
	int new_iy = p3_new->GetYaxis()->FindBin(y_center);
	int new_iz = p3_new->GetZaxis()->FindBin(z_center);

	// Get the new bin index
	int bin_new = p3_new->GetBin(new_ix, new_iy, new_iz);

	// profile keeps track of sumw, sumwt, sumwt2, sumw2
	// sumw=fArray, sumwt=fBinEntries.fArray, 
	// sumwt2 = fBinSumw2.fArray, sumw2 = fSum2.fArray
	// GetBinContent = sumwt/sumw
	// https://root-forum.cern.ch/t/copy-entries-of-tprofile/11828
	
	// Add content and entries from the original bin to the new one
	(*p3_new)[bin_new] = (*p3)[bin_orig] + (*p3_new)[bin_new]; // sumwt
	(*p3_new->GetSumw2())[bin_new] = (*p3->GetSumw2())[bin_orig] +
	  (*p3_new->GetSumw2())[bin_new]; // sumwt2
	p3_new->SetBinEntries(bin_new, p3->GetBinEntries(bin_orig) + 
			      p3_new->GetBinEntries(bin_new)); // sumw
	
	// Copy (if needed) bin sum of weight square
	if (p3->GetBinSumw2()->fN > bin_orig) {
	  p3_new->Sumw2();
	  (*p3_new->GetBinSumw2())[bin_new] = (*p3->GetBinSumw2())[bin_orig] +
	    (*p3_new->GetBinSumw2())[bin_new]; // sum2
	}

	// Accumulate overall profile entries
	p3_new->SetEntries(p3_new->GetEntries() + p3->GetEntries());
      } // for iz
    } // for iy
  } // for ix
} // rebinProfile3DCustom
#endif


