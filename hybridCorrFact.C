// Purpose: hybridize two sets of CorrFact.root,
//          1) first (iov) to give the depth-independent C(ieta) in narrow bins
//          2) second to give relative depth-dependence in wide bins
// => C_depth(ieta) = C_DI,IOV(ieta) * C_DD,ref(|ieta|)/C_DI,ref(|ieta|)
#ifndef __hybridCorrFact_C__
#define __hybridCorrFact_C__
#include "TFile.h"
#include "TH1D.h"

#include <string>
#include <iostream>
#include <cmath>

void hybridCorrFact(string file_out = "rootfiles/CorrFact_hybrid_lxplus_v8_24CDEF.root",
		    string file_iov = "rootfiles/CorrFact_lxplus_v8_24CDEF.root",
		    //string file_ref = "rootfiles/CorrFact_wide_lxplus_v8_24CDEF.root"
		    string file_ref = "rootfiles/CorrFact_abs_lxplus_v8_24CDEF.root"
		    ) {

  cout << "Running hybridCorrFact(\""<<file_out<<",\""<<file_iov<<"\")..."
       <<",\""<<file_ref<<"\")..." << endl << flush;
  
  TDirectory *curdir = gDirectory;

  TFile *fiov = new TFile(Form("%s",file_iov.c_str()),"READ");
  assert(fiov && !fiov->IsZombie());
  TFile *fref = new TFile(Form("%s",file_ref.c_str()),"READ");
  assert(fref && !fref->IsZombie());

  TFile *fout = new TFile(Form("%s",file_out.c_str()),"RECREATE");
  assert(fout && !fout->IsZombie());

  // Depth-independent scales (different eta bin widths, presumably)
  TH1D *hdi_iov = (TH1D*)fiov->Get("h_di"); assert(hdi_iov);
  hdi_iov->SetName("hdi_iov");
  TH1D *hdi_ref = (TH1D*)fref->Get("h_di"); assert(hdi_ref);
  hdi_ref->SetName("hdi_ref");
  
  fout->cd();
  TH1D *hdi_out = (TH1D*)hdi_iov->Clone("h_di");
  hdi_out->Write("h_di");
  
  // Process depths
  for (int id = 1; id != 8; ++id) {
    TH1D *hdd_iov = (TH1D*)fiov->Get(Form("hf_dd_%d",id)); assert(hdd_iov);
    hdd_iov->SetName(Form("hdd_iov_%d",id));
    TH1D *hdd_ref = (TH1D*)fref->Get(Form("hf_dd_%d",id)); assert(hdd_ref);
    hdd_ref->SetName(Form("hdd_ref_%d",id));

    fout->cd();
    TH1D *hdd_out = (TH1D*)hdd_iov->Clone(Form("hf_dd_%d",id));
    hdd_out->Clear();

    for (int i = 1; i != hdd_out->GetNbinsX()+1; ++i) {
      double ix = hdd_out->GetBinCenter(i);
      int i2 = hdi_ref->FindBin(fabs(ix));
      int i3 = hdd_ref->FindBin(fabs(ix));
      if (hdi_ref->GetBinContent(i2)!=0) {
	hdd_out->SetBinContent(i, hdi_iov->GetBinContent(i) *
			       hdd_ref->GetBinContent(i3) /
			       hdi_ref->GetBinContent(i2));
	double edi = hdi_iov->GetBinError(i);
	double edd = hdd_ref->GetBinError(i3);
	double err = sqrt(pow(edi,2)+pow(edd,2));
	hdd_out->SetBinError(i, err);
      } // non-zero bin
    } // for i

    fout->cd();
    hdd_out->Write(Form("hf_dd_%d",id));
  } // for id

  fout->Close();
  fiov->Close();
  fref->Close();
  
  cout << "Done with hybridCorrFact(\""<<file_out<<",\""<<file_iov<<"\")..."
       <<",\""<<file_ref<<"\")." << endl << flush;
  cout << "Output written to " << file_out << "." << endl << flush;
} // hybridCorrFact
#endif
