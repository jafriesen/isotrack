// Purpose: solve HCAL corrections per depth and plot them
#ifndef __drawIsoTrack_C__
#define __drawIsoTrack_C__
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TVectorD.h"
#include "TMatrixD.h"

#include "tdrstyle_mod22.C"

#include <iostream>

bool isValidEta(int ieta, int depth) {
  int absi = abs(ieta);
  if (depth==0 && absi<=29) return true;
  if (depth==1 && (absi<=16 || (absi>=19 && absi<=29))) return true;
  if (depth==2 && absi<=29) return true;
  if (depth==3 && absi<=29) return true;
  if (depth==4 && (absi<=16 || (absi>=18 && absi<=28))) return true;
  if (depth==5 && (absi>=18 && absi<=28)) return true;
  if (depth==6 && (absi>=19 && absi<=28)) return true;
  if (depth==7 && (absi>=26 && absi<=28)) return true;
  if (depth==8 && (absi<=15)) return true;
  if (depth==9 && (absi>=29)) return true;
  if (depth==10 && (absi>=29)) return true;
  return false;
}

void drawIsoTracks(string mode, string era, string version);

void drawIsoTrack(string era="24F", string version="vX", bool doEvenOdd=true) {

  //drawIsoTracks("_even_wide", era, version);
  //drawIsoTracks("_odd_wide", era, version);
  //drawIsoTracks("_wide", era, version);

  if (era=="24CDEFGHI") {
    if (doEvenOdd) {
      drawIsoTracks("_even_abs", era, version);
      drawIsoTracks("_odd_abs", era, version);
    }
    drawIsoTracks("_abs", era, version);
  }

  if (doEvenOdd) {
    drawIsoTracks("_even", era, version);
    drawIsoTracks("_odd", era, version);
  }
  drawIsoTracks("", era, version);

} // drawIsoTrack

void drawIsoTracks(string mode, string era, string version) {

  cout << "Running drawIsoTrack(\""<<mode<<",\""<<era
       << ",\""<<version<<"\")..."<<endl<<flush;
  
  TDirectory *curdir = gDirectory;
  setTDRStyle();

  const char *ce = era.c_str();
  const char *cv = version.c_str();
  //TFile *f = new TFile(Form("IsoTrack%s.root",ce),"READ");
  TFile *f = new TFile(Form("IsoTrack_%s_%s.root",cv,ce),"READ");
  //TFile *f = new TFile("IsoTrack_wide.root","READ");
  //TFile *f = new TFile("IsoTrack_lxplus_v4.root","READ");
  //TFile *f = new TFile("IsoTrack_lxplus_v3_odd.root","READ");
  assert(f && !f->IsZombie());
  curdir->cd();

  cout << "Read inputs from " << f->GetName() << endl << flush;
  
  // Some quick plotting for sanity checks
  //TH2D *h2 = (TH2D*)f->Get("h2c"); assert(h2);
  //h2->Draw("COLZ");
  //h2->GetYaxis()->SetRangeUser(0,2);
  //double ntrk = h2->Integral();

  const char *cm = mode.c_str();
  TH2D *h2prevcorr = (TH2D*)f->Get("h2prevcorr"); // only since v26
  TProfile *p = (TProfile*)f->Get(Form("pc%s",cm));
  TProfile2D *p2 = (TProfile2D*)f->Get(Form("p2c%s",cm));
  TProfile3D *p3 = (TProfile3D*)f->Get(Form("p3c%s",cm));
  assert(p);   if (!p)  exit(1);
  assert(p2);  if (!p2) exit(2);
  assert(p3);  if (!p3) exit(3);

  cout << "Inverting depth-independent corrections" << endl << flush;
  
  TH1D *h = p->ProjectionX(Form("h%s",cm));
  // Reverse to get correction factor
  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    if (h->GetBinContent(i)!=0) {
      double k = h->GetBinContent(i);
      h->SetBinContent(i, 1./k);
      h->SetBinError(i, h->GetBinError(i)/k);
    }
  }

  // Calculate corrections per depth
  cout << "Calculating depth-dependent corrections" << endl << flush;
    
  const int ndepth = 10;
  //int color[ndepth] =
  //{kBlack, kRed, kGreen+2, kBlue, kMagenta+1,
  // kOrange+1,kCyan+2,kGray+2, 1,1};
  int color[ndepth] =
    {kBlack, kBlue, kOrange+1, kGreen+1, kRed,
     kYellow+1,kOrange+3,kGray+2, 1,1};
  int marker[ndepth] =
    {1, kOpenTriangleDown, kOpenSquare, kOpenCircle, kOpenTriangleUp,
     kOpenStar,kOpenDiamond,kOpenCross, 1,1};
  //int size[ndepth]
  vector<TH1D*> vhd(10,0);
  vector<TH1D*> vhdf(10,0);
  vector<TH1D*> vhrms(10,0);
  vector<TH1D*> vhmean(10,0);
  for (int ieta = 1; ieta != h->GetNbinsX()+1; ++ieta) {
    if (h->GetBinContent(ieta)!=0) {
      for (int i = 0; i != ndepth; ++i) { // loop depths (include ECAL: i==0)

	if (vhd[i]==0) {
	  vhd[i] = p->ProjectionX(Form("hd_ieta%d_%d%s",ieta,i,cm));
	  vhd[i]->Reset();
	  vhdf[i] = p->ProjectionX(Form("hdf_ieta%d_%d%s",ieta,i,cm));
	  vhdf[i]->Reset();
	  vhrms[i] = p->ProjectionX(Form("hrms_ieta%d_%d%s",ieta,i,cm));
	  vhrms[i]->Reset();
	  vhmean[i] = p->ProjectionX(Form("hmean_ieta%d_%d%s",ieta,i,cm));
	  vhmean[i]->Reset();
	}
	TH1D *hd = vhd[i];

	// Optimization problem to solve to get the ideal correction per layer:
	// minimize sum_i C_i^2*sigma_i^2 given sum_i C_i*mu_i = p,track
	// using the method of Lagrange multipliers
	// => C_i = p,track * ( (mu_i/sigma_i^2) / (sum_j mu_j^2/sigma_j^2) )
	//
	// Above neglects the covariances arising from shower fluctuations
	// vec_c = 1/(vec_mu_T inv_mat_Sigma vec_mu) * inv_mat_Sigma vec_mu
	
	double mi = p2->GetBinContent(ieta, i+1); // Mean
	p2->SetErrorOption("s"); // return RMS
	double si = p2->GetBinError(ieta, i+1); // RMS
	p2->SetErrorOption("e"); // return error on mean
	double ei = p2->GetBinError(ieta, i+1); // error
	double mi1si2 = (si!=0 ? mi/(si*si) : 0.);
	double ei1si2 = (si!=0 ? ei/(si*si) : 0.);

	double mj2sj2(0);
	double vmj[ndepth], vsj2[ndepth], vej[ndepth];
	for (int j = 1; j != ndepth; ++j) { // loop over depths again (exclude ECAL: j!=0)
	  double mj = p2->GetBinContent(ieta, j+1); // Mean
	  p2->SetErrorOption("s"); // return RMS
	  double sj = p2->GetBinError(ieta, j+1); // RMS
	  p2->SetErrorOption("e"); // return error on mean
	  double ej = p2->GetBinError(ieta, j+1); // error
	  mj2sj2 += (sj!=0 ? mj*mj/(sj*sj) : 0.);

	  vmj[j] = mj;
	  vsj2[j] = sj*sj;
	  vej[j] = ej;
	}
	
	// Full error propagation
	// df_j = dC_i/dmu_j
	// var = sum_i dfi*dfj*epsi*epsj
	double df[ndepth];
	double mi2 = mi*mi;
	double si2 = si*si;
	double si4 = si2*si2;
	double err2(0);
	for (int j = 1; j != ndepth; ++j) {
	  double mj = vmj[j];
	  double sj2 = vsj2[j];
	  double ej = vej[j];
	  if (si2==0||sj2==0||mj2sj2==0) df[j] = 0;
	  else if (i==j) df[j] = (1./si2)/mj2sj2 - (2*mi2/si4)/(mj2sj2*mj2sj2);
	  else           df[j] = -(mi/si2)*(2*mj/sj2)/(mj2sj2*mj2sj2);

	  // ignoring off-diagonal covariances here
	  err2 += df[j]*df[j]*ej*ej;
	} // for depth2
	
	if (mj2sj2!=0) {
	  double corr = mi1si2 / mj2sj2;
	  hd->SetBinContent(ieta, corr);
	  //double err = ei1si2 / mj2sj2; // earlier approximation
	  double err = sqrt(fabs(err2)); // full calculation
	  hd->SetBinError(ieta, err);

	  vhrms[i]->SetBinContent(ieta, (mi!=0 ? si/mi * 0.5: 0.));
	  vhmean[i]->SetBinContent(ieta, mi);
	}
	if (i==0 && p->GetBinContent(ieta)!=0) {
	  double mu = p->GetBinContent(ieta);
	  vhmean[0]->SetBinContent(ieta, mu/4.);
	  p->SetErrorOption("s"); // return RMS
	  vhrms[0]->SetBinContent(ieta, p->GetBinError(ieta)/mu);
	  p->SetErrorOption("e"); // return error
	}

      } // for depth1

      /////////////////////////////////////////////////////
      // Account for covariances between depths properly //
      /////////////////////////////////////////////////////
      //double vc[ndepth];
      double vm[ndepth];
      double ve[ndepth];
      double as[ndepth][ndepth];
      double ae[ndepth][ndepth];

      // Save means vec_mu and it's error
      for (int i = 1; i != p2->GetNbinsY()+1; ++i) {
      	vm[i-1] = p2->GetBinContent(ieta, i);
	p2->SetErrorOption("e");
	ve[i-1] = p2->GetBinError(ieta, i);
      }
	   
      // Covariance matrix mat_Sigma for improved method
      // (Draw this with drawCovariance.C)
      p3->GetXaxis()->SetRange(ieta, ieta);
      TProfile2D *p2p = (TProfile2D*)p3->Project3DProfile("yz");
      p2p->SetName(Form("p3yz_%d%s",ieta,cm));
      TH2D *h2 = p2p->ProjectionXY(); h2->Reset();
      for (int i = 1; i != p2p->GetNbinsX()+1; ++i) {
	for (int j = 1; j != p2p->GetNbinsY()+1; ++j) {
	  double xy = p2p->GetBinContent(i, j);
	  double xx = p2p->GetBinContent(i, i);
	  double yy = p2p->GetBinContent(j, j);
	  if (xx!=0 && yy!=0)
	    h2->SetBinContent(i, j, xy/sqrt(xx*yy));
	  else
	    h2->SetBinContent(i, j, 0.);

	  // Var[x] = 1/N * sum_i (xi - mu_x)^2 = E[x^2]-E[x]^2
	  // Cov[x,y] = 1/N * sum_i (xi - mu_x) * (yi - mu_y)
	  //          = E[xy] - E[x]*E[y]
	  p2->SetErrorOption("s");
	  double Ex = p2->GetBinContent(ieta, i);
	  double Ey = p2->GetBinContent(ieta, j);
	  double Exy = p2p->GetBinContent(i, j);
	  as[i-1][j-1] = Exy - Ex*Ey;
	} // for j
      } // for i

      // Set vector length to active depths and exclude ECAL
      int ndepth2 = 4; //
      int jeta = p2->GetXaxis()->GetBinCenter(ieta);
      /*
      // Depths for original case (and Sunanda) with wide cone
      if (fabs(jeta)==15) ndepth2 = 5;
      if (fabs(jeta)>=16) ndepth2 = 6;
      if (fabs(jeta)>=23) ndepth2 = 7;
      */
      // Depths for new 3x5
      if (fabs(jeta)==17) ndepth2 = 5;
      if (fabs(jeta)>=18) ndepth2 = 6;
      //if (fabs(jeta)>=18) ndepth2 = 6; // just duplicate line?
      if (fabs(jeta)>=25) ndepth2 = 7;

      TVectorD vec_c(ndepth2);
      TVectorD vec_err(ndepth2);
      TVectorD vec_mu(ndepth2);
      TMatrixD mat_Sigma(ndepth2, ndepth2);
      for (int i = 0; i != ndepth2; ++i) {
	vec_mu[i] = vm[i+1]; // exclude ECAL (+1)
	for (int j = 0; j != ndepth2; ++j) {
	  mat_Sigma[i][j] = as[i+1][j+1]; // exclude ECAL (+1)
	  //if (i!=j) mat_Sigma[i][j] = 0; // diagonal case as sanity check
	} // for j
      } // for i

      // Core linear algebra calculations
      TMatrixD inv_mat_Sigma = mat_Sigma;
      inv_mat_Sigma.Invert(); // In-place inversion
      // Check if inversion succeeded
      if (inv_mat_Sigma.IsValid()) {
	// Calculate corrections using:
	// vec_c = 1/(vec_mu_T inv_mat_Sigma vec_mu) * inv_mat_Sigma vec_mu
	double muT_invSigma_mu = vec_mu * (inv_mat_Sigma * vec_mu);
	TVectorD invSigma_mu = inv_mat_Sigma * vec_mu;
	vec_c = (1. / muT_invSigma_mu) * invSigma_mu;

	// Calculate the uncertainties for the corrections
	// Jacobian is covInverse / muT_invSigma_mu
	TMatrixD jacobian = inv_mat_Sigma;//covInverse;
	jacobian *= (1.0 / muT_invSigma_mu);

	// Propagate the uncertainties: vc_err = sqrt(J * covariances * J^T)
	TMatrixD tempMatrix = jacobian * mat_Sigma * jacobian.T();
	for (int i = 0; i != ndepth2; ++i) {
	  vec_err[i] = std::sqrt(tempMatrix(i, i)) / muT_invSigma_mu;
	}

	// Save corrections and uncertainties
	for (int i = 0; i != ndepth2; ++i) {
	  TH1D *hd = vhdf[i+1];
	  double prevcorr(1);
	  if (h2prevcorr) {
	    prevcorr = h2prevcorr->GetBinContent(ieta, i+2);
	    if (prevcorr==0) prevcorr = 1; // safety for empty h2prevcorr
	  }
	  hd->SetBinContent(ieta, vec_c[i] * prevcorr);
	  hd->SetBinError(ieta, vec_err[i] * prevcorr);
	}

      }
      else {
	cout << "Matrix inversion failed for ieta = " << ieta << endl;
      }

    } // non-empty bin
  } // for ieta

  cout << "Making plots of corrections factors" << endl << flush;

  TH1D *h_1 = tdrHist(Form("h_1%s",cm),"Correction Factor",0.15,2.7,
		      "i#eta",-29,29);
  lumi_136TeV = Form("20%s EGamma",ce);
  extraText = "Private";
  TCanvas *c1 = tdrCanvas(Form("c1%s",cm),h_1,8,11,kRectangular);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kBlue);
  //l->DrawLine(-25,1,25,1);
  l->DrawLine(-29,1,29,1);
  l->SetLineStyle(kDotted);
  l->SetLineColor(kGray+1);
  double y1 = h_1->GetMinimum();
  double y2 = h_1->GetMaximum();
  l->DrawLine(-15.5,y1,-15.5,y2);
  l->DrawLine(+15.5,y1,+15.5,y2);
  l->SetLineColor(kGray+1);
  l->DrawLine(-17.5,y1,-17.5,y2);
  l->DrawLine(+17.5,y1,+17.5,y2);

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);

  TLegend *leg = tdrLeg(0.40,0.90-5*0.05,0.65,0.90);
  TLegend *leg2 = tdrLeg(0.40,0.15,0.65,0.15+3*0.05);
  
  leg->AddEntry(h,"Depth independent","PLE");
  for (int i = 1; i != 8; ++i) {
    TH1D *hd = vhdf[i]; // vhd[i]
    // Compare to free depths: https://indico.cern.ch/event/1424547/contributions/5991794/attachments/2872957/5030614/IsoTrackN156.pdf?#page=25
    tdrDraw(hd,"Pz",marker[i],color[i], kSolid,-1,kNone,0, 0.8);
    if (i<5) leg->AddEntry(hd,Form("DD depth %d",i),"PLE");
    else     leg2->AddEntry(hd,Form("DD depth %d (HE)",i),"PLE");
  }
  
  // Compare to single depth: https://indico.cern.ch/event/1460156/contributions/6147826/attachments/2935195/5155251/IsoTrackN163.pdf?#page=9
  tdrDraw(h,"Pz",kFullCircle,kBlack, kSolid,-1,kNone,0, 0.8);

  c1->SaveAs(Form("pdf/IsoTrack_CorrFact%s_%s_%s.pdf",cm,cv,ce));


  // Loook at means and RMS per depth vs total
  TH1D *h_2 = tdrHist(Form("h_2%s",cm),"Mean and RMS/Mean",0.,1.5,
		      "i#eta",-29.5,29.5);
  lumi_136TeV = Form("20%s EGamma Depth Dependent",ce);
  extraText = "Private";
  TCanvas *c2 = tdrCanvas(Form("c2%s",cm),h_2,8,11,kRectangular);

  l->SetLineStyle(kDashed);
  l->SetLineColor(kBlue);
  l->DrawLine(-25,1,25,1);

  for (int i = 0; i != 5; ++i) {
    tdrDraw(vhmean[i],"Pz",kOpenSquare, color[i], kSolid,-1,kNone,0, 0.8);
    tdrDraw(vhrms[i],"Pz",kFullCircle, color[i], kSolid,-1,kNone,0, 0.8);
  }

  c2->SaveAs(Form("pdf/IsoTrack_RMSoMean%s_%s_%s.pdf",cm,cv,ce));

  
  // Store results for e.g. comparing even-odd cases or different eras
  TFile *fout = new TFile(Form("rootfiles/CorrFact%s_%s_%s.root",cm,cv,ce),
			  "RECREATE");
  cout << "Storing results to " << fout->GetName() << endl << flush;
  
  h->Write("h_di");
  for (int i = 0; i != ndepth; ++i) {
    vhd[i]->Write(Form("h_dd_%d",i));
  }
  for (int i = 0; i != ndepth; ++i) {
    vhdf[i]->Write(Form("hf_dd_%d",i));
  }
  fout->Close();

  // Store results also in jecsys3/minitools/HcalDepthsFromIsoTrack.C style
  if (true) {
    
    TFile *fout2 = new TFile(Form("rootfiles/HcalDepthsFromIsoTrackOrig_%s_%s.root",cv,ce),
			     "RECREATE");
    cout << "Storing results to " << fout2->GetName() << endl << flush;
    
    int ietas[] = { // ieta
      -41, -40, -39, -38, -37, -36, -35, -34, -33, -32, -31, -30, -29, -28, -27, -26, -25, -24, -23, -22, -21, -20, -19, -18, -17, -16, -15, -14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41
    };
    const int n_ieta = sizeof(ietas)/sizeof(ietas[0]);
    //assert(n_eta==n_ieta); // n_eta is undefined in this repo (leftover check); n_ieta is 82 by construction
    double vieta[n_ieta+1];
    //for (int i = 0; i != n_ieta+2; ++i) { vieta[i] = ietas[i]; }
    for (int i = 0; i != n_ieta/2; ++i) { vieta[i] = ietas[i]-0.5; }
    vieta[n_ieta/2] = 0.;
    for (int i = n_ieta/2; i != n_ieta+1; ++i) { vieta[i+1] = ietas[i]+0.5; }

    // Depths (0=ECAL, 1-7=HCAL, 8=HO, 9=hfLong, 10=hfShort)
    int num_depths = 9;
    std::string depth_str[11] = {"0", "1", "2", "3", "4", "5", "6", "7", "8"};// "hfLong", "hfShort"};
    const char* label[11] = {"ECAL","Depth 1","Depth 2","Depth 3","Depth 4","5 (HE)","6 (HE)","7 (HE)","HO","HF Long","HF Short"};
    int colors[11] = {kCyan+1, kBlue, kOrange+1, kGreen+1, kRed, kYellow+1, kOrange+3, kGray+2, kGray+1, kBlue, kRed};
    int markers[11] = {20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30};
    
    for (int d = 0; d != num_depths; ++d) {
      
      TH1D *h = vhdf[d];
      TH1D *hnew = new TH1D(Form("h_depth_%d",d),Form(";Track i#eta;Correction factor (%s)",label[d]),n_ieta,vieta);
      for (int i = 1; i != hnew->GetNbinsX()+1; ++i) {
	double ieta = hnew->GetBinCenter(i);
	int j = h->GetXaxis()->FindBin(ieta);
	if (j>=1 && j<=h->GetNbinsX() && isValidEta(ieta,d)) {
	  
	  hnew->SetBinContent(i, h->GetBinContent(j));
	  hnew->SetBinError(i, h->GetBinError(j));
	}
      } // for i

      hnew->SetLineColor(colors[d]);
      hnew->SetMarkerColor(colors[d]);
      hnew->SetMarkerStyle(markers[d]);
      hnew->Write(Form("h_depth_%d",d));
    }
    fout2->Close();
  }
  
  cout << "File closed, finished drawIsoTracks(\""<<mode<<",\""
       <<era<<"\",\""<<version<<"\").\n"<<endl<<flush;
} // drawIsoTracks(mode)
#endif
