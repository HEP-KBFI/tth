#include "THStack.h"
#include "norm.h"

void plot_HN(){
  TString hn  = "HN";
  TString hna = "HNa";
  TString hnb = "HNb";
  TString hnc = "HNc";
  TString hnd = "HNd";

  int reBin = 2; int fillStyle = 3001;
  const char *hNameX = "XNAME";
  const char *hNameY = "n events";
  gROOT->ProcessLine(".L func.C");
  func();

  const int nf = 20;
  TFile* f[nf];
  f[0]   = new TFile("../hist/tth125.root"  ,"READ");
  f[1]   = new TFile("../hist/dyjets_m10to50.root"  ,"READ");
  f[2]   = new TFile("../hist/dyjets_m50.root"      ,"READ");
  f[3]   = new TFile("../hist/wjets.root" ,"READ");
  f[4]   = new TFile("../hist/ttjets.root"  ,"READ");
  f[5]   = new TFile("../hist/ttwjets.root" ,"READ");
  f[6]   = new TFile("../hist/ttzjets.root" ,"READ");
  f[7]   = new TFile("../hist/t_t.root"     ,"READ");
  f[8]   = new TFile("../hist/t_tw.root"    ,"READ");
  f[9]   = new TFile("../hist/t_s.root"     ,"READ");
  f[10]  = new TFile("../hist/tbar_t.root"  ,"READ");
  f[11]  = new TFile("../hist/tbar_tw.root" ,"READ");
  f[12]  = new TFile("../hist/tbar_s.root"  ,"READ");
  f[13]  = new TFile("../hist/ww.root"      ,"READ");
  f[14]  = new TFile("../hist/wz.root"      ,"READ");
  f[15]  = new TFile("../hist/zz.root"      ,"READ");
  f[16]  = new TFile("../hist/singleMuRun2012A.root" ,"READ");
  f[17]  = new TFile("../hist/singleMuRun2012B.root" ,"READ");
  f[18]  = new TFile("../hist/singleMuRun2012C.root" ,"READ");
  f[19]  = new TFile("../hist/singleMuRun2012D.root" ,"READ");

  TH1D *ha[nf], *hb[nf], *hc[nf], *hd[nf]; 
  for( int i = 0; i <nf; ++i){
	  //signal region
	  ha[i] = (TH1D*)f[i]->Get(hna); ha[i]->Rebin(reBin); ha[i]->Scale(sh[i]); ha[i]->SetFillStyle(fillStyle); ha[i]->Sumw2();
	  //qcd region 
	  hb[i] = (TH1D*)f[i]->Get(hnb); hb[i]->Rebin(reBin); hb[i]->Scale(sh[i]); hb[i]->SetFillStyle(fillStyle); hb[i]->Sumw2();
	  hc[i] = (TH1D*)f[i]->Get(hnc); hc[i]->Rebin(reBin); hc[i]->Scale(sh[i]); hc[i]->SetFillStyle(fillStyle); hc[i]->Sumw2();
	  hd[i] = (TH1D*)f[i]->Get(hnd); hd[i]->Rebin(reBin); hd[i]->Scale(sh[i]); hd[i]->SetFillStyle(fillStyle); hd[i]->Sumw2();
  }

  //gather sample
  TH1D *htth    = (TH1D*)ha[0]->Clone();   htth->Scale(10.);
  TH1D *hdy     = (TH1D*)ha[1]->Clone();   hdy->Add(ha[2]);
  TH1D *hwjets  = (TH1D*)ha[3]->Clone();
  TH1D *httjets = (TH1D*)ha[4]->Clone();
  TH1D *httV    = (TH1D*)ha[5]->Clone();   httV->Add(ha[6]); 
  TH1D *ht      = (TH1D*)ha[7]->Clone();   for( int i=8; i <=12; ++i) ht->Add(ha[i]);         
  TH1D *hVV     = (TH1D*)ha[13]->Clone();  hVV->Add(ha[14]); hVV->Add(ha[15]); 
  TH1D *hdata   = (TH1D*)ha[16]->Clone();  for( int i=17; i <=19; ++i) hdata->Add(ha[i]);  

  //substract EW bkg from data to get Mult. bkg
  TH1D *hqa     = (TH1D*)hdata->Clone();                                               for( int i=0; i <=15; ++i) hqa->Add(ha[i],-1);
  TH1D *hqb     = (TH1D*)hb[16]->Clone(); for( int i=17; i <=19; ++i) hqb->Add(hb[i]); for( int i=0; i <=15; ++i) hqb->Add(hb[i],-1); 
  TH1D *hqc     = (TH1D*)hc[16]->Clone(); for( int i=17; i <=19; ++i) hqc->Add(hc[i]); for( int i=0; i <=15; ++i) hqc->Add(hc[i],-1);
  TH1D *hqd     = (TH1D*)hd[16]->Clone(); for( int i=17; i <=19; ++i) hqd->Add(hd[i]); for( int i=0; i <=15; ++i) hqd->Add(hd[i],-1); 

  //qcd normalisation !!do not work!!
  //double wq = hqb->Integral()/hqc->Integral(); 
  //hqd->Scale(wq);
  TH1D *hqcd     = (TH1D*)hqb->Clone(); hqcd->Add(hqc);

  htth->SetLineColor(ktth); htth->SetLineWidth(2); hdy->SetFillColor(kdy); hwjets->SetFillColor(kwjets); httjets->SetFillColor(kttjets);
  httV->SetFillColor(kttV); ht->SetFillColor(kt); hVV->SetFillColor(kVV); hdata->SetMarkerStyle(kFullCircle); hqcd->SetFillColor(kqcd);

  THStack *hs = new THStack("hs","stacked histograms"); 
  hs->Add(httV); hs->Add(ht); hs->Add(hVV); hs->Add(httjets); hs->Add(hdy); hs->Add(hwjets); hs->Add(hqcd);

  TCanvas *c = new TCanvas("c", hn, 700, 700); SetCanvas (c, false); //log

  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);  //x1,y1, x2, y2
  pad1->Draw();
  pad1->cd();

  //draw data 1st to set axis
  SetAxesTitle( hdata, hNameX, hNameY, hs->GetMaximum()+100 );  
  hdata->Draw("hist p");
  hs->Draw("hist same");
  hdata->Draw("hist same&p");
  htth->Draw("hist same&l");

  drawLegend( hdata, htth, hdy, hwjets, httjets, httV, ht, hVV, hqcd );
 
  //draw ratio plot
  c->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);  //x1,y1, x2, y2
  pad1->SetTopMargin(0.05);
  pad2->Draw();
  pad2->cd();

  TH1D *hmc = (TH1D*)httV->Clone(); hmc->Add(ht); hmc->Add(hVV); hmc->Add(httjets); hmc->Add(hdy); hmc->Add(hwjets); hmc->Add(hqcd);
  hmc->Divide(hdata);
  hmc->SetMarkerStyle(20);
  hmc->SetMarkerSize(0.7);
  hmc->SetXTitle(hn);
  hmc->SetYTitle("mc/Data");
  hmc->Draw("p&e");

  TString dir = "plot/"; TString fileName = "HN.png"; dir += fileName; c->Print(dir); 

}


