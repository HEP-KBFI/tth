void ratio() {

  TFile *f = new TFile("skim_TTH_HToTauTau_M-125.root");
  //f.ls();

  TH1F * h0 = (TH1F*)f->Get("/demo/h0");
  TH1F * h1 = (TH1F*)f->Get("/demo/h1"); 
  //TH1F * h2 = (TH1F*)f->Get("/demo/h2");
  //TH1F * h3 = (TH1F*)f->Get("/demo/h3"); 

  TLegend *leg=new TLegend(0.5,0.68,0.8,0.99);
  leg->SetTextSize(0.02);

  TCanvas *c1 = new TCanvas("c1","hadronic tau eff.",600,700);
  c1->cd();

  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);  //x1,y1, x2, y2

  //pad1->SetBottomMargin(0);
  pad1->Draw();
  pad1->cd();


  h0->SetStats(0);
  h1->SetStats(0);
  //h2->SetStats(0);
  //h3->SetStats(0);
  h0 ->SetTitle(0);
  h1 ->SetTitle(0);
  //h2 ->SetTitle(0);
  //h3 ->SetTitle(0);
  h0->SetXTitle("tau pt");
  h1->SetXTitle("tau pt");
  //h2->SetXTitle("tau pt");
  //h3->SetXTitle("tau pt");

  h0->SetLineColor(1);
  h0->DrawCopy("");
  h1->SetLineColor(2);
  h1->DrawCopy("same");
  //h2->SetLineColor(3);
  //h2->DrawCopy("same");
  //h3->SetLineColor(4);
  //h3->DrawCopy("same");
  leg->AddEntry(h0, "RecoMatchGen: pt>20, |n|<0.3","lp");
  leg->AddEntry(h1,"LooseIsoDBCorr3Hits ","lp");
  //leg->AddEntry(h2,"MediumIsoDBCorr3Hits","lp");
  //leg->AddEntry(h3,"TightIsoDBCorr3Hits","lp");
  leg->Draw("same");

  c1->cd();

  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
  pad2->SetTopMargin(0);
  pad2->Draw();
  pad2->cd();

  h0->Sumw2();
  h1->Sumw2();
  //h2->Sumw2();
  //h3->Sumw2();
  
  h1->Divide(h0);
  //h2->Divide(h0);
  //h3->Divide(h0);

  //h1->SetTitle("Bin by Bin Ratio of pt (h2/h1)");
  h1->Draw("ep");
  //h2->Draw("same");
  //h3->Draw("same");
}
