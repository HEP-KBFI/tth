// style for plots for CMS - to be addapted to official style

void func(){
  //use in a ROOT macro like this://
  //gROOT->ProcessLine(".L ~/your_path/CMSStyle.C");
  //CMSStyle();

  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetLabelFont(102,"xyz");
  gStyle->SetTitleFont(102);
  gStyle->SetTitleFont(102,"xyz");
  gStyle->SetStatFont(102);
  gStyle->SetOptFit(111);
  gStyle->SetOptTitle(0);
  gStyle->SetPaintTextFormat(".3f");

  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetCanvasColor(kWhite);  
  gStyle->SetOptStat(0); //don't show statistics box
  gStyle->SetTitleSize(0.05, "xyz"); 
  gStyle->SetLegendBorderSize(1);

  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadRightMargin(0.05);

  gROOT->ForceStyle();

  TH1::SetDefaultSumw2(kTRUE);	
}

/***********************************************************************/
/*   Useful functions:                                                 */
/*       set axes title for a histogram                                */
/*                                                                     */
/***********************************************************************/
/*
void SetAxesTitle(TH1 *hist, const char* titlex, const char* titley, int range1, int range2, double offset=1.3 ){
  hist->GetXaxis()->SetTitle(titlex);
  hist->GetYaxis()->SetTitle(titley);
  hist->GetXaxis()->SetRangeUser(range1, range2);
  hist->GetYaxis()->SetTitleOffset(offset);
}
*/
void SetAxesTitle(TH1D *h, const char* tx, const char* ty, double max, double os=1.1 ){
  h->GetXaxis()->SetTitle(tx);
  h->GetYaxis()->SetTitle(ty);
  h->SetMaximum(max);
  h->GetYaxis()->SetTitleOffset(os);
}

void SetCanvas( TCanvas *c , bool log=false, double lM=0.13, double rM=0.05, double tM=0.04) { // left, right, top margin
c->SetLeftMargin(lM);
c->SetRightMargin(rM);
c->SetTopMargin(tM);
if(log) c->SetLogy();
}

void drawLegend( TH1D *hd, TH1D *htth, TH1D *hdy, TH1D *hwjets, TH1D *httjets, TH1D *httV, TH1D *ht, TH1D *hVV, TH1D *hqcd ){
  //TLegend *legend = new TLegend( 0.5, 0.53, 0.95, 0.97, NULL,"brNDC");
  TLegend *legend = new TLegend( 0.8, 0.55, 1., 1., NULL,"brNDC");
  legend->SetTextAlign(12); legend->SetTextSize(0.03); legend->SetTextFont(102);
  //legend->SetHeader("CMS work in progress");

  TString c[9];
  c[0]   = "Data 18.5 fb^{-1}";
  c[1]   = "tth125*10";
  c[2]   = "#gamma*/Z^{0}#rightarrowl^{+}l^{-}";
  c[3]   = "wjets#rightarrowl#nu";
  c[4]   = "t#bar{t}jets";
  c[5]   = "t#bar{t}+V";
  c[6]   = "t,#bar{t}";
  c[7]   = "VV";
  c[8]   = "qcd";

  TLegendEntry* entry;
  entry = legend->AddEntry( hd,     c[0], "lp");
  entry = legend->AddEntry( htth,   c[1], "l"); 
  entry = legend->AddEntry( hdy,    c[2], "f"); 
  entry = legend->AddEntry( hwjets, c[3], "f"); 
  entry = legend->AddEntry( httjets,c[4], "f"); 
  entry = legend->AddEntry( httV,   c[5], "f"); 
  entry = legend->AddEntry( ht,     c[6], "f"); 
  entry = legend->AddEntry( hVV,    c[7], "f"); 
  entry = legend->AddEntry( hqcd,   c[8], "f"); 

  legend->Draw();
}



/***********************************************************************/
/*                                                                     */
/*       draw legend for 1 TGraph or TGraphErrors                      */
/*                                                                     */
/***********************************************************************/
void DrawLegend(TGraph *histo1,
		const Char_t *header,
		const Char_t *label1){
  TLegend *legend = new TLegend(0.553482,0.653319,0.950072,0.948683,
				NULL,"brNDC");
  legend->SetTextAlign(22);/*this is for centering the header*/
                           /*(see root documentation, TAttText constructor)*/
  legend->SetTextFont(102);/*font for the header: 102 -courier bold*/
  legend->SetHeader(header);  /*leg. header (name of histogram or sth. else)*/
  legend->SetTextFont(102);/*font for other legend's entries*/
  
  TLegendEntry* entry1 = legend->AddEntry(histo1,label1,"LP");  
  entry1->SetTextColor(histo1->GetLineColor());
  
  legend->SetFillColor(kWhite);
  legend->Draw(); 
}

/***********************************************************************/
/*                                                                     */
/*       draw legend for 1 histo                                       */
/*                                                                     */
/***********************************************************************/
void DrawLegend(TH1 *histo1,
		const Char_t *header,
		const Char_t *label1){
  TLegend *legend = new TLegend(0.537356,0.601695,0.899425,0.900424,
				NULL,"brNDC");
  legend->SetTextAlign(22);
  legend->SetTextFont(102);
  if (header != "")  legend->SetHeader(header); 
  legend->SetTextFont(102);

  TLegendEntry* entry1 = legend->AddEntry(histo1,label1,"LPF");  
  entry1->SetTextColor(histo1->GetLineColor());
 
  legend->SetFillColor(kWhite);
  legend->Draw();
  legend->SetBorderSize(1);
}

/***********************************************************************/
/*                                                                     */
/*       draw legend for 2 histos on the same plot                     */
/*                                                                     */
/***********************************************************************/
void Draw2Legend(TH1 *histo1, 
		 TH1 *histo2,
		 const Char_t *label1, 
		 const Char_t *label2,
		 const Char_t *header=""){

  Float_t max1 = histo1->GetMaximum();
  Float_t max2 = histo2->GetMaximum();
  if (max1 >= max2) histo1->SetMaximum(max1 * 1.1);
  else histo1->SetMaximum(max2 * 1.1);

  TLegend *legend = new TLegend(0.553482,0.653319,0.950072,0.948683,
				NULL,"brNDC");
  legend->SetTextAlign(22);
  legend->SetTextSize(0.1);
  legend->SetTextSize(0.06);
  legend->SetTextFont(102);
  
  TLegendEntry* entry1 = legend->AddEntry(histo1,label1,"LPF");  
  entry1->SetTextColor(histo1->GetLineColor());
  
  TLegendEntry* entry2 = legend->AddEntry(histo2,label2,"LPF");
  entry2->SetTextColor(histo2->GetLineColor());

  if (header != "") legend->SetHeader(header); 

  legend->SetFillColor(kWhite); 
  legend->Draw(); 
}
/***********************************************************************/
/*                                                                     */
/*       draw legend for 2 histos on the same plot                     */
/*                                                                     */
/***********************************************************************/
void Draw2LegendUL(TH1 *histo1, 
		   TH1 *histo2,
		   const Char_t *label1, 
		   const Char_t *label2,
		   const Char_t *header=""){

  Float_t max1 = histo1->GetMaximum();
  Float_t max2 = histo2->GetMaximum();
  if (max1 >= max2) histo1->SetMaximum(max1 * 1.1);
  else histo1->SetMaximum(max2 * 1.1);

    TLegend *legend = new TLegend(0.170738,0.619174,0.529933,0.919315,NULL,"brNDC");
  legend->SetTextAlign(22);
  legend->SetTextSize(0.1);
  legend->SetTextSize(0.06);
  legend->SetTextFont(102);
  
  TLegendEntry* entry1 = legend->AddEntry(histo1,label1,"LPF");  
  entry1->SetTextColor(histo1->GetMarkerColor());
  
  TLegendEntry* entry2 = legend->AddEntry(histo2,label2,"LPF");
  entry2->SetTextColor(histo2->GetMarkerColor());

  if (header != "") legend->SetHeader(header); 

  legend->SetFillColor(kWhite); 
  legend->Draw(); 
}


/***********************************************************************/
/*                                                                     */
/*       draw legend for 2 graphs on the same plot                     */
/*                                                                     */
/***********************************************************************/
void Draw2Legend(TGraph *histo1, 
		 TGraph *histo2,
		 const Char_t *label1, 
		 const Char_t *label2,
		 const Char_t *header="")
{
  TLegend *legend = new TLegend(0.574832,0.645657,0.946001,0.93697,NULL,"brNDC");
  legend->SetTextAlign(22);
  legend->SetTextFont(102);
  legend->SetTextSize(0.1);
  if (header != "") legend->SetHeader(header); 
  legend->SetTextSize(0.04);

  legend->SetTextFont(102);
  
  TLegendEntry* entry1 = legend->AddEntry(histo1,label1,"LP");  
  entry1->SetTextColor(histo1->GetLineColor());

  TLegendEntry* entry2 = legend->AddEntry(histo2,label2,"LP");
  entry2->SetTextColor(histo2->GetLineColor());

  legend->SetFillColor(kWhite); 
  legend->Draw(); 
}

/***********************************************************************/
/*                                                                     */
/*       draw a text box, for example "CALICE preliminary"             */
/*                                                                     */
/***********************************************************************/
void DrawText(TString string)
{
  TPaveText *pt = new TPaveText(0.1494042,0.9099783,0.4793767,0.9598698,"tbNDC");
  pt->SetTextSize(0.03);
  pt->SetFillColor(0);
  pt->SetBorderSize(1.);
  pt->SetTextAlign(22); //center
  pt->AddText(string);
  pt->Draw();

  pt->SetTextSize(0.05);
  pt->SetTextFont(72);
 
}


/***********************************************************************/
/*                                                                     */
/*       print parameters of a Gaussian fit                            */
/*                                                                     */
/***********************************************************************/
void PrintFitGausParameters(TF1* fitfctn)
{
  TPaveText *pt = new TPaveText(0.383985,0.612734,0.866147,0.842713,"tbNDC");
  pt->SetTextSize(0.04);
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  pt->SetBorderSize(1);

  TText *pte;
  /*char array to hold the legend's entries:*/
  char display[50];

  /*extract the fit parameters*/
  Double_t p0 = fitfctn->GetParameter(0);
  Double_t e0 = fitfctn->GetParError(0);
  sprintf(display, "const = %5.2f #\pm %3.2f", p0, e0);
  pte = pt->AddText(display);

  Double_t p1 = fitfctn->GetParameter(1);
  Double_t e1 = fitfctn->GetParError(1);
  sprintf(display, "#mu = %5.4f #\pm %3.2f", p1, e1);
  pte = pt->AddText(display);

  Double_t p2 = fitfctn->GetParameter(2);
  Double_t e2 = fitfctn->GetParError(2);
  sprintf(display, "#sigma = %5.4f #\pm %3.2f", p2, e2);
  pte = pt->AddText(display);

  pt->SetTextFont(102);
  pt->Draw();

}
