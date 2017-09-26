void ClosureTest(TString histname, TString njet){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();
  //TString dir = "/gatbawi/palgongsan/chuh/susy/susy170427/added/";
  TString dir = "/gatbawi/palgongsan/chuh/susy/susy170801/";
  TFile* file1 = TFile::Open(dir+"Multijet.root");
  TFile* file2 = TFile::Open(dir+"TTST.root");
  TFile* file3 = TFile::Open(dir+"WJet.root");
  TFile* file4 = TFile::Open(dir+"ZJet.root");
  TFile* file5 = TFile::Open(dir+"Multiboson.root");
  TFile* file6 = TFile::Open(dir+"DYToLL.root");
  TFile* file7 = TFile::Open(dir+"data.root");

  //TString histname = "MR";
  //TString histname = "R2";
  //TString njet		 = "";
  //TString njet		 = "_nj35";
  //TString njet		 = "_nj6";

  TString histname_S = histname+"_S"+njet;
  TString histname_Q = histname+"_Q"+njet;
  TString histname_T = histname+"_T"+njet;
  TString histname_W = histname+"_W"+njet;
  TString histname_s = histname+"_s"+njet;
  TString histname_q = histname+"_q"+njet;

  TH1F*   Multijet_S = (TH1F*)file1->Get(histname_S);
  TH1F*       TTST_S = (TH1F*)file2->Get(histname_S);
  TH1F*       WJet_S = (TH1F*)file3->Get(histname_S);
  TH1F*       ZJet_S = (TH1F*)file4->Get(histname_S);
  TH1F* Multiboson_S = (TH1F*)file5->Get(histname_S);
  TH1F*     DYToLL_S = (TH1F*)file6->Get(histname_S);
  TH1F*       data_S = (TH1F*)file7->Get(histname_S);
  TH1F*   Multijet_Q = (TH1F*)file1->Get(histname_Q);
  TH1F*       TTST_Q = (TH1F*)file2->Get(histname_Q);
  TH1F*       WJet_Q = (TH1F*)file3->Get(histname_Q);
  TH1F*       ZJet_Q = (TH1F*)file4->Get(histname_Q);
  TH1F* Multiboson_Q = (TH1F*)file5->Get(histname_Q);
  TH1F*     DYToLL_Q = (TH1F*)file6->Get(histname_Q);
  TH1F*       data_Q = (TH1F*)file7->Get(histname_Q);
  TH1F*   Multijet_T = (TH1F*)file1->Get(histname_T);
  TH1F*       TTST_T = (TH1F*)file2->Get(histname_T);
  TH1F*       WJet_T = (TH1F*)file3->Get(histname_T);
  TH1F*       ZJet_T = (TH1F*)file4->Get(histname_T);
  TH1F* Multiboson_T = (TH1F*)file5->Get(histname_T);
  TH1F*     DYToLL_T = (TH1F*)file6->Get(histname_T);
  TH1F*       data_T = (TH1F*)file7->Get(histname_T);
  TH1F*   Multijet_W = (TH1F*)file1->Get(histname_W);
  TH1F*       TTST_W = (TH1F*)file2->Get(histname_W);
  TH1F*       WJet_W = (TH1F*)file3->Get(histname_W);
  TH1F*       ZJet_W = (TH1F*)file4->Get(histname_W);
  TH1F* Multiboson_W = (TH1F*)file5->Get(histname_W);
  TH1F*     DYToLL_W = (TH1F*)file6->Get(histname_W);
  TH1F*       data_W = (TH1F*)file7->Get(histname_W);
  TH1F*   Multijet_s = (TH1F*)file1->Get(histname_s);
  TH1F*       TTST_s = (TH1F*)file2->Get(histname_s);
  TH1F*       WJet_s = (TH1F*)file3->Get(histname_s);
  TH1F*       ZJet_s = (TH1F*)file4->Get(histname_s);
  TH1F* Multiboson_s = (TH1F*)file5->Get(histname_s);
  TH1F*     DYToLL_s = (TH1F*)file6->Get(histname_s);
  TH1F*       data_s = (TH1F*)file7->Get(histname_s);
  TH1F*   Multijet_q = (TH1F*)file1->Get(histname_q);
  TH1F*       TTST_q = (TH1F*)file2->Get(histname_q);
  TH1F*       WJet_q = (TH1F*)file3->Get(histname_q);
  TH1F*       ZJet_q = (TH1F*)file4->Get(histname_q);
  TH1F* Multiboson_q = (TH1F*)file5->Get(histname_q);
  TH1F*     DYToLL_q = (TH1F*)file6->Get(histname_q);
  TH1F*       data_q = (TH1F*)file7->Get(histname_q);


  Multijet_S->SetName("Multijet_S");
  TTST_S->SetName("TTST_S");
  WJet_S->SetName("WJet_S");
  ZJet_S->SetName("ZJet_S");
  Multiboson_S->SetName("Multiboson_S");
  DYToLL_S->SetName("DYToLL_S");
  data_S->SetName("data_S");
  Multijet_Q->SetName("Multijet_Q");
  TTST_Q->SetName("TTST_Q");
  WJet_Q->SetName("WJet_Q");
  ZJet_Q->SetName("ZJet_Q");
  Multiboson_Q->SetName("Multiboson_Q");
  DYToLL_Q->SetName("DYToLL_Q");
  data_Q->SetName("data_Q");
  Multijet_T->SetName("Multijet_T");
  TTST_T->SetName("TTST_T");
  WJet_T->SetName("WJet_T");
  ZJet_T->SetName("ZJet_T");
  Multiboson_T->SetName("Multiboson_T");
  DYToLL_T->SetName("DYToLL_T");
  data_T->SetName("data_T");
  Multijet_W->SetName("Multijet_W");
  TTST_W->SetName("TTST_W");
  WJet_W->SetName("WJet_W");
  ZJet_W->SetName("ZJet_W");
  Multiboson_W->SetName("Multiboson_W");
  DYToLL_W->SetName("DYToLL_W");
  data_W->SetName("data_W");
  Multijet_s->SetName("Multijet_s");
  TTST_s->SetName("TTST_s");
  WJet_s->SetName("WJet_s");
  ZJet_s->SetName("ZJet_s");
  Multiboson_s->SetName("Multiboson_s");
  DYToLL_s->SetName("DYToLL_s");
  data_s->SetName("data_s");
  Multijet_q->SetName("Multijet_q");
  TTST_q->SetName("TTST_q");
  WJet_q->SetName("WJet_q");
  ZJet_q->SetName("ZJet_q");
  Multiboson_q->SetName("Multiboson_q");
  DYToLL_q->SetName("DYToLL_q");
  data_q->SetName("data_q");

  TH1F* NT_QCD  = (TH1F*)Multijet_S->Clone();
  TH1F* NSp_QCD = (TH1F*)Multijet_S->Clone();
  TH1F* NSp_T   = (TH1F*)Multijet_S->Clone();
  TH1F* NSp_W   = (TH1F*)Multijet_S->Clone();
  TH1F* NSp     = (TH1F*)Multijet_S->Clone();
  TH1F* NQp     = (TH1F*)Multijet_S->Clone();
  TH1F* NQp1    = (TH1F*)Multijet_S->Clone();
  TH1F* NQp2    = (TH1F*)Multijet_S->Clone();

  NT_QCD = (TH1F*)data_Q->Clone();
  NT_QCD->Add(ZJet_Q,-1);
  NT_QCD->Add(WJet_Q,-1);
  NT_QCD->Add(TTST_Q,-1);
  NT_QCD->Add(Multiboson_Q,-1);
  NT_QCD->Add(DYToLL_Q,-1);
  //NT_QCD->Multiply(Multijet_T);
  //NT_QCD->Divide(Multijet_Q);
  NT_QCD->Scale(Multijet_T->Integral()/Multijet_Q->Integral());

  NSp_QCD = (TH1F*)data_Q->Clone();
  NSp_QCD->Add(ZJet_Q,-1);
  NSp_QCD->Add(WJet_Q,-1);
  NSp_QCD->Add(TTST_Q,-1);
  NSp_QCD->Add(Multiboson_Q,-1);
  NSp_QCD->Add(DYToLL_Q,-1);
  //NSp_QCD->Multiply(Multijet_s);
  //NSp_QCD->Divide(Multijet_Q);
  NSp_QCD->Scale(Multijet_s->Integral()/Multijet_Q->Integral());

  NSp_W = (TH1F*)data_W->Clone();
  NSp_W->Add(ZJet_W,-1);
  NSp_W->Add(Multijet_W,-1);
  NSp_W->Add(TTST_W,-1);
  NSp_W->Add(Multiboson_W,-1);
  NSp_W->Add(DYToLL_W,-1);
  //NSp_W->Multiply(WJet_s);
  //NSp_W->Divide(WJet_W);
  NSp_W->Scale(WJet_s->Integral()/WJet_W->Integral());

  NSp_T = (TH1F*)data_T->Clone();
  NSp_T->Add(ZJet_T,-1);
  NSp_T->Add(WJet_T,-1);
  NSp_T->Add(Multijet_T,-1);
  NSp_T->Add(Multiboson_T,-1);
  NSp_T->Add(DYToLL_T,-1);
  NSp_T->Add(NT_QCD,-1);
  //NSp_T->Multiply(TTST_s);
  //NSp_T->Divide(TTST_T);
  NSp_T->Scale(TTST_s->Integral()/TTST_T->Integral());

  NSp = (TH1F*)NSp_QCD->Clone();
  NSp->Add(NSp_W);
  NSp->Add(NSp_T);
  NSp->Add(ZJet_s);
  NSp->Add(Multiboson_s);
  NSp->Add(DYToLL_s);

  NQp = (TH1F*)data_Q->Clone();
  NQp1 = (TH1F*)Multijet_q->Clone();
  NQp1->Add(TTST_q);
  NQp1->Add(WJet_q);
  NQp1->Add(ZJet_q);
  NQp1->Add(Multiboson_q);
  NQp1->Add(DYToLL_q);
  NQp2 = (TH1F*)Multijet_Q->Clone();
  NQp2->Add(TTST_Q);
  NQp2->Add(WJet_Q);
  NQp2->Add(ZJet_Q);
  NQp2->Add(Multiboson_Q);
  NQp2->Add(DYToLL_Q);
  //NQp->Multiply(NQp1);
  //NQp->Divide(NQp2);
  NQp->Scale(NQp1->Integral()/NQp2->Integral());
/*
  for(int i=1;i<=NT_QCD->GetNbinsX();i++){
    //cout << NT_QCD->GetBinContent(i) << endl;
    cout << i << " bin" << endl;
    cout << "Multijet_S' : " << Multijet_s->GetBinContent(i) << endl;
    cout << "Multijet_Q  : " << Multijet_Q->GetBinContent(i) << endl;
    cout << "Sp_Q'       : " << NSp_QCD->GetBinContent(i) << endl;
    cout << "Wjet_S'     : " << WJet_s->GetBinContent(i) << endl;
    cout << "Sp_W        : " << NSp_W->GetBinContent(i) << endl;
    cout << "TTST_S'     : " << TTST_s->GetBinContent(i) << endl;
    cout << "TTST_T      : " << TTST_T->GetBinContent(i) << endl;
    cout << "Sp_T        : " << NSp_T->GetBinContent(i) << endl;
    cout << "Total       : " << NSp_T->GetBinContent(i)+NSp_QCD->GetBinContent(i)+NSp_W->GetBinContent(i)+ZJet_s->GetBinContent(i)+Multiboson_s->GetBinContent(i)+DYToLL_s->GetBinContent(i) << endl;
    //cout << NSp_T->GetBinContent(i)+NSp_QCD->GetBinContent(i)+NSp_W->GetBinContent(i) << endl;
    //cout << NSp->GetBinContent(i) << endl;
  }
*/

  THStack *hs = new THStack("hs","");
  //NSp_W->SetLineColorAlpha(kGreen,30);
  NSp_W->SetLineColor(418);
  NSp_W->SetLineWidth(4);
  hs->Add(NSp_W);
  //NSp_T->SetLineColorAlpha(28,30);
  NSp_T->SetLineColor(633);
  NSp_T->SetLineWidth(4);
  hs->Add(NSp_T);
  //NSp_QCD->SetLineColorAlpha(kMagenta,30);
  NSp_QCD->SetLineColor(619);
  NSp_QCD->SetLineWidth(4);
  hs->Add(NSp_QCD);

  TCanvas* comSp = new TCanvas("comSp", "", 800, 800);
  comSp->Divide(1,2);
  TPad *canvas_up = (TPad*)comSp->GetListOfPrimitives()->FindObject("comSp_1");
  TPad *canvas_dw = (TPad*)comSp->GetListOfPrimitives()->FindObject("comSp_2");

  double up_height     = 0.7; // please tune so that the upper figures size will meet your requirement
  double dw_correction = 1.210;//1.390 // please tune so that the smaller canvas size will work in your environment
  double dw_height    = (1. - up_height) * dw_correction;

  canvas_up->SetPad(0., 1 - up_height, 1., 1.);
  canvas_dw->SetPad(0., 0., 1., dw_height);
  canvas_up->SetLogy();
  canvas_up->SetFrameFillColor(0);
  canvas_up->SetFillColor(0);
  canvas_dw->SetFillColor(0);
  canvas_dw->SetFrameFillColor(0);
  canvas_up->SetLeftMargin(0.1);
  canvas_dw->SetTopMargin(0.03);//0.03
  canvas_dw->SetLeftMargin(0.1);
  canvas_dw->SetBottomMargin(0.3);

  canvas_up->cd();
  data_s->SetMarkerColor(kBlack);
  data_s->SetLineColor(kBlack);
  data_s->SetMarkerStyle(20);
  data_s->SetMarkerSize(1.3);
  NSp->SetFillColor(kTeal);
  NSp->SetLineColor(kTeal);
  NSp->SetMarkerColor(kTeal);
  NSp->SetFillStyle(3003);
  NSp->SetMaximum(1e4);
  NSp->SetMinimum(1e-1);
  NSp->GetYaxis()->SetTitle("Events/bin");
  NSp->Draw("E2");
  hs->Draw("HISTsame");
  data_s->Draw("EPsame");

  float xmin = NSp->GetXaxis()->GetBinLowEdge(data_s->GetXaxis()->GetFirst());
  float xmax = NSp->GetXaxis()->GetBinUpEdge(data_s->GetXaxis()->GetLast());

  string text = "CMS #scale[0.7]{#font[52]{Work in progress 2016}}";
  TLatex* cms_lat = new TLatex(xmin, 1.8e4, text.c_str());
  cms_lat->SetTextSize(0.07);
  cms_lat->SetLineWidth(2);
  cms_lat->Draw();
  text = "#scale[0.7]{W ana, 35.9 fb^{-1} (13 TeV)}";
  TLatex* era_lat = new TLatex(xmax,2e4, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.06);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();

  leg = new TLegend(0.51,0.62,0.9,0.9);
  if(TString(njet).Contains("nj35")) leg->SetHeader("S' Closure test, 3 #leq jet < 5");
  else if(TString(njet).Contains("nj35")) leg->SetHeader("S' Closure test, 6 #leq jet");
  else leg->SetHeader("S' Closure test");

  char legentry[10];

  leg->SetNColumns(2);
  leg->SetTextSize(0.03);
  sprintf(legentry, "%.2f",NSp_QCD->Integral());
  leg->AddEntry(NSp_QCD,  "Multijet estimate", "l");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.2f",NSp_T->Integral());
  leg->AddEntry(NSp_T,  "t#bar{t} + single top estimate", "l");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.2f",NSp_W->Integral());
  leg->AddEntry(NSp_W,  "W(#rightarrowl#nu) estimate", "l");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.2f",data_s->Integral());
  leg->AddEntry(data_s, "data", "lp");
  leg->AddEntry((TObject*)0,legentry,"");
  leg->Draw();
  
  TH1F* compareSp = (TH1F*)NSp->Clone(); 
  TH1F* compareSp1 = (TH1F*)NSp->Clone(); 
  compareSp->Divide(data_s);

  canvas_dw->cd();
  gPad->SetGridy();
  compareSp->SetLineColor(kBlack);
  compareSp->SetMarkerColor(kBlack);
  compareSp->SetMarkerStyle(20);
  compareSp->SetMarkerSize(1.3);
  compareSp->GetYaxis()->SetRangeUser(0.0,2.0);
  compareSp->SetTitle(0);
  compareSp->SetStats(0);
  compareSp->GetYaxis()->SetLabelSize(0.09);
  compareSp->GetXaxis()->SetLabelSize(0.09);
  compareSp->GetYaxis()->SetNdivisions(502);
  compareSp->GetXaxis()->SetTickLength(0.07);
  compareSp->GetYaxis()->SetTickLength(0.03);
  compareSp->GetYaxis()->SetTitle("#frac{Data}{Estimate}");
  compareSp->GetXaxis()->SetTitleOffset(0.81);
  compareSp->GetXaxis()->SetTitleSize(0.12);
  compareSp->GetYaxis()->SetTitleOffset(0.35);
  compareSp->GetYaxis()->SetTitleSize(0.10);//dividend->SetTitleSize(0.12,"Y");
  for(int i=1;i<=NSp->GetNbinsX();i++){
    compareSp1->SetBinContent(i,1);
    if(NSp->GetBinContent(i)==0) compareSp1->SetBinError(i,0);
    else compareSp1->SetBinError(i,NSp->GetBinError(i)/NSp->GetBinContent(i));
  }
  compareSp1->SetFillColor(13);
  compareSp1->SetFillStyle(3235);
  compareSp1->SetMarkerStyle(1);

  compareSp->Draw("EP");
  compareSp1->Draw("E2same");

  cout << "S' region drawn" << endl;

  TCanvas* comQp = new TCanvas("comQp", "", 800, 800);
  comQp->Divide(1,2);
  TPad *canvsa_up = (TPad*)comQp->GetListOfPrimitives()->FindObject("comQp_1");
  TPad *canvsa_dw = (TPad*)comQp->GetListOfPrimitives()->FindObject("comQp_2");

  up_height     = 0.7; // please tune so that the upper figures size will meet your requirement
  dw_correction = 1.210;//1.390 // please tune so that the smaller canvsa size will work in your environment
  dw_height    = (1. - up_height) * dw_correction;

  canvsa_up->SetPad(0., 1 - up_height, 1., 1.);
  canvsa_dw->SetPad(0., 0., 1., dw_height);
  canvsa_up->SetLogy();
  canvsa_up->SetFrameFillColor(0);
  canvsa_up->SetFillColor(0);
  canvsa_dw->SetFillColor(0);
  canvsa_dw->SetFrameFillColor(0);
  canvsa_up->SetLeftMargin(0.1);
  canvsa_dw->SetTopMargin(0.03);//0.03
  canvsa_dw->SetLeftMargin(0.1);
  canvsa_dw->SetBottomMargin(0.3);

  canvsa_up->cd();
  data_q->SetMarkerColor(kBlack);
  data_q->SetLineColor(kBlack);
  data_q->SetMarkerStyle(20);
  data_q->SetMarkerSize(1.3);
  NQp->SetFillColor(kTeal);
  NQp->SetLineColor(kTeal);
  NQp->SetMarkerColor(kTeal);
  NQp->SetFillStyle(3003);
  NQp->SetMaximum(1e4);
  NQp->SetMinimum(1e-1);
  NQp->GetYaxis()->SetTitle("Events/bin");
  NQp->Draw("E2");
  data_q->Draw("EPsame");

  xmin = NQp->GetXaxis()->GetBinLowEdge(data_s->GetXaxis()->GetFirst());
  xmax = NQp->GetXaxis()->GetBinUpEdge(data_s->GetXaxis()->GetLast());

  text = "CMS #scale[0.7]{#font[52]{Work in progress 2016}}";
  cms_lat = new TLatex(xmin, 1.8e4, text.c_str());
  cms_lat->SetTextSize(0.07);
  cms_lat->SetLineWidth(2);
  cms_lat->Draw();
  text = "#scale[0.7]{W ana, 35.9 fb^{-1} (13 TeV)}";
  era_lat = new TLatex(xmax,2e4, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.06);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();

  leg = new TLegend(0.51,0.70,0.9,0.9);
  if(TString(njet).Contains("nj35")) leg->SetHeader("Q' Closure test, 3 #leq jet < 5");
  else if(TString(njet).Contains("nj35")) leg->SetHeader("Q' Closure test, 6 #leq jet");
  else leg->SetHeader("Q' Closure test");

  leg->SetNColumns(2);
  leg->SetTextSize(0.03);
  sprintf(legentry, "%.2f",NQp->Integral());
  leg->AddEntry(NQp,  "Total background estimate", "LPE");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.2f",data_q->Integral());
  leg->AddEntry(data_q, "data", "lp");
  leg->AddEntry((TObject*)0,legentry,"");
  leg->Draw();
  
  TH1F* compareQp = (TH1F*)NQp->Clone(); 
  TH1F* compareQp1 = (TH1F*)NQp->Clone(); 
  compareQp->Divide(data_q);

  canvsa_dw->cd();
  gPad->SetGridy();
  compareQp->SetLineColor(kBlack);
  compareQp->SetMarkerColor(kBlack);
  compareQp->SetMarkerStyle(20);
  compareQp->SetMarkerSize(1.3);
  compareQp->GetYaxis()->SetRangeUser(0.0,2.0);
  compareQp->SetTitle(0);
  compareQp->SetStats(0);
  compareQp->GetYaxis()->SetLabelSize(0.09);
  compareQp->GetXaxis()->SetLabelSize(0.09);
  compareQp->GetYaxis()->SetNdivisions(502);
  compareQp->GetXaxis()->SetTickLength(0.07);
  compareQp->GetYaxis()->SetTickLength(0.03);
  compareQp->GetYaxis()->SetTitle("#frac{Data}{Estimate}");
  compareQp->GetXaxis()->SetTitleOffset(0.81);
  compareQp->GetXaxis()->SetTitleSize(0.12);
  compareQp->GetYaxis()->SetTitleOffset(0.35);
  compareQp->GetYaxis()->SetTitleSize(0.10);//dividend->SetTitleSize(0.12,"Y");
  for(int i=1;i<=NQp->GetNbinsX();i++){
    compareQp1->SetBinContent(i,1);
    if(NQp->GetBinContent(i)==0) compareQp1->SetBinError(i,0);
    else compareQp1->SetBinError(i,NQp->GetBinError(i)/NQp->GetBinContent(i));
  }
  compareQp1->SetFillColor(13);
  compareQp1->SetFillStyle(3235);
  compareQp1->SetMarkerStyle(1);

  compareQp->Draw("EP");
  compareQp1->Draw("E2same");

  cout << "Q' region drawn" << endl;

  TString Dir = "./plot2/";
  TString type1 = "_s.png";
  TString type2 = "_q.png";
  comSp->SaveAs(Dir+histname+njet+type1); 
  comQp->SaveAs(Dir+histname+njet+type2); 

}

void BGest(){
  ClosureTest("MR", "");
  ClosureTest("MR", "_nj35");
  ClosureTest("MR", "_nj6");
  ClosureTest("R2", "");
  ClosureTest("R2", "_nj35");
  ClosureTest("R2", "_nj6");
  ClosureTest("MR_topana", "");
  ClosureTest("R2_topana", "");
}
