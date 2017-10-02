//void DrawScaleFactor(TString histname, TString sample){
TH1D* DrawScaleFactor(TString histname, TString sample){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();
  TString dir = "/uscms/homes/c/chuh/work/susy170922/";
  TFile* file6 = TFile::Open(dir+"TTJets_madgraphMLM_FastSim.root");
  TFile* file7 = TFile::Open(dir+"TTJets_madgraphMLM-pythia8.root");

  TString histname_W  = histname+"_Gen"+sample+"_"+sample;
  TString histname_Wo = histname+"_Gen"+sample+"_no_"+sample;

  TH1D* FullSim_W  = (TH1D*)file7->Get(histname_W);
  TH1D* FullSim_Wo = (TH1D*)file7->Get(histname_Wo);
  TH1D* FastSim_W  = (TH1D*)file6->Get(histname_W);
  TH1D* FastSim_Wo = (TH1D*)file6->Get(histname_Wo);

  TH1D* WtagFakeRate_FullSim  = (TH1D*)FullSim_W->Clone();
  TH1D* WtagFakeRate_FastSim  = (TH1D*)FastSim_W->Clone();

  WtagFakeRate_FullSim->Divide(FullSim_Wo);
  WtagFakeRate_FastSim->Divide(FastSim_Wo);

  TH1D* WtagFakeSF = (TH1D*)WtagFakeRate_FullSim->Clone();
  WtagFakeSF->Divide(WtagFakeRate_FastSim);

  TString name  = WtagFakeSF->GetName();
  TString title = WtagFakeSF->GetTitle();
  Int_t nx    = WtagFakeSF->GetNbinsX()+1;
  Double_t x1 = WtagFakeSF->GetBinLowEdge(1);
  Double_t bw = WtagFakeSF->GetBinWidth(nx);
  Double_t x2 = WtagFakeSF->GetBinLowEdge(nx)+bw;

  TH1D* hWtagFakeSF = new TH1D(name, title, nx, x1, x2);
  for (Int_t i=1; i<=nx; i++) {
      hWtagFakeSF->SetBinContent(hWtagFakeSF->FindBin(hWtagFakeSF->GetBinCenter(i)), WtagFakeSF->GetBinContent(i));
      hWtagFakeSF->SetBinError(  hWtagFakeSF->FindBin(hWtagFakeSF->GetBinCenter(i)), WtagFakeSF->GetBinError(i));
   }


  TCanvas* cWtagFakeRate_FullSim = new TCanvas("cWtagFakeRate_FullSim", "", 800, 800);
  WtagFakeRate_FullSim->SetMarkerColor(kBlack);
  WtagFakeRate_FullSim->SetLineColor(kBlack);
  WtagFakeRate_FullSim->SetMarkerStyle(20);
  WtagFakeRate_FullSim->SetMarkerSize(1.3);
  WtagFakeRate_FullSim->SetMaximum(0.5);
  WtagFakeRate_FullSim->SetMinimum(0);
  //WtagFakeRate_FullSim->GetYaxis()->SetTitle("W tagging fake rate");
  WtagFakeRate_FullSim->Draw("EP");

  float xmin = WtagFakeRate_FullSim->GetXaxis()->GetBinLowEdge(WtagFakeRate_FullSim->GetXaxis()->GetFirst());
  float xmax = WtagFakeRate_FullSim->GetXaxis()->GetBinUpEdge(WtagFakeRate_FullSim->GetXaxis()->GetLast());

  string text = "CMS #scale[0.7]{#font[52]{Work in progress 2016}}";
  TLatex* cms_lat = new TLatex(xmin, 0.51, text.c_str());
  cms_lat->SetTextSize(0.04);
  cms_lat->SetLineWidth(2);
  cms_lat->Draw();
  text = "#scale[0.7]{35.9 fb^{-1} (13 TeV)}";
  TLatex* era_lat = new TLatex(xmax,0.52, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.03);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();

  leg = new TLegend(0.51,0.82,0.9,0.9);
  leg->SetHeader("FullSim");
  leg->SetTextSize(0.03);
  leg->AddEntry(WtagFakeRate_FullSim,  "stat only", "ep");
  leg->Draw();
  
  TCanvas* cWtagFakeRate_FastSim = new TCanvas("cWtagFakeRate_FastSim", "", 800, 800);
  WtagFakeRate_FastSim->SetMarkerColor(kBlack);
  WtagFakeRate_FastSim->SetLineColor(kBlack);
  WtagFakeRate_FastSim->SetMarkerStyle(20);
  WtagFakeRate_FastSim->SetMarkerSize(1.3);
  WtagFakeRate_FastSim->SetMaximum(0.5);
  WtagFakeRate_FastSim->SetMinimum(0);
  //WtagFakeRate_FastSim->GetYaxis()->SetTitle("W tagging fake rate");
  WtagFakeRate_FastSim->Draw("EP");

  xmin = WtagFakeRate_FullSim->GetXaxis()->GetBinLowEdge(WtagFakeRate_FastSim->GetXaxis()->GetFirst());
  xmax = WtagFakeRate_FullSim->GetXaxis()->GetBinUpEdge(WtagFakeRate_FastSim->GetXaxis()->GetLast());

  text = "CMS #scale[0.7]{#font[52]{Work in progress 2016}}";
  cms_lat = new TLatex(xmin, 0.51, text.c_str());
  cms_lat->SetTextSize(0.04);
  cms_lat->SetLineWidth(2);
  cms_lat->Draw();
  text = "#scale[0.7]{35.9 fb^{-1} (13 TeV)}";
  era_lat = new TLatex(xmax,0.52, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.03);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();

  leg = new TLegend(0.51,0.82,0.9,0.9);
  leg->SetHeader("FastSim");
  leg->SetTextSize(0.03);
  leg->AddEntry(WtagFakeRate_FastSim,  "stat only", "ep");
  leg->Draw();
  
  TCanvas* cWtagFakeSF = new TCanvas("cWtagFakeSF", "", 800, 800);
  WtagFakeSF->SetMarkerColor(kBlack);
  WtagFakeSF->SetLineColor(kBlack);
  WtagFakeSF->SetMarkerStyle(20);
  WtagFakeSF->SetMarkerSize(1.3);
  WtagFakeSF->SetMaximum(1.5);
  WtagFakeSF->SetMinimum(0.5);
  //WtagFakeSF->GetYaxis()->SetTitle("W tagging fake rate");
  WtagFakeSF->Draw("EP");

  xmin = WtagFakeRate_FullSim->GetXaxis()->GetBinLowEdge(WtagFakeSF->GetXaxis()->GetFirst());
  xmax = WtagFakeRate_FullSim->GetXaxis()->GetBinUpEdge(WtagFakeSF->GetXaxis()->GetLast());

  text = "CMS #scale[0.7]{#font[52]{Work in progress 2016}}";
  cms_lat = new TLatex(xmin, 1.53, text.c_str());
  cms_lat->SetTextSize(0.04);
  cms_lat->SetLineWidth(2);
  cms_lat->Draw();
  text = "#scale[0.7]{35.9 fb^{-1} (13 TeV)}";
  era_lat = new TLatex(xmax,1.54, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.03);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();

  leg = new TLegend(0.51,0.82,0.9,0.9);
  //leg->SetHeader("W boson tag fake rate");
  leg->SetTextSize(0.03);
  leg->AddEntry(WtagFakeSF,  "stat only", "ep");
  leg->Draw();
  

  TCanvas* comQp = new TCanvas("comQp", "", 800, 800);
  comQp->Divide(1,2);
  TPad *canvsa_up = (TPad*)comQp->GetListOfPrimitives()->FindObject("comQp_1");
  TPad *canvsa_dw = (TPad*)comQp->GetListOfPrimitives()->FindObject("comQp_2");

  up_height     = 0.7; // please tune so that the upper figures size will meet your requirement
  dw_correction = 1.210;//1.390 // please tune so that the smaller canvsa size will work in your environment
  dw_height    = (1. - up_height) * dw_correction;

  canvsa_up->SetPad(0., 1 - up_height, 1., 1.);
  canvsa_dw->SetPad(0., 0., 1., dw_height);
  //canvsa_up->SetLogy();
  canvsa_up->SetFrameFillColor(0);
  canvsa_up->SetFillColor(0);
  canvsa_dw->SetFillColor(0);
  canvsa_dw->SetFrameFillColor(0);
  canvsa_up->SetLeftMargin(0.1);
  canvsa_dw->SetTopMargin(0.03);//0.03
  canvsa_dw->SetLeftMargin(0.1);
  canvsa_dw->SetBottomMargin(0.3);

  canvsa_up->cd();
  FullSim_W->SetLineColor(kBlack);
  FastSim_W->SetLineColor(kRed);
  FullSim_W->SetMaximum(1);
  FullSim_W->SetMinimum(0);
  FullSim_W->GetYaxis()->SetTitle("A.U.");
  FullSim_W->Scale(1./FullSim_W->Integral());
  FastSim_W->Scale(1./FastSim_W->Integral());
  FullSim_W->Draw("HIST");
  FastSim_W->Draw("HISTsame");

  xmin = FullSim_W->GetXaxis()->GetBinLowEdge(FullSim_W->GetXaxis()->GetFirst());
  xmax = FullSim_W->GetXaxis()->GetBinUpEdge(FullSim_W->GetXaxis()->GetLast());

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

  leg = new TLegend(0.51,0.82,0.9,0.9);

  leg->SetTextSize(0.03);
  leg->AddEntry(FullSim_W, "FullSim", "L");
  leg->AddEntry(FastSim_W, "FastSim", "l");
  leg->Draw();
  
  TH1D* compareQp = (TH1D*)FullSim_W->Clone(); 
  TH1D* compareQp1 = (TH1D*)FullSim_W->Clone(); 
  compareQp->Divide(FastSim_W);

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
  compareQp->GetYaxis()->SetTitle("#frac{FullSim}{FastSim}");
  compareQp->GetXaxis()->SetTitleOffset(0.81);
  compareQp->GetXaxis()->SetTitleSize(0.12);
  compareQp->GetYaxis()->SetTitleOffset(0.35);
  compareQp->GetYaxis()->SetTitleSize(0.10);//dividend->SetTitleSize(0.12,"Y");
  for(int i=1;i<=FullSim_W->GetNbinsX();i++){
    compareQp1->SetBinContent(i,1);
    if(FullSim_W->GetBinContent(i)==0) compareQp1->SetBinError(i,0);
    else compareQp1->SetBinError(i,FullSim_W->GetBinError(i)/FullSim_W->GetBinContent(i));
  }
  compareQp1->SetFillColor(13);
  compareQp1->SetFillStyle(3235);
  compareQp1->SetMarkerStyle(1);

  compareQp->Draw("EP");
  compareQp1->Draw("E2same");


  TCanvas* comSp = new TCanvas("comSp", "", 800, 800);
  comSp->Divide(1,2);
  TPad *canvas_up = (TPad*)comSp->GetListOfPrimitives()->FindObject("comSp_1");
  TPad *canvas_dw = (TPad*)comSp->GetListOfPrimitives()->FindObject("comSp_2");

  up_height     = 0.7; // please tune so that the upper figures size will meet your requirement
  dw_correction = 1.210;//1.390 // please tune so that the smaller canvas size will work in your environment
  dw_height    = (1. - up_height) * dw_correction;

  canvas_up->SetPad(0., 1 - up_height, 1., 1.);
  canvas_dw->SetPad(0., 0., 1., dw_height);
  //canvas_up->SetLogy();
  canvas_up->SetFrameFillColor(0);
  canvas_up->SetFillColor(0);
  canvas_dw->SetFillColor(0);
  canvas_dw->SetFrameFillColor(0);
  canvas_up->SetLeftMargin(0.1);
  canvas_dw->SetTopMargin(0.03);//0.03
  canvas_dw->SetLeftMargin(0.1);
  canvas_dw->SetBottomMargin(0.3);

  canvas_up->cd();
  FullSim_Wo->SetLineColor(kBlack);
  FastSim_Wo->SetLineColor(kRed);
  FullSim_Wo->SetMaximum(1);
  FullSim_Wo->SetMinimum(0);
  FullSim_Wo->GetYaxis()->SetTitle("A.U.");
  FullSim_Wo->Scale(1./FullSim_Wo->Integral());
  FastSim_Wo->Scale(1./FastSim_Wo->Integral());
  FullSim_Wo->Draw("HIST");
  FastSim_Wo->Draw("HISTsame");

  xmin = FullSim_Wo->GetXaxis()->GetBinLowEdge(FullSim_Wo->GetXaxis()->GetFirst());
  xmax = FullSim_Wo->GetXaxis()->GetBinUpEdge(FullSim_Wo->GetXaxis()->GetLast());

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

  leg = new TLegend(0.51,0.82,0.9,0.9);

  leg->SetTextSize(0.03);
  leg->AddEntry(FullSim_Wo, "FullSim", "L");
  leg->AddEntry(FastSim_Wo, "FastSim", "l");
  leg->Draw();
  
  TH1D* compareSp = (TH1D*)FullSim_Wo->Clone(); 
  TH1D* compareSp1 = (TH1D*)FullSim_Wo->Clone(); 
  compareSp->Divide(FastSim_Wo);

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
  compareSp->GetYaxis()->SetTitle("#frac{FullSim}{FastSim}");
  compareSp->GetXaxis()->SetTitleOffset(0.81);
  compareSp->GetXaxis()->SetTitleSize(0.12);
  compareSp->GetYaxis()->SetTitleOffset(0.35);
  compareSp->GetYaxis()->SetTitleSize(0.10);//dividend->SetTitleSize(0.12,"Y");
  for(int i=1;i<=FullSim_Wo->GetNbinsX();i++){
    compareSp1->SetBinContent(i,1);
    if(FullSim_Wo->GetBinContent(i)==0) compareSp1->SetBinError(i,0);
    else compareSp1->SetBinError(i,FullSim_Wo->GetBinError(i)/FullSim_Wo->GetBinContent(i));
  }
  compareSp1->SetFillColor(13);
  compareSp1->SetFillStyle(3235);
  compareSp1->SetMarkerStyle(1);

  compareSp->Draw("EP");
  compareSp1->Draw("E2same");


  TString Dir = "./plot/";
  TString type1 = "_FullSim.png";
  TString type2 = "_FastSim.png";
  TString type3 = "_SF.png";
  TString type4 = "_wSelection.png";
  TString type5 = "_woSelection.png";
  cWtagFakeRate_FullSim->SaveAs(Dir+histname_W+type1); 
  cWtagFakeRate_FastSim->SaveAs(Dir+histname_W+type2); 
  cWtagFakeSF->SaveAs(Dir+histname_W+type3); 
  comQp->SaveAs(Dir+histname_W+type4);
  comSp->SaveAs(Dir+histname_W+type5);

  return hWtagFakeSF;
}

void FullFastSimScaleFactor(){
  DrawScaleFactor("tau1", "W");
  DrawScaleFactor("tau2", "W");
  DrawScaleFactor("tau3", "W");
  DrawScaleFactor("tau21", "W");
  DrawScaleFactor("tau32", "W");
  DrawScaleFactor("tau1", "Top");
  DrawScaleFactor("tau2", "Top");
  DrawScaleFactor("tau3", "Top");
  DrawScaleFactor("tau21", "Top");
  DrawScaleFactor("tau32", "Top");
  DrawScaleFactor("SubjetBTag", "Top");
  DrawScaleFactor("AK8Jet1Pt", "W");
  DrawScaleFactor("AK8Jet1Pt", "Top");
/*
  TH1D* hFullFastSimWTagSF = DrawScaleFactor("AK8Jet1Pt_diffBin", "W");
  TH1D* hFullFastSimTopTagSF = DrawScaleFactor("AK8Jet1Pt_diffBin", "Top");
  hFullFastSimWTagSF->SetNameTitle("hFullFastSimWTagSF", "hFullFastSimWTagSF");
  hFullFastSimTopTagSF->SetNameTitle("hFullFastSimTopTagSF", "hFullFastSimTopTagSF");
  hFullFastSimWTagSF->GetYaxis()->UnZoom();
  hFullFastSimTopTagSF->GetYaxis()->UnZoom();

  TFile* output = new TFile("FullFastSimTagSF.root","recreate");
  hFullFastSimWTagSF->Write();
  hFullFastSimTopTagSF->Write();
*/
}
