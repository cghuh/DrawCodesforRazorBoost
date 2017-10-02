//void DrawScaleFactor(TString histname, TString sample){
pair<TH1D*, TH1D*> DrawScaleFactor(TString histname, TString sample){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();
  TString dir = "/uscms/homes/c/chuh/work/susy170922/";
  TFile* file6 = TFile::Open(dir+"TTJets_madgraphMLM_FastSim.root");
  TFile* file7 = TFile::Open(dir+"TTJets_madgraphMLM-pythia8.root");

  TString histname_W  = histname+"_Gen"+sample+"_"+sample;
  TString histname_Wo = histname+"_Gen"+sample+"_no_"+sample;

  TH2D* FullSim_W  = (TH2D*)file7->Get(histname_W);
  TH2D* FullSim_Wo = (TH2D*)file7->Get(histname_Wo);
  TH2D* FastSim_W  = (TH2D*)file6->Get(histname_W);
  TH2D* FastSim_Wo = (TH2D*)file6->Get(histname_Wo);
  TH1D* htool	   = (TH1D*)file7->Get("AK8Jet1Pt_GenW_W");

  TH1F* htemp1 = (TH1F*)htool->Clone();
  TH1F* htemp2 = (TH1F*)htool->Clone();
  TH1F* htemp3 = (TH1F*)htool->Clone();
  TH1F* htemp4 = (TH1F*)htool->Clone();
  TH1F* htemp5 = (TH1F*)htool->Clone();
  TH1F* htemp6 = (TH1F*)htool->Clone();
  TH1F* htemp7 = (TH1F*)htool->Clone();
  TH1F* htemp8 = (TH1F*)htool->Clone();

  for(Int_t i=1;i<=htool->GetNbinsX();i++){
    htemp1->SetBinContent(i,FullSim_W->GetBinContent(i,1));
    htemp1->SetBinError(i,FullSim_W->GetBinError(i,1));
    htemp2->SetBinContent(i,FullSim_W->GetBinContent(i,2));
    htemp2->SetBinError(i,FullSim_W->GetBinError(i,2));
    htemp3->SetBinContent(i,FullSim_Wo->GetBinContent(i,1));
    htemp3->SetBinError(i,FullSim_Wo->GetBinError(i,1));
    htemp4->SetBinContent(i,FullSim_Wo->GetBinContent(i,2));
    htemp4->SetBinError(i,FullSim_Wo->GetBinError(i,2));
    htemp5->SetBinContent(i,FastSim_W->GetBinContent(i,1));
    htemp5->SetBinError(i,FastSim_W->GetBinError(i,1));
    htemp6->SetBinContent(i,FastSim_W->GetBinContent(i,2));
    htemp6->SetBinError(i,FastSim_W->GetBinError(i,2));
    htemp7->SetBinContent(i,FastSim_Wo->GetBinContent(i,1));
    htemp7->SetBinError(i,FastSim_Wo->GetBinError(i,1));
    htemp8->SetBinContent(i,FastSim_Wo->GetBinContent(i,2));
    htemp8->SetBinError(i,FastSim_Wo->GetBinError(i,2));
  }


  TH1D* WtagFakeRate_EB_FullSim  = (TH1D*)htemp1->Clone();
  TH1D* WtagFakeRate_EE_FullSim  = (TH1D*)htemp2->Clone();
  TH1D* WtagFakeRate_EB_FastSim  = (TH1D*)htemp5->Clone();
  TH1D* WtagFakeRate_EE_FastSim  = (TH1D*)htemp6->Clone();

  WtagFakeRate_EB_FullSim->Divide(htemp3);
  WtagFakeRate_EE_FullSim->Divide(htemp4);
  WtagFakeRate_EB_FastSim->Divide(htemp7);
  WtagFakeRate_EE_FastSim->Divide(htemp8);

  TH1D* WtagFakeSF_EB = (TH1D*)WtagFakeRate_EB_FullSim->Clone();
  WtagFakeSF_EB->Divide(WtagFakeRate_EB_FastSim);
  WtagFakeSF_EB->SetName("b"+histname);
  TH1D* WtagFakeSF_EE = (TH1D*)WtagFakeRate_EE_FullSim->Clone();
  WtagFakeSF_EE->Divide(WtagFakeRate_EE_FastSim);
  WtagFakeSF_EE->SetName("b"+histname);

  TCanvas* c1 = new TCanvas("c1", "", 800, 800);
  WtagFakeSF_EB->SetMarkerColor(kBlack);
  WtagFakeSF_EB->SetLineColor(kBlack);
  WtagFakeRate_EB_FastSim->SetMarkerColor(kRed);
  WtagFakeRate_EB_FastSim->SetLineColor(kRed);
  WtagFakeRate_EB_FullSim->SetMarkerColor(kAzure);
  WtagFakeRate_EB_FullSim->SetLineColor(kAzure);
  WtagFakeSF_EE->SetMarkerColor(kCyan);
  WtagFakeSF_EE->SetLineColor(kCyan);
  WtagFakeRate_EE_FastSim->SetMarkerColor(kMagenta);
  WtagFakeRate_EE_FastSim->SetLineColor(kMagenta);
  WtagFakeRate_EE_FullSim->SetMarkerColor(kViolet);
  WtagFakeRate_EE_FullSim->SetLineColor(kViolet);
  WtagFakeSF_EB->SetMarkerStyle(20);
  WtagFakeSF_EB->SetMarkerSize(1.3);
  WtagFakeSF_EE->SetMarkerStyle(20);
  WtagFakeSF_EE->SetMarkerSize(1.3);
  WtagFakeSF_EB->SetMaximum(1.4);
  WtagFakeSF_EB->SetMinimum(0.0);
  WtagFakeSF_EB->Draw("EP");
  WtagFakeRate_EB_FullSim->Draw("EPsame");
  WtagFakeRate_EB_FastSim->Draw("EPsame");
  WtagFakeSF_EE->Draw("EPsame");
  WtagFakeRate_EE_FullSim->Draw("EPsame");
  WtagFakeRate_EE_FastSim->Draw("EPsame");

  float xmin = WtagFakeRate_EB_FullSim->GetXaxis()->GetBinLowEdge(WtagFakeSF_EB->GetXaxis()->GetFirst());
  float xmax = WtagFakeRate_EB_FullSim->GetXaxis()->GetBinUpEdge(WtagFakeSF_EB->GetXaxis()->GetLast());
  string text = "CMS #scale[0.7]{#font[52]{Work in progress 2016}}";
  TLatex* cms_lat = new TLatex(xmin, 1.43, text.c_str());
  cms_lat->SetTextSize(0.04);
  cms_lat->SetLineWidth(2);
  cms_lat->Draw();
  text = "#scale[0.7]{35.9 fb^{-1} (13 TeV)}";
  TLatex* era_lat = new TLatex(xmax,1.44, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.03);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();

  leg = new TLegend(0.48,0.20,0.9,0.42);
  leg->SetTextSize(0.03);
  if(TString(sample).Contains("W"))		leg->SetHeader("FullSim / FastSim W boson tag");
  else if(TString(sample).Contains("Top"))	leg->SetHeader("FullSim / FastSim Top quark tag");
  leg->SetTextSize(0.02);
  leg->AddEntry(WtagFakeSF_EB,            "Scale Factor, Barrel, stat only", "ep");
  leg->AddEntry(WtagFakeRate_EB_FullSim,  "#epsilon_{FullSim}, Barrel, stat only", "ep");
  leg->AddEntry(WtagFakeRate_EB_FastSim,  "#epsilon_{FastSim}, Barrel, stat only", "ep");
  leg->AddEntry(WtagFakeSF_EE,            "Scale Factor, Endcap, stat only", "ep");
  leg->AddEntry(WtagFakeRate_EE_FullSim,  "#epsilon_{FullSim}, Endcap, stat only", "ep");
  leg->AddEntry(WtagFakeRate_EE_FastSim,  "#epsilon_{FastSim}, Endcap, stat only", "ep");
  leg->Draw();

  TString Dir = "./plot/";
  TString type4 = ".png";
  c1->SaveAs(Dir+histname_W+type4);

  return make_pair(WtagFakeSF_EB, WtagFakeSF_EE);
}

void FullFastSimScaleFactor(){

  pair<TH1D*, TH1D*> p = DrawScaleFactor("AK8Jet1Pt_eta", "W");
  TH1D* hbWFF = p.first;
  hbWFF->SetName("bWFF");
  TH1D* heWFF = p.second;
  heWFF->SetName("eWFF");
  p = DrawScaleFactor("AK8Jet1Pt_eta", "Top");
  TH1D* hbTopFF = p.first;
  hbTopFF->SetName("bTopFF");
  TH1D* heTopFF = p.second;
  heTopFF->SetName("eTopFF");

  //hFullFastSimWTagSF->SetNameTitle("hFullFastSimWTagSF", "hFullFastSimWTagSF");
  //hFullFastSimTopTagSF->SetNameTitle("hFullFastSimTopTagSF", "hFullFastSimTopTagSF");
  //hFullFastSimWTagSF->GetYaxis()->UnZoom();
  //hFullFastSimTopTagSF->GetYaxis()->UnZoom();

  TFile* output = new TFile("FullFastSimTagSF.root","recreate");
  hbWFF->Write();
  heWFF->Write();
  hbTopFF->Write();
  heTopFF->Write();
}
