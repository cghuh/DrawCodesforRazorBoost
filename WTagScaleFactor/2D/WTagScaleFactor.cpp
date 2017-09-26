//void DrawScaleFactor(TString histname){
TH1F* DrawScaleFactor(TString histname){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();
  TString dir = "/gatbawi/palgongsan/chuh/susy/susy170922/";
  TFile* file6 = TFile::Open(dir+"background.root");
  TFile* file7 = TFile::Open(dir+"data.root");

  TString histname_W = "AK8Jet1Pt_Eta_"+histname+"_fakerate";
  TString histname_Wo = "AK8Jet1Pt_Eta_no_"+histname+"_fakerate";
  TString tool = "AK8Jet1Pt_"+histname+"_fakerate";

  TH2D*       data_W  = (TH2D*)file7->Get(histname_W);
  TH2D*       data_Wo = (TH2D*)file7->Get(histname_Wo);
  TH2D* background_W  = (TH2D*)file6->Get(histname_W);
  TH2D* background_Wo = (TH2D*)file6->Get(histname_Wo);
  TH1D*       	htool	= (TH1D*)file7->Get(tool);


  TH1F* htemp1 = (TH1F*)htool->Clone();
  TH1F* htemp2 = (TH1F*)htool->Clone();
  TH1F* htemp3 = (TH1F*)htool->Clone();
  TH1F* htemp4 = (TH1F*)htool->Clone();
  TH1F* htemp5 = (TH1F*)htool->Clone();
  TH1F* htemp6 = (TH1F*)htool->Clone();
  TH1F* htemp7 = (TH1F*)htool->Clone();
  TH1F* htemp8 = (TH1F*)htool->Clone();

  for(Int_t i=1;i<=htool->GetNbinsX();i++){
    htemp1->SetBinContent(i,data_W->GetBinContent(i,1));
    htemp1->SetBinError(i,data_W->GetBinError(i,1));
    htemp2->SetBinContent(i,data_W->GetBinContent(i,2));
    htemp2->SetBinError(i,data_W->GetBinError(i,2));
    htemp3->SetBinContent(i,data_Wo->GetBinContent(i,1));
    htemp3->SetBinError(i,data_Wo->GetBinError(i,1));
    htemp4->SetBinContent(i,data_Wo->GetBinContent(i,2));
    htemp4->SetBinError(i,data_Wo->GetBinError(i,2));
    htemp5->SetBinContent(i,background_W->GetBinContent(i,1));
    htemp5->SetBinError(i,background_W->GetBinError(i,1));
    htemp6->SetBinContent(i,background_W->GetBinContent(i,2));
    htemp6->SetBinError(i,background_W->GetBinError(i,2));
    htemp7->SetBinContent(i,background_Wo->GetBinContent(i,1));
    htemp7->SetBinError(i,background_Wo->GetBinError(i,1));
    htemp8->SetBinContent(i,background_Wo->GetBinContent(i,2));
    htemp8->SetBinError(i,background_Wo->GetBinError(i,2));
  }

  TH1F* WtagFakeRate_Barrel_data  = (TH1F*)htemp1->Clone();
  TH1F* WtagFakeRate_Endcap_data  = (TH1F*)htemp2->Clone();
  TH1F* WtagFakeRate_Barrel_background  = (TH1F*)htemp5->Clone();
  TH1F* WtagFakeRate_Endcap_background  = (TH1F*)htemp6->Clone();

  WtagFakeRate_Barrel_data->Divide(htemp3);
  WtagFakeRate_Endcap_data->Divide(htemp4);
  WtagFakeRate_Barrel_background->Divide(htemp7);
  WtagFakeRate_Endcap_background->Divide(htemp8);

  TH1F* WtagFakeSF_Barrel = (TH1F*)WtagFakeRate_Barrel_data->Clone();
  WtagFakeSF_Barrel->Divide(WtagFakeRate_Barrel_background);
  TH1F* WtagFakeSF_Endcap = (TH1F*)WtagFakeRate_Endcap_data->Clone();
  WtagFakeSF_Endcap->Divide(WtagFakeRate_Endcap_background);

/*
  TString    name  = WtagFakeRate_Barrel_data->GetName();
  TString    title = WtagFakeRate_Barrel_data->GetTitle();
  Int_t    nx    = WtagFakeRate_Barrel_data->GetNbinsX()+1;
  Double_t x1 = WtagFakeRate_Barrel_data->GetBinLowEdge(1);
  Double_t bw = WtagFakeRate_Barrel_data->GetBinWidth(nx);
  Double_t x2 = WtagFakeRate_Barrel_data->GetBinLowEdge(nx)+bw;

  TH1F* hWtagFakeRate_Barrel_data = new TH1F(name, title, nx, x1, x2);
  for (Int_t i=1; i<=nx; i++) {
      hWtagFakeRate_Barrel_data->Fill(hWtagFakeRate_Barrel_data->GetBinCenter(i), WtagFakeRate_Barrel_data->GetBinContent(i));
  }

  name  = WtagFakeRate_Barrel_background->GetName();
  title = WtagFakeRate_Barrel_background->GetTitle();
  nx    = WtagFakeRate_Barrel_background->GetNbinsX()+1;
  x1 = WtagFakeRate_Barrel_background->GetBinLowEdge(1);
  bw = WtagFakeRate_Barrel_background->GetBinWidth(nx);
  x2 = WtagFakeRate_Barrel_background->GetBinLowEdge(nx)+bw;

  TH1F* hWtagFakeRate_Barrel_background = new TH1F(name, title, nx, x1, x2);
  for (Int_t i=1; i<=nx; i++) {
      hWtagFakeRate_Barrel_background->Fill(hWtagFakeRate_Barrel_background->GetBinCenter(i), WtagFakeRate_Barrel_background->GetBinContent(i));
  }

  name  = WtagFakeSF_Barrel->GetName();
  title = WtagFakeSF_Barrel->GetTitle();
  nx    = WtagFakeSF_Barrel->GetNbinsX()+1;
  //x1 = WtagFakeSF_Barrel->GetBinLowEdge(1);
  //bw = WtagFakeSF_Barrel->GetBinWidth(nx);
  //x2 = WtagFakeSF_Barrel->GetBinLowEdge(nx)+bw;
  Double_t *xbins= new Double_t[nx+1];
  for (UInt_t i=0;i<nx;i++) xbins[i]=*(WtagFakeSF_Barrel->GetXaxis()->GetXbins()->GetArray()+i);
  xbins[nx]=xbins[nx-1]+WtagFakeSF_Barrel->GetBinWidth(nx);

  TH1F* hWtagFakeSF_Barrel = new TH1F(name, title, nx, xbins);
  for (Int_t i=1; i<=nx; i++) {
      hWtagFakeSF_Barrel->SetBinContent(hWtagFakeSF_Barrel->FindBin(hWtagFakeSF_Barrel->GetBinCenter(i)), WtagFakeSF_Barrel->GetBinContent(i));
      hWtagFakeSF_Barrel->SetBinError(  hWtagFakeSF_Barrel->FindBin(hWtagFakeSF_Barrel->GetBinCenter(i)), WtagFakeSF_Barrel->GetBinError(i));
   }
*/

  TCanvas* cWtagFakeRate_Barrel_data = new TCanvas("cWtagFakeRate_Barrel_data", "", 800, 800);
  WtagFakeRate_Barrel_data->SetMarkerColor(kBlack);
  WtagFakeRate_Barrel_data->SetLineColor(kBlack);
  WtagFakeRate_Barrel_data->SetMarkerStyle(20);
  WtagFakeRate_Barrel_data->SetMarkerSize(1.3);
  WtagFakeRate_Barrel_data->SetMaximum(0.2);
  if(TString(histname_W).Contains("_Top")) WtagFakeRate_Barrel_data->SetMaximum(0.01);
  WtagFakeRate_Barrel_data->SetMinimum(0);
  //WtagFakeRate_Barrel_data->GetYaxis()->SetTitle("W tagging efficiency");
  WtagFakeRate_Barrel_data->Draw("EP");

  float xmin = WtagFakeRate_Barrel_data->GetXaxis()->GetBinLowEdge(WtagFakeRate_Barrel_data->GetXaxis()->GetFirst());
  float xmax = WtagFakeRate_Barrel_data->GetXaxis()->GetBinUpEdge(WtagFakeRate_Barrel_data->GetXaxis()->GetLast());

  string text = "CMS #scale[0.7]{#font[52]{Work in progress 2016}}";
  TLatex* cms_lat = new TLatex(xmin, 0.21, text.c_str());
  cms_lat->SetTextSize(0.04);
  cms_lat->SetLineWidth(2);
  cms_lat->Draw();
  text = "#scale[0.7]{35.9 fb^{-1} (13 TeV)}";
  TLatex* era_lat = new TLatex(xmax,0.22, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.03);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();

  leg = new TLegend(0.51,0.77,0.9,0.9);
  if(TString(histname).Contains("mW")) leg->SetHeader("W boson mass-tag efficiency");
  else if(TString(histname).Contains("aW")) leg->SetHeader("W boson anti-tag efficiency");
  else if(TString(histname).Contains("W"))leg->SetHeader("W boson tag efficiency");
  else if(TString(histname).Contains("mTop")) 			leg->SetHeader("Top quark mass-tag efficiency");
  else if(TString(histname).Contains("aTop")) leg->SetHeader("Top quark anti-tag efficiency");
  else if(TString(histname).Contains("Top"))  leg->SetHeader("Top quark tag efficiency");
  leg->SetTextSize(0.03);
  leg->AddEntry(WtagFakeRate_Barrel_data,  "stat only", "ep");
  leg->Draw();
  
  TCanvas* cWtagFakeRate_Barrel_background = new TCanvas("cWtagFakeRate_Barrel_background", "", 800, 800);
  WtagFakeRate_Barrel_background->SetMarkerColor(kBlack);
  WtagFakeRate_Barrel_background->SetLineColor(kBlack);
  WtagFakeRate_Barrel_background->SetMarkerStyle(20);
  WtagFakeRate_Barrel_background->SetMarkerSize(1.3);
  WtagFakeRate_Barrel_background->SetMaximum(0.2);
  if(TString(histname_W).Contains("_Top")) WtagFakeRate_Barrel_background->SetMaximum(0.01);
  WtagFakeRate_Barrel_background->SetMinimum(0);
  //WtagFakeRate_Barrel_background->GetYaxis()->SetTitle("W tagging efficiency");
  WtagFakeRate_Barrel_background->Draw("EP");

  xmin = WtagFakeRate_Barrel_data->GetXaxis()->GetBinLowEdge(WtagFakeRate_Barrel_background->GetXaxis()->GetFirst());
  xmax = WtagFakeRate_Barrel_data->GetXaxis()->GetBinUpEdge(WtagFakeRate_Barrel_background->GetXaxis()->GetLast());

  text = "CMS #scale[0.7]{#font[52]{Work in progress 2016}}";
  cms_lat = new TLatex(xmin, 0.21, text.c_str());
  cms_lat->SetTextSize(0.04);
  cms_lat->SetLineWidth(2);
  cms_lat->Draw();
  text = "#scale[0.7]{35.9 fb^{-1} (13 TeV)}";
  era_lat = new TLatex(xmax,0.22, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.03);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();

  leg = new TLegend(0.51,0.77,0.9,0.9);
  if(TString(histname).Contains("mW")) 			leg->SetHeader("W boson mass-tag efficiency");
  else if(TString(histname).Contains("aW")) leg->SetHeader("W boson anti-tag efficiency");
  else if(TString(histname).Contains("W"))  leg->SetHeader("W boson tag efficiency");
  else if(TString(histname).Contains("mTop")) 			leg->SetHeader("Top quark mass-tag efficiency");
  else if(TString(histname).Contains("aTop")) leg->SetHeader("Top quark anti-tag efficiency");
  else if(TString(histname).Contains("Top"))  leg->SetHeader("Top quark tag efficiency");
  leg->SetTextSize(0.03);
  leg->AddEntry(WtagFakeRate_Barrel_background,  "stat only", "ep");
  leg->Draw();
  
  TCanvas* cWtagFakeSF_Barrel = new TCanvas("cWtagFakeSF_Barrel", "", 800, 800);
  WtagFakeSF_Barrel->SetMarkerColor(kBlack);
  WtagFakeSF_Barrel->SetLineColor(kBlack);
  WtagFakeSF_Barrel->SetMarkerStyle(20);
  WtagFakeSF_Barrel->SetMarkerSize(1.3);
  WtagFakeSF_Barrel->SetMaximum(1.4);
  WtagFakeSF_Barrel->SetMinimum(0.5);
  WtagFakeSF_Barrel->Draw("EP");

  xmin = WtagFakeRate_Barrel_data->GetXaxis()->GetBinLowEdge(WtagFakeSF_Barrel->GetXaxis()->GetFirst());
  xmax = WtagFakeRate_Barrel_data->GetXaxis()->GetBinUpEdge(WtagFakeSF_Barrel->GetXaxis()->GetLast());

  text = "CMS #scale[0.7]{#font[52]{Work in progress 2016}}";
  cms_lat = new TLatex(xmin, 1.43, text.c_str());
  cms_lat->SetTextSize(0.04);
  cms_lat->SetLineWidth(2);
  cms_lat->Draw();
  text = "#scale[0.7]{35.9 fb^{-1} (13 TeV)}";
  era_lat = new TLatex(xmax,1.44, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.03);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();

  leg = new TLegend(0.51,0.77,0.9,0.9);
  if(TString(histname).Contains("mW")) 			leg->SetHeader("W boson mass-tag scale factor");
  else if(TString(histname).Contains("aW")) leg->SetHeader("W boson anti-tag scale factor");
  else if(TString(histname).Contains("W"))  leg->SetHeader("W boson tag scale factor");
  else if(TString(histname).Contains("mTop")) 			leg->SetHeader("Top quark mass-tag scale factor");
  else if(TString(histname).Contains("aTop")) leg->SetHeader("Top quark anti-tag scale factor");
  else if(TString(histname).Contains("Top"))  leg->SetHeader("Top quark tag scale factor");
  leg->SetTextSize(0.03);
  leg->AddEntry(WtagFakeSF_Barrel,  						"Scale Factor, stat only", "ep");
  leg->Draw();

  TString Dir = "./plot/";
  TString type1 = "_data.png";
  TString type2 = "_background.png";
  TString type3 = "_SF.png";
  TString type4 = "_barrel.png";
  cWtagFakeRate_Barrel_data->SaveAs(Dir+histname_W+type1); 
  cWtagFakeRate_Barrel_background->SaveAs(Dir+histname_W+type2); 
  cWtagFakeSF_Barrel->SaveAs(Dir+histname_W+type3); 


  TCanvas* c1 = new TCanvas("c1", "", 800, 800);
  WtagFakeSF_Barrel->SetMarkerColor(kBlack);
  WtagFakeSF_Barrel->SetLineColor(kBlack);
  WtagFakeRate_Barrel_background->SetMarkerColor(kRed);
  WtagFakeRate_Barrel_background->SetLineColor(kRed);
  WtagFakeRate_Barrel_data->SetMarkerColor(kAzure);
  WtagFakeRate_Barrel_data->SetLineColor(kAzure);
  WtagFakeRate_Barrel_data->SetMarkerStyle(20);
  WtagFakeRate_Barrel_data->SetMarkerSize(1.3);
  WtagFakeRate_Barrel_background->SetMarkerStyle(20);
  WtagFakeRate_Barrel_background->SetMarkerSize(1.3);
  WtagFakeSF_Barrel->SetMarkerStyle(20);
  WtagFakeSF_Barrel->SetMarkerSize(1.3);
  WtagFakeSF_Barrel->SetMaximum(1.4);
  WtagFakeSF_Barrel->SetMinimum(0.0);
  WtagFakeSF_Barrel->Draw("EP");
  WtagFakeRate_Barrel_data->Draw("EPsame");
  WtagFakeRate_Barrel_background->Draw("EPsame");

  xmin = WtagFakeRate_Barrel_data->GetXaxis()->GetBinLowEdge(WtagFakeSF_Barrel->GetXaxis()->GetFirst());
  xmax = WtagFakeRate_Barrel_data->GetXaxis()->GetBinUpEdge(WtagFakeSF_Barrel->GetXaxis()->GetLast());

  text = "CMS #scale[0.7]{#font[52]{Work in progress 2016}}";
  cms_lat = new TLatex(xmin, 1.43, text.c_str());
  cms_lat->SetTextSize(0.04);
  cms_lat->SetLineWidth(2);
  cms_lat->Draw();
  text = "#scale[0.7]{35.9 fb^{-1} (13 TeV)}";
  era_lat = new TLatex(xmax,1.44, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.03);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();

  leg = new TLegend(0.51,0.77,0.9,0.9);
  if(TString(histname).Contains("mW")) 			leg->SetHeader("W boson mass-tag");
  else if(TString(histname).Contains("aW")) leg->SetHeader("W boson anti-tag");
  else if(TString(histname).Contains("W"))  leg->SetHeader("W boson tag");
  else if(TString(histname).Contains("mTop")) 			leg->SetHeader("Top quark mass-tag");
  else if(TString(histname).Contains("aTop")) leg->SetHeader("Top quark anti-tag");
  else if(TString(histname).Contains("Top"))  leg->SetHeader("Top quark tag");
  leg->SetTextSize(0.03);
  leg->AddEntry(WtagFakeSF_Barrel,  							"Scale Factor, stat only", "ep");
  leg->AddEntry(WtagFakeRate_Barrel_data,  			"#epsilon_{data}, stat only", "ep");
  leg->AddEntry(WtagFakeRate_Barrel_background,  "#epsilon_{background}, stat only", "ep");
  leg->Draw();

  c1->SaveAs(Dir+histname_W+type4); 
  
  c1 = new TCanvas("c1", "", 800, 800);
  WtagFakeSF_Endcap->SetMarkerColor(kBlack);
  WtagFakeSF_Endcap->SetLineColor(kBlack);
  WtagFakeRate_Endcap_background->SetMarkerColor(kRed);
  WtagFakeRate_Endcap_background->SetLineColor(kRed);
  WtagFakeRate_Endcap_data->SetMarkerColor(kAzure);
  WtagFakeRate_Endcap_data->SetLineColor(kAzure);
  WtagFakeRate_Endcap_data->SetMarkerStyle(20);
  WtagFakeRate_Endcap_data->SetMarkerSize(1.3);
  WtagFakeRate_Endcap_background->SetMarkerStyle(20);
  WtagFakeRate_Endcap_background->SetMarkerSize(1.3);
  WtagFakeSF_Endcap->SetMarkerStyle(20);
  WtagFakeSF_Endcap->SetMarkerSize(1.3);
  WtagFakeSF_Endcap->SetMaximum(1.4);
  WtagFakeSF_Endcap->SetMinimum(0.0);
  WtagFakeSF_Endcap->Draw("EP");
  WtagFakeRate_Endcap_data->Draw("EPsame");
  WtagFakeRate_Endcap_background->Draw("EPsame");

  xmin = WtagFakeRate_Endcap_data->GetXaxis()->GetBinLowEdge(WtagFakeSF_Endcap->GetXaxis()->GetFirst());
  xmax = WtagFakeRate_Endcap_data->GetXaxis()->GetBinUpEdge(WtagFakeSF_Endcap->GetXaxis()->GetLast());

  text = "CMS #scale[0.7]{#font[52]{Work in progress 2016}}";
  cms_lat = new TLatex(xmin, 1.43, text.c_str());
  cms_lat->SetTextSize(0.04);
  cms_lat->SetLineWidth(2);
  cms_lat->Draw();
  text = "#scale[0.7]{35.9 fb^{-1} (13 TeV)}";
  era_lat = new TLatex(xmax,1.44, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.03);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();

  leg = new TLegend(0.51,0.77,0.9,0.9);
  if(TString(histname).Contains("mW")) 			leg->SetHeader("W boson mass-tag");
  else if(TString(histname).Contains("aW")) leg->SetHeader("W boson anti-tag");
  else if(TString(histname).Contains("W"))  leg->SetHeader("W boson tag");
  else if(TString(histname).Contains("mTop")) 			leg->SetHeader("Top quark mass-tag");
  else if(TString(histname).Contains("aTop")) leg->SetHeader("Top quark anti-tag");
  else if(TString(histname).Contains("Top"))  leg->SetHeader("Top quark tag");
  leg->SetTextSize(0.03);
  leg->AddEntry(WtagFakeSF_Endcap,  							"Scale Factor, stat only", "ep");
  leg->AddEntry(WtagFakeRate_Endcap_data,  			"#epsilon_{data}, stat only", "ep");
  leg->AddEntry(WtagFakeRate_Endcap_background,  "#epsilon_{background}, stat only", "ep");
  leg->Draw();

  type4 = "_endcap.png";
  c1->SaveAs(Dir+histname_W+type4); 
  
  c1 = new TCanvas("c1", "", 800, 800);
  WtagFakeSF_Barrel->SetMarkerColor(kBlack);
  WtagFakeSF_Barrel->SetLineColor(kBlack);
  WtagFakeRate_Barrel_background->SetMarkerColor(kRed);
  WtagFakeRate_Barrel_background->SetLineColor(kRed);
  WtagFakeRate_Barrel_data->SetMarkerColor(kAzure);
  WtagFakeRate_Barrel_data->SetLineColor(kAzure);
  WtagFakeSF_Endcap->SetMarkerColor(kCyan);
  WtagFakeSF_Endcap->SetLineColor(kCyan);
  WtagFakeRate_Endcap_background->SetMarkerColor(kMagenta);
  WtagFakeRate_Endcap_background->SetLineColor(kMagenta);
  WtagFakeRate_Endcap_data->SetMarkerColor(kViolet);
  WtagFakeRate_Endcap_data->SetLineColor(kViolet);
  WtagFakeSF_Barrel->SetMarkerStyle(20);
  WtagFakeSF_Barrel->SetMarkerSize(1.3);
  WtagFakeSF_Endcap->SetMarkerStyle(20);
  WtagFakeSF_Endcap->SetMarkerSize(1.3);
  WtagFakeSF_Barrel->SetMaximum(1.4);
  WtagFakeSF_Barrel->SetMinimum(0.0);
  WtagFakeSF_Barrel->Draw("EP");
  WtagFakeRate_Barrel_data->Draw("EPsame");
  WtagFakeRate_Barrel_background->Draw("EPsame");
  WtagFakeSF_Endcap->Draw("EPsame");
  WtagFakeRate_Endcap_data->Draw("EPsame");
  WtagFakeRate_Endcap_background->Draw("EPsame");

  xmin = WtagFakeRate_Barrel_data->GetXaxis()->GetBinLowEdge(WtagFakeSF_Barrel->GetXaxis()->GetFirst());
  xmax = WtagFakeRate_Barrel_data->GetXaxis()->GetBinUpEdge(WtagFakeSF_Barrel->GetXaxis()->GetLast());

  text = "CMS #scale[0.7]{#font[52]{Work in progress 2016}}";
  cms_lat = new TLatex(xmin, 1.43, text.c_str());
  cms_lat->SetTextSize(0.04);
  cms_lat->SetLineWidth(2);
  cms_lat->Draw();
  text = "#scale[0.7]{35.9 fb^{-1} (13 TeV)}";
  era_lat = new TLatex(xmax,1.44, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.03);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();

  leg = new TLegend(0.31,0.20,0.9,0.45);
  if(TString(histname).Contains("mW")) 			leg->SetHeader("W boson mass-tag");
  else if(TString(histname).Contains("aW")) leg->SetHeader("W boson anti-tag");
  else if(TString(histname).Contains("W"))  leg->SetHeader("W boson tag");
  else if(TString(histname).Contains("mTop")) 			leg->SetHeader("Top quark mass-tag");
  else if(TString(histname).Contains("aTop")) leg->SetHeader("Top quark anti-tag");
  else if(TString(histname).Contains("Top"))  leg->SetHeader("Top quark tag");
  leg->SetTextSize(0.03);
  leg->AddEntry(WtagFakeSF_Barrel,  							"Scale Factor, #eta < 1.5, stat only", "ep");
  leg->AddEntry(WtagFakeRate_Barrel_data,  			"#epsilon_{data}, #eta < 1.5, stat only", "ep");
  leg->AddEntry(WtagFakeRate_Barrel_background,  "#epsilon_{background}, #eta < 1.5, stat only", "ep");
  leg->AddEntry(WtagFakeSF_Endcap,  							"Scale Factor, 1.5 < #eta < 2.4, stat only", "ep");
  leg->AddEntry(WtagFakeRate_Endcap_data,  			"#epsilon_{data}, 1.5 < #eta < 2.4, stat only", "ep");
  leg->AddEntry(WtagFakeRate_Endcap_background,  "#epsilon_{background}, 1.5 < #eta < 2.4, stat only", "ep");
  leg->Draw();

  type4 = ".png";
  c1->SaveAs(Dir+histname_W+type4); 
  
  //TFile* output = new TFile("W_Top_TagScaleFactor.root", "update");
  //WtagFakeSF_Barrel->SetNameTitle(histname, histname);
  //WtagFakeSF_Barrel->Write(); 
   //return hWtagFakeSF_Barrel;
   //return htool;
   return 0;

}

void WTagScaleFactor(){
  TH1F* hW = DrawScaleFactor("W");
  //hW->SetNameTitle("W", "W");
  TH1F* haW = DrawScaleFactor("aW");
  //haW->SetNameTitle("aW", "aW");
  TH1F* hmW = DrawScaleFactor("mW");
  //hmW->SetNameTitle("mW", "mW");
  TH1F* hTop = DrawScaleFactor("Top");
  //hTop->SetNameTitle("Top", "Top");
  TH1F* haTop = DrawScaleFactor("aTop");
  //haTop->SetNameTitle("aTop", "aTop");
  TH1F* hmTop = DrawScaleFactor("mTop");
  //hmTop->SetNameTitle("mTop", "mTop");
  //hW->GetYaxis()->UnZoom();
  //haW->GetYaxis()->UnZoom();
  //hmW->GetYaxis()->UnZoom();
  //hTop->GetYaxis()->UnZoom();
  //haTop->GetYaxis()->UnZoom();
  //hmTop->GetYaxis()->UnZoom();
//
 // TFile* output = new TFile("WTopTagSF.root", "recreate");
 // hW->Write(); 
 // haW->Write(); 
 // hmW->Write(); 
 // hTop->Write(); 
 // haTop->Write(); 
 // hmTop->Write(); 
}
