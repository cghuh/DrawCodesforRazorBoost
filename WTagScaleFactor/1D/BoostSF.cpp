TH1F* BoostSF(TString histname){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();
  //TString dir = "/gatbawi/palgongsan/chuh/susy/susy170818/";
  TString dir = "/uscms/homes/c/chuh/work/susy170922/;
  TFile* file6 = TFile::Open(dir+"background.root");
  TFile* file7 = TFile::Open(dir+"data.root");

  TString histname_W = "AK8Jet1Pt_"+histname+"_fakerate";
  TString histname_Wo = "AK8Jet1Pt_no_"+histname+"_fakerate";
  cout << histname_W << endl;

  TH1F*       data_W  = (TH1F*)file7->Get(histname_W);
  TH1F*       data_Wo = (TH1F*)file7->Get(histname_Wo);
  TH1F* background_W  = (TH1F*)file6->Get(histname_W);
  TH1F* background_Wo = (TH1F*)file6->Get(histname_Wo);

  TH1F* WtagFakeRate_data  = (TH1F*)data_W->Clone();
  TH1F* WtagFakeRate_background  = (TH1F*)background_W->Clone();

  WtagFakeRate_data->Divide(data_Wo);
  WtagFakeRate_background->Divide(background_Wo);

  TH1F* WtagFakeSF = (TH1F*)WtagFakeRate_data->Clone();
  WtagFakeSF->Divide(WtagFakeRate_background);

  TString    name  = WtagFakeRate_data->GetName();
  TString    title = WtagFakeRate_data->GetTitle();
  Int_t    nx    = WtagFakeRate_data->GetNbinsX()+1;
  Double_t x1 = WtagFakeRate_data->GetBinLowEdge(1);
  Double_t bw = WtagFakeRate_data->GetBinWidth(nx);
  Double_t x2 = WtagFakeRate_data->GetBinLowEdge(nx)+bw;

  TH1F* hWtagFakeRate_data = new TH1F(name, title, nx, x1, x2);
  for (Int_t i=1; i<=nx; i++) {
      hWtagFakeRate_data->Fill(hWtagFakeRate_data->GetBinCenter(i), WtagFakeRate_data->GetBinContent(i));
  }

  name  = WtagFakeRate_background->GetName();
  title = WtagFakeRate_background->GetTitle();
  nx    = WtagFakeRate_background->GetNbinsX()+1;
  x1 = WtagFakeRate_background->GetBinLowEdge(1);
  bw = WtagFakeRate_background->GetBinWidth(nx);
  x2 = WtagFakeRate_background->GetBinLowEdge(nx)+bw;

  TH1F* hWtagFakeRate_background = new TH1F(name, title, nx, x1, x2);
  for (Int_t i=1; i<=nx; i++) {
      hWtagFakeRate_background->Fill(hWtagFakeRate_background->GetBinCenter(i), WtagFakeRate_background->GetBinContent(i));
  }

  name  = WtagFakeSF->GetName();
  title = WtagFakeSF->GetTitle();
  nx    = WtagFakeSF->GetNbinsX()+1;
  //x1 = WtagFakeSF->GetBinLowEdge(1);
  //bw = WtagFakeSF->GetBinWidth(nx);
  //x2 = WtagFakeSF->GetBinLowEdge(nx)+bw;
  Double_t *xbins= new Double_t[nx+1];
  for (UInt_t i=0;i<nx;i++) xbins[i]=*(WtagFakeSF->GetXaxis()->GetXbins()->GetArray()+i);
  xbins[nx]=xbins[nx-1]+WtagFakeSF->GetBinWidth(nx);

  TH1F* hWtagFakeSF = new TH1F(name, title, nx, xbins);
  for (Int_t i=1; i<=nx; i++) {
      hWtagFakeSF->SetBinContent(hWtagFakeSF->FindBin(hWtagFakeSF->GetBinCenter(i)), WtagFakeSF->GetBinContent(i));
      hWtagFakeSF->SetBinError(  hWtagFakeSF->FindBin(hWtagFakeSF->GetBinCenter(i)), WtagFakeSF->GetBinError(i));
   }

  TCanvas* cWtagFakeRate_data = new TCanvas("cWtagFakeRate_data", "", 800, 800);
  WtagFakeRate_data->SetMarkerColor(kBlack);
  WtagFakeRate_data->SetLineColor(kBlack);
  WtagFakeRate_data->SetMarkerStyle(20);
  WtagFakeRate_data->SetMarkerSize(1.3);
  WtagFakeRate_data->SetMaximum(0.2);
  if(TString(histname_W).Contains("_Top")) WtagFakeRate_data->SetMaximum(0.01);
  WtagFakeRate_data->SetMinimum(0);
  //WtagFakeRate_data->GetYaxis()->SetTitle("W tagging efficiency");
  WtagFakeRate_data->Draw("EP");

  float xmin = WtagFakeRate_data->GetXaxis()->GetBinLowEdge(WtagFakeRate_data->GetXaxis()->GetFirst());
  float xmax = WtagFakeRate_data->GetXaxis()->GetBinUpEdge(WtagFakeRate_data->GetXaxis()->GetLast());

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
  leg->AddEntry(WtagFakeRate_data,  "stat only", "ep");
  leg->Draw();
  
  TCanvas* cWtagFakeRate_background = new TCanvas("cWtagFakeRate_background", "", 800, 800);
  WtagFakeRate_background->SetMarkerColor(kBlack);
  WtagFakeRate_background->SetLineColor(kBlack);
  WtagFakeRate_background->SetMarkerStyle(20);
  WtagFakeRate_background->SetMarkerSize(1.3);
  WtagFakeRate_background->SetMaximum(0.2);
  if(TString(histname_W).Contains("_Top")) WtagFakeRate_background->SetMaximum(0.01);
  WtagFakeRate_background->SetMinimum(0);
  //WtagFakeRate_background->GetYaxis()->SetTitle("W tagging efficiency");
  WtagFakeRate_background->Draw("EP");

  xmin = WtagFakeRate_data->GetXaxis()->GetBinLowEdge(WtagFakeRate_background->GetXaxis()->GetFirst());
  xmax = WtagFakeRate_data->GetXaxis()->GetBinUpEdge(WtagFakeRate_background->GetXaxis()->GetLast());

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
  leg->AddEntry(WtagFakeRate_background,  "stat only", "ep");
  leg->Draw();
  
  TCanvas* cWtagFakeSF = new TCanvas("cWtagFakeSF", "", 800, 800);
  WtagFakeSF->SetMarkerColor(kBlack);
  WtagFakeSF->SetLineColor(kBlack);
  WtagFakeSF->SetMarkerStyle(20);
  WtagFakeSF->SetMarkerSize(1.3);
  WtagFakeSF->SetMaximum(1.4);
  WtagFakeSF->SetMinimum(0.5);
  WtagFakeSF->Draw("EP");

  xmin = WtagFakeRate_data->GetXaxis()->GetBinLowEdge(WtagFakeSF->GetXaxis()->GetFirst());
  xmax = WtagFakeRate_data->GetXaxis()->GetBinUpEdge(WtagFakeSF->GetXaxis()->GetLast());

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
  leg->AddEntry(WtagFakeSF,  						"Scale Factor, stat only", "ep");
  leg->Draw();

  TString Dir = "./plot/";
  TString type1 = "_data.png";
  TString type2 = "_background.png";
  TString type3 = "_SF.png";
  TString type4 = ".png";
  cWtagFakeRate_data->SaveAs(Dir+histname_W+type1); 
  cWtagFakeRate_background->SaveAs(Dir+histname_W+type2); 
  cWtagFakeSF->SaveAs(Dir+histname_W+type3); 


  TCanvas* c1 = new TCanvas("c1", "", 800, 800);
  WtagFakeSF->SetMarkerColor(kBlack);
  WtagFakeSF->SetLineColor(kBlack);
  WtagFakeRate_background->SetMarkerColor(kRed);
  WtagFakeRate_background->SetLineColor(kRed);
  WtagFakeRate_data->SetMarkerColor(kAzure);
  WtagFakeRate_data->SetLineColor(kAzure);
  WtagFakeSF->SetMarkerStyle(20);
  WtagFakeSF->SetMarkerSize(1.3);
  WtagFakeSF->SetMaximum(1.4);
  WtagFakeSF->SetMinimum(0.0);
  WtagFakeSF->Draw("EP");
  WtagFakeRate_data->Draw("EPsame");
  WtagFakeRate_background->Draw("EPsame");

  xmin = WtagFakeRate_data->GetXaxis()->GetBinLowEdge(WtagFakeSF->GetXaxis()->GetFirst());
  xmax = WtagFakeRate_data->GetXaxis()->GetBinUpEdge(WtagFakeSF->GetXaxis()->GetLast());

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
  leg->AddEntry(WtagFakeSF,  							"Scale Factor, stat only", "ep");
  leg->AddEntry(WtagFakeRate_data,  			"#epsilon_{data}, stat only", "ep");
  leg->AddEntry(WtagFakeRate_background,  "#epsilon_{background}, stat only", "ep");
  leg->Draw();

  c1->SaveAs(Dir+histname_W+type4); 
  
  //TFile* output = new TFile("W_Top_TagScaleFactor.root", "update");
  //WtagFakeSF->SetNameTitle(histname, histname);
  //WtagFakeSF->Write(); 
   return hWtagFakeSF;

}

void WTagScaleFactor(){
  TH1F* hW = DrawScaleFactor("W");
  hW->SetNameTitle("W", "W");
  TH1F* haW = DrawScaleFactor("aW");
  haW->SetNameTitle("aW", "aW");
  TH1F* hmW = DrawScaleFactor("mW");
  hmW->SetNameTitle("mW", "mW");
  TH1F* hTop = DrawScaleFactor("Top");
  hTop->SetNameTitle("Top", "Top");
  TH1F* haTop = DrawScaleFactor("aTop");
  haTop->SetNameTitle("aTop", "aTop");
  TH1F* hmTop = DrawScaleFactor("mTop");
  hmTop->SetNameTitle("mTop", "mTop");
  hW->GetYaxis()->UnZoom();
  haW->GetYaxis()->UnZoom();
  hmW->GetYaxis()->UnZoom();
  hTop->GetYaxis()->UnZoom();
  haTop->GetYaxis()->UnZoom();
  hmTop->GetYaxis()->UnZoom();

  TFile* output = new TFile("WTopTagSF.root", "recreate");
  hW->Write(); 
  haW->Write(); 
  hmW->Write(); 
  hTop->Write(); 
  haTop->Write(); 
  hmTop->Write(); 
}
