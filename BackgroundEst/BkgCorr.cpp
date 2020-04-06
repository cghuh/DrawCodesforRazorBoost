TH1D* Correction(TString sample, TString region){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();
  TString dir = "/eos/cms/store/user/chuh/RazorBoost/191127/added/";
  TFile* file1 = TFile::Open(dir+"Multijet.root");
  TFile* file2 = TFile::Open(dir+"TTST.root");
  TFile* file3 = TFile::Open(dir+"WJet.root");
  TFile* file4 = TFile::Open(dir+"ZJet.root");
  TFile* file5 = TFile::Open(dir+"DYToLL.root");
  TFile* file6 = TFile::Open(dir+"Multiboson+TTX.root");
  TFile* file7 = TFile::Open(dir+"GJet.root");
  TFile* file0 = TFile::Open(dir+"data.root");

  TString histname;
  if(sample == "TTST")      histname = "R2MR_T"+region;
  else if(sample == "WJet") histname = "R2MR_W"+region;
  else                      histname = "R2MR_Q"+region;

  TH1D* bkg[8];
  bkg[0] = (TH1D*)file0->Get(histname);
  bkg[1] = (TH1D*)file1->Get(histname);
  bkg[2] = (TH1D*)file2->Get(histname);
  bkg[3] = (TH1D*)file3->Get(histname);
  bkg[4] = (TH1D*)file4->Get(histname);
  bkg[5] = (TH1D*)file5->Get(histname);
  bkg[6] = (TH1D*)file6->Get(histname);
  bkg[7] = (TH1D*)file7->Get(histname);

  TH1D* CF = (TH1D*)bkg[0]->Clone();
  for(int i=1;i<8;i++){
    if     (sample == "Multijet" && i==1) continue;
    else if(sample == "TTST"     && i==2) continue;
    else if(sample == "WJet"     && i==3) continue;
    //cout << i << endl;
    CF->Add(bkg[i], -1);
  }
  if     (sample == "Multijet") CF->Divide(bkg[1]);
  else if(sample == "TTST"    ) CF->Divide(bkg[2]);
  else if(sample == "WJet"    ) CF->Divide(bkg[3]);

  TCanvas* c1 = new TCanvas("c1","",700,700);
  
  c1->SetGridx();
  c1->SetGridy();
  CF->SetMaximum(2);
  CF->SetMinimum(0);
  CF->GetYaxis()->SetTitle("Correction Factor");
  CF->GetYaxis()->SetTitleOffset(0.9);
  CF->Draw("EP");

  float xmin = CF->GetXaxis()->GetBinLowEdge(CF->GetXaxis()->GetFirst());
  float xmax = CF->GetXaxis()->GetBinUpEdge( CF->GetXaxis()->GetLast());
 
  string text = "CMS #scale[0.7]{#font[52]{Work in progress 2017}}";
  TLatex* cms_lat = new TLatex(xmin, 2.05, text.c_str());
  cms_lat->SetTextSize(0.04);
  cms_lat->SetLineWidth(2);
  cms_lat->Draw();
  text = "#scale[0.7]{41.5 fb^{-1} (13 TeV)}";
  TLatex* era_lat = new TLatex(xmax,2.1, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.03);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();

  c1->SaveAs("plot/CF_"+sample+region+".png");
  return CF;
}

TH1D* NjetCorrection(TString sample, TString region){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();
  TString dir = "/eos/cms/store/user/chuh/RazorBoost/191210/added/";
  TFile* file1 = TFile::Open(dir+"Multijet.root");
  TFile* file2 = TFile::Open(dir+"TTST.root");
  TFile* file3 = TFile::Open(dir+"WJet.root");
  TFile* file4 = TFile::Open(dir+"ZJet.root");
  TFile* file5 = TFile::Open(dir+"DYToLL.root");
  TFile* file6 = TFile::Open(dir+"Multiboson+TTX.root");
  TFile* file7 = TFile::Open(dir+"GJet.root");
  TFile* file0 = TFile::Open(dir+"data.root");

  TString histname;
  if(sample == "TTST")      histname = "njet_T"+region;
  else if(sample == "WJet") histname = "njet_W"+region;
  else                      histname = "njet_Q"+region;

  TH1D* bkg[8];
  bkg[0] = (TH1D*)file0->Get(histname);
  bkg[1] = (TH1D*)file1->Get(histname);
  bkg[2] = (TH1D*)file2->Get(histname);
  bkg[3] = (TH1D*)file3->Get(histname);
  bkg[4] = (TH1D*)file4->Get(histname);
  bkg[5] = (TH1D*)file5->Get(histname);
  bkg[6] = (TH1D*)file6->Get(histname);
  bkg[7] = (TH1D*)file7->Get(histname);

  TH1D* CF = (TH1D*)bkg[0]->Clone();
  for(int i=1;i<8;i++){
    if     (sample == "Multijet" && i==1) continue;
    else if(sample == "TTST"     && i==2) continue;
    else if(sample == "WJet"     && i==3) continue;
    //cout << i << endl;
    CF->Add(bkg[i], -1);
  }
  if     (sample == "Multijet") CF->Divide(bkg[1]);
  else if(sample == "TTST"    ) CF->Divide(bkg[2]);
  else if(sample == "WJet"    ) CF->Divide(bkg[3]);

  TCanvas* c1 = new TCanvas("c1","",700,700);
  
  c1->SetGridx();
  c1->SetGridy();
  CF->SetMaximum(2);
  CF->SetMinimum(0);
  CF->GetYaxis()->SetTitle("Correction Factor");
  CF->GetYaxis()->SetTitleOffset(0.9);
  CF->Draw("EP");

  float xmin = CF->GetXaxis()->GetBinLowEdge(CF->GetXaxis()->GetFirst());
  float xmax = CF->GetXaxis()->GetBinUpEdge( CF->GetXaxis()->GetLast());
 
  string text = "CMS #scale[0.7]{#font[52]{Work in progress 2017}}";
  TLatex* cms_lat = new TLatex(xmin, 2.05, text.c_str());
  cms_lat->SetTextSize(0.04);
  cms_lat->SetLineWidth(2);
  cms_lat->Draw();
  text = "#scale[0.7]{41.5 fb^{-1} (13 TeV)}";
  TLatex* era_lat = new TLatex(xmax,2.1, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.03);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();

  c1->SaveAs("plot/CF_Njet_"+sample+region+".png");
  return CF;
}

void BkgCorr(){
  TH1D* h1[10];
  h1[0] = Correction("TTST", "");
  h1[1] = Correction("TTST", "_nj45");
  h1[2] = Correction("TTST", "_nj6");
  h1[3] = Correction("WJet", "");
  h1[4] = Correction("WJet", "_nj45");
  h1[5] = Correction("WJet", "_nj6");
  h1[6] = Correction("Multijet", "");
  h1[7] = Correction("Multijet", "_nj45");
  h1[8] = Correction("Multijet", "_nj6");
  TFile* output = new TFile("CF_T.root", "recreate");
  for(int i=0;i<3;i++) h1[i]->Write();
  output = new TFile("CF_W.root", "recreate");
  for(int i=3;i<6;i++) h1[i]->Write();
  output = new TFile("CF_Q.root", "recreate");
  for(int i=6;i<9;i++) h1[i]->Write();
  TH1D* h2[3];
  h2[0] = NjetCorrection("TTST", "");
  h2[1] = NjetCorrection("WJet", "");
  h2[2] = NjetCorrection("Multijet", "");
  output = new TFile("CF_Njet.root", "recreate");
  for(int i=0;i<3;i++) h2[i]->Write();
}
