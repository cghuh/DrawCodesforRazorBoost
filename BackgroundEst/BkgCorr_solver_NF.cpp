TH1D* Correction(TString period, TString region){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();
  TString dir = "/Users/huhchanggi/temp/200423/";
  TFile* file0 = TFile::Open(dir+"run_2020_04_28.root");

  TString histname[3][9];
  file0->cd("Counts_vs_MRR2Bins/Syst_vs_MRR2Bins");
  histname[0][0] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Data_JetHTMET_"+period+"_CR_QCD17_"+region;
  histname[0][1] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Multijet_JetHTMET_"+period+"_CR_QCD17_"+region;
  histname[0][2] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Top_JetHTMET_"+period+"_CR_QCD17_"+region;
  histname[0][3] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/TT_powheg_pythia8_JetHTMET_"+period+"_CR_QCD17_"+region;
  histname[0][4] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/WToLNu_JetHTMET_"+period+"_CR_QCD17_"+region;
  histname[0][5] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/ZToNuNu_JetHTMET_"+period+"_CR_QCD17_"+region;
  histname[0][6] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Multiboson_JetHTMET_"+period+"_CR_QCD17_"+region;
  histname[0][7] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/GJets_JetHTMET_"+period+"_CR_QCD17_"+region;
  histname[0][8] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/DYToLL_JetHTMET_"+period+"_CR_QCD17_"+region;
  histname[1][0] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Data_JetHTMET_"+period+"_CR_Top17_"+region;
  histname[1][1] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Multijet_JetHTMET_"+period+"_CR_Top17_"+region;
  histname[1][2] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Top_JetHTMET_"+period+"_CR_Top17_"+region;
  histname[1][3] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/TT_powheg_pythia8_JetHTMET_"+period+"_CR_Top17_"+region;
  histname[1][4] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/WToLNu_JetHTMET_"+period+"_CR_Top17_"+region;
  histname[1][5] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/ZToNuNu_JetHTMET_"+period+"_CR_Top17_"+region;
  histname[1][6] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Multiboson_JetHTMET_"+period+"_CR_Top17_"+region;
  histname[1][7] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/GJets_JetHTMET_"+period+"_CR_Top17_"+region;
  histname[1][8] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/DYToLL_JetHTMET_"+period+"_CR_Top17_"+region;
  histname[2][0] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Data_JetHTMET_"+period+"_CR_W17_"+region;
  histname[2][1] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Multijet_JetHTMET_"+period+"_CR_W17_"+region;
  histname[2][2] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Top_JetHTMET_"+period+"_CR_W17_"+region;
  histname[2][3] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/TT_powheg_pythia8_JetHTMET_"+period+"_CR_W17_"+region;
  histname[2][4] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/WToLNu_JetHTMET_"+period+"_CR_W17_"+region;
  histname[2][5] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/ZToNuNu_JetHTMET_"+period+"_CR_W17_"+region;
  histname[2][6] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Multiboson_JetHTMET_"+period+"_CR_W17_"+region;
  histname[2][7] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/GJets_JetHTMET_"+period+"_CR_W17_"+region;
  histname[2][8] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/DYToLL_JetHTMET_"+period+"_CR_W17_"+region;

  TH2D* bkg2D[3][9];
  TH1D* bkg[3][9];
  TH1D* top[3];
  TH1D* CF[3];

  for(int i=0; i<3; i++){
    for(int j=0; j<9; j++) bkg2D[i][j] = (TH2D*)file0->Get(histname[i][j]);
  }
  int binsize = bkg2D[0][0]->GetNbinsX();
  double bin[binsize+1];
  string name;

  for(int i=1;i<=bkg2D[0][0]->GetNbinsX()+1;i++) bin[i-1] = bkg2D[0][0]->GetXaxis()->GetBinLowEdge(i);

  for(int i=0; i<3; i++) {
    for(int j=0; j<9; j++) {
      name = "bkg_"+to_string(i)+to_string(j);
      bkg[i][j] = new TH1D(name.c_str(), "", binsize, bin);
      //bkg[i][j] = (TH1D*)bkg2D[0][0]->ProjectionX("temp");
      for(int k=1;k<=bkg2D[0][0]->GetNbinsX();k++) {
        if(bkg2D[i][j] == NULL) {
          bkg[i][j]->SetBinContent(k, 0);
          bkg[i][j]->SetBinError(k, 0);
        }
        else {
          bkg[i][j]->SetBinContent(k, bkg2D[i][j]->GetBinContent(k,1));
          bkg[i][j]->SetBinError(k, bkg2D[i][j]->GetBinError(k,1));
        }
      }
    }
    bkg[i][2]->Add(bkg[i][3]);
    CF[i] = (TH1D*)bkg[i][0]->Clone();
    for(int j=5; j<9; j++){
      CF[i]->Add(bkg[i][j], -1);
    }
  }

  TMatrix D(3,3);
  TMatrix Q(3,3);
  TMatrix T(3,3);
  TMatrix W(3,3);
  double bnR2MR[5] = {0.,200.,400,600.,3000.};
  TH1D* CF_Q = new TH1D("CF_Q", "", binsize,bin);
  TH1D* CF_T = new TH1D("CF_T", "", binsize,bin);
  TH1D* CF_W = new TH1D("CF_W", "", binsize,bin);

  for(int i=1; i<= CF_Q->GetNbinsX();i++){
    for(int j=0; j<3;j++){
      D[j][0] = bkg[j][1]->GetBinContent(i);
      D[j][1] = bkg[j][2]->GetBinContent(i);
      D[j][2] = bkg[j][4]->GetBinContent(i);
      Q[j][0] = CF[j]->GetBinContent(i);
      Q[j][1] = bkg[j][2]->GetBinContent(i);
      Q[j][2] = bkg[j][4]->GetBinContent(i);
      T[j][0] = bkg[j][1]->GetBinContent(i);
      T[j][1] = CF[j]->GetBinContent(i);
      T[j][2] = bkg[j][4]->GetBinContent(i);
      W[j][0] = bkg[j][1]->GetBinContent(i);
      W[j][1] = bkg[j][2]->GetBinContent(i);
      W[j][2] = CF[j]->GetBinContent(i);
    }
    CF_Q->SetBinContent(i, Q.Determinant()/D.Determinant());
    CF_T->SetBinContent(i, T.Determinant()/D.Determinant());
    CF_W->SetBinContent(i, W.Determinant()/D.Determinant());
  }

  TCanvas* c1 = new TCanvas("c1","",700,700);
  
  c1->SetGridx();
  c1->SetGridy();
  CF_Q->SetMaximum(2);
  CF_Q->SetMinimum(-1);
  CF_Q->GetYaxis()->SetTitle("Correction Factor");
  CF_Q->GetYaxis()->SetTitleOffset(0.9);
  CF_Q->SetMarkerStyle(21);
  CF_T->SetMarkerStyle(21);
  CF_W->SetMarkerStyle(21);
  CF_T->SetMarkerColor(kRed);
  CF_W->SetMarkerColor(kBlue);
  CF_Q->Draw("P");
  CF_T->Draw("Psame");
  CF_W->Draw("Psame");

  auto leg = new TLegend(0.45,0.77,0.9,0.9);
  leg->SetHeader("boost jet final state");
  if(TString(region).Contains("1boost"))     leg->SetHeader("1 boost jet final state");
  else if(TString(region).Contains("2boost")) leg->SetHeader("#geq 2 boost jet final state");
  else leg->SetHeader("boost jet final state");
  leg->AddEntry(CF_Q,  "Multijet CF", "p");
  leg->AddEntry(CF_T,  "Top(TT+ST) CF", "p");
  leg->AddEntry(CF_W,  "W(#rightarrowl#nu)+jet CF", "p");
  leg->Draw();


  float xmin = CF_Q->GetXaxis()->GetBinLowEdge(CF_Q->GetXaxis()->GetFirst());
  float xmax = CF_Q->GetXaxis()->GetBinUpEdge( CF_Q->GetXaxis()->GetLast());
 
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

  c1->SaveAs("plot/"+period+"CF_"+region+".png");
  return CF_Q;
}

TH1D* NjetCorrection(TString period, TString region){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();
  TString dir = "/Users/huhchanggi/temp/200423/";
  TFile* file0 = TFile::Open(dir+"run_2020_04_28.root");

  TString histname[3][9];
  file0->cd("Counts_vs_NJetBins/Syst_vs_NJetBins");
  histname[0][0] = "/Counts_vs_NJetBins/Syst_vs_NJetBins/Data_JetHTMET_"+period+"_CR_QCD17_"+region;
  histname[0][1] = "/Counts_vs_NJetBins/Syst_vs_NJetBins/Multijet_JetHTMET_"+period+"_CR_QCD17_"+region;
  histname[0][2] = "/Counts_vs_NJetBins/Syst_vs_NJetBins/Top_JetHTMET_"+period+"_CR_QCD17_"+region;
  histname[0][3] = "/Counts_vs_NJetBins/Syst_vs_NJetBins/TT_powheg_pythia8_JetHTMET_"+period+"_CR_QCD17_"+region;
  histname[0][4] = "/Counts_vs_NJetBins/Syst_vs_NJetBins/WToLNu_JetHTMET_"+period+"_CR_QCD17_"+region;
  histname[0][5] = "/Counts_vs_NJetBins/Syst_vs_NJetBins/ZToNuNu_JetHTMET_"+period+"_CR_QCD17_"+region;
  histname[0][6] = "/Counts_vs_NJetBins/Syst_vs_NJetBins/Multiboson_JetHTMET_"+period+"_CR_QCD17_"+region;
  histname[0][7] = "/Counts_vs_NJetBins/Syst_vs_NJetBins/GJets_JetHTMET_"+period+"_CR_QCD17_"+region;
  histname[0][8] = "/Counts_vs_NJetBins/Syst_vs_NJetBins/DYToLL_JetHTMET_"+period+"_CR_QCD17_"+region;
  histname[1][0] = "/Counts_vs_NJetBins/Syst_vs_NJetBins/Data_JetHTMET_"+period+"_CR_Top17_"+region;
  histname[1][1] = "/Counts_vs_NJetBins/Syst_vs_NJetBins/Multijet_JetHTMET_"+period+"_CR_Top17_"+region;
  histname[1][2] = "/Counts_vs_NJetBins/Syst_vs_NJetBins/Top_JetHTMET_"+period+"_CR_Top17_"+region;
  histname[1][3] = "/Counts_vs_NJetBins/Syst_vs_NJetBins/TT_powheg_pythia8_JetHTMET_"+period+"_CR_Top17_"+region;
  histname[1][4] = "/Counts_vs_NJetBins/Syst_vs_NJetBins/WToLNu_JetHTMET_"+period+"_CR_Top17_"+region;
  histname[1][5] = "/Counts_vs_NJetBins/Syst_vs_NJetBins/ZToNuNu_JetHTMET_"+period+"_CR_Top17_"+region;
  histname[1][6] = "/Counts_vs_NJetBins/Syst_vs_NJetBins/Multiboson_JetHTMET_"+period+"_CR_Top17_"+region;
  histname[1][7] = "/Counts_vs_NJetBins/Syst_vs_NJetBins/GJets_JetHTMET_"+period+"_CR_Top17_"+region;
  histname[1][8] = "/Counts_vs_NJetBins/Syst_vs_NJetBins/DYToLL_JetHTMET_"+period+"_CR_Top17_"+region;
  histname[2][0] = "/Counts_vs_NJetBins/Syst_vs_NJetBins/Data_JetHTMET_"+period+"_CR_W17_"+region;
  histname[2][1] = "/Counts_vs_NJetBins/Syst_vs_NJetBins/Multijet_JetHTMET_"+period+"_CR_W17_"+region;
  histname[2][2] = "/Counts_vs_NJetBins/Syst_vs_NJetBins/Top_JetHTMET_"+period+"_CR_W17_"+region;
  histname[2][3] = "/Counts_vs_NJetBins/Syst_vs_NJetBins/TT_powheg_pythia8_JetHTMET_"+period+"_CR_W17_"+region;
  histname[2][4] = "/Counts_vs_NJetBins/Syst_vs_NJetBins/WToLNu_JetHTMET_"+period+"_CR_W17_"+region;
  histname[2][5] = "/Counts_vs_NJetBins/Syst_vs_NJetBins/ZToNuNu_JetHTMET_"+period+"_CR_W17_"+region;
  histname[2][6] = "/Counts_vs_NJetBins/Syst_vs_NJetBins/Multiboson_JetHTMET_"+period+"_CR_W17_"+region;
  histname[2][7] = "/Counts_vs_NJetBins/Syst_vs_NJetBins/GJets_JetHTMET_"+period+"_CR_W17_"+region;
  histname[2][8] = "/Counts_vs_NJetBins/Syst_vs_NJetBins/DYToLL_JetHTMET_"+period+"_CR_W17_"+region;

  TH2D* bkg2D[3][9];
  TH1D* bkg[3][9];
  TH1D* top[3];
  TH1D* CF[3];

  for(int i=0; i<3; i++){
    for(int j=0; j<9; j++) bkg2D[i][j] = (TH2D*)file0->Get(histname[i][j]);
  }
  int binsize = bkg2D[0][0]->GetNbinsX();
  double bin[binsize+1];
  string name;

  for(int i=1;i<=bkg2D[0][0]->GetNbinsX()+1;i++) bin[i-1] = bkg2D[0][0]->GetXaxis()->GetBinLowEdge(i);

  for(int i=0; i<3; i++) {
    for(int j=0; j<9; j++) {
      name = "bkg_"+to_string(i)+to_string(j);
      bkg[i][j] = new TH1D(name.c_str(), "", binsize, bin);
      for(int k=1;k<=bkg2D[0][0]->GetNbinsX();k++) {
        if(bkg2D[i][j] == NULL) {
          cout << "NULL " << i << ", " << j << endl;
          bkg[i][j]->SetBinContent(k, 0);
          bkg[i][j]->SetBinError(k, 0);
        }
        else {
          bkg[i][j]->SetBinContent(k, bkg2D[i][j]->GetBinContent(k,1));
          bkg[i][j]->SetBinError(k, bkg2D[i][j]->GetBinError(k,1));
        }
      }
    }
    bkg[i][2]->Add(bkg[i][3]);
    CF[i] = (TH1D*)bkg[i][0]->Clone();
    for(int j=5; j<9; j++){
      CF[i]->Add(bkg[i][j], -1);
    }
  }

  TMatrix D(3,3);
  TMatrix Q(3,3);
  TMatrix T(3,3);
  TMatrix W(3,3);
  double bnR2MR[5] = {0.,200.,400,600.,3000.};
  TH1D* CF_Q = new TH1D("CF_Q", "", binsize,bin);
  TH1D* CF_T = new TH1D("CF_T", "", binsize,bin);
  TH1D* CF_W = new TH1D("CF_W", "", binsize,bin);

  for(int i=1; i<= CF_Q->GetNbinsX();i++){
    for(int j=0; j<3;j++){
      D[j][0] = bkg[j][1]->GetBinContent(i);
      D[j][1] = bkg[j][2]->GetBinContent(i);
      D[j][2] = bkg[j][4]->GetBinContent(i);
      Q[j][0] = CF[j]->GetBinContent(i);
      Q[j][1] = bkg[j][2]->GetBinContent(i);
      Q[j][2] = bkg[j][4]->GetBinContent(i);
      T[j][0] = bkg[j][1]->GetBinContent(i);
      T[j][1] = CF[j]->GetBinContent(i);
      T[j][2] = bkg[j][4]->GetBinContent(i);
      W[j][0] = bkg[j][1]->GetBinContent(i);
      W[j][1] = bkg[j][2]->GetBinContent(i);
      W[j][2] = CF[j]->GetBinContent(i);
    }
    CF_Q->SetBinContent(i, Q.Determinant()/D.Determinant());
    CF_T->SetBinContent(i, T.Determinant()/D.Determinant());
    CF_W->SetBinContent(i, W.Determinant()/D.Determinant());
  }

  TCanvas* c1 = new TCanvas("c1","",700,700);
  
  c1->SetGridx();
  c1->SetGridy();
  CF_Q->SetMaximum(2);
  CF_Q->SetMinimum(-1);
  CF_Q->GetYaxis()->SetTitle("Correction Factor");
  CF_Q->GetYaxis()->SetTitleOffset(0.9);
  CF_Q->SetMarkerStyle(21);
  CF_T->SetMarkerStyle(21);
  CF_W->SetMarkerStyle(21);
  CF_T->SetMarkerColor(kRed);
  CF_W->SetMarkerColor(kBlue);
  CF_Q->Draw("P");
  CF_T->Draw("Psame");
  CF_W->Draw("Psame");

  auto leg = new TLegend(0.45,0.77,0.9,0.9);
  leg->SetHeader("boost jet final state");
  if(TString(region).Contains("1boost"))     leg->SetHeader("1 boost jet final state");
  else if(TString(region).Contains("2boost")) leg->SetHeader("#geq 2 boost jet final state");
  else leg->SetHeader("boost jet final state");
  leg->AddEntry(CF_Q,  "Multijet CF", "p");
  leg->AddEntry(CF_T,  "Top(TT+ST) CF", "p");
  leg->AddEntry(CF_W,  "W(#rightarrowl#nu)+jet CF", "p");
  leg->Draw();


  float xmin = CF_Q->GetXaxis()->GetBinLowEdge(CF_Q->GetXaxis()->GetFirst());
  float xmax = CF_Q->GetXaxis()->GetBinUpEdge( CF_Q->GetXaxis()->GetLast());
 
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

  c1->SaveAs("plot/"+period+"NJet_CF_"+region+".png");
  return CF_Q;

}

void BkgCorr_solver_NF(){
  TH1D* h1[10];
  h1[1] = Correction("2016", "1Boost");
  h1[2] = Correction("2016", "2Boost");
  h1[1] = Correction("2017", "1Boost");
  h1[2] = Correction("2017", "2Boost");
  h1[1] = Correction("2018", "1Boost");
  h1[2] = Correction("2018", "2Boost");
  h1[1] = NjetCorrection("2016", "1Boost");
  h1[2] = NjetCorrection("2016", "2Boost");
  h1[1] = NjetCorrection("2017", "1Boost");
  h1[2] = NjetCorrection("2017", "2Boost");
  h1[1] = NjetCorrection("2018", "1Boost");
  h1[2] = NjetCorrection("2018", "2Boost");
/*
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
*/
}
