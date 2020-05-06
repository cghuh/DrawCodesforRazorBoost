tuple<TH1D*, TH1D*, TH1D*> Correction(TString period, TString region, TString obj){
  delete gROOT->FindObject("c1");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();
  TString dir = "/Users/huhchanggi/temp/200423/";
  TFile* file0 = TFile::Open(dir+"run_2020_04_28.root");

  TString histname[3][9];
  file0->cd("Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins");
  histname[0][0] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Data_JetHTMET_"+period+"_CR_QCD17_"+region;
  histname[0][1] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Multijet_JetHTMET_"+period+"_CR_QCD17_"+region;
  histname[0][2] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Top_JetHTMET_"+period+"_CR_QCD17_"+region;
  histname[0][3] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/TT_powheg_pythia8_JetHTMET_"+period+"_CR_QCD17_"+region;
  histname[0][4] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/WToLNu_JetHTMET_"+period+"_CR_QCD17_"+region;
  histname[0][5] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/ZToNuNu_JetHTMET_"+period+"_CR_QCD17_"+region;
  histname[0][6] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Multiboson_JetHTMET_"+period+"_CR_QCD17_"+region;
  histname[0][7] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/GJets_JetHTMET_"+period+"_CR_QCD17_"+region;
  histname[0][8] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/DYToLL_JetHTMET_"+period+"_CR_QCD17_"+region;
  histname[1][0] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Data_JetHTMET_"+period+"_CR_Top17_"+region;
  histname[1][1] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Multijet_JetHTMET_"+period+"_CR_Top17_"+region;
  histname[1][2] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Top_JetHTMET_"+period+"_CR_Top17_"+region;
  histname[1][3] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/TT_powheg_pythia8_JetHTMET_"+period+"_CR_Top17_"+region;
  histname[1][4] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/WToLNu_JetHTMET_"+period+"_CR_Top17_"+region;
  histname[1][5] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/ZToNuNu_JetHTMET_"+period+"_CR_Top17_"+region;
  histname[1][6] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Multiboson_JetHTMET_"+period+"_CR_Top17_"+region;
  histname[1][7] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/GJets_JetHTMET_"+period+"_CR_Top17_"+region;
  histname[1][8] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/DYToLL_JetHTMET_"+period+"_CR_Top17_"+region;
  histname[2][0] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Data_JetHTMET_"+period+"_CR_W17_"+region;
  histname[2][1] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Multijet_JetHTMET_"+period+"_CR_W17_"+region;
  histname[2][2] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Top_JetHTMET_"+period+"_CR_W17_"+region;
  histname[2][3] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/TT_powheg_pythia8_JetHTMET_"+period+"_CR_W17_"+region;
  histname[2][4] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/WToLNu_JetHTMET_"+period+"_CR_W17_"+region;
  histname[2][5] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/ZToNuNu_JetHTMET_"+period+"_CR_W17_"+region;
  histname[2][6] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Multiboson_JetHTMET_"+period+"_CR_W17_"+region;
  histname[2][7] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/GJets_JetHTMET_"+period+"_CR_W17_"+region;
  histname[2][8] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/DYToLL_JetHTMET_"+period+"_CR_W17_"+region;

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
  string name_Q = "CF_Q_"+(string)period+"_"+(string)region;
  string name_T = "CF_T_"+(string)period+"_"+(string)region;
  string name_W = "CF_W_"+(string)period+"_"+(string)region;
  TH1D* CF_Q = new TH1D(name_Q.c_str(), "", binsize,bin);
  TH1D* CF_T = new TH1D(name_T.c_str(), "", binsize,bin);
  TH1D* CF_W = new TH1D(name_W.c_str(), "", binsize,bin);

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
  CF_Q->SetMinimum(0);
  CF_Q->GetYaxis()->SetTitle("Correction Factor");
  CF_Q->GetXaxis()->SetTitle("M_{R} #times R^{2} [GeV]");
  if(TString(obj).Contains("NJet")) CF_Q->GetXaxis()->SetTitle("N_{Jets}");
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
  if(TString(region).Contains("1Boost"))     leg->SetHeader("1 boost jet final state");
  else if(TString(region).Contains("2Boost")) leg->SetHeader("#geq 2 boost jet final state");
  else leg->SetHeader("boost jet final state");
  leg->AddEntry(CF_Q,  "Multijet CF", "p");
  leg->AddEntry(CF_T,  "Top(TT+ST) CF", "p");
  leg->AddEntry(CF_W,  "W(#rightarrowl#nu)+jet CF", "p");
  leg->Draw();


  float xmin = CF_Q->GetXaxis()->GetBinLowEdge(CF_Q->GetXaxis()->GetFirst());
  float xmax = CF_Q->GetXaxis()->GetBinUpEdge( CF_Q->GetXaxis()->GetLast());
 
  string text;
  if(TString(region).Contains("1Boost"))      text = "CMS #scale[0.7]{#font[52]{Work in progress 2016}}";
  else if(TString(region).Contains("2Boost")) text = "CMS #scale[0.7]{#font[52]{Work in progress 2017}}";
  else                                        text = "CMS #scale[0.7]{#font[52]{Work in progress 2018}}";
  TLatex* cms_lat = new TLatex(xmin, 2.05, text.c_str());
  cms_lat->SetTextSize(0.04);
  cms_lat->SetLineWidth(2);
  cms_lat->Draw();
  if(TString(region).Contains("1Boost"))      text = "#scale[0.7]{35.9 fb^{-1} (13 TeV)}";
  else if(TString(region).Contains("2Boost")) text = "#scale[0.7]{41.5 fb^{-1} (13 TeV)}";
  else                                        text = "#scale[0.7]{59.7 fb^{-1} (13 TeV)}";
  TLatex* era_lat = new TLatex(xmax,2.1, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.03);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();

  c1->SaveAs("plot/"+period+obj+"CF_"+region+".png");
  return make_tuple(CF_Q, CF_T, CF_W);
}

tuple<TH1D*, TH1D*> LCorrection(TString period, TString region, TString obj){
  delete gROOT->FindObject("c1");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();
  TString dir = "/Users/huhchanggi/temp/200423/";
  TFile* file0 = TFile::Open(dir+"run_2020_04_28.root");

  TString histname[2][9];
  file0->cd("Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins");
  histname[0][0] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Data_JetHTMET_"+period+"_CR_L17_"+region;
  histname[0][1] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Multijet_JetHTMET_"+period+"_CR_L17_"+region;
  histname[0][2] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Top_JetHTMET_"+period+"_CR_L17_"+region;
  histname[0][3] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/TT_powheg_pythia8_JetHTMET_"+period+"_CR_L17_"+region;
  histname[0][4] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/WToLNu_JetHTMET_"+period+"_CR_L17_"+region;
  histname[0][5] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/ZToNuNu_JetHTMET_"+period+"_CR_L17_"+region;
  histname[0][6] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Multiboson_JetHTMET_"+period+"_CR_L17_"+region;
  histname[0][7] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/GJets_JetHTMET_"+period+"_CR_L17_"+region;
  histname[0][8] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/DYToLL_JetHTMET_"+period+"_CR_L17_"+region;
  histname[1][0] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Data_JetHTMET_"+period+"_CR_LTop17_"+region;
  histname[1][1] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Multijet_JetHTMET_"+period+"_CR_LTop17_"+region;
  histname[1][2] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Top_JetHTMET_"+period+"_CR_LTop17_"+region;
  histname[1][3] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/TT_powheg_pythia8_JetHTMET_"+period+"_CR_LTop17_"+region;
  histname[1][4] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/WToLNu_JetHTMET_"+period+"_CR_LTop17_"+region;
  histname[1][5] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/ZToNuNu_JetHTMET_"+period+"_CR_LTop17_"+region;
  histname[1][6] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Multiboson_JetHTMET_"+period+"_CR_LTop17_"+region;
  histname[1][7] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/GJets_JetHTMET_"+period+"_CR_LTop17_"+region;
  histname[1][8] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/DYToLL_JetHTMET_"+period+"_CR_LTop17_"+region;

  TH2D* bkg2D[2][9];
  TH1D* bkg[2][9];
  TH1D* top[2];
  TH1D* CF[2];

  for(int i=0; i<2; i++){
    for(int j=0; j<9; j++) bkg2D[i][j] = (TH2D*)file0->Get(histname[i][j]);
  }
  int binsize = bkg2D[0][0]->GetNbinsX();
  double bin[binsize+1];
  string name;

  for(int i=1;i<=bkg2D[0][0]->GetNbinsX()+1;i++) bin[i-1] = bkg2D[0][0]->GetXaxis()->GetBinLowEdge(i);

  for(int i=0; i<2; i++) {
    for(int j=0; j<9; j++) {
      name = "bkg_"+to_string(i)+to_string(j);
      bkg[i][j] = new TH1D(name.c_str(), "", binsize, bin);
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

  TMatrix D(2,2);
  TMatrix L(2,2);
  TMatrix T(2,2);
  string name_L = "CF_L_"+(string)period+"_"+(string)region;
  string name_T = "CF_T_"+(string)period+"_"+(string)region;
  TH1D* CF_L = new TH1D(name_L.c_str(), "", binsize,bin);
  TH1D* CF_T = new TH1D(name_T.c_str(), "", binsize,bin);

  for(int i=1; i<= CF_L->GetNbinsX();i++){
    for(int j=0; j<2;j++){
      D[j][0] = bkg[j][4]->GetBinContent(i);
      D[j][1] = bkg[j][2]->GetBinContent(i);
      L[j][0] = CF[j]->GetBinContent(i);
      L[j][1] = bkg[j][2]->GetBinContent(i);
      T[j][0] = bkg[j][4]->GetBinContent(i);
      T[j][1] = CF[j]->GetBinContent(i);
    }
    CF_L->SetBinContent(i, L.Determinant()/D.Determinant());
    CF_T->SetBinContent(i, T.Determinant()/D.Determinant());
  }

  TCanvas* c1 = new TCanvas("c1","",700,700);
  
  c1->SetGridx();
  c1->SetGridy();
  CF_L->SetMaximum(2);
  CF_L->SetMinimum(0);
  CF_L->GetYaxis()->SetTitle("Correction Factor");
  CF_L->GetXaxis()->SetTitle("M_{R} #times R^{2} [GeV]");
  if(TString(obj).Contains("NJet")) CF_L->GetXaxis()->SetTitle("N_{Jets}");
  CF_L->GetYaxis()->SetTitleOffset(0.9);
  CF_L->SetMarkerStyle(21);
  CF_T->SetMarkerStyle(21);
  CF_T->SetMarkerColor(kRed);
  CF_L->Draw("P");
  CF_T->Draw("Psame");

  auto leg = new TLegend(0.45,0.77,0.9,0.9);
  if(TString(region).Contains("1Boost"))     leg->SetHeader("1 boost jet with Lep+MET final state");
  else if(TString(region).Contains("2Boost")) leg->SetHeader("#geq 2 boost jet with Lep+MET final state");
  else leg->SetHeader("boost jet with Lep+MET final state");
  leg->AddEntry(CF_L,  "W(#rightarrowl#nu)+jet CF", "p");
  leg->AddEntry(CF_T,  "Top(TT+ST) CF", "p");
  leg->Draw();


  float xmin = CF_L->GetXaxis()->GetBinLowEdge(CF_L->GetXaxis()->GetFirst());
  float xmax = CF_L->GetXaxis()->GetBinUpEdge( CF_L->GetXaxis()->GetLast());
 
  string text;
  if(TString(region).Contains("1Boost"))      text = "CMS #scale[0.7]{#font[52]{Work in progress 2016}}";
  else if(TString(region).Contains("2Boost")) text = "CMS #scale[0.7]{#font[52]{Work in progress 2017}}";
  else                                        text = "CMS #scale[0.7]{#font[52]{Work in progress 2018}}";
  TLatex* cms_lat = new TLatex(xmin, 2.05, text.c_str());
  cms_lat->SetTextSize(0.04);
  cms_lat->SetLineWidth(2);
  cms_lat->Draw();
  if(TString(region).Contains("1Boost"))      text = "#scale[0.7]{35.9 fb^{-1} (13 TeV)}";
  else if(TString(region).Contains("2Boost")) text = "#scale[0.7]{41.5 fb^{-1} (13 TeV)}";
  else                                        text = "#scale[0.7]{59.7 fb^{-1} (13 TeV)}";
  TLatex* era_lat = new TLatex(xmax,2.1, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.03);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();

  c1->SaveAs("plot/"+period+obj+"CF_"+region+".png");
  return make_tuple(CF_L, CF_T);
}

void BkgCorr_solver_NF(){
  tuple<TH1D*, TH1D*, TH1D*> CFs;
  tuple<TH1D*, TH1D*> LCFs;
  TH1D* h1[5][6][2];
  CFs = Correction("2016", "1Boost", "MRR2");
  h1[0][0][0] = get<0>(CFs); h1[1][0][0] = get<1>(CFs); h1[2][0][0] = get<2>(CFs);
  CFs = Correction("2016", "2Boost", "MRR2");
  h1[0][1][0] = get<0>(CFs); h1[1][1][0] = get<1>(CFs); h1[2][1][0] = get<2>(CFs);
  CFs = Correction("2017", "1Boost", "MRR2");
  h1[0][2][0] = get<0>(CFs); h1[1][2][0] = get<1>(CFs); h1[2][2][0] = get<2>(CFs);
  CFs = Correction("2017", "2Boost", "MRR2");
  h1[0][3][0] = get<0>(CFs); h1[1][3][0] = get<1>(CFs); h1[2][3][0] = get<2>(CFs);
  CFs = Correction("2018", "1Boost", "MRR2");
  h1[0][4][0] = get<0>(CFs); h1[1][4][0] = get<1>(CFs); h1[2][4][0] = get<2>(CFs);
  CFs = Correction("2018", "2Boost", "MRR2");
  h1[0][5][0] = get<0>(CFs); h1[1][5][0] = get<1>(CFs); h1[2][5][0] = get<2>(CFs);

  LCFs = LCorrection("2016", "1Boost", "MRR2");
  h1[3][0][0] = get<0>(LCFs); h1[4][0][0] = get<1>(LCFs);
  LCFs = LCorrection("2016", "2Boost", "MRR2");
  h1[3][1][0] = get<0>(LCFs); h1[4][1][0] = get<1>(LCFs);
  LCFs = LCorrection("2017", "1Boost", "MRR2");
  h1[3][2][0] = get<0>(LCFs); h1[4][2][0] = get<1>(LCFs);
  LCFs = LCorrection("2017", "2Boost", "MRR2");
  h1[3][3][0] = get<0>(LCFs); h1[4][3][0] = get<1>(LCFs);
  LCFs = LCorrection("2018", "1Boost", "MRR2");
  h1[3][4][0] = get<0>(LCFs); h1[4][4][0] = get<1>(LCFs);
  LCFs = LCorrection("2018", "2Boost", "MRR2");
  h1[3][5][0] = get<0>(LCFs); h1[4][5][0] = get<1>(LCFs);

  CFs = Correction("2016", "1Boost", "NJet");
  h1[0][0][1] = get<0>(CFs); h1[1][0][1] = get<1>(CFs); h1[2][0][1] = get<2>(CFs);
  CFs = Correction("2016", "2Boost", "NJet");
  h1[0][1][1] = get<0>(CFs); h1[1][1][1] = get<1>(CFs); h1[2][1][1] = get<2>(CFs);
  CFs = Correction("2017", "1Boost", "NJet");
  h1[0][2][1] = get<0>(CFs); h1[1][2][1] = get<1>(CFs); h1[2][2][1] = get<2>(CFs);
  CFs = Correction("2017", "2Boost", "NJet");
  h1[0][3][1] = get<0>(CFs); h1[1][3][1] = get<1>(CFs); h1[2][3][1] = get<2>(CFs);
  CFs = Correction("2018", "1Boost", "NJet");
  h1[0][4][1] = get<0>(CFs); h1[1][4][1] = get<1>(CFs); h1[2][4][1] = get<2>(CFs);
  CFs = Correction("2018", "2Boost", "NJet");
  h1[0][5][1] = get<0>(CFs); h1[1][5][1] = get<1>(CFs); h1[2][5][1] = get<2>(CFs);

  LCFs = LCorrection("2016", "1Boost", "NJet");
  h1[3][0][1] = get<0>(LCFs); h1[4][0][1] = get<1>(LCFs);
  LCFs = LCorrection("2016", "2Boost", "NJet");
  h1[3][1][1] = get<0>(LCFs); h1[4][1][1] = get<1>(LCFs);
  LCFs = LCorrection("2017", "1Boost", "NJet");
  h1[3][2][1] = get<0>(LCFs); h1[4][2][1] = get<1>(LCFs);
  LCFs = LCorrection("2017", "2Boost", "NJet");
  h1[3][3][1] = get<0>(LCFs); h1[4][3][1] = get<1>(LCFs);
  LCFs = LCorrection("2018", "1Boost", "NJet");
  h1[3][4][1] = get<0>(LCFs); h1[4][4][1] = get<1>(LCFs);
  LCFs = LCorrection("2018", "2Boost", "NJet");
  h1[3][5][1] = get<0>(LCFs); h1[4][5][1] = get<1>(LCFs);

  TFile* output = new TFile("CFs.root", "recreate");
  for(int i=0;i<5;i++) {
     for(int j=0;j<6;j++) h1[i][j][0]->Write();
  }
  output = new TFile("NJet_CFs.root", "recreate");
  for(int i=0;i<5;i++) {
     for(int j=0;j<6;j++) h1[i][j][1]->Write();
  }
}
