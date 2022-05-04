void DrawHist(TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4, TH1D* h5, TH1D* h6, TH1D* h7, TH1D* h8, string region) {
  delete gROOT->FindObject("c1");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();

  h1->Add(h2, -1);
  h3->Add(h4, -1);
  h5->Add(h6, -1);
  h7->Add(h8, -1);
  
  h1->Divide(h3);
  h5->Divide(h7);

  TCanvas* c1 = new TCanvas("c1","",700,700);

  c1->SetGridx();
  c1->SetGridy();
  h1->SetMaximum(1e-1);
  h1->SetMinimum(1e-3);
  c1->GetFrame()->SetBorderSize(12);
  h1->SetMarkerStyle(21);
  h5->SetMarkerStyle(21);
  h1->SetLineColor(kRed);
  h5->SetLineColor(kBlue);
  h1->SetMarkerColor(kRed);
  h5->SetMarkerColor(kBlue);
  h1->Draw("EP");
  h5->Draw("EPsame");

  auto leg = new TLegend(0.48,0.72,0.9,0.9);
  leg->SetHeader(TString(region)+", Higgs Correction Factors");
  leg->AddEntry(h1,  "CF_{Multijet}", "p");
  leg->AddEntry(h5,  "CF_{TT}", "p");
  leg->Draw();

  c1->SaveAs("higgs/"+TString(region)+"_fakerate.png");

}

void Correction(TString period, TString sample){
  cout << "period : " << period << endl;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();
  TString dir = "/Users/huhchanggi/temp/";
  TFile* file0 = TFile::Open(dir+sample);

  TString histname[9][9];
  file0->cd("Counts_vs_MRR2Bins/Syst_vs_MRR2Bins");
  histname[0][0] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Data_"+period+"_CR_QCD17_1Boost";
  histname[0][1] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Multijet_"+period+"_CR_QCD17_1Boost";
  histname[0][2] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Top_"+period+"_CR_QCD17_1Boost";
  histname[0][3] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/TT_powheg_pythia8_"+period+"_CR_QCD17_1Boost";
  histname[0][4] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/WToLNu_"+period+"_CR_QCD17_1Boost";
  histname[0][5] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/ZToNuNu_"+period+"_CR_QCD17_1Boost";
  histname[0][6] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Multiboson_"+period+"_CR_QCD17_1Boost";
  histname[0][7] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/GJets_"+period+"_CR_QCD17_1Boost";
  histname[0][8] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/DYToLL_"+period+"_CR_QCD17_1Boost";
  histname[1][0] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Data_"+period+"_CR_Top17_1Boost";
  histname[1][1] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Multijet_"+period+"_CR_Top17_1Boost";
  histname[1][2] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Top_"+period+"_CR_Top17_1Boost";
  histname[1][3] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/TT_powheg_pythia8_"+period+"_CR_Top17_1Boost";
  histname[1][4] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/WToLNu_"+period+"_CR_Top17_1Boost";
  histname[1][5] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/ZToNuNu_"+period+"_CR_Top17_1Boost";
  histname[1][6] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Multiboson_"+period+"_CR_Top17_1Boost";
  histname[1][7] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/GJets_"+period+"_CR_Top17_1Boost";
  histname[1][8] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/DYToLL_"+period+"_CR_Top17_1Boost";
  histname[2][0] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Data_"+period+"_CR_W17_1Boost";
  histname[2][1] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Multijet_"+period+"_CR_W17_1Boost";
  histname[2][2] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Top_"+period+"_CR_W17_1Boost";
  histname[2][3] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/TT_powheg_pythia8_"+period+"_CR_W17_1Boost";
  histname[2][4] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/WToLNu_"+period+"_CR_W17_1Boost";
  histname[2][5] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/ZToNuNu_"+period+"_CR_W17_1Boost";
  histname[2][6] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Multiboson_"+period+"_CR_W17_1Boost";
  histname[2][7] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/GJets_"+period+"_CR_W17_1Boost";
  histname[2][8] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/DYToLL_"+period+"_CR_W17_1Boost";

  histname[3][0] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Data_"+period+"_CR_QCD17_2Boost";
  histname[3][1] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Multijet_"+period+"_CR_QCD17_2Boost";
  histname[3][2] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Top_"+period+"_CR_QCD17_2Boost";
  histname[3][3] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/TT_powheg_pythia8_"+period+"_CR_QCD17_2Boost";
  histname[3][4] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/WToLNu_"+period+"_CR_QCD17_2Boost";
  histname[3][5] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/ZToNuNu_"+period+"_CR_QCD17_2Boost";
  histname[3][6] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Multiboson_"+period+"_CR_QCD17_2Boost";
  histname[3][7] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/GJets_"+period+"_CR_QCD17_2Boost";
  histname[3][8] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/DYToLL_"+period+"_CR_QCD17_2Boost";
  histname[4][0] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Data_"+period+"_CR_Top17_2Boost";
  histname[4][1] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Multijet_"+period+"_CR_Top17_2Boost";
  histname[4][2] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Top_"+period+"_CR_Top17_2Boost";
  histname[4][3] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/TT_powheg_pythia8_"+period+"_CR_Top17_2Boost";
  histname[4][4] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/WToLNu_"+period+"_CR_Top17_2Boost";
  histname[4][5] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/ZToNuNu_"+period+"_CR_Top17_2Boost";
  histname[4][6] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Multiboson_"+period+"_CR_Top17_2Boost";
  histname[4][7] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/GJets_"+period+"_CR_Top17_2Boost";
  histname[4][8] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/DYToLL_"+period+"_CR_Top17_2Boost";
  histname[5][0] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Data_"+period+"_CR_W17_2Boost";
  histname[5][1] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Multijet_"+period+"_CR_W17_2Boost";
  histname[5][2] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Top_"+period+"_CR_W17_2Boost";
  histname[5][3] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/TT_powheg_pythia8_"+period+"_CR_W17_2Boost";
  histname[5][4] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/WToLNu_"+period+"_CR_W17_2Boost";
  histname[5][5] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/ZToNuNu_"+period+"_CR_W17_2Boost";
  histname[5][6] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Multiboson_"+period+"_CR_W17_2Boost";
  histname[5][7] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/GJets_"+period+"_CR_W17_2Boost";
  histname[5][8] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/DYToLL_"+period+"_CR_W17_2Boost";

  histname[6][0] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Data_"+period+"_CR_QCD16_H";
  histname[6][1] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Multijet_"+period+"_CR_QCD16_H";
  histname[6][2] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Top_"+period+"_CR_QCD16_H";
  histname[6][3] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/TT_powheg_pythia8_"+period+"_CR_QCD16_H";
  histname[6][4] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/WToLNu_"+period+"_CR_QCD16_H";
  histname[6][5] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/ZToNuNu_"+period+"_CR_QCD16_H";
  histname[6][6] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Multiboson_"+period+"_CR_QCD16_H";
  histname[6][7] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/GJets_"+period+"_CR_QCD16_H";
  histname[6][8] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/DYToLL_"+period+"_CR_QCD16_H";
  histname[7][0] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Data_"+period+"_CR_Top16_H";
  histname[7][1] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Multijet_"+period+"_CR_Top16_H";
  histname[7][2] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Top_"+period+"_CR_Top16_H";
  histname[7][3] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/TT_powheg_pythia8_"+period+"_CR_Top16_H";
  histname[7][4] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/WToLNu_"+period+"_CR_Top16_H";
  histname[7][5] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/ZToNuNu_"+period+"_CR_Top16_H";
  histname[7][6] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Multiboson_"+period+"_CR_Top16_H";
  histname[7][7] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/GJets_"+period+"_CR_Top16_H";
  histname[7][8] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/DYToLL_"+period+"_CR_Top16_H";
  histname[8][0] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Data_"+period+"_CR_W16_H";
  histname[8][1] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Multijet_"+period+"_CR_W16_H";
  histname[8][2] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Top_"+period+"_CR_W16_H";
  histname[8][3] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/TT_powheg_pythia8_"+period+"_CR_W16_H";
  histname[8][4] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/WToLNu_"+period+"_CR_W16_H";
  histname[8][5] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/ZToNuNu_"+period+"_CR_W16_H";
  histname[8][6] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/Multiboson_"+period+"_CR_W16_H";
  histname[8][7] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/GJets_"+period+"_CR_W16_H";
  histname[8][8] = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/DYToLL_"+period+"_CR_W16_H";

  TH2D* bkg2D[9][9];
  TH1D* temp[9][9];
  TH1D* bkg[9];
  TH1D* data[9];

  TH1D* data_LCR[3];
  TH1D* bkg_LCR[3];
  TH1D* data_HCR[3];
  TH1D* bkg_HCR[3];

  for(int i=0; i<9; i++){
    for(int j=0; j<9; j++) bkg2D[i][j] = (TH2D*)file0->Get(histname[i][j]);
  }
  int binsize = bkg2D[0][0]->GetNbinsX();
  double bin[binsize+1];
  string name;

  for(int i=1;i<=binsize+1;i++) bin[i-1] = bkg2D[0][0]->GetXaxis()->GetBinLowEdge(i);

  for(int i=0; i<9; i++) {
    delete gROOT->FindObject("remove");
    bkg[i] = new TH1D("remove", "", 1,0,1);
    for(int j=0; j<9; j++) {
      name = "bkg_"+to_string(i)+to_string(j);
      temp[i][j] = new TH1D(name.c_str(), "", binsize, bin);
      for(int k=1;k<=bkg2D[0][0]->GetNbinsX();k++) {
        if(bkg2D[i][j] == NULL) {
          temp[i][j]->SetBinContent(k, 0);
          temp[i][j]->SetBinError(k, 0);
        }
        else {
          temp[i][j]->SetBinContent(k, bkg2D[i][j]->GetBinContent(k,1));
          temp[i][j]->SetBinError(k, bkg2D[i][j]->GetBinError(k,1));
        }
      }
      if      (j==0) data[i] = (TH1D*)temp[i][j]->Clone();
      if( i%3==0) {
      //     if (j==1) bkg[i]  = (TH1D*)temp[i][j]->Clone();
      //else if (j>=4) bkg[i]->Add(temp[i][j]);
           if (j==2) bkg[i]  = (TH1D*)temp[i][j]->Clone();
      else if (j==3) bkg[i]->Add(temp[i][j]);
      }
      else if( i%3==1) {
      //     if (j==2) bkg[i]  = (TH1D*)temp[i][j]->Clone();
      //else if (j==3) bkg[i]->Add(temp[i][j]);
           if (j==1) bkg[i]  = (TH1D*)temp[i][j]->Clone();
      else if (j>=4) bkg[i]->Add(temp[i][j]);
      }
    }
  }

  data_LCR[0] = (TH1D*)data[0]->Clone();
  data_LCR[0]->Add(data[3]);
  data_LCR[1] = (TH1D*)data[1]->Clone();
  data_LCR[1]->Add(data[4]);
  //data_LCR[2] = (TH1D*)data[2]->Clone();
  //data_LCR[2]->Add(data[5]);

  data_HCR[0] = (TH1D*)data[6]->Clone();
  data_HCR[1] = (TH1D*)data[7]->Clone();
  //data_HCR[2] = (TH1D*)data[8]->Clone();

  bkg_LCR[0] = (TH1D*)bkg[0]->Clone();
  bkg_LCR[0]->Add(bkg[3]);
  bkg_LCR[1] = (TH1D*)bkg[1]->Clone();
  bkg_LCR[1]->Add(bkg[4]);
  //bkg_LCR[2] = (TH1D*)bkg[2]->Clone();
  //bkg_LCR[2]->Add(bkg[5]);

  bkg_HCR[0] = (TH1D*)bkg[6]->Clone();
  bkg_HCR[1] = (TH1D*)bkg[7]->Clone();
  //bkg_HCR[2] = (TH1D*)bkg[8]->Clone();

  name = period;
  DrawHist(data_HCR[0], bkg_HCR[0], data_LCR[0], bkg_LCR[0], data_HCR[1], bkg_HCR[1], data_LCR[1], bkg_LCR[1], name);

}

void Higgs(){
  Correction("2016", "210604/run_2021_06_01.root");
  Correction("2017", "210604/run_2021_06_01.root");
  Correction("2018", "210604/run_2021_06_01.root");
}
