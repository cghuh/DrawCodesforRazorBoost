void DrawHist(TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4, TH1D* h5, TH1D* h6, string region) {
  delete gROOT->FindObject("c1");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();

  h1->Divide(h3);
  h2->Divide(h4);

  TCanvas* c1 = new TCanvas("c1","",700,700);

  c1->SetGridx();
  c1->SetGridy();
  if(TString(region).Contains("T")) {
    h1->SetMaximum(1e-1);
    h1->SetMinimum(1e-3);
  }
  else {
    h1->SetMaximum(2e-3);
    h1->SetMinimum(1e-4);
  }
  c1->GetFrame()->SetBorderSize(12);
  h1->SetMarkerStyle(21);
  h2->SetMarkerStyle(21);
  h5->SetMarkerStyle(21);
  h6->SetMarkerStyle(21);
  h1->SetLineColor(kRed);
  h2->SetLineColor(kBlue);
  h5->SetLineColor(kMagenta);
  h6->SetLineColor(kBlack);
  h1->SetMarkerColor(kRed);
  h2->SetMarkerColor(kBlue);
  h5->SetMarkerColor(kMagenta);
  h6->SetMarkerColor(kBlack);
  h1->Draw("EP");
  h2->Draw("EPsame");
  h5->Draw("EPsame");
  h6->Draw("EPsame");

  auto leg = new TLegend(0.48,0.72,0.9,0.9);
  if(TString(region).Contains("Q"))      leg->SetHeader("Multijet dominant CRH/CR");
  else if(TString(region).Contains("T")) leg->SetHeader("Top dominant CRH/CR");
  else                                   leg->SetHeader("Wjet dominant CRH/CR");
  leg->AddEntry(h1,  "CRH/CR, data", "p");
  leg->AddEntry(h2,  "CRH/CR, MC", "p");
  leg->AddEntry(h5,  "CRH/CR, data", "p");
  leg->AddEntry(h6,  "CRH/CR, MC", "p");
  leg->Draw();

  c1->SaveAs("higgs/"+TString(region)+".png");

}

void Correction(TString period, TString sample){
  cout << "period : " << period << endl;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();
  TString dir = "/Users/huhchanggi/temp/200827/";
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
    //delete gROOT->FindObject("remove");
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
      //else if (j==1) bkg[i]  = (TH1D*)temp[i][j]->Clone();
      //else if (j>=2) bkg[i]->Add(temp[i][j]);
      if( i==0 || i==6) {
             if (j==1) bkg[i]  = (TH1D*)temp[i][j]->Clone();
        else if (j>=4) bkg[i]->Add(temp[i][j]);
      }
      else if( i==1 || i==7) {
             if (j==2) bkg[i]  = (TH1D*)temp[i][j]->Clone();
        else if (j==3) bkg[i]->Add(temp[i][j]);
      }
    }
  }

  data_LCR[0] = (TH1D*)data[0]->Clone();
  //data_LCR[0]->Add(data[3]);
  data_LCR[1] = (TH1D*)data[1]->Clone();
  //data_LCR[1]->Add(data[4]);
  data_LCR[2] = (TH1D*)data[2]->Clone();
  //data_LCR[2]->Add(data[5]);

  data_HCR[0] = (TH1D*)data[6]->Clone();
  data_HCR[1] = (TH1D*)data[7]->Clone();
  data_HCR[2] = (TH1D*)data[8]->Clone();

  bkg_LCR[0] = (TH1D*)bkg[0]->Clone();
  //bkg_LCR[0]->Add(bkg[3]);
  bkg_LCR[1] = (TH1D*)bkg[1]->Clone();
  //bkg_LCR[1]->Add(bkg[4]);
  bkg_LCR[2] = (TH1D*)bkg[2]->Clone();
  //bkg_LCR[2]->Add(bkg[5]);

  bkg_HCR[0] = (TH1D*)bkg[6]->Clone();
  bkg_HCR[1] = (TH1D*)bkg[7]->Clone();
  bkg_HCR[2] = (TH1D*)bkg[8]->Clone();

  TH1D* h1[4][binsize+1];

  for(int i=1;i<=binsize;i++){
    delete gROOT->FindObject("h1");
    delete gROOT->FindObject("h2");
    delete gROOT->FindObject("h3");
    delete gROOT->FindObject("h4");
    h1[0][i] = new TH1D("h1", "", 1,0,3000);
    h1[0][i]->SetBinContent(1,data[6]->GetBinContent(i));
    h1[0][i]->SetBinError(1,data[6]->GetBinError(i));
    h1[1][i] = new TH1D("h2", "", 1,0,3000);
    h1[1][i]->SetBinContent(1,bkg[6]->GetBinContent(i));
    h1[1][i]->SetBinError(1,bkg[6]->GetBinError(i));
    h1[2][i] = new TH1D("h3", "", 1,0,3000);
    h1[2][i]->SetBinContent(1,data_LCR[0]->GetBinContent(i));
    h1[2][i]->SetBinError(1,data_LCR[0]->GetBinError(i));
    h1[3][i] = new TH1D("h4", "", 1,0,3000);
    h1[3][i]->SetBinContent(1,bkg_LCR[0]->GetBinContent(i));
    h1[3][i]->SetBinError(1,bkg_LCR[0]->GetBinError(i));
    if(i==1) {
      h1[0][0] = (TH1D*)h1[0][1]->Clone();
      h1[1][0] = (TH1D*)h1[1][1]->Clone();
      h1[2][0] = (TH1D*)h1[2][1]->Clone();
      h1[3][0] = (TH1D*)h1[3][1]->Clone();
      h1[0][0]->SetName("h00");
      h1[1][0]->SetName("h10");
      h1[2][0]->SetName("h20");
      h1[3][0]->SetName("h30");
    }
    else {
      h1[0][0]->Add(h1[0][i]);
      h1[1][0]->Add(h1[1][i]);
      h1[2][0]->Add(h1[2][i]);
      h1[3][0]->Add(h1[3][i]);
    }
  }
  cout << "CF_Q_data : " << h1[0][0]->GetBinContent(1) << " +-" << h1[0][0]->GetBinError(1) << endl;
  cout << "CF_Q_bkg  : " << h1[1][0]->GetBinContent(1) << " +-" << h1[1][0]->GetBinError(1) << endl;
  cout << "CF_Q_data : " << h1[2][0]->GetBinContent(1) << " +-" << h1[2][0]->GetBinError(1) << endl;
  cout << "CF_Q_bkg  : " << h1[3][0]->GetBinContent(1) << " +-" << h1[3][0]->GetBinError(1) << endl;
  h1[0][0]->Divide(h1[2][0]);
  h1[1][0]->Divide(h1[3][0]);
  cout << "CF_Q_data : " << h1[0][0]->GetBinContent(1) << " +-" << h1[0][0]->GetBinError(1) << endl;
  cout << "CF_Q_bkg  : " << h1[1][0]->GetBinContent(1) << " +-" << h1[1][0]->GetBinError(1) << endl;

  TCanvas *c2 = new TCanvas("c2", "", 700, 700);
  c2->SetGridx();
  c2->SetGridy();
  c2->GetFrame()->SetBorderSize(12);
  h1[0][0]->SetMaximum(2e-3);
  h1[0][0]->SetMinimum(1e-4);
  h1[0][0]->SetMarkerStyle(22);
  h1[1][0]->SetMarkerStyle(22);
  h1[0][0]->SetMarkerColor(kBlack);
  h1[0][0]->SetLineColor(kBlack);
  h1[1][0]->SetMarkerColor(kRed);
  h1[1][0]->SetLineColor(kRed);
  h1[0][0]->Draw("EP");
  h1[1][0]->Draw("EPsame");
  c2->SaveAs("higgs/"+period+"_fakerate_Q.png");

  name = period+"_Q";
  DrawHist(data[6], bkg[6], data_LCR[0], bkg_LCR[0], h1[0][0], h1[1][0], name);


  for(int i=1;i<=binsize;i++){
    delete gROOT->FindObject("h1");
    delete gROOT->FindObject("h2");
    delete gROOT->FindObject("h3");
    delete gROOT->FindObject("h4");
    h1[0][i] = new TH1D("h1", "", 1,0,3000);
    h1[0][i]->SetBinContent(1,data[7]->GetBinContent(i));
    h1[0][i]->SetBinError(1,data[7]->GetBinError(i));
    h1[1][i] = new TH1D("h2", "", 1,0,3000);
    h1[1][i]->SetBinContent(1,bkg[7]->GetBinContent(i));
    h1[1][i]->SetBinError(1,bkg[7]->GetBinError(i));
    h1[2][i] = new TH1D("h3", "", 1,0,3000);
    h1[2][i]->SetBinContent(1,data_LCR[1]->GetBinContent(i));
    h1[2][i]->SetBinError(1,data_LCR[1]->GetBinError(i));
    h1[3][i] = new TH1D("h4", "", 1,0,3000);
    h1[3][i]->SetBinContent(1,bkg_LCR[1]->GetBinContent(i));
    h1[3][i]->SetBinError(1,bkg_LCR[1]->GetBinError(i));
    if(i==1) {
      h1[0][0] = (TH1D*)h1[0][i]->Clone();
      h1[1][0] = (TH1D*)h1[1][i]->Clone();
      h1[2][0] = (TH1D*)h1[2][i]->Clone();
      h1[3][0] = (TH1D*)h1[3][i]->Clone();
      h1[0][0]->SetName("h00");
      h1[1][0]->SetName("h10");
      h1[2][0]->SetName("h20");
      h1[3][0]->SetName("h30");
    }
    else {
      h1[0][0]->Add(h1[0][i]);
      h1[1][0]->Add(h1[1][i]);
      h1[2][0]->Add(h1[2][i]);
      h1[3][0]->Add(h1[3][i]);
    }
  }
  cout << "CF_T_data : " << h1[0][0]->GetBinContent(1) << " +-" << h1[0][0]->GetBinError(1) << endl;
  cout << "CF_T_bkg  : " << h1[1][0]->GetBinContent(1) << " +-" << h1[1][0]->GetBinError(1) << endl;
  cout << "CF_T_data : " << h1[2][0]->GetBinContent(1) << " +-" << h1[2][0]->GetBinError(1) << endl;
  cout << "CF_T_bkg  : " << h1[3][0]->GetBinContent(1) << " +-" << h1[3][0]->GetBinError(1) << endl;
  h1[0][0]->Divide(h1[2][0]);
  h1[1][0]->Divide(h1[3][0]);
  cout << "CF_T_data : " << h1[0][0]->GetBinContent(1) << " +-" << h1[0][0]->GetBinError(1) << endl;
  cout << "CF_T_bkg  : " << h1[1][0]->GetBinContent(1) << " +-" << h1[1][0]->GetBinError(1) << endl;

  h1[0][0]->SetMaximum(1e-1);
  h1[0][0]->SetMinimum(1e-3);
  h1[0][0]->SetMarkerStyle(22);
  h1[1][0]->SetMarkerStyle(22);
  h1[0][0]->SetMarkerColor(kBlack);
  h1[0][0]->SetLineColor(kBlack);
  h1[1][0]->SetMarkerColor(kRed);
  h1[1][0]->SetLineColor(kRed);
  h1[0][0]->Draw("EP");
  h1[1][0]->Draw("EPsame");
  c2->SaveAs("higgs/"+period+"_fakerate_T.png");

  name = period+"_T";
  DrawHist(data[7], bkg[7], data_LCR[1], bkg_LCR[1], h1[0][0], h1[1][0], name);
  //name = period+"_W";
  //DrawHist(data_HCR[2], bkg_HCR[2], data_LCR[2], bkg_LCR[2], name);
  //DrawHist(data[8], bkg[8], data_LCR[2], bkg_LCR[2], name);

}

void Higgs(){
  Correction("2016", "run_2020_08_30.root");
  Correction("2017", "run_2020_08_30.root");
  Correction("2018", "run_2020_08_30.root");
}
