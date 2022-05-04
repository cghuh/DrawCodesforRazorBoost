tuple<TH1D*, TH1D*, TH1D*, TH1D*, TH1D*, TH1D*> CalcErrors(TH1D* d1[3], TH1D* b1[3][9], string sample) {
  cout << sample << endl;
  TRandom3 r(0);
  gRandom->SetSeed(time(0));
  TFile* output;
  if(sample=="MRR220171Boost") output = new TFile("CF_Error.root", "recreate");
  else output = new TFile("CF_Error.root", "update");
  TH1D* data[3];
  TH1D* QB[3];
  TH1D* TB[3];
  TH1D* WB[3];
  for(int i=0; i<3; i++) {
    data[i]  = (TH1D*)d1[i]->Clone();
    QB[i]    = (TH1D*)b1[i][1]->Clone();
    TB[i]    = (TH1D*)b1[i][2]->Clone();
    WB[i]    = (TH1D*)b1[i][4]->Clone();
  }
  TMatrix D(3,3);
  TMatrix Q(3,3);
  TMatrix T(3,3);
  TMatrix W(3,3);
  string name;
  TH1D* QError[data[0]->GetNbinsX()];
  TH1D* TError[data[0]->GetNbinsX()];
  TH1D* WError[data[0]->GetNbinsX()];
  double nom_Q, nom_T, nom_W;
  vector<double> err_Q, err_T, err_W;
  TF1* gaus_D;
  TF1* gaus_Q;
  TF1* gaus_T;
  TF1* gaus_W;
  double ran_Q, ran_T, ran_W, ran_D, temp;
  double bin_up, bin_dw;
  int num=1e6;
  TH1D* Err_Up_Q = (TH1D*)d1[0]->Clone();
  TH1D* Err_Up_T = (TH1D*)d1[0]->Clone();
  TH1D* Err_Up_W = (TH1D*)d1[0]->Clone();
  TH1D* Err_Dw_Q = (TH1D*)d1[0]->Clone();
  TH1D* Err_Dw_T = (TH1D*)d1[0]->Clone();
  TH1D* Err_Dw_W = (TH1D*)d1[0]->Clone();
  TCanvas* c1, *c2, *c3;
  TLine* line;

  for(int i=1; i<= data[0]->GetNbinsX();i++) {
    name = "QError_"+sample+to_string(i);
    QError[i-1] = new TH1D(name.c_str(), "", 8000, -2., 2.);
    c1 = new TCanvas(name.c_str(), "", 700, 700);
    name = "TError_"+sample+to_string(i);
    TError[i-1] = new TH1D(name.c_str(), "", 8000, -2., 2.);
    c2 = new TCanvas(name.c_str(), "", 700, 700);
    name = "WError_"+sample+to_string(i);
    WError[i-1] = new TH1D(name.c_str(), "", 8000, -2., 2.);
    c3 = new TCanvas(name.c_str(), "", 700, 700);
    if(i > 1 && TString(sample).Contains("MRR2")) {
      QError[i-1]->Rebin(10);
      TError[i-1]->Rebin(10);
      WError[i-1]->Rebin(10);
    }
    for(int m=0; m<3; m++) { // region
      D[m][0] = QB[m]->GetBinContent(i);
      D[m][1] = TB[m]->GetBinContent(i);
      D[m][2] = WB[m]->GetBinContent(i);
      Q[m][0] = data[m]->GetBinContent(i);
      Q[m][1] = TB[m]->GetBinContent(i);
      Q[m][2] = WB[m]->GetBinContent(i);
      T[m][0] = QB[m]->GetBinContent(i);
      T[m][1] = data[m]->GetBinContent(i);
      T[m][2] = WB[m]->GetBinContent(i);
      W[m][0] = QB[m]->GetBinContent(i);
      W[m][1] = TB[m]->GetBinContent(i);
      W[m][2] = data[m]->GetBinContent(i);
    }
    nom_Q = Q.Determinant()/D.Determinant();
    nom_T = T.Determinant()/D.Determinant();
    nom_W = W.Determinant()/D.Determinant();
    for(int j=0; j<num; j++) {
      if(j%200000 == 0) cout << j << ", " << (double)j/1e6*100 << "[%] finsihed" << endl;
      for(int k=0; k<3; k++) { // region
        ran_D = gRandom->Gaus(data[k]->GetBinContent(i), data[k]->GetBinError(i));
        ran_Q = gRandom->Gaus(QB[k]->GetBinContent(i), QB[k]->GetBinError(i));
        ran_T = gRandom->Gaus(TB[k]->GetBinContent(i), TB[k]->GetBinError(i));
        ran_W = gRandom->Gaus(WB[k]->GetBinContent(i), WB[k]->GetBinError(i));
        D[k][0] = ran_Q;
        D[k][1] = ran_T;
        D[k][2] = ran_W;
        Q[k][0] = ran_D;
        Q[k][1] = ran_T;
        Q[k][2] = ran_W;
        T[k][0] = ran_Q;
        T[k][1] = ran_D;
        T[k][2] = ran_W;
        W[k][0] = ran_Q;
        W[k][1] = ran_T;
        W[k][2] = ran_D;
      }
      QError[i-1]->Fill(nom_Q-(Q.Determinant()/D.Determinant()));
      TError[i-1]->Fill(nom_T-(T.Determinant()/D.Determinant()));
      WError[i-1]->Fill(nom_W-(W.Determinant()/D.Determinant()));
      err_Q.push_back(nom_Q-Q.Determinant()/D.Determinant());
      err_T.push_back(nom_T-T.Determinant()/D.Determinant());
      err_W.push_back(nom_W-W.Determinant()/D.Determinant());
    }
    //QError[i-1]->Write();
    //TError[i-1]->Write();
    //WError[i-1]->Write();

    temp = 0;
    bin_up = 0; bin_dw = 0;
    for(int bin=QError[i-1]->GetMaximumBin();bin<=QError[i-1]->GetNbinsX();bin++) {
      temp += QError[i-1]->GetBinContent(bin);
      temp += QError[i-1]->GetBinContent(2*QError[i-1]->GetMaximumBin()-bin-1);
      if(temp/num > 0.68) continue;
      else {
        bin_up = QError[i-1]->GetBinCenter(bin); 
        bin_dw = QError[i-1]->GetBinCenter(2.*QError[i-1]->GetMaximumBin()-bin-1);
      }
    }
    cout << "Q : " << bin_up << ", " << bin_dw << endl;
    Err_Up_Q->SetBinContent(i, bin_up);
    Err_Up_Q->SetBinError(i, 0);
    Err_Dw_Q->SetBinContent(i, bin_dw);
    Err_Dw_Q->SetBinError(i, 0);

    if(TString(sample).Contains("MRR2")) {
      QError[i-1]->GetXaxis()->SetRangeUser(-0.5,0.5);
      TError[i-1]->GetXaxis()->SetRangeUser(-0.5,0.5);
      WError[i-1]->GetXaxis()->SetRangeUser(-0.5,0.5);
      if(i==1) {
        QError[i-1]->GetXaxis()->SetRangeUser(-0.1,0.1);
        TError[i-1]->GetXaxis()->SetRangeUser(-0.1,0.1);
        WError[i-1]->GetXaxis()->SetRangeUser(-0.1,0.1);
      }
      else if (i > 2) {
        QError[i-1]->GetXaxis()->SetRangeUser(-1.0,1.0);
        TError[i-1]->GetXaxis()->SetRangeUser(-1.0,1.0);
      }
    } else {
      QError[i-1]->GetXaxis()->SetRangeUser(-1.0,1.0);
      TError[i-1]->GetXaxis()->SetRangeUser(-1.0,1.0);
      WError[i-1]->GetXaxis()->SetRangeUser(-1.0,1.0);
    }

    c1->cd();
    QError[i-1]->Draw();
    line = new TLine(bin_up,0,bin_up,1.05*QError[i-1]->GetMaximum());
    line->SetLineColor(kRed);
    line->Draw("same");
    line = new TLine(bin_dw,0,bin_dw,1.05*QError[i-1]->GetMaximum());
    line->SetLineColor(kRed);
    line->Draw("same");
    c1->Write();

    temp = 0;
    bin_up = 0; bin_dw = 0;
    for(int bin=TError[i-1]->GetMaximumBin();bin<=TError[i-1]->GetNbinsX();bin++) {
      temp += TError[i-1]->GetBinContent(bin);
      temp += TError[i-1]->GetBinContent(2*TError[i-1]->GetMaximumBin()-bin-1);
      if(temp/num > 0.68) continue;
      else {
        bin_up = TError[i-1]->GetBinCenter(bin); 
        bin_dw = TError[i-1]->GetBinCenter(2.*TError[i-1]->GetMaximumBin()-bin-1);
      }
    }
    cout << "T : " << bin_up << ", " << bin_dw << endl;
    Err_Up_T->SetBinContent(i, bin_up);
    Err_Up_T->SetBinError(i, 0);
    Err_Dw_T->SetBinContent(i, bin_dw);
    Err_Dw_T->SetBinError(i, 0);

    c2->cd();
    TError[i-1]->Draw();
    line = new TLine(bin_up,0,bin_up,1.05*TError[i-1]->GetMaximum());
    line->SetLineColor(kRed);
    line->Draw("same");
    line = new TLine(bin_dw,0,bin_dw,1.05*TError[i-1]->GetMaximum());
    line->SetLineColor(kRed);
    line->Draw("same");
    c2->Write();

    temp = 0;
    bin_up = 0; bin_dw = 0;
    for(int bin=WError[i-1]->GetMaximumBin();bin<=WError[i-1]->GetNbinsX();bin++) {
      temp += WError[i-1]->GetBinContent(bin);
      temp += WError[i-1]->GetBinContent(2*WError[i-1]->GetMaximumBin()-bin-1);
      if(temp/num > 0.68) continue;
      else {
        bin_up = WError[i-1]->GetBinCenter(bin); 
        bin_dw = WError[i-1]->GetBinCenter(2.*WError[i-1]->GetMaximumBin()-bin-1);
      }
    }
    cout << "W : " << bin_up << ", " << bin_dw << endl;
    Err_Up_W->SetBinContent(i, bin_up);
    Err_Up_W->SetBinError(i, 0);
    Err_Dw_W->SetBinContent(i, bin_dw);
    Err_Dw_W->SetBinError(i, 0);

    c3->cd();
    WError[i-1]->Draw();
    line = new TLine(bin_up,0,bin_up,1.05*WError[i-1]->GetMaximum());
    line->SetLineColor(kRed);
    line->Draw("same");
    line = new TLine(bin_dw,0,bin_dw,1.05*WError[i-1]->GetMaximum());
    line->SetLineColor(kRed);
    line->Draw("same");
    c3->Write();
  }
  output->Close();
  return make_tuple(Err_Up_Q, Err_Up_T, Err_Up_W, Err_Dw_Q, Err_Dw_T, Err_Dw_W);
}

tuple<TH1D*, TH1D*, TH1D*, TH1D*, TH1D*, TH1D*> Errors(TH1D* d1[3], TH1D* d2[3], TH1D* d3[3], TH1D* b1[3][9], TH1D* b2[3][9], TH1D* b3[3][9]) {
  TFile* output = new TFile("CF_Error.root", "recreate");
  TH1D* data[3][3];
  TH1D* QB[3][3];
  TH1D* TB[3][3];
  TH1D* WB[3][3];
  for(int i=0; i<3; i++) {
    data[i][0] = (TH1D*)d1[i]->Clone();
    data[i][1] = (TH1D*)d2[i]->Clone();
    data[i][2] = (TH1D*)d3[i]->Clone();
    QB[i][0]    = (TH1D*)b1[i][1]->Clone();
    QB[i][1]    = (TH1D*)b2[i][1]->Clone();
    QB[i][2]    = (TH1D*)b3[i][1]->Clone();
    TB[i][0]    = (TH1D*)b1[i][2]->Clone();
    TB[i][1]    = (TH1D*)b2[i][2]->Clone();
    TB[i][2]    = (TH1D*)b3[i][2]->Clone();
    WB[i][0]    = (TH1D*)b1[i][4]->Clone();
    WB[i][1]    = (TH1D*)b2[i][4]->Clone();
    WB[i][2]    = (TH1D*)b3[i][4]->Clone();
  }
  TMatrix D(3,3);
  TMatrix Q(3,3);
  TMatrix T(3,3);
  TMatrix W(3,3);
  string name;
  TH1D* QError[data[0][0]->GetNbinsX()];
  TH1D* TError[data[0][0]->GetNbinsX()];
  TH1D* WError[data[0][0]->GetNbinsX()];
  double nom_Q, nom_T, nom_W;
  vector<double> err_Q, err_T, err_W;
  TH1D* Err_Up_Q = (TH1D*)d1[0]->Clone();
  TH1D* Err_Up_T = (TH1D*)d1[0]->Clone();
  TH1D* Err_Up_W = (TH1D*)d1[0]->Clone();
  TH1D* Err_Dw_Q = (TH1D*)d1[0]->Clone();
  TH1D* Err_Dw_T = (TH1D*)d1[0]->Clone();
  TH1D* Err_Dw_W = (TH1D*)d1[0]->Clone();

  for(int i=1; i<= data[0][0]->GetNbinsX();i++) {
    name = "QError_"+to_string(i);
    QError[i-1] = new TH1D(name.c_str(), "", 50, -0.5, 0.5);
    name = "TError_"+to_string(i);
    TError[i-1] = new TH1D(name.c_str(), "", 50, -0.5, 0.5);
    name = "WError_"+to_string(i);
    WError[i-1] = new TH1D(name.c_str(), "", 50, -0.5, 0.5);
    for(int j=0; j<3; j++) { // data variation
      for(int k=0; k<3; k++) { // Multijet variation
        for(int l=0; l<3; l++) { // TT+ST variation
          for(int n=0; n<3; n++) { // WJet variation
            for(int m=0; m<3; m++) { // region
              D[m][0] = QB[m][k]->GetBinContent(i);
              D[m][1] = TB[m][l]->GetBinContent(i);
              D[m][2] = WB[m][n]->GetBinContent(i);
              Q[m][0] = data[m][j]->GetBinContent(i);
              Q[m][1] = TB[m][l]->GetBinContent(i);
              Q[m][2] = WB[m][n]->GetBinContent(i);
              T[m][0] = QB[m][k]->GetBinContent(i);
              T[m][1] = data[m][j]->GetBinContent(i);
              T[m][2] = WB[m][n]->GetBinContent(i);
              W[m][0] = QB[m][k]->GetBinContent(i);
              W[m][1] = TB[m][l]->GetBinContent(i);
              W[m][2] = data[m][j]->GetBinContent(i);
            }
            if(j==0 && k==0 && l == 0 && n == 0) {
              nom_Q = Q.Determinant()/D.Determinant();
              nom_T = T.Determinant()/D.Determinant();
              nom_W = W.Determinant()/D.Determinant();
            }
            QError[i-1]->Fill(nom_Q-Q.Determinant()/D.Determinant());
            TError[i-1]->Fill(nom_T-T.Determinant()/D.Determinant());
            WError[i-1]->Fill(nom_W-W.Determinant()/D.Determinant());
            err_Q.push_back(nom_Q-Q.Determinant()/D.Determinant());
            err_T.push_back(nom_T-T.Determinant()/D.Determinant());
            err_W.push_back(nom_W-W.Determinant()/D.Determinant());
          }
        }
      }
    }
    QError[i-1]->Write();
    TError[i-1]->Write();
    WError[i-1]->Write();
    Err_Up_Q->SetBinContent(i, *max_element(err_Q.begin(), err_Q.end()));
    Err_Up_Q->SetBinError(i, 0);
    Err_Up_T->SetBinContent(i, *max_element(err_T.begin(), err_T.end()));
    Err_Up_T->SetBinError(i, 0);
    Err_Up_W->SetBinContent(i, *max_element(err_W.begin(), err_W.end()));
    Err_Up_W->SetBinError(i, 0);
    Err_Dw_Q->SetBinContent(i, *min_element(err_Q.begin(), err_Q.end()));
    Err_Dw_Q->SetBinError(i, 0);
    Err_Dw_T->SetBinContent(i, *min_element(err_T.begin(), err_T.end()));
    Err_Dw_T->SetBinError(i, 0);
    Err_Dw_W->SetBinContent(i, *min_element(err_W.begin(), err_W.end()));
    Err_Dw_W->SetBinError(i, 0);
  }
  Err_Up_Q->Write();
  Err_Up_T->Write();
  Err_Up_W->Write();
  Err_Dw_Q->Write();
  Err_Dw_T->Write();
  Err_Dw_W->Write();
  output->Close();
  return make_tuple(Err_Up_Q, Err_Up_T, Err_Up_W, Err_Dw_Q, Err_Dw_T, Err_Dw_W);
}

tuple<TGraphAsymmErrors*, TGraphAsymmErrors*, TGraphAsymmErrors*> Correction(TString period, TString region, TString obj, TString sample){
  cout << "period : " << period << ", region : " << region << ", obj : " << obj << endl;
  delete gROOT->FindObject("c1");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();
  TString dir = "/Users/huhchanggi/temp/200702/";
  TFile* file0 = TFile::Open(dir+sample);

  TString histname[3][9];
  file0->cd("Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins");
  histname[0][0] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Data_"+period+"_CR_QCD17_"+region;
  histname[0][1] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Multijet_"+period+"_CR_QCD17_"+region;
  histname[0][2] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Top_"+period+"_CR_QCD17_"+region;
  histname[0][3] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/TT_powheg_pythia8_"+period+"_CR_QCD17_"+region;
  histname[0][4] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/WToLNu_"+period+"_CR_QCD17_"+region;
  histname[0][5] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/ZToNuNu_"+period+"_CR_QCD17_"+region;
  histname[0][6] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Multiboson_"+period+"_CR_QCD17_"+region;
  histname[0][7] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/GJets_"+period+"_CR_QCD17_"+region;
  histname[0][8] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/DYToLL_"+period+"_CR_QCD17_"+region;
  histname[1][0] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Data_"+period+"_CR_Top17_"+region;
  histname[1][1] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Multijet_"+period+"_CR_Top17_"+region;
  histname[1][2] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Top_"+period+"_CR_Top17_"+region;
  histname[1][3] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/TT_powheg_pythia8_"+period+"_CR_Top17_"+region;
  histname[1][4] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/WToLNu_"+period+"_CR_Top17_"+region;
  histname[1][5] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/ZToNuNu_"+period+"_CR_Top17_"+region;
  histname[1][6] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Multiboson_"+period+"_CR_Top17_"+region;
  histname[1][7] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/GJets_"+period+"_CR_Top17_"+region;
  histname[1][8] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/DYToLL_"+period+"_CR_Top17_"+region;
  histname[2][0] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Data_"+period+"_CR_W17_"+region;
  histname[2][1] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Multijet_"+period+"_CR_W17_"+region;
  histname[2][2] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Top_"+period+"_CR_W17_"+region;
  histname[2][3] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/TT_powheg_pythia8_"+period+"_CR_W17_"+region;
  histname[2][4] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/WToLNu_"+period+"_CR_W17_"+region;
  histname[2][5] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/ZToNuNu_"+period+"_CR_W17_"+region;
  histname[2][6] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Multiboson_"+period+"_CR_W17_"+region;
  histname[2][7] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/GJets_"+period+"_CR_W17_"+region;
  histname[2][8] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/DYToLL_"+period+"_CR_W17_"+region;

  TH2D* bkg2D[3][9];
  TH1D* bkg[3][9];
  TH1D* bkg_up[3][9];
  TH1D* bkg_dw[3][9];
  TH1D* top[3];
  TH1D* CF[3];
  TH1D* CF_up[3];
  TH1D* CF_dw[3];

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
      name = "bkg_up_"+to_string(i)+to_string(j);
      bkg_up[i][j] = new TH1D(name.c_str(), "", binsize, bin);
      name = "bkg_dw_"+to_string(i)+to_string(j);
      bkg_dw[i][j] = new TH1D(name.c_str(), "", binsize, bin);
      for(int k=1;k<=bkg2D[0][0]->GetNbinsX();k++) {
        if(bkg2D[i][j] == NULL) {
          bkg[i][j]->SetBinContent(k, 0);
          bkg[i][j]->SetBinError(k, 0);
        }
        else {
          bkg[i][j]->SetBinContent(k, bkg2D[i][j]->GetBinContent(k,1));
          bkg[i][j]->SetBinError(k, bkg2D[i][j]->GetBinError(k,1));
          bkg_up[i][j]->SetBinContent(k, bkg2D[i][j]->GetBinContent(k,1)+bkg2D[i][j]->GetBinError(k,1));
          bkg_up[i][j]->SetBinError(k, bkg2D[i][j]->GetBinError(k,1));
          bkg_dw[i][j]->SetBinContent(k, bkg2D[i][j]->GetBinContent(k,1)-bkg2D[i][j]->GetBinError(k,1));
          bkg_dw[i][j]->SetBinError(k, bkg2D[i][j]->GetBinError(k,1));
        }
      }
    }
    bkg[i][2]->Add(bkg[i][3]);
    bkg_up[i][2]->Add(bkg_up[i][3]);
    bkg_dw[i][2]->Add(bkg_dw[i][3]);
    CF[i] = (TH1D*)bkg[i][0]->Clone();
    CF_up[i] = (TH1D*)bkg_up[i][0]->Clone();
    CF_dw[i] = (TH1D*)bkg_dw[i][0]->Clone();
    for(int j=5; j<9; j++){
      CF[i]->Add(bkg[i][j], -1);
      CF_up[i]->Add(bkg_up[i][j], -1);
      CF_dw[i]->Add(bkg_dw[i][j], -1);
    }
  }

  TMatrix D(3,3);
  TMatrix Q(3,3);
  TMatrix T(3,3);
  TMatrix W(3,3);
  TMatrix D_up(3,3);
  TMatrix Q_up(3,3);
  TMatrix T_up(3,3);
  TMatrix W_up(3,3);
  TMatrix D_dw(3,3);
  TMatrix Q_dw(3,3);
  TMatrix T_dw(3,3);
  TMatrix W_dw(3,3);

  double cf_Q[binsize];
  double cf_Q_up[binsize];
  double cf_Q_dw[binsize];
  double cf_T[binsize];
  double cf_T_up[binsize];
  double cf_T_dw[binsize];
  double cf_W[binsize];
  double cf_W_up[binsize];
  double cf_W_dw[binsize];
  double binerr_x[binsize];
  double binx[binsize];

  double qbinerr_x[binsize-1];
  double qbinx[binsize-1];

  tuple<TH1D*, TH1D*, TH1D*, TH1D*, TH1D*, TH1D*> Err_CFs;
  TH1D* temp[6];
  int flag = 0;

  for(int i=1; i<= binsize;i++){
    for(int j=0; j<3;j++) {
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
      D_up[j][0] = bkg_up[j][1]->GetBinContent(i);
      D_up[j][1] = bkg_up[j][2]->GetBinContent(i);
      D_up[j][2] = bkg_up[j][4]->GetBinContent(i);
      Q_up[j][0] = CF_up[j]->GetBinContent(i);
      Q_up[j][1] = bkg_up[j][2]->GetBinContent(i);
      Q_up[j][2] = bkg_up[j][4]->GetBinContent(i);
      T_up[j][0] = bkg_up[j][1]->GetBinContent(i);
      T_up[j][1] = CF_up[j]->GetBinContent(i);
      T_up[j][2] = bkg_up[j][4]->GetBinContent(i);
      W_up[j][0] = bkg_up[j][1]->GetBinContent(i);
      W_up[j][1] = bkg_up[j][2]->GetBinContent(i);
      W_up[j][2] = CF_up[j]->GetBinContent(i);
      D_dw[j][0] = bkg_dw[j][1]->GetBinContent(i);
      D_dw[j][1] = bkg_dw[j][2]->GetBinContent(i);
      D_dw[j][2] = bkg_dw[j][4]->GetBinContent(i);
      Q_dw[j][0] = CF_dw[j]->GetBinContent(i);
      Q_dw[j][1] = bkg_dw[j][2]->GetBinContent(i);
      Q_dw[j][2] = bkg_dw[j][4]->GetBinContent(i);
      T_dw[j][0] = bkg_dw[j][1]->GetBinContent(i);
      T_dw[j][1] = CF_dw[j]->GetBinContent(i);
      T_dw[j][2] = bkg_dw[j][4]->GetBinContent(i);
      W_dw[j][0] = bkg_dw[j][1]->GetBinContent(i);
      W_dw[j][1] = bkg_dw[j][2]->GetBinContent(i);
      W_dw[j][2] = CF_dw[j]->GetBinContent(i);
      binerr_x[i-1] = CF[j]->GetXaxis()->GetBinWidth(i)/2;
      binx[i-1] = CF[j]->GetXaxis()->GetBinCenter(i);
    }

    cf_Q[i-1] = Q.Determinant()/D.Determinant();
    cf_T[i-1] = T.Determinant()/D.Determinant();
    cf_W[i-1] = W.Determinant()/D.Determinant();

    if(cf_Q[i-1] < 0) {
      cout << "CF is smaller than 0" << endl;
      for(int j=0; j<3;j++) {
        D[j][0] = bkg[j][1]->GetBinContent(i) + bkg[j][1]->GetBinContent(i-1);
        D[j][1] = bkg[j][2]->GetBinContent(i) + bkg[j][2]->GetBinContent(i-1);
        D[j][2] = bkg[j][4]->GetBinContent(i) + bkg[j][4]->GetBinContent(i-1);
        Q[j][0] = CF[j]->GetBinContent(i) + CF[j]->GetBinContent(i-1);
        Q[j][1] = bkg[j][2]->GetBinContent(i) + bkg[j][2]->GetBinContent(i-1);
        Q[j][2] = bkg[j][4]->GetBinContent(i) + bkg[j][4]->GetBinContent(i-1);
        T[j][0] = bkg[j][1]->GetBinContent(i) + bkg[j][1]->GetBinContent(i-1);
        T[j][1] = CF[j]->GetBinContent(i) + CF[j]->GetBinContent(i-1);
        T[j][2] = bkg[j][4]->GetBinContent(i) + bkg[j][4]->GetBinContent(i-1);
        W[j][0] = bkg[j][1]->GetBinContent(i) + bkg[j][1]->GetBinContent(i-1);
        W[j][1] = bkg[j][2]->GetBinContent(i) + bkg[j][2]->GetBinContent(i-1);
        W[j][2] = CF[j]->GetBinContent(i) + CF[j]->GetBinContent(i-1);
      }
      qbinerr_x[0] = binerr_x[0];
      qbinerr_x[1] = binerr_x[1];
      qbinx[0] = binx[0];
      qbinx[1] = binx[1];
      qbinerr_x[3] = (CF[0]->GetXaxis()->GetBinWidth(i) + CF[0]->GetXaxis()->GetBinWidth(i-1))/2.;
      qbinx[3] = (CF[0]->GetXaxis()->GetBinCenter(i-1) - CF[0]->GetXaxis()->GetBinWidth(i-1)/2 + CF[0]->GetXaxis()->GetBinCenter(i) + CF[0]->GetXaxis()->GetBinWidth(i)/2 )/2.;
      cf_Q[i-2] = Q.Determinant()/D.Determinant();
      flag = 1;
    }
  }

  Err_CFs = CalcErrors(CF, bkg, (string)obj+(string)period+(string)region);

  temp[0] = get<0>(Err_CFs);
  temp[1] = get<1>(Err_CFs);
  temp[2] = get<2>(Err_CFs);
  temp[3] = get<3>(Err_CFs);
  temp[4] = get<4>(Err_CFs);
  temp[5] = get<5>(Err_CFs);

  for(int i=1; i<= binsize;i++){
    cf_Q_up[i-1] = abs(temp[0]->GetBinContent(i));
    cf_T_up[i-1] = abs(temp[1]->GetBinContent(i));
    cf_W_up[i-1] = abs(temp[2]->GetBinContent(i));
    cf_Q_dw[i-1] = abs(temp[3]->GetBinContent(i));
    cf_T_dw[i-1] = abs(temp[4]->GetBinContent(i));
    cf_W_dw[i-1] = abs(temp[5]->GetBinContent(i));
    cout << cf_Q_up[i-1] << "," << cf_Q_dw[i-1] << ", " << cf_T_up[i-1] << "," << cf_T_dw[i-1] << ", " << cf_W_up[i-1] << "," << cf_W_dw[i-1] << endl;
  }

  auto CF_Q = new TGraphAsymmErrors(binsize, binx, cf_Q, binerr_x, binerr_x, cf_Q_up, cf_Q_dw);
  if(flag == 1) CF_Q = new TGraphAsymmErrors(binsize-1, qbinx, cf_Q, qbinerr_x, qbinerr_x, cf_Q_up, cf_Q_dw);
  auto CF_T = new TGraphAsymmErrors(binsize, binx, cf_T, binerr_x, binerr_x, cf_T_up, cf_T_dw);
  auto CF_W = new TGraphAsymmErrors(binsize, binx, cf_W, binerr_x, binerr_x, cf_W_up, cf_W_dw);

  string name_Q = "CF_Q_"+(string)period+"_"+(string)region;
  string name_T = "CF_T_"+(string)period+"_"+(string)region;
  string name_W = "CF_W_"+(string)period+"_"+(string)region;
  CF_Q->SetName(name_Q.c_str());
  CF_T->SetName(name_T.c_str());
  CF_W->SetName(name_W.c_str());


  TCanvas* c1 = new TCanvas("c1","",700,700);

  c1->SetGridx();
  c1->SetGridy();
  //TH1F *hr = c1->DrawFrame(0,3000,0.,2.0);
  //hr->SetMaximum(2);
  //hr->SetMinimum(0);
  //hr->GetYaxis()->SetTitle("Correction Factor");
  //hr->GetXaxis()->SetTitle("M_{R} #times R^{2} [GeV]");
  //if(TString(obj).Contains("NJet")) hr->GetXaxis()->SetTitle("N_{Jets}");
  //hr->GetYaxis()->SetTitleOffset(0.9);
  c1->GetFrame()->SetBorderSize(12);
  //CF_Q->SetMaximum(2);
  //CF_Q->SetMinimum(0);
  //CF_Q->GetYaxis()->SetTitle("Correction Factor");
  //CF_Q->GetXaxis()->SetTitle("M_{R} #times R^{2} [GeV]");
  //if(TString(obj).Contains("NJet")) CF_Q->GetXaxis()->SetTitle("N_{Jets}");
  //CF_Q->GetYaxis()->SetTitleOffset(0.9);
  CF_Q->SetMarkerStyle(21);
  CF_T->SetMarkerStyle(21);
  CF_W->SetMarkerStyle(21);
  CF_T->SetMarkerColor(kRed);
  CF_W->SetMarkerColor(kBlue);
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(CF_Q);
  mg->Add(CF_T);
  mg->Add(CF_W);
  mg->Draw("AP");

  auto leg = new TLegend(0.45,0.77,0.9,0.9);
  if(TString(region).Contains("1Boost"))     leg->SetHeader("1 boost jet final state");
  else if(TString(region).Contains("2Boost")) leg->SetHeader("#geq 2 boost jet final state");
  else leg->SetHeader("boost jet final state");
  leg->AddEntry(CF_Q,  "Multijet CF", "p");
  leg->AddEntry(CF_T,  "Top(TT+ST) CF", "p");
  leg->AddEntry(CF_W,  "W(#rightarrowl#nu)+jet CF", "p");
  leg->Draw();


  /*
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
   */
  c1->SaveAs("plot/"+period+obj+"CF_"+region+".png");
  return make_tuple(CF_Q, CF_T, CF_W);
}

tuple<TGraphAsymmErrors*, TGraphAsymmErrors*> LCorrection(TString period, TString region, TString obj, TString sample){
  delete gROOT->FindObject("c1");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();
  TString dir = "/Users/huhchanggi/temp/200702/";
  TFile* file0 = TFile::Open(dir+sample);

  TString histname[2][9];
  file0->cd("Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins");
  histname[0][0] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Data_"+period+"_CR_L17_"+region;
  histname[0][1] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Multijet_"+period+"_CR_L17_"+region;
  histname[0][2] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Top_"+period+"_CR_L17_"+region;
  histname[0][3] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/TT_powheg_pythia8_"+period+"_CR_L17_"+region;
  histname[0][4] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/WToLNu_"+period+"_CR_L17_"+region;
  histname[0][5] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/ZToNuNu_"+period+"_CR_L17_"+region;
  histname[0][6] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Multiboson_"+period+"_CR_L17_"+region;
  histname[0][7] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/GJets_"+period+"_CR_L17_"+region;
  histname[0][8] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/DYToLL_"+period+"_CR_L17_"+region;
  histname[1][0] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Data_"+period+"_CR_LTop17_"+region;
  histname[1][1] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Multijet_"+period+"_CR_LTop17_"+region;
  histname[1][2] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Top_"+period+"_CR_LTop17_"+region;
  histname[1][3] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/TT_powheg_pythia8_"+period+"_CR_LTop17_"+region;
  histname[1][4] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/WToLNu_"+period+"_CR_LTop17_"+region;
  histname[1][5] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/ZToNuNu_"+period+"_CR_LTop17_"+region;
  histname[1][6] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/Multiboson_"+period+"_CR_LTop17_"+region;
  histname[1][7] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/GJets_"+period+"_CR_LTop17_"+region;
  histname[1][8] = "/Counts_vs_"+obj+"Bins/Syst_vs_"+obj+"Bins/DYToLL_"+period+"_CR_LTop17_"+region;

  TH2D* bkg2D[2][9];
  TH1D* bkg[2][9];
  TH1D* bkg_up[2][9];
  TH1D* bkg_dw[2][9];
  TH1D* top[2];
  TH1D* CF[2];
  TH1D* CF_up[3];
  TH1D* CF_dw[3];

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
      name = "bkg_up_"+to_string(i)+to_string(j);
      bkg_up[i][j] = new TH1D(name.c_str(), "", binsize, bin);
      name = "bkg_dw_"+to_string(i)+to_string(j);
      bkg_dw[i][j] = new TH1D(name.c_str(), "", binsize, bin);
      for(int k=1;k<=bkg2D[0][0]->GetNbinsX();k++) {
        if(bkg2D[i][j] == NULL) {
          bkg[i][j]->SetBinContent(k, 0);
          bkg[i][j]->SetBinError(k, 0);
        }
        else {
          bkg[i][j]->SetBinContent(k, bkg2D[i][j]->GetBinContent(k,1));
          bkg[i][j]->SetBinError(k, bkg2D[i][j]->GetBinError(k,1));
          bkg_up[i][j]->SetBinContent(k, bkg2D[i][j]->GetBinContent(k,1)+bkg2D[i][j]->GetBinError(k,1));
          bkg_up[i][j]->SetBinError(k, bkg2D[i][j]->GetBinError(k,1));
          bkg_dw[i][j]->SetBinContent(k, bkg2D[i][j]->GetBinContent(k,1)-bkg2D[i][j]->GetBinError(k,1));
          bkg_dw[i][j]->SetBinError(k, bkg2D[i][j]->GetBinError(k,1));
        }
      }
    }
    bkg[i][2]->Add(bkg[i][3]);
    bkg_up[i][2]->Add(bkg_up[i][3]);
    bkg_dw[i][2]->Add(bkg_dw[i][3]);
    CF[i] = (TH1D*)bkg[i][0]->Clone();
    CF[i]->Add(bkg[i][1], -1);
    CF_up[i] = (TH1D*)bkg_up[i][0]->Clone();
    CF_up[i]->Add(bkg_up[i][1], -1);
    CF_dw[i] = (TH1D*)bkg_dw[i][0]->Clone();
    CF_dw[i]->Add(bkg_dw[i][1], -1);
    for(int j=5; j<9; j++){
      CF[i]->Add(bkg[i][j], -1);
      CF_up[i]->Add(bkg_up[i][j], -1);
      CF_dw[i]->Add(bkg_dw[i][j], -1);
    }
  }

  TMatrix D(2,2);
  TMatrix L(2,2);
  TMatrix T(2,2);
  TMatrix D_up(2,2);
  TMatrix L_up(2,2);
  TMatrix T_up(2,2);
  TMatrix D_dw(2,2);
  TMatrix L_dw(2,2);
  TMatrix T_dw(2,2);
  /*
     string name_L = "CF_L_"+(string)period+"_"+(string)region;
     string name_T = "CF_LT_"+(string)period+"_"+(string)region;
     TH1D* CF_L = new TH1D(name_L.c_str(), "", binsize,bin);
     TH1D* CF_T = new TH1D(name_T.c_str(), "", binsize,bin);
   */

  double cf_L[binsize];
  double cf_L_up[binsize];
  double cf_L_dw[binsize];
  double cf_T[binsize];
  double cf_T_up[binsize];
  double cf_T_dw[binsize];
  double binerr_x[binsize];
  double binx[binsize];

  for(int i=1; i<= binsize;i++){
    for(int j=0; j<2;j++){
      D[j][0] = bkg[j][4]->GetBinContent(i);
      D[j][1] = bkg[j][2]->GetBinContent(i);
      L[j][0] = CF[j]->GetBinContent(i);
      L[j][1] = bkg[j][2]->GetBinContent(i);
      T[j][0] = bkg[j][4]->GetBinContent(i);
      T[j][1] = CF[j]->GetBinContent(i);
      D_up[j][0] = bkg_up[j][4]->GetBinContent(i);
      D_up[j][1] = bkg_up[j][2]->GetBinContent(i);
      L_up[j][0] = CF_up[j]->GetBinContent(i);
      L_up[j][1] = bkg_up[j][2]->GetBinContent(i);
      T_up[j][0] = bkg_up[j][4]->GetBinContent(i);
      T_up[j][1] = CF_up[j]->GetBinContent(i);
      D_dw[j][0] = bkg_dw[j][4]->GetBinContent(i);
      D_dw[j][1] = bkg_dw[j][2]->GetBinContent(i);
      L_dw[j][0] = CF_dw[j]->GetBinContent(i);
      L_dw[j][1] = bkg_dw[j][2]->GetBinContent(i);
      T_dw[j][0] = bkg_dw[j][4]->GetBinContent(i);
      T_dw[j][1] = CF_dw[j]->GetBinContent(i);
      binerr_x[i-1] = CF[j]->GetXaxis()->GetBinWidth(i)/2;
      binx[i-1] = CF[j]->GetXaxis()->GetBinCenter(i);
    }
    //CF_L->SetBinContent(i, L.Determinant()/D.Determinant());
    //CF_T->SetBinContent(i, T.Determinant()/D.Determinant());
    cf_L[i-1] = L.Determinant()/D.Determinant();
    cf_T[i-1] = T.Determinant()/D.Determinant();
    cf_L_up[i-1] = abs((L.Determinant()/D.Determinant())-(L_up.Determinant()/D_up.Determinant()));
    cf_T_up[i-1] = abs((T.Determinant()/D.Determinant())-(T_up.Determinant()/D_up.Determinant()));
    cf_L_dw[i-1] = abs((L.Determinant()/D.Determinant())-(L_dw.Determinant()/D_dw.Determinant()));
    cf_T_dw[i-1] = abs((T.Determinant()/D.Determinant())-(T_dw.Determinant()/D_dw.Determinant()));
  }

  auto CF_L = new TGraphAsymmErrors(binsize, binx, cf_L, binerr_x, binerr_x, cf_L_up, cf_L_dw);
  auto CF_T = new TGraphAsymmErrors(binsize, binx, cf_T, binerr_x, binerr_x, cf_T_up, cf_T_dw);

  string name_L = "CF_L_"+(string)period+"_"+(string)region;
  string name_T = "CF_LT_"+(string)period+"_"+(string)region;
  CF_L->SetName(name_L.c_str());
  CF_T->SetName(name_T.c_str());

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
  //CF_L->Draw("P");
  //CF_T->Draw("Psame");
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(CF_L);
  mg->Add(CF_T);
  mg->Draw("AP");


  auto leg = new TLegend(0.45,0.77,0.9,0.9);
  if(TString(region).Contains("1Boost"))     leg->SetHeader("1 boost jet with Lep+MET final state");
  else if(TString(region).Contains("2Boost")) leg->SetHeader("#geq 2 boost jet with Lep+MET final state");
  else leg->SetHeader("boost jet with Lep+MET final state");
  leg->AddEntry(CF_L,  "W(#rightarrowl#nu)+jet CF", "p");
  leg->AddEntry(CF_T,  "Top(TT+ST) CF", "p");
  leg->Draw();

  /*
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
   */

  c1->SaveAs("plot/"+period+obj+"_L_CF_"+region+".png");
  return make_tuple(CF_L, CF_T);
}

void BkgCorr_solver_NF(){
  tuple<TGraphAsymmErrors*, TGraphAsymmErrors*, TGraphAsymmErrors*> CFs;
  tuple<TGraphAsymmErrors*, TGraphAsymmErrors*> LCFs;
  TGraphAsymmErrors* h1[5][6][2];
  //CFs = Correction("2016", "1Boost", "MRR2", "run_2020_07_07_QTW.root");
  //h1[0][0][0] = get<0>(CFs); h1[1][0][0] = get<1>(CFs); h1[2][0][0] = get<2>(CFs);
  //CFs = Correction("2016", "2Boost", "MRR2", "run_2020_07_07_QTW.root");
  //h1[0][1][0] = get<0>(CFs); h1[1][1][0] = get<1>(CFs); h1[2][1][0] = get<2>(CFs);
  //CFs = Correction("2017", "1Boost", "MRR2", "run_2020_07_03_QTW.root");
  //h1[0][2][0] = get<0>(CFs); h1[1][2][0] = get<1>(CFs); h1[2][2][0] = get<2>(CFs);
  CFs = Correction("2017", "2Boost", "MRR2", "run_2020_07_03_QTW.root");
  h1[0][3][0] = get<0>(CFs); h1[1][3][0] = get<1>(CFs); h1[2][3][0] = get<2>(CFs);
  CFs = Correction("2018", "1Boost", "MRR2", "run_2020_07_03_QTW.root");
  h1[0][4][0] = get<0>(CFs); h1[1][4][0] = get<1>(CFs); h1[2][4][0] = get<2>(CFs);
  CFs = Correction("2018", "2Boost", "MRR2", "run_2020_07_03_QTW.root");
  h1[0][5][0] = get<0>(CFs); h1[1][5][0] = get<1>(CFs); h1[2][5][0] = get<2>(CFs);

  LCFs = LCorrection("2016", "1Boost", "MRR2","run_2020_07_06_L.root");
  h1[3][0][0] = get<0>(LCFs); h1[4][0][0] = get<1>(LCFs);
  LCFs = LCorrection("2016", "2Boost", "MRR2","run_2020_07_06_L.root");
  h1[3][1][0] = get<0>(LCFs); h1[4][1][0] = get<1>(LCFs);
  LCFs = LCorrection("2017", "1Boost", "MRR2","run_2020_07_01_L.root");
  h1[3][2][0] = get<0>(LCFs); h1[4][2][0] = get<1>(LCFs);
  LCFs = LCorrection("2017", "2Boost", "MRR2","run_2020_07_01_L.root");
  h1[3][3][0] = get<0>(LCFs); h1[4][3][0] = get<1>(LCFs);
  LCFs = LCorrection("2018", "1Boost", "MRR2","run_2020_07_01_L.root");
  h1[3][4][0] = get<0>(LCFs); h1[4][4][0] = get<1>(LCFs);
  LCFs = LCorrection("2018", "2Boost", "MRR2","run_2020_07_01_L.root");
  h1[3][5][0] = get<0>(LCFs); h1[4][5][0] = get<1>(LCFs);

  CFs = Correction("2016", "1Boost", "NJet", "run_2020_07_08_QTW.root");
  h1[0][0][1] = get<0>(CFs); h1[1][0][1] = get<1>(CFs); h1[2][0][1] = get<2>(CFs);
  CFs = Correction("2016", "2Boost", "NJet", "run_2020_07_08_QTW.root");
  h1[0][1][1] = get<0>(CFs); h1[1][1][1] = get<1>(CFs); h1[2][1][1] = get<2>(CFs);
  CFs = Correction("2017", "1Boost", "NJet", "run_2020_07_04_QTW.root");
  h1[0][2][1] = get<0>(CFs); h1[1][2][1] = get<1>(CFs); h1[2][2][1] = get<2>(CFs);
  CFs = Correction("2017", "2Boost", "NJet", "run_2020_07_04_QTW.root");
  h1[0][3][1] = get<0>(CFs); h1[1][3][1] = get<1>(CFs); h1[2][3][1] = get<2>(CFs);
  CFs = Correction("2018", "1Boost", "NJet", "run_2020_07_04_QTW.root");
  h1[0][4][1] = get<0>(CFs); h1[1][4][1] = get<1>(CFs); h1[2][4][1] = get<2>(CFs);
  CFs = Correction("2018", "2Boost", "NJet", "run_2020_07_04_QTW.root");
  h1[0][5][1] = get<0>(CFs); h1[1][5][1] = get<1>(CFs); h1[2][5][1] = get<2>(CFs);

  LCFs = LCorrection("2016", "1Boost", "NJet","run_2020_07_07_L.root");
  h1[3][0][1] = get<0>(LCFs); h1[4][0][1] = get<1>(LCFs);
  LCFs = LCorrection("2016", "2Boost", "NJet","run_2020_07_07_L.root");
  h1[3][1][1] = get<0>(LCFs); h1[4][1][1] = get<1>(LCFs);
  LCFs = LCorrection("2017", "1Boost", "NJet","run_2020_07_02_L.root");
  h1[3][2][1] = get<0>(LCFs); h1[4][2][1] = get<1>(LCFs);
  LCFs = LCorrection("2017", "2Boost", "NJet","run_2020_07_02_L.root");
  h1[3][3][1] = get<0>(LCFs); h1[4][3][1] = get<1>(LCFs);
  LCFs = LCorrection("2018", "1Boost", "NJet","run_2020_07_02_L.root");
  h1[3][4][1] = get<0>(LCFs); h1[4][4][1] = get<1>(LCFs);
  LCFs = LCorrection("2018", "2Boost", "NJet","run_2020_07_02_L.root");
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
