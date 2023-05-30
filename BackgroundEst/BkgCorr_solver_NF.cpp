tuple<vector<double>, vector<double>, vector<double>, vector<double>, vector<double>, vector<double>> CalcErrors(TH1D* d1[3], TH1D* b1[3][9], string sample) {
  cout << sample << endl;
  TRandom3 r(0);
  gRandom->SetSeed(time(0));
  TFile* output;
  TString name1, name2;
  if(sample=="MRR2Bin20161Boost") output = new TFile("CF_Error.root", "recreate");
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
  std::vector<double> temp1, temp2, temp3, temp4, temp5, temp6;

  for(int i=1; i<= data[0]->GetNbinsX();i++) {
    if(TString(sample).Contains("MRR2") && TString(sample).Contains("APV")) num = 1;
    name = "QError_"+sample+to_string(i);
    QError[i-1] = new TH1D(name.c_str(), "", 4000, -2., 2.);
    c1 = new TCanvas(name.c_str(), "", 900, 900);
    name = "TError_"+sample+to_string(i);
    TError[i-1] = new TH1D(name.c_str(), "", 4000, -2., 2.);
    c2 = new TCanvas(name.c_str(), "", 900, 900);
    name = "WError_"+sample+to_string(i);
    WError[i-1] = new TH1D(name.c_str(), "", 4000, -2., 2.);
    c3 = new TCanvas(name.c_str(), "", 900, 900);
    if(i == 1 && TString(sample).Contains("MRR2")) {
      QError[i-1]->Rebin(1);
    } else if(i > 3 && TString(sample).Contains("MRR2")) {
      QError[i-1]->Rebin(5);
      TError[i-1]->Rebin(5);
      WError[i-1]->Rebin(5);
    } else { 
      QError[i-1]->Rebin(2);
      TError[i-1]->Rebin(2);
      WError[i-1]->Rebin(2);
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
      if(j%200000 == 0) cout << j << ", " << (double)j/num*100 << "[%] finsihed" << endl;
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
      if(!Q.Determinant() || !T.Determinant() || !W.Determinant()) {err_Q.push_back(0); err_T.push_back(0); err_W.push_back(0);continue;}
      QError[i-1]->Fill(nom_Q-(Q.Determinant()/D.Determinant()));
      TError[i-1]->Fill(nom_T-(T.Determinant()/D.Determinant()));
      WError[i-1]->Fill(nom_W-(W.Determinant()/D.Determinant()));
      err_Q.push_back(nom_Q-Q.Determinant()/D.Determinant());
      err_T.push_back(nom_T-T.Determinant()/D.Determinant());
      err_W.push_back(nom_W-W.Determinant()/D.Determinant());
    }

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
    temp1.push_back(bin_up);
    temp4.push_back(bin_dw);

    if(TString(sample).Contains("MRR2Bin")) {
      QError[i-1]->GetXaxis()->SetRangeUser(-1.0,1.0);
      TError[i-1]->GetXaxis()->SetRangeUser(-1.0,1.0);
      WError[i-1]->GetXaxis()->SetRangeUser(-1.0,1.0);
      if(i==1) {
        QError[i-1]->GetXaxis()->SetRangeUser(-0.05,0.05);
        TError[i-1]->GetXaxis()->SetRangeUser(-0.2,0.2);
        WError[i-1]->GetXaxis()->SetRangeUser(-0.2,0.2);
      } else if(i==2) {
        QError[i-1]->GetXaxis()->SetRangeUser(-0.05,0.05);
        TError[i-1]->GetXaxis()->SetRangeUser(-0.2,0.2);
        WError[i-1]->GetXaxis()->SetRangeUser(-0.2,0.2);
      } else if(i==3) {
        QError[i-1]->GetXaxis()->SetRangeUser(-0.5,0.5);
        TError[i-1]->GetXaxis()->SetRangeUser(-0.2,0.2);
        WError[i-1]->GetXaxis()->SetRangeUser(-0.5,0.5);
      } else if(i==4) {
        WError[i-1]->GetXaxis()->SetRangeUser(-0.5,0.5);
      } else if(i>=5) {
        QError[i-1]->GetXaxis()->SetRangeUser(-2.0,2.0);
      }
    } else {
      QError[i-1]->GetXaxis()->SetRangeUser(-0.2,0.2);
      TError[i-1]->GetXaxis()->SetRangeUser(-1.0,1.0);
      WError[i-1]->GetXaxis()->SetRangeUser(-0.4,0.4);
      if(i>2) {
        TError[i-1]->GetXaxis()->SetRangeUser(-0.2,0.2);
      }
    }

    c1->cd();
    QError[i-1]->GetXaxis()->SetTitle("CF-Iteration");
    QError[i-1]->Draw();

    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.03);
    latex->SetTextAlign(31); // align right
    latex->DrawLatex(0.90, 0.93, Form("138 fb^{-1} (#sqrt{s} = 13 TeV)"));
    latex->SetTextAlign(11); // align left
    latex->DrawLatex(0.10,0.93,"CMS Work in Progress");

    auto leg = new TLegend(0.60,0.80,0.9,0.88);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.02);
    if(TString(sample).Contains("2016APV")) name1 = "2016APV";
    else if(TString(sample).Contains("2016")) name1 = "2016";
    else if(TString(sample).Contains("2017")) name1 = "2017";
    else if(TString(sample).Contains("2018")) name1 = "2018";
    else name1 = "Run2";
    if(TString(sample).Contains("1Boost")) name2 = "1 Boost jet";
    else if(TString(sample).Contains("2Boost")) name2 = "2 Boost jet";
    leg->SetHeader(name1+", CF_{QCD}"+name2+" "+(TString)(TString)to_string(i)+" bin");
    leg->Draw();

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
    temp2.push_back(bin_up);
    temp5.push_back(bin_dw);

    c2->cd();
    TError[i-1]->GetXaxis()->SetTitle("CF-Iteration");
    TError[i-1]->Draw();

    latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.03);
    latex->SetTextAlign(31); // align right
    latex->DrawLatex(0.90, 0.93, Form("138 fb^{-1} (#sqrt{s} = 13 TeV)"));
    latex->SetTextAlign(11); // align left
    latex->DrawLatex(0.10,0.93,"CMS Work in Progress");

    leg = new TLegend(0.60,0.80,0.9,0.88);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.02);
    if(TString(sample).Contains("2016APV")) name1 = "2016APV";
    else if(TString(sample).Contains("2016")) name1 = "2016";
    else if(TString(sample).Contains("2017")) name1 = "2017";
    else if(TString(sample).Contains("2018")) name1 = "2018";
    else name1 = "Run2";
    if(TString(sample).Contains("1Boost")) name2 = "1 Boost jet";
    else if(TString(sample).Contains("2Boost")) name2 = "2 Boost jet";
    leg->SetHeader(name1+", CF_{Top}"+name2+" "+(TString)to_string(i)+" bin");
    leg->Draw();

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
    temp3.push_back(bin_up);
    temp6.push_back(bin_dw);

    c3->cd();
    WError[i-1]->GetXaxis()->SetTitle("CF-Iteration");
    WError[i-1]->Draw();

    latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.03);
    latex->SetTextAlign(31); // align right
    latex->DrawLatex(0.90, 0.93, Form("138 fb^{-1} (#sqrt{s} = 13 TeV)"));
    latex->SetTextAlign(11); // align left
    latex->DrawLatex(0.10,0.93,"CMS Work in Progress");

    leg = new TLegend(0.60,0.80,0.9,0.88);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.02);
    if(TString(sample).Contains("2016APV")) name1 = "2016APV";
    else if(TString(sample).Contains("2016")) name1 = "2016";
    else if(TString(sample).Contains("2017")) name1 = "2017";
    else if(TString(sample).Contains("2018")) name1 = "2018";
    else name1 = "Run2";
    if(TString(sample).Contains("1Boost")) name2 = "1 Boost jet";
    else if(TString(sample).Contains("2Boost")) name2 = "2 Boost jet";
    leg->SetHeader(name1+", CF_{Wjets}"+name2+" "+(TString)to_string(i)+" bin");
    leg->Draw();
    line = new TLine(bin_up,0,bin_up,1.05*WError[i-1]->GetMaximum());
    line->SetLineColor(kRed);
    line->Draw("same");
    line = new TLine(bin_dw,0,bin_dw,1.05*WError[i-1]->GetMaximum());
    line->SetLineColor(kRed);
    line->Draw("same");
    c3->Write();
  }
  output->Close();
  return make_tuple(temp1, temp2, temp3, temp4, temp5, temp6);
}

tuple<vector<double>, vector<double>, vector<double>, vector<double>> CalcLErrors(TH1D* d1[2], TH1D* b1[2][9], string sample) {
  cout << sample << endl;
  TRandom3 r(0);
  gRandom->SetSeed(time(0));
  TFile* output;
  TString name1, name2;
  output = new TFile("CF_Error.root", "update");
  TH1D* data[2];
  TH1D* LB[2];
  TH1D* TB[2];
  for(int i=0; i<2; i++) {
    data[i]  = (TH1D*)d1[i]->Clone();
    LB[i]    = (TH1D*)b1[i][4]->Clone();
    TB[i]    = (TH1D*)b1[i][2]->Clone();
  }
  TMatrix D(2,2);
  TMatrix L(2,2);
  TMatrix T(2,2);
  string name;
  TH1D* LError[data[0]->GetNbinsX()];
  TH1D* TError[data[0]->GetNbinsX()];
  double nom_L, nom_T;
  vector<double> err_L, err_T;
  TF1* gaus_D;
  TF1* gaus_L;
  TF1* gaus_T;
  double ran_L, ran_T, ran_D, temp;
  double bin_up, bin_dw;
  int num=1e6;
  TH1D* Err_Up_L = (TH1D*)d1[0]->Clone();
  TH1D* Err_Up_T = (TH1D*)d1[0]->Clone();
  TH1D* Err_Dw_L = (TH1D*)d1[0]->Clone();
  TH1D* Err_Dw_T = (TH1D*)d1[0]->Clone();
  TCanvas* c1, *c2, *c3;
  TLine* line;
  std::vector<double> temp1, temp2, temp3, temp4, temp5, temp6;

  for(int i=1; i<= data[0]->GetNbinsX();i++) {
    name = "LError_"+sample+to_string(i);
    if(TString(sample).Contains("Iso")) name = "NonIso_WError_"+sample+to_string(i);
    LError[i-1] = new TH1D(name.c_str(), "", 4000, -4., 4.);
    c1 = new TCanvas(name.c_str(), "", 900, 900);
    name = "LTError_"+sample+to_string(i);
    if(TString(sample).Contains("Iso")) name = "NonIso_TError_"+sample+to_string(i);
    TError[i-1] = new TH1D(name.c_str(), "", 4000, -4., 4.);
    c2 = new TCanvas(name.c_str(), "", 900, 900);
    if(i > 1 && TString(sample).Contains("MRR2")) {
      LError[i-1]->Rebin(5);
      TError[i-1]->Rebin(5);
    }
    for(int m=0; m<2; m++) { // region
      D[m][0] = LB[m]->GetBinContent(i);
      D[m][1] = TB[m]->GetBinContent(i);
      L[m][0] = data[m]->GetBinContent(i);
      L[m][1] = TB[m]->GetBinContent(i);
      T[m][0] = LB[m]->GetBinContent(i);
      T[m][1] = data[m]->GetBinContent(i);
    }
    nom_L = L.Determinant()/D.Determinant();
    nom_T = T.Determinant()/D.Determinant();
    for(int j=0; j<num; j++) {
      if(j%200000 == 0) cout << j << ", " << (double)j/num*100 << "[%] finsihed" << endl;
      for(int k=0; k<2; k++) { // region
        ran_D = gRandom->Gaus(data[k]->GetBinContent(i), data[k]->GetBinError(i));
        ran_L = gRandom->Gaus(LB[k]->GetBinContent(i), LB[k]->GetBinError(i));
        ran_T = gRandom->Gaus(TB[k]->GetBinContent(i), TB[k]->GetBinError(i));
        D[k][0] = ran_L;
        D[k][1] = ran_T;
        L[k][0] = ran_D;
        L[k][1] = ran_T;
        T[k][0] = ran_L;
        T[k][1] = ran_D;
      }
      LError[i-1]->Fill(nom_L-(L.Determinant()/D.Determinant()));
      TError[i-1]->Fill(nom_T-(T.Determinant()/D.Determinant()));
      err_L.push_back(nom_L-L.Determinant()/D.Determinant());
      err_T.push_back(nom_T-T.Determinant()/D.Determinant());
    }

    temp = 0;
    bin_up = 0; bin_dw = 0;
    for(int bin=LError[i-1]->GetMaximumBin();bin<=LError[i-1]->GetNbinsX();bin++) {
      temp += LError[i-1]->GetBinContent(bin);
      temp += LError[i-1]->GetBinContent(2*LError[i-1]->GetMaximumBin()-bin-1);
      if(temp/num > 0.68) continue;
      else {
        bin_up = LError[i-1]->GetBinCenter(bin); 
        bin_dw = LError[i-1]->GetBinCenter(2.*LError[i-1]->GetMaximumBin()-bin-1);
      }
    }
    cout << "L : " << bin_up << ", " << bin_dw << endl;
    Err_Up_L->SetBinContent(i, bin_up);
    Err_Up_L->SetBinError(i, 0);
    Err_Dw_L->SetBinContent(i, bin_dw);
    Err_Dw_L->SetBinError(i, 0);
    temp1.push_back(bin_up);
    temp3.push_back(bin_dw);

    if(TString(sample).Contains("MRR2")) {
      LError[i-1]->GetXaxis()->SetRangeUser(-1.0,1.0);
      TError[i-1]->GetXaxis()->SetRangeUser(-0.5,0.5);
      if(i==1) {
        LError[i-1]->GetXaxis()->SetRangeUser(-0.3,0.3);
        TError[i-1]->GetXaxis()->SetRangeUser(-0.2,0.2);
      } else if(i==2) {
        LError[i-1]->GetXaxis()->SetRangeUser(-0.3,0.3);
        TError[i-1]->GetXaxis()->SetRangeUser(-0.2,0.2);
      } else if(i==3) {
        LError[i-1]->GetXaxis()->SetRangeUser(-0.5,0.5);
        TError[i-1]->GetXaxis()->SetRangeUser(-0.2,0.2);
      } else if(i==4) {
        LError[i-1]->GetXaxis()->SetRangeUser(-0.5,0.5);
        TError[i-1]->GetXaxis()->SetRangeUser(-0.2,0.2);
      }
    } else {
      LError[i-1]->GetXaxis()->SetRangeUser(-0.5,0.5);
      TError[i-1]->GetXaxis()->SetRangeUser(-0.4,0.4);
      if(i<2) {
        TError[i-1]->GetXaxis()->SetRangeUser(-1.0,1.0);
      } else if(i>3) {
        TError[i-1]->GetXaxis()->SetRangeUser(-0.2,0.2);
      }
    }

    c1->cd();
    LError[i-1]->GetXaxis()->SetTitle("CF-Iteration");
    LError[i-1]->Draw();

    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.03);
    latex->SetTextAlign(31); // align right
    latex->DrawLatex(0.90, 0.93, Form("138 fb^{-1} (#sqrt{s} = 13 TeV)"));
    latex->SetTextAlign(11); // align left
    latex->DrawLatex(0.10,0.93,"CMS Work in Progress");

    auto leg = new TLegend(0.60,0.80,0.9,0.88);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.02);
    if(TString(sample).Contains("2016APV")) name1 = "2016APV";
    else if(TString(sample).Contains("2016")) name1 = "2016";
    else if(TString(sample).Contains("2017")) name1 = "2017";
    else if(TString(sample).Contains("2018")) name1 = "2018";
    else name1 = "Run2";
    if(TString(sample).Contains("1Boost")) name2 = "1 Boost jet";
    else if(TString(sample).Contains("2Boost")) name2 = "2 Boost jet";
    leg->SetHeader(name1+", CF_{Lep+MET, Wjets}"+name2+" "+(TString)to_string(i)+" bin");
    if(TString(sample).Contains("NonIso")) leg->SetHeader(name1+", CF_{NonIso, W}"+name2+" "+(TString)to_string(i)+" bin");
    leg->Draw();

    line = new TLine(bin_up,0,bin_up,1.05*LError[i-1]->GetMaximum());
    line->SetLineColor(kRed);
    line->Draw("same");
    line = new TLine(bin_dw,0,bin_dw,1.05*LError[i-1]->GetMaximum());
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
    temp2.push_back(bin_up);
    temp4.push_back(bin_dw);

    c2->cd();
    TError[i-1]->GetXaxis()->SetTitle("CF-Iteration");
    TError[i-1]->Draw();

    latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.03);
    latex->SetTextAlign(31); // align right
    latex->DrawLatex(0.90, 0.93, Form("138 fb^{-1} (#sqrt{s} = 13 TeV)"));
    latex->SetTextAlign(11); // align left
    latex->DrawLatex(0.10,0.93,"CMS Work in Progress");

    leg = new TLegend(0.60,0.80,0.9,0.88);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.02);
    if(TString(sample).Contains("2016APV")) name1 = "2016APV";
    else if(TString(sample).Contains("2016")) name1 = "2016";
    else if(TString(sample).Contains("2017")) name1 = "2017";
    else if(TString(sample).Contains("2018")) name1 = "2018";
    else name1 = "Run2";
    if(TString(sample).Contains("1Boost")) name2 = "1 Boost jet";
    else if(TString(sample).Contains("2Boost")) name2 = "2 Boost jet";
    leg->SetHeader(name1+", CF_{Lep+MET, Top}"+name2+" "+(TString)to_string(i)+" bin");
    if(TString(sample).Contains("NonIso")) leg->SetHeader(name1+", CF_{NonIso, Top}"+name2+" "+(TString)to_string(i)+" bin");
    leg->Draw();

    line = new TLine(bin_up,0,bin_up,1.05*TError[i-1]->GetMaximum());
    line->SetLineColor(kRed);
    line->Draw("same");
    line = new TLine(bin_dw,0,bin_dw,1.05*TError[i-1]->GetMaximum());
    line->SetLineColor(kRed);
    line->Draw("same");
    c2->Write();
  }
  output->Close();
  return make_tuple(temp1, temp2, temp3, temp4);
}

tuple<TGraphAsymmErrors*, TGraphAsymmErrors*, TGraphAsymmErrors*> Correction(TString period, TString region, TString obj, TString sample){
  cout << "period : " << period << ", region : " << region << ", obj : " << obj << endl;
  delete gROOT->FindObject("c1");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();
  TString dir = "/Users/chuh/Dropbox/Analysis/razor/";
  TFile* file0 = TFile::Open(dir+sample);

  TString histname[3][9];
  TString path = "/Counts_vs_"+obj+"/Syst_vs_"+obj;
  if(TString(obj).Contains("NJet")||TString(period).Contains("run2")) path = "/"+obj;
  histname[0][0] = path+"/Data_"+period+"_CR_QCD17_"+region;
  histname[0][1] = path+"/Multijet_"+period+"_CR_QCD17_"+region;
  histname[0][2] = path+"/Top_"+period+"_CR_QCD17_"+region;
  histname[0][3] = path+"/TT_powheg_pythia8_"+period+"_CR_QCD17_"+region;
  histname[0][4] = path+"/WToLNu_"+period+"_CR_QCD17_"+region;
  histname[0][5] = path+"/ZToNuNu_"+period+"_CR_QCD17_"+region;
  histname[0][6] = path+"/Multiboson_"+period+"_CR_QCD17_"+region;
  histname[0][7] = path+"/GJets_"+period+"_CR_QCD17_"+region;
  histname[0][8] = path+"/DYToLL_"+period+"_CR_QCD17_"+region;
  histname[1][0] = path+"/Data_"+period+"_CR_Top17_"+region;
  histname[1][1] = path+"/Multijet_"+period+"_CR_Top17_"+region;
  histname[1][2] = path+"/Top_"+period+"_CR_Top17_"+region;
  histname[1][3] = path+"/TT_powheg_pythia8_"+period+"_CR_Top17_"+region;
  histname[1][4] = path+"/WToLNu_"+period+"_CR_Top17_"+region;
  histname[1][5] = path+"/ZToNuNu_"+period+"_CR_Top17_"+region;
  histname[1][6] = path+"/Multiboson_"+period+"_CR_Top17_"+region;
  histname[1][7] = path+"/GJets_"+period+"_CR_Top17_"+region;
  histname[1][8] = path+"/DYToLL_"+period+"_CR_Top17_"+region;
  histname[2][0] = path+"/Data_"+period+"_CR_W17_"+region;
  histname[2][1] = path+"/Multijet_"+period+"_CR_W17_"+region;
  histname[2][2] = path+"/Top_"+period+"_CR_W17_"+region;
  histname[2][3] = path+"/TT_powheg_pythia8_"+period+"_CR_W17_"+region;
  histname[2][4] = path+"/WToLNu_"+period+"_CR_W17_"+region;
  histname[2][5] = path+"/ZToNuNu_"+period+"_CR_W17_"+region;
  histname[2][6] = path+"/Multiboson_"+period+"_CR_W17_"+region;
  histname[2][7] = path+"/GJets_"+period+"_CR_W17_"+region;
  histname[2][8] = path+"/DYToLL_"+period+"_CR_W17_"+region;

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
      name = "bkg_"+(TString)to_string(i)+(TString)to_string(j);
      bkg[i][j] = new TH1D(name.c_str(), "", binsize, bin);
      name = "bkg_up_"+(TString)to_string(i)+(TString)to_string(j);
      bkg_up[i][j] = new TH1D(name.c_str(), "", binsize, bin);
      name = "bkg_dw_"+(TString)to_string(i)+(TString)to_string(j);
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

  int flag = 0;
  double qbinerr_x[binsize-1];
  double qbinx[binsize-1];

  tuple<vector<double>, vector<double>, vector<double>, vector<double>, vector<double>, vector<double>> Err_CFs;
  TH1D* temp;
  vector<double> temp1;
  vector<double> temp2;
  vector<double> temp3;
  vector<double> temp4;
  vector<double> temp5;
  vector<double> temp6;

  Err_CFs = CalcErrors(CF, bkg, (string)obj+(string)period+(string)region);
  temp1 = get<0>(Err_CFs);
  temp2 = get<1>(Err_CFs);
  temp3 = get<2>(Err_CFs);
  temp4 = get<3>(Err_CFs);
  temp5 = get<4>(Err_CFs);
  temp6 = get<5>(Err_CFs);

  for(int i=1; i<= binsize;i++){

    cf_Q_up[i-1] = abs(temp1.at(i-1));
    cf_T_up[i-1] = abs(temp2.at(i-1));
    cf_W_up[i-1] = abs(temp3.at(i-1));
    cf_Q_dw[i-1] = abs(temp4.at(i-1));
    cf_T_dw[i-1] = abs(temp5.at(i-1));
    cf_W_dw[i-1] = abs(temp6.at(i-1));

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
      binerr_x[i-1] = CF[j]->GetXaxis()->GetBinWidth(i)/2;
      binx[i-1] = CF[j]->GetXaxis()->GetBinCenter(i);
    }

    cf_Q[i-1] = Q.Determinant()/D.Determinant();
    cf_T[i-1] = T.Determinant()/D.Determinant();
    cf_W[i-1] = W.Determinant()/D.Determinant();

    if(cf_Q[i-1] + temp4.at(i-1) < 0) cf_Q_dw[i-1] = cf_Q[i-1];
    if(cf_T[i-1] + temp5.at(i-1) < 0) cf_T_dw[i-1] = cf_T[i-1];
    if(cf_W[i-1] + temp6.at(i-1) < 0) cf_W_dw[i-1] = cf_W[i-1];


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
      qbinerr_x[2] = binerr_x[2];
      qbinerr_x[3] = binerr_x[3];
      qbinx[0] = binx[0];
      qbinx[1] = binx[1];
      qbinx[2] = binx[2];
      qbinx[3] = binx[3];
      qbinerr_x[4] = (CF[0]->GetXaxis()->GetBinWidth(5) + CF[0]->GetXaxis()->GetBinWidth(6))/2.;
      qbinx[4] = (CF[0]->GetXaxis()->GetBinCenter(5) - CF[0]->GetXaxis()->GetBinWidth(5)/2 + CF[0]->GetXaxis()->GetBinCenter(6) + CF[0]->GetXaxis()->GetBinWidth(6)/2 )/2.;
      cf_Q[4] = Q.Determinant()/D.Determinant();
      cf_Q_dw[4] = abs(temp4.at(4));
      flag = 1;
    }
  }
  auto CF_Q = new TGraphAsymmErrors(binsize, binx, cf_Q, binerr_x, binerr_x, cf_Q_dw, cf_Q_up);
  if(flag == 1) CF_Q = new TGraphAsymmErrors(binsize-1, qbinx, cf_Q, qbinerr_x, qbinerr_x, cf_Q_dw, cf_Q_up);
  auto CF_T = new TGraphAsymmErrors(binsize, binx, cf_T, binerr_x, binerr_x, cf_T_dw, cf_T_up);
  auto CF_W = new TGraphAsymmErrors(binsize, binx, cf_W, binerr_x, binerr_x, cf_W_dw, cf_W_up);

  string name_Q = "CF_Q_"+(string)period+"_"+(string)region;
  string name_T = "CF_T_"+(string)period+"_"+(string)region;
  string name_W = "CF_W_"+(string)period+"_"+(string)region;
  CF_Q->SetName(name_Q.c_str());
  CF_T->SetName(name_T.c_str());
  CF_W->SetName(name_W.c_str());


  TCanvas* c1 = new TCanvas("c1","",700,700);

  c1->SetGridx();
  c1->SetGridy();
  c1->GetFrame()->SetBorderSize(12);
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

  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.03);
  latex->SetTextAlign(31); // align right
  latex->DrawLatex(0.90, 0.93, Form("138 fb^{-1} (#sqrt{s} = 13 TeV)"));
  latex->SetTextAlign(11); // align left
  latex->DrawLatex(0.10,0.93,"CMS Work in Progress");

  auto leg = new TLegend(0.45,0.77,0.9,0.9);
  if(TString(region).Contains("1Boost"))     leg->SetHeader("1 boost jet final state");
  else if(TString(region).Contains("2Boost")) leg->SetHeader("#geq 2 boost jet final state");
  else leg->SetHeader("boost jet final state");
  leg->AddEntry(CF_Q,  "Multijet CF", "p");
  leg->AddEntry(CF_T,  "Top(TT+ST) CF", "p");
  leg->AddEntry(CF_W,  "W(#rightarrowl#nu)+jet CF", "p");
  leg->Draw();

  c1->SaveAs("plot/"+period+obj+"CF_"+region+".png");
  return make_tuple(CF_Q, CF_T, CF_W);
}

tuple<TGraphAsymmErrors*, TGraphAsymmErrors*> LCorrection(TString period, TString region, TString obj, TString sample){
  delete gROOT->FindObject("c1");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();
  TString dir = "/Users/chuh/Dropbox/Analysis/razor/";
  TFile* file0 = TFile::Open(dir+sample);

  TString histname[2][9];
  TString path = "/Counts_vs_"+obj+"/Syst_vs_"+obj;
  if(TString(obj).Contains("NJet")||TString(period).Contains("run2")) path = "/"+obj;
  histname[0][0] = path+"/Data_"+period+"_CR_L17_"+region;
  histname[0][1] = path+"/Multijet_"+period+"_CR_L17_"+region;
  histname[0][2] = path+"/Top_"+period+"_CR_L17_"+region;
  histname[0][3] = path+"/TT_powheg_pythia8_"+period+"_CR_L17_"+region;
  histname[0][4] = path+"/WToLNu_"+period+"_CR_L17_"+region;
  histname[0][5] = path+"/ZToNuNu_"+period+"_CR_L17_"+region;
  histname[0][6] = path+"/Multiboson_"+period+"_CR_L17_"+region;
  histname[0][7] = path+"/GJets_"+period+"_CR_L17_"+region;
  histname[0][8] = path+"/DYToLL_"+period+"_CR_L17_"+region;
  histname[1][0] = path+"/Data_"+period+"_CR_LTop17_"+region;
  histname[1][1] = path+"/Multijet_"+period+"_CR_LTop17_"+region;
  histname[1][2] = path+"/Top_"+period+"_CR_LTop17_"+region;
  histname[1][3] = path+"/TT_powheg_pythia8_"+period+"_CR_LTop17_"+region;
  histname[1][4] = path+"/WToLNu_"+period+"_CR_LTop17_"+region;
  histname[1][5] = path+"/ZToNuNu_"+period+"_CR_LTop17_"+region;
  histname[1][6] = path+"/Multiboson_"+period+"_CR_LTop17_"+region;
  histname[1][7] = path+"/GJets_"+period+"_CR_LTop17_"+region;
  histname[1][8] = path+"/DYToLL_"+period+"_CR_LTop17_"+region;

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
      name = "bkg_"+(TString)to_string(i)+(TString)to_string(j);
      bkg[i][j] = new TH1D(name.c_str(), "", binsize, bin);
      name = "bkg_up_"+(TString)to_string(i)+(TString)to_string(j);
      bkg_up[i][j] = new TH1D(name.c_str(), "", binsize, bin);
      name = "bkg_dw_"+(TString)to_string(i)+(TString)to_string(j);
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

  double cf_L[binsize];
  double cf_L_up[binsize];
  double cf_L_dw[binsize];
  double cf_T[binsize];
  double cf_T_up[binsize];
  double cf_T_dw[binsize];
  double cf_t[binsize-1];
  double cf_t_up[binsize-1];
  double cf_t_dw[binsize-1];
  double binerr_x[binsize];
  double binx[binsize];

  int flag = 0;
  double tbinerr_x[binsize-1];
  double tbinx[binsize-1];


  tuple<vector<double>, vector<double>, vector<double>, vector<double>> Err_CFs;
  TH1D* temp;
  vector<double> temp1;
  vector<double> temp2;
  vector<double> temp3;
  vector<double> temp4;

  Err_CFs = CalcLErrors(CF, bkg, (string)obj+(string)period+(string)region);
  temp1 = get<0>(Err_CFs);
  temp2 = get<1>(Err_CFs);
  temp3 = get<2>(Err_CFs);
  temp4 = get<3>(Err_CFs);

  for(int i=1; i<= binsize;i++){
    cf_L_up[i-1] = abs(temp1.at(i-1));
    cf_T_up[i-1] = abs(temp2.at(i-1));
    cf_L_dw[i-1] = abs(temp3.at(i-1));
    cf_T_dw[i-1] = abs(temp4.at(i-1));

    for(int j=0; j<2;j++){
      D[j][0] = bkg[j][4]->GetBinContent(i);
      D[j][1] = bkg[j][2]->GetBinContent(i);
      L[j][0] = CF[j]->GetBinContent(i);
      L[j][1] = bkg[j][2]->GetBinContent(i);
      T[j][0] = bkg[j][4]->GetBinContent(i);
      T[j][1] = CF[j]->GetBinContent(i);
      binerr_x[i-1] = CF[j]->GetXaxis()->GetBinWidth(i)/2;
      binx[i-1] = CF[j]->GetXaxis()->GetBinCenter(i);
    }
    cf_L[i-1] = L.Determinant()/D.Determinant();
    cf_T[i-1] = T.Determinant()/D.Determinant();
    if(flag == 2 && i > 2) {
      cf_t[i-2] = T.Determinant()/D.Determinant();
      cf_t_up[i-2] = abs(temp2.at(i-1));
      cf_t_dw[i-2] = abs(temp2.at(i-1));
    }

    if(cf_L[i-1] + temp3.at(i-1) < 0) cf_L_dw[i-1] = cf_L[i-1];
    if(cf_T[i-1] + temp4.at(i-1) < 0) cf_T_dw[i-1] = cf_T[i-1];

    if(cf_T[i-1] < 0) {
      cout << "CF is smaller than 0" << endl;
      if(TString(obj).Contains("NJet")) {
        for(int j=0; j<2;j++) {
          D[j][0] = bkg[j][4]->GetBinContent(i) + bkg[j][4]->GetBinContent(i+1);
          D[j][1] = bkg[j][2]->GetBinContent(i) + bkg[j][2]->GetBinContent(i+1);
          L[j][0] = CF[j]->GetBinContent(i) + CF[j]->GetBinContent(i+1);
          L[j][1] = bkg[j][2]->GetBinContent(i) + bkg[j][2]->GetBinContent(i+1);
          T[j][0] = bkg[j][4]->GetBinContent(i) + bkg[j][4]->GetBinContent(i+1);
          T[j][1] = CF[j]->GetBinContent(i) + CF[j]->GetBinContent(i+1);
        }
        tbinerr_x[0] = (CF[0]->GetXaxis()->GetBinWidth(1) + CF[0]->GetXaxis()->GetBinWidth(2))/2.;
        tbinx[0] = (CF[0]->GetXaxis()->GetBinCenter(1) - CF[0]->GetXaxis()->GetBinWidth(1)/2 + CF[0]->GetXaxis()->GetBinCenter(2) + CF[0]->GetXaxis()->GetBinWidth(2)/2 )/2.;
        tbinerr_x[1] = CF[0]->GetXaxis()->GetBinWidth(3)/2;
        tbinerr_x[2] = CF[0]->GetXaxis()->GetBinWidth(4)/2;
        tbinerr_x[3] = CF[0]->GetXaxis()->GetBinWidth(5)/2;
        tbinx[1] = CF[0]->GetXaxis()->GetBinCenter(3);
        tbinx[2] = CF[0]->GetXaxis()->GetBinCenter(4);
        tbinx[3] = CF[0]->GetXaxis()->GetBinCenter(5);
        cf_t[i-1] = T.Determinant()/D.Determinant();
        cf_t_up[i-1] = abs(temp2.at(1));
        cf_t_dw[i-1] = abs(temp4.at(1));
        flag = 2;
      } else {
        for(int j=0; j<2;j++) {
          D[j][0] = bkg[j][4]->GetBinContent(i) + bkg[j][4]->GetBinContent(i-1);
          D[j][1] = bkg[j][2]->GetBinContent(i) + bkg[j][2]->GetBinContent(i-1);
          L[j][0] = CF[j]->GetBinContent(i) + CF[j]->GetBinContent(i-1);
          L[j][1] = bkg[j][2]->GetBinContent(i) + bkg[j][2]->GetBinContent(i-1);
          T[j][0] = bkg[j][4]->GetBinContent(i) + bkg[j][4]->GetBinContent(i-1);
          T[j][1] = CF[j]->GetBinContent(i) + CF[j]->GetBinContent(i-1);
        }
        tbinerr_x[0] = binerr_x[0];
        tbinerr_x[1] = binerr_x[1];
        tbinerr_x[2] = binerr_x[2];
        tbinerr_x[3] = binerr_x[3];
        tbinx[0] = binx[0];
        tbinx[1] = binx[1];
        tbinx[2] = binx[2];
        tbinx[3] = binx[3];
        tbinerr_x[4] = (CF[0]->GetXaxis()->GetBinWidth(i) + CF[0]->GetXaxis()->GetBinWidth(i-1))/2.;
        tbinx[4] = (CF[0]->GetXaxis()->GetBinCenter(i-1) - CF[0]->GetXaxis()->GetBinWidth(i-1)/2 + CF[0]->GetXaxis()->GetBinCenter(i) + CF[0]->GetXaxis()->GetBinWidth(i)/2 )/2.;
        cf_T[i-2] = T.Determinant()/D.Determinant();
        flag = 1;
      }
    }
  }
  //for(int i=0;i<binsize-1;i++) cout << i << ", " << tbinx[i] << ", " << tbinerr_x[i] << ", " << cf_t[i] << endl;

  auto CF_L = new TGraphAsymmErrors(binsize, binx, cf_L, binerr_x, binerr_x, cf_L_dw, cf_L_up);
  auto CF_T = new TGraphAsymmErrors(binsize, binx, cf_T, binerr_x, binerr_x, cf_T_dw, cf_T_up);
  if(flag == 1) CF_T = new TGraphAsymmErrors(binsize-1, tbinx, cf_T, tbinerr_x, tbinerr_x, cf_T_dw, cf_T_up);
  if(flag == 2) CF_T = new TGraphAsymmErrors(binsize-1, tbinx, cf_t, tbinerr_x, tbinerr_x, cf_t_dw, cf_T_up);

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
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(CF_L);
  mg->Add(CF_T);
  mg->Draw("AP");

  auto latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.03);
  latex->SetTextAlign(31); // align right
  latex->DrawLatex(0.90, 0.93, Form("138 fb^{-1} (#sqrt{s} = 13 TeV)"));
  latex->SetTextAlign(11); // align left
  latex->DrawLatex(0.10,0.93,"CMS Work in Progress");

  auto leg = new TLegend(0.45,0.77,0.9,0.9);
  if(TString(region).Contains("1Boost"))     leg->SetHeader("1 boost jet with Lep+MET final state");
  else if(TString(region).Contains("2Boost")) leg->SetHeader("#geq 2 boost jet with Lep+MET final state");
  else leg->SetHeader("boost jet with Lep+MET final state");
  leg->AddEntry(CF_L,  "W(#rightarrowl#nu)+jet CF", "p");
  leg->AddEntry(CF_T,  "Top(TT+ST) CF", "p");
  leg->Draw();

  c1->SaveAs("plot/"+period+obj+"_L_CF_"+region+".png");
  return make_tuple(CF_L, CF_T);
}

tuple<TGraphAsymmErrors*, TGraphAsymmErrors*> NonIsoCorrection(TString period, TString region, TString obj, TString sample){
  delete gROOT->FindObject("c1");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();
  TString dir = "/Users/chuh/Dropbox/Analysis/razor/";
  TFile* file0 = TFile::Open(dir+sample);
  TString iso="NonIso";

  TString histname[2][9];
  TString path = "/Counts_vs_"+obj+"/Syst_vs_"+obj;
  if(TString(obj).Contains("NJet")||TString(period).Contains("run2")) path = "/"+obj;
  histname[0][0] = path+"/Data_"+period+"_CR_NonIso_0b_RMTdPhi_"+region;
  histname[0][1] = path+"/Multijet_"+period+"_CR_NonIso_0b_RMTdPhi_"+region;
  histname[0][2] = path+"/Top_"+period+"_CR_NonIso_0b_RMTdPhi_"+region;
  histname[0][3] = path+"/TT_powheg_pythia8_"+period+"_CR_NonIso_0b_RMTdPhi_"+region;
  histname[0][4] = path+"/WToLNu_"+period+"_CR_NonIso_0b_RMTdPhi_"+region;
  histname[0][5] = path+"/ZToNuNu_"+period+"_CR_NonIso_0b_RMTdPhi_"+region;
  histname[0][6] = path+"/Multiboson_"+period+"_CR_NonIso_0b_RMTdPhi_"+region;
  histname[0][7] = path+"/GJets_"+period+"_CR_NonIso_0b_RMTdPhi_"+region;
  histname[0][8] = path+"/DYToLL_"+period+"_CR_NonIso_0b_RMTdPhi_"+region;
  histname[1][0] = path+"/Data_"+period+"_CR_NonIso_b_RMTdPhi_"+region;
  histname[1][1] = path+"/Multijet_"+period+"_CR_NonIso_b_RMTdPhi_"+region;
  histname[1][2] = path+"/Top_"+period+"_CR_NonIso_b_RMTdPhi_"+region;
  histname[1][3] = path+"/TT_powheg_pythia8_"+period+"_CR_NonIso_b_RMTdPhi_"+region;
  histname[1][4] = path+"/WToLNu_"+period+"_CR_NonIso_b_RMTdPhi_"+region;
  histname[1][5] = path+"/ZToNuNu_"+period+"_CR_NonIso_b_RMTdPhi_"+region;
  histname[1][6] = path+"/Multiboson_"+period+"_CR_NonIso_b_RMTdPhi_"+region;
  histname[1][7] = path+"/GJets_"+period+"_CR_NonIso_b_RMTdPhi_"+region;
  histname[1][8] = path+"/DYToLL_"+period+"_CR_NonIso_b_RMTdPhi_"+region;
  cout << histname[0][0] << endl;

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
      name = "bkg_"+(TString)to_string(i)+(TString)to_string(j);
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
    CF[i]->Add(bkg[i][1], -1);
    for(int j=5; j<9; j++){
      CF[i]->Add(bkg[i][j], -1);
    }
  }

  TMatrix D(2,2);
  TMatrix W(2,2);
  TMatrix T(2,2);
  TMatrix D_up(2,2);
  TMatrix W_up(2,2);
  TMatrix T_up(2,2);
  TMatrix D_dw(2,2);
  TMatrix W_dw(2,2);
  TMatrix T_dw(2,2);

  double cf_W[binsize];
  double cf_W_up[binsize];
  double cf_W_dw[binsize];
  double cf_T[binsize];
  double cf_T_up[binsize];
  double cf_T_dw[binsize];
  double cf_t[binsize-1];
  double cf_t_up[binsize-1];
  double cf_t_dw[binsize-1];
  double binerr_x[binsize];
  double binx[binsize];

  int flag = 0;
  double tbinerr_x[binsize-1];
  double tbinx[binsize-1];


  tuple<vector<double>, vector<double>, vector<double>, vector<double>> Err_CFs;
  TH1D* temp;
  vector<double> temp1;
  vector<double> temp2;
  vector<double> temp3;
  vector<double> temp4;

  Err_CFs = CalcLErrors(CF, bkg, (string)iso+(string)obj+(string)period+(string)region);
  temp1 = get<0>(Err_CFs);
  temp2 = get<1>(Err_CFs);
  temp3 = get<2>(Err_CFs);
  temp4 = get<3>(Err_CFs);

  for(int i=1; i<= binsize;i++){
    cf_W_up[i-1] = abs(temp1.at(i-1));
    cf_T_up[i-1] = abs(temp2.at(i-1));
    cf_W_dw[i-1] = abs(temp3.at(i-1));
    cf_T_dw[i-1] = abs(temp4.at(i-1));

    for(int j=0; j<2;j++){
      D[j][0] = bkg[j][4]->GetBinContent(i);
      D[j][1] = bkg[j][2]->GetBinContent(i);
      W[j][0] = CF[j]->GetBinContent(i);
      W[j][1] = bkg[j][2]->GetBinContent(i);
      T[j][0] = bkg[j][4]->GetBinContent(i);
      T[j][1] = CF[j]->GetBinContent(i);
      binerr_x[i-1] = CF[j]->GetXaxis()->GetBinWidth(i)/2;
      binx[i-1] = CF[j]->GetXaxis()->GetBinCenter(i);
    }
    cf_W[i-1] = W.Determinant()/D.Determinant();
    cf_T[i-1] = T.Determinant()/D.Determinant();
    if(flag == 2 && i > 2) {
      cf_t[i-2] = T.Determinant()/D.Determinant();
      cf_t_up[i-2] = abs(temp2.at(i-1));
      cf_t_dw[i-2] = abs(temp2.at(i-1));
    }

    if(cf_W[i-1] + temp3.at(i-1) < 0) cf_W_dw[i-1] = cf_W[i-1];
    if(cf_T[i-1] + temp4.at(i-1) < 0) cf_T_dw[i-1] = cf_T[i-1];

    if(cf_T[i-1] < 0) {
      cout << "CF is smaller than 0" << endl;
      if(TString(obj).Contains("NJet")) {
        for(int j=0; j<2;j++) {
          D[j][0] = bkg[j][4]->GetBinContent(i) + bkg[j][4]->GetBinContent(i+1);
          D[j][1] = bkg[j][2]->GetBinContent(i) + bkg[j][2]->GetBinContent(i+1);
          W[j][0] = CF[j]->GetBinContent(i) + CF[j]->GetBinContent(i+1);
          W[j][1] = bkg[j][2]->GetBinContent(i) + bkg[j][2]->GetBinContent(i+1);
          T[j][0] = bkg[j][4]->GetBinContent(i) + bkg[j][4]->GetBinContent(i+1);
          T[j][1] = CF[j]->GetBinContent(i) + CF[j]->GetBinContent(i+1);
        }
        tbinerr_x[0] = (CF[0]->GetXaxis()->GetBinWidth(1) + CF[0]->GetXaxis()->GetBinWidth(2))/2.;
        tbinx[0] = (CF[0]->GetXaxis()->GetBinCenter(1) - CF[0]->GetXaxis()->GetBinWidth(1)/2 + CF[0]->GetXaxis()->GetBinCenter(2) + CF[0]->GetXaxis()->GetBinWidth(2)/2 )/2.;
        tbinerr_x[1] = CF[0]->GetXaxis()->GetBinWidth(3)/2;
        tbinerr_x[2] = CF[0]->GetXaxis()->GetBinWidth(4)/2;
        tbinerr_x[3] = CF[0]->GetXaxis()->GetBinWidth(5)/2;
        tbinx[1] = CF[0]->GetXaxis()->GetBinCenter(3);
        tbinx[2] = CF[0]->GetXaxis()->GetBinCenter(4);
        tbinx[3] = CF[0]->GetXaxis()->GetBinCenter(5);
        cf_t[i-1] = T.Determinant()/D.Determinant();
        cf_t_up[i-1] = abs(temp2.at(1));
        cf_t_dw[i-1] = abs(temp4.at(1));
        flag = 2;
      } else {
        for(int j=0; j<2;j++) {
          D[j][0] = bkg[j][4]->GetBinContent(i) + bkg[j][4]->GetBinContent(i-1);
          D[j][1] = bkg[j][2]->GetBinContent(i) + bkg[j][2]->GetBinContent(i-1);
          W[j][0] = CF[j]->GetBinContent(i) + CF[j]->GetBinContent(i-1);
          W[j][1] = bkg[j][2]->GetBinContent(i) + bkg[j][2]->GetBinContent(i-1);
          T[j][0] = bkg[j][4]->GetBinContent(i) + bkg[j][4]->GetBinContent(i-1);
          T[j][1] = CF[j]->GetBinContent(i) + CF[j]->GetBinContent(i-1);
        }
        tbinerr_x[0] = binerr_x[0];
        tbinerr_x[1] = binerr_x[1];
        tbinerr_x[2] = binerr_x[2];
        tbinerr_x[3] = binerr_x[3];
        tbinx[0] = binx[0];
        tbinx[1] = binx[1];
        tbinx[2] = binx[2];
        tbinx[3] = binx[3];
        tbinerr_x[4] = (CF[0]->GetXaxis()->GetBinWidth(i) + CF[0]->GetXaxis()->GetBinWidth(i-1))/2.;
        tbinx[4] = (CF[0]->GetXaxis()->GetBinCenter(i-1) - CF[0]->GetXaxis()->GetBinWidth(i-1)/2 + CF[0]->GetXaxis()->GetBinCenter(i) + CF[0]->GetXaxis()->GetBinWidth(i)/2 )/2.;
        cf_T[i-2] = T.Determinant()/D.Determinant();
        flag = 1;
      }
    }
  }
  //for(int i=0;i<binsize-1;i++) cout << i << ", " << tbinx[i] << ", " << tbinerr_x[i] << ", " << cf_t[i] << endl;

  auto CF_W = new TGraphAsymmErrors(binsize, binx, cf_W, binerr_x, binerr_x, cf_W_dw, cf_W_up);
  auto CF_T = new TGraphAsymmErrors(binsize, binx, cf_T, binerr_x, binerr_x, cf_T_dw, cf_T_up);
  if(flag == 1) CF_T = new TGraphAsymmErrors(binsize-1, tbinx, cf_T, tbinerr_x, tbinerr_x, cf_T_dw, cf_T_up);
  if(flag == 2) CF_T = new TGraphAsymmErrors(binsize-1, tbinx, cf_t, tbinerr_x, tbinerr_x, cf_t_dw, cf_T_up);

  string name_W = "CF_NonIso_W_"+(string)period+"_"+(string)region;
  string name_T = "CF_NonIso_T_"+(string)period+"_"+(string)region;
  CF_W->SetName(name_W.c_str());
  CF_T->SetName(name_T.c_str());

  TCanvas* c1 = new TCanvas("c1","",700,700);

  c1->SetGridx();
  c1->SetGridy();
  CF_W->SetMaximum(2);
  CF_W->SetMinimum(0);
  CF_W->GetYaxis()->SetTitle("Correction Factor");
  CF_W->GetXaxis()->SetTitle("M_{R} #times R^{2} [GeV]");
  if(TString(obj).Contains("NJet")) CF_W->GetXaxis()->SetTitle("N_{Jets}");
  CF_W->GetYaxis()->SetTitleOffset(0.9);
  CF_W->SetMarkerStyle(21);
  CF_T->SetMarkerStyle(21);
  CF_T->SetMarkerColor(kRed);
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(CF_W);
  mg->Add(CF_T);
  mg->Draw("AP");

  auto latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.03);
  latex->SetTextAlign(31); // align right
  latex->DrawLatex(0.90, 0.93, Form("138 fb^{-1} (#sqrt{s} = 13 TeV)"));
  latex->SetTextAlign(11); // align left
  latex->DrawLatex(0.10,0.93,"CMS Work in Progress");

  auto leg = new TLegend(0.45,0.77,0.9,0.9);
  if(TString(region).Contains("1Boost"))     leg->SetHeader("1 boost jet with Non-isolated lepton final state");
  else if(TString(region).Contains("2Boost")) leg->SetHeader("#geq 2 boost jet with Non-isolated lepton final state");
  leg->AddEntry(CF_W,  "W(#rightarrowl#nu)+jet CF", "p");
  leg->AddEntry(CF_T,  "Top(TT+ST) CF", "p");
  leg->Draw();

  c1->SaveAs("plot/"+period+obj+"_NonIso_CF_"+region+".png");
  return make_tuple(CF_W, CF_T);
}

void DrawCF(TString Title, TGraphAsymmErrors* g2, TGraphAsymmErrors* g3, TGraphAsymmErrors* g4, TGraphAsymmErrors* g5, TGraphAsymmErrors* g1){
  delete gROOT->FindObject("c1");
  TCanvas* c1 = new TCanvas("c1","",700,700);
  c1->SetGridx();
  c1->SetGridy();
  c1->GetFrame()->SetBorderSize(12);
  g1->SetMarkerStyle(21);
  g2->SetMarkerStyle(21);
  g3->SetMarkerStyle(21);
  g4->SetMarkerStyle(21);
  g5->SetMarkerStyle(21);
  g1->SetMarkerColor(kBlack);
  g2->SetMarkerColor(kRed);
  g3->SetMarkerColor(kBlue);
  g4->SetMarkerColor(kMagenta);
  g5->SetMarkerColor(kCyan);
  g1->SetLineColor(kBlack);
  g2->SetLineColor(kRed);
  g3->SetLineColor(kBlue);
  g4->SetLineColor(kMagenta);
  g5->SetLineColor(kCyan);
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(g1);
  mg->Add(g2);
  mg->Add(g3);
  mg->Add(g4);
  mg->Add(g5);
  if(TString(Title).Contains("MRR2")) mg->GetXaxis()->SetTitle("M_{R} #times R^{2}");
  if(TString(Title).Contains("NJet")) mg->GetXaxis()->SetTitle("N_{AK4 jet}");
  mg->Draw("AP");

  auto latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.03);
  latex->SetTextAlign(31); // align right
  latex->DrawLatex(0.90, 0.93, Form("138 fb^{-1} (#sqrt{s} = 13 TeV)"));
  latex->SetTextAlign(11); // align left
  latex->DrawLatex(0.10,0.93,"CMS Work in Progress");

  auto leg = new TLegend(0.55,0.57,0.9,0.9);
  leg->SetHeader(Title);
  leg->AddEntry(g1, "run2", "p");
  leg->AddEntry(g2, "2016", "p");
  leg->AddEntry(g3, "2016APV", "p");
  leg->AddEntry(g4, "2017", "p");
  leg->AddEntry(g5, "2018", "p");
  leg->Draw();

  Title.ReplaceAll(", ", "_");
  c1->SaveAs("plot/CF_"+Title+".png");
}

void DrawCFs(TGraphAsymmErrors* h1[7][10][2]){

  DrawCF("Q, MRR2, 1Boost", h1[0][0][0], h1[0][2][0], h1[0][4][0], h1[0][6][0], h1[0][8][0]);
  DrawCF("T, MRR2, 1Boost", h1[1][0][0], h1[1][2][0], h1[1][4][0], h1[1][6][0], h1[1][8][0]);
  DrawCF("W, MRR2, 1Boost", h1[2][0][0], h1[2][2][0], h1[2][4][0], h1[2][6][0], h1[2][8][0]);
  DrawCF("Q, MRR2, 2Boost", h1[0][1][0], h1[0][3][0], h1[0][5][0], h1[0][7][0], h1[0][9][0]);
  DrawCF("T, MRR2, 2Boost", h1[1][1][0], h1[1][3][0], h1[1][5][0], h1[1][7][0], h1[1][9][0]);
  DrawCF("W, MRR2, 2Boost", h1[2][1][0], h1[2][3][0], h1[2][5][0], h1[2][7][0], h1[2][9][0]);

  DrawCF("L, MRR2, 1Boost" , h1[3][0][0], h1[3][2][0], h1[3][4][0], h1[3][6][0], h1[3][8][0]);
  DrawCF("LT, MRR2, 1Boost", h1[4][0][0], h1[4][2][0], h1[4][4][0], h1[4][6][0], h1[4][9][0]);
  DrawCF("L, MRR2, 2Boost" , h1[3][1][0], h1[3][3][0], h1[3][5][0], h1[3][7][0], h1[3][8][0]);
  DrawCF("LT, MRR2, 2Boost", h1[4][1][0], h1[4][3][0], h1[4][5][0], h1[4][7][0], h1[4][9][0]);

  DrawCF("W, MRR2, NonIso, 1Boost", h1[5][0][0], h1[5][2][0], h1[5][4][0], h1[5][6][0], h1[5][8][0]);
  DrawCF("W, MRR2, NonIso, 2Boost", h1[5][1][0], h1[5][3][0], h1[5][5][0], h1[5][7][0], h1[5][9][0]);
  DrawCF("T, MRR2, NonIso, 1Boost", h1[6][0][0], h1[6][2][0], h1[6][4][0], h1[6][6][0], h1[6][8][0]);
  DrawCF("T, MRR2, NonIso, 2Boost", h1[6][1][0], h1[6][3][0], h1[6][5][0], h1[6][7][0], h1[6][9][0]);

  DrawCF("Q, NJet, 1Boost", h1[0][0][1], h1[0][2][1], h1[0][4][1], h1[0][6][1], h1[0][8][1]);
  DrawCF("T, NJet, 1Boost", h1[1][0][1], h1[1][2][1], h1[1][4][1], h1[1][6][1], h1[1][8][1]);
  DrawCF("W, NJet, 1Boost", h1[2][0][1], h1[2][2][1], h1[2][4][1], h1[2][6][1], h1[2][8][1]);
  DrawCF("Q, NJet, 2Boost", h1[0][1][1], h1[0][3][1], h1[0][5][1], h1[0][7][1], h1[0][9][1]);
  DrawCF("T, NJet, 2Boost", h1[1][1][1], h1[1][3][1], h1[1][5][1], h1[1][7][1], h1[1][9][1]);
  DrawCF("W, NJet, 2Boost", h1[2][1][1], h1[2][3][1], h1[2][5][1], h1[2][7][1], h1[2][9][1]);

  DrawCF("L, NJet, 1Boost" , h1[3][0][1], h1[3][2][1], h1[3][4][1], h1[3][6][1], h1[3][8][1]);
  DrawCF("LT, NJet, 1Boost", h1[4][0][1], h1[4][2][1], h1[4][4][1], h1[4][6][1], h1[4][9][1]);
  DrawCF("L, NJet, 2Boost" , h1[3][1][1], h1[3][3][1], h1[3][5][1], h1[3][7][1], h1[3][8][1]);
  DrawCF("LT, NJet, 2Boost", h1[4][1][1], h1[4][3][1], h1[4][5][1], h1[4][7][1], h1[4][9][1]);

  DrawCF("W, NJet, NonIso, 1Boost", h1[5][0][1], h1[5][2][1], h1[5][4][1], h1[5][6][1], h1[5][8][1]);
  DrawCF("W, NJet, NonIso, 2Boost", h1[5][1][1], h1[5][3][1], h1[5][5][1], h1[5][7][1], h1[5][9][1]);
  DrawCF("T, NJet, NonIso, 1Boost", h1[6][0][1], h1[6][2][1], h1[6][4][1], h1[6][6][1], h1[6][8][1]);
  DrawCF("T, NJet, NonIso, 2Boost", h1[6][1][1], h1[6][3][1], h1[6][5][1], h1[6][7][1], h1[6][9][1]);

}

void BkgCorr_solver_NF(){
  tuple<TGraphAsymmErrors*, TGraphAsymmErrors*, TGraphAsymmErrors*> CFs;
  tuple<TGraphAsymmErrors*, TGraphAsymmErrors*> LCFs;
  tuple<TGraphAsymmErrors*, TGraphAsymmErrors*> NonIsoCFs;
  TGraphAsymmErrors* h1[7][10][2];
  TFile* output;

  CFs = Correction("2016", "1Boost", "MRR2Bin", "230131/run_2023_02_08.root");
  h1[0][0][0] = get<0>(CFs); h1[1][0][0] = get<1>(CFs); h1[2][0][0] = get<2>(CFs);
  CFs = Correction("2016", "2Boost", "MRR2Bin", "230131/run_2023_02_08.root");
  h1[0][1][0] = get<0>(CFs); h1[1][1][0] = get<1>(CFs); h1[2][1][0] = get<2>(CFs);
  CFs = Correction("2016APV", "1Boost", "MRR2Bin", "230131/run_2023_02_08.root");
  h1[0][2][0] = get<0>(CFs); h1[1][2][0] = get<1>(CFs); h1[2][2][0] = get<2>(CFs);
  CFs = Correction("2016APV", "2Boost", "MRR2Bin", "230131/run_2023_02_08.root");
  h1[0][3][0] = get<0>(CFs); h1[1][3][0] = get<1>(CFs); h1[2][3][0] = get<2>(CFs);
  CFs = Correction("2017", "1Boost", "MRR2Bin", "230131/run_2023_02_08.root");
  h1[0][4][0] = get<0>(CFs); h1[1][4][0] = get<1>(CFs); h1[2][4][0] = get<2>(CFs);
  CFs = Correction("2017", "2Boost", "MRR2Bin", "230131/run_2023_02_08.root");
  h1[0][5][0] = get<0>(CFs); h1[1][5][0] = get<1>(CFs); h1[2][5][0] = get<2>(CFs);
  CFs = Correction("2018", "1Boost", "MRR2Bin", "230131/run_2023_02_08.root");
  h1[0][6][0] = get<0>(CFs); h1[1][6][0] = get<1>(CFs); h1[2][6][0] = get<2>(CFs);
  CFs = Correction("2018", "2Boost", "MRR2Bin", "230131/run_2023_02_08.root");
  h1[0][7][0] = get<0>(CFs); h1[1][7][0] = get<1>(CFs); h1[2][7][0] = get<2>(CFs);
  CFs = Correction("run2", "1Boost", "MRR2Bin", "230131/run_2023_02_08.root");
  h1[0][8][0] = get<0>(CFs); h1[1][8][0] = get<1>(CFs); h1[2][8][0] = get<2>(CFs);
  CFs = Correction("run2", "2Boost", "MRR2Bin", "230131/run_2023_02_08.root");
  h1[0][9][0] = get<0>(CFs); h1[1][9][0] = get<1>(CFs); h1[2][9][0] = get<2>(CFs);

  LCFs = LCorrection("2016", "1Boost", "MRR21vlBin","230131/run_2023_02_08.root");
  h1[3][0][0] = get<0>(LCFs); h1[4][0][0] = get<1>(LCFs);
  LCFs = LCorrection("2016", "2Boost", "MRR21vlBin","230131/run_2023_02_08.root");
  h1[3][1][0] = get<0>(LCFs); h1[4][1][0] = get<1>(LCFs);
  LCFs = LCorrection("2016APV", "1Boost", "MRR21vlBin","230131/run_2023_02_08.root");
  h1[3][2][0] = get<0>(LCFs); h1[4][2][0] = get<1>(LCFs);
  LCFs = LCorrection("2016APV", "2Boost", "MRR21vlBin","230131/run_2023_02_08.root");
  h1[3][3][0] = get<0>(LCFs); h1[4][3][0] = get<1>(LCFs);
  LCFs = LCorrection("2017", "1Boost", "MRR21vlBin","230131/run_2023_02_08.root");
  h1[3][4][0] = get<0>(LCFs); h1[4][4][0] = get<1>(LCFs);
  LCFs = LCorrection("2017", "2Boost", "MRR21vlBin","230131/run_2023_02_08.root");
  h1[3][5][0] = get<0>(LCFs); h1[4][5][0] = get<1>(LCFs);
  LCFs = LCorrection("2018", "1Boost", "MRR21vlBin","230131/run_2023_02_08.root");
  h1[3][6][0] = get<0>(LCFs); h1[4][6][0] = get<1>(LCFs);
  LCFs = LCorrection("2018", "2Boost", "MRR21vlBin","230131/run_2023_02_08.root");
  h1[3][7][0] = get<0>(LCFs); h1[4][7][0] = get<1>(LCFs);
  LCFs = LCorrection("run2", "1Boost", "MRR21vlBin","230131/run_2023_02_08.root");
  h1[3][8][0] = get<0>(LCFs); h1[4][8][0] = get<1>(LCFs);
  LCFs = LCorrection("run2", "2Boost", "MRR21vlBin","230131/run_2023_02_08.root");
  h1[3][9][0] = get<0>(LCFs); h1[4][9][0] = get<1>(LCFs);


  NonIsoCFs = NonIsoCorrection("2016", "1Boost", "MRR2Bin","230131/run_2023_02_11.root");
  h1[5][0][0] = get<0>(NonIsoCFs); h1[6][0][0] = get<1>(NonIsoCFs);
  NonIsoCFs = NonIsoCorrection("2016", "2Boost", "MRR2Bin","230131/run_2023_02_11.root");
  h1[5][1][0] = get<0>(NonIsoCFs); h1[6][1][0] = get<1>(NonIsoCFs);
  NonIsoCFs = NonIsoCorrection("2016APV", "1Boost", "MRR2Bin","230131/run_2023_02_11.root");
  h1[5][2][0] = get<0>(NonIsoCFs); h1[6][2][0] = get<1>(NonIsoCFs);
  NonIsoCFs = NonIsoCorrection("2016APV", "2Boost", "MRR2Bin","230131/run_2023_02_11.root");
  h1[5][3][0] = get<0>(NonIsoCFs); h1[6][3][0] = get<1>(NonIsoCFs);
  NonIsoCFs = NonIsoCorrection("2017", "1Boost", "MRR2Bin","230131/run_2023_02_11.root");
  h1[5][4][0] = get<0>(NonIsoCFs); h1[6][4][0] = get<1>(NonIsoCFs);
  NonIsoCFs = NonIsoCorrection("2017", "2Boost", "MRR2Bin","230131/run_2023_02_11.root");
  h1[5][5][0] = get<0>(NonIsoCFs); h1[6][5][0] = get<1>(NonIsoCFs);
  NonIsoCFs = NonIsoCorrection("2018", "1Boost", "MRR2Bin","230131/run_2023_02_11.root");
  h1[5][6][0] = get<0>(NonIsoCFs); h1[6][6][0] = get<1>(NonIsoCFs);
  NonIsoCFs = NonIsoCorrection("2018", "2Boost", "MRR2Bin","230131/run_2023_02_11.root");
  h1[5][7][0] = get<0>(NonIsoCFs); h1[6][7][0] = get<1>(NonIsoCFs);
  NonIsoCFs = NonIsoCorrection("run2", "1Boost", "MRR2Bin","230131/run_2023_02_11.root");
  h1[5][8][0] = get<0>(NonIsoCFs); h1[6][8][0] = get<1>(NonIsoCFs);
  NonIsoCFs = NonIsoCorrection("run2", "2Boost", "MRR2Bin","230131/run_2023_02_11.root");
  h1[5][9][0] = get<0>(NonIsoCFs); h1[6][9][0] = get<1>(NonIsoCFs);

  output = new TFile("CFs.root", "recreate");
  for(int i=0;i<10;i++) {
    //if(i==2 || i==3) continue;
    for(int j=0;j<7;j++) h1[j][i][0]->Write();
  }

  CFs = Correction("2016", "1Boost", "NJetBins", "230131/run_2023_02_10.root");
  h1[0][0][1] = get<0>(CFs); h1[1][0][1] = get<1>(CFs); h1[2][0][1] = get<2>(CFs);
  CFs = Correction("2016", "2Boost", "NJetBins", "230131/run_2023_02_10.root");
  h1[0][1][1] = get<0>(CFs); h1[1][1][1] = get<1>(CFs); h1[2][1][1] = get<2>(CFs);
  CFs = Correction("2016APV", "1Boost", "NJetBins", "230131/run_2023_02_10.root");
  h1[0][2][1] = get<0>(CFs); h1[1][2][1] = get<1>(CFs); h1[2][2][1] = get<2>(CFs);
  CFs = Correction("2016APV", "2Boost", "NJetBins", "230131/run_2023_02_10.root");
  h1[0][3][1] = get<0>(CFs); h1[1][3][1] = get<1>(CFs); h1[2][3][1] = get<2>(CFs);
  CFs = Correction("2017", "1Boost", "NJetBins", "230131/run_2023_02_10.root");
  h1[0][4][1] = get<0>(CFs); h1[1][4][1] = get<1>(CFs); h1[2][4][1] = get<2>(CFs);
  CFs = Correction("2017", "2Boost", "NJetBins", "230131/run_2023_02_10.root");
  h1[0][5][1] = get<0>(CFs); h1[1][5][1] = get<1>(CFs); h1[2][5][1] = get<2>(CFs);
  CFs = Correction("2018", "1Boost", "NJetBins", "230131/run_2023_02_10.root");
  h1[0][6][1] = get<0>(CFs); h1[1][6][1] = get<1>(CFs); h1[2][6][1] = get<2>(CFs);
  CFs = Correction("2018", "2Boost", "NJetBins", "230131/run_2023_02_10.root");
  h1[0][7][1] = get<0>(CFs); h1[1][7][1] = get<1>(CFs); h1[2][7][1] = get<2>(CFs);
  CFs = Correction("run2", "1Boost", "NJetBins", "230131/run_2023_02_10.root");
  h1[0][8][1] = get<0>(CFs); h1[1][8][1] = get<1>(CFs); h1[2][8][1] = get<2>(CFs);
  CFs = Correction("run2", "2Boost", "NJetBins", "230131/run_2023_02_10.root");
  h1[0][9][1] = get<0>(CFs); h1[1][9][1] = get<1>(CFs); h1[2][9][1] = get<2>(CFs);

  LCFs = LCorrection("2016", "1Boost", "NJetBins","230131/run_2023_02_10.root");
  h1[3][0][1] = get<0>(LCFs); h1[4][0][1] = get<1>(LCFs);
  LCFs = LCorrection("2016", "2Boost", "NJetBins","230131/run_2023_02_10.root");
  h1[3][1][1] = get<0>(LCFs); h1[4][1][1] = get<1>(LCFs);
  LCFs = LCorrection("2016APV", "1Boost", "NJetBins","230131/run_2023_02_10.root");
  h1[3][2][1] = get<0>(LCFs); h1[4][2][1] = get<1>(LCFs);
  LCFs = LCorrection("2016APV", "2Boost", "NJetBins","230131/run_2023_02_10.root");
  h1[3][3][1] = get<0>(LCFs); h1[4][3][1] = get<1>(LCFs);
  LCFs = LCorrection("2017", "1Boost", "NJetBins","230131/run_2023_02_10.root");
  h1[3][4][1] = get<0>(LCFs); h1[4][4][1] = get<1>(LCFs);
  LCFs = LCorrection("2017", "2Boost", "NJetBins","230131/run_2023_02_10.root");
  h1[3][5][1] = get<0>(LCFs); h1[4][5][1] = get<1>(LCFs);
  LCFs = LCorrection("2018", "1Boost", "NJetBins","230131/run_2023_02_10.root");
  h1[3][6][1] = get<0>(LCFs); h1[4][6][1] = get<1>(LCFs);
  LCFs = LCorrection("2018", "2Boost", "NJetBins","230131/run_2023_02_10.root");
  h1[3][7][1] = get<0>(LCFs); h1[4][7][1] = get<1>(LCFs);
  LCFs = LCorrection("run2", "1Boost", "NJetBins","230131/run_2023_02_10.root");
  h1[3][8][1] = get<0>(LCFs); h1[4][8][1] = get<1>(LCFs);
  LCFs = LCorrection("run2", "2Boost", "NJetBins","230131/run_2023_02_10.root");
  h1[3][9][1] = get<0>(LCFs); h1[4][9][1] = get<1>(LCFs);

  NonIsoCFs = NonIsoCorrection("2016", "1Boost", "NJetBins","230131/run_2023_02_12.root");
  h1[5][0][1] = get<0>(NonIsoCFs); h1[6][0][1] = get<1>(NonIsoCFs);
  NonIsoCFs = NonIsoCorrection("2016", "2Boost", "NJetBins","230131/run_2023_02_12.root");
  h1[5][1][1] = get<0>(NonIsoCFs); h1[6][1][1] = get<1>(NonIsoCFs);
  NonIsoCFs = NonIsoCorrection("2016APV", "1Boost", "NJetBins","230131/run_2023_02_12.root");
  h1[5][2][1] = get<0>(NonIsoCFs); h1[6][2][1] = get<1>(NonIsoCFs);
  NonIsoCFs = NonIsoCorrection("2016APV", "2Boost", "NJetBins","230131/run_2023_02_12.root");
  h1[5][3][1] = get<0>(NonIsoCFs); h1[6][3][1] = get<1>(NonIsoCFs);
  NonIsoCFs = NonIsoCorrection("2017", "1Boost", "NJetBins","230131/run_2023_02_12.root");
  h1[5][4][1] = get<0>(NonIsoCFs); h1[6][4][1] = get<1>(NonIsoCFs);
  NonIsoCFs = NonIsoCorrection("2017", "2Boost", "NJetBins","230131/run_2023_02_12.root");
  h1[5][5][1] = get<0>(NonIsoCFs); h1[6][5][1] = get<1>(NonIsoCFs);
  NonIsoCFs = NonIsoCorrection("2018", "1Boost", "NJetBins","230131/run_2023_02_12.root");
  h1[5][6][1] = get<0>(NonIsoCFs); h1[6][6][1] = get<1>(NonIsoCFs);
  NonIsoCFs = NonIsoCorrection("2018", "2Boost", "NJetBins","230131/run_2023_02_12.root");
  h1[5][7][1] = get<0>(NonIsoCFs); h1[6][7][1] = get<1>(NonIsoCFs);
  NonIsoCFs = NonIsoCorrection("run2", "1Boost", "NJetBins","230131/run_2023_02_12.root");
  h1[5][8][1] = get<0>(NonIsoCFs); h1[6][8][1] = get<1>(NonIsoCFs);
  NonIsoCFs = NonIsoCorrection("run2", "2Boost", "NJetBins","230131/run_2023_02_12.root");
  h1[5][9][1] = get<0>(NonIsoCFs); h1[6][9][1] = get<1>(NonIsoCFs);


  output = new TFile("NJet_CFs.root", "recreate");
  for(int i=0;i<10;i++) {
    //if(i==2 || i==3) continue;
    for(int j=0;j<7;j++) h1[j][i][1]->Write();
  }

  DrawCFs(h1);

}
