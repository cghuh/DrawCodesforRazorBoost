void CalcComparison(TString period="2016", TString region="Val_Signal"){
  cout << period << ", " << region << endl;
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TString locate = "/Users/huhchanggi/temp/221110/run_2022_11_21";
  //TString locate = "/Users/huhchanggi/temp/211021/run_2021_10_25";
  TFile* file1 = TFile::Open(locate+".root");

  TString histname[9];
  TString path = "/Counts_vs_MRR2/Syst_vs_MRR2/";
  //TString path = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/";
  histname[0] = path+"DYToLL_"+period+"_"+region;
  histname[1] = path+"GJets_"+period+"_"+region;
  histname[2] = path+"Higgs_"+period+"_"+region;
  histname[3] = path+"Multiboson_"+period+"_"+region;
  histname[4] = path+"Multijet_"+period+"_"+region;
  histname[5] = path+"Top_"+period+"_"+region;
  histname[6] = path+"TT_powheg_pythia8_"+period+"_"+region;
  histname[7] = path+"WToLNu_"+period+"_"+region;
  histname[8] = path+"ZToNuNu_"+period+"_"+region;
  //histname[9] = path+"Data_"+period+"_"+region;

  TH2D* bkg2D[9];
  for(int i=0; i<9;i++) bkg2D[i] = (TH2D*)file1->Get(histname[i]);

  const int binsizex = bkg2D[0]->GetNbinsX();
  const int binsizey = bkg2D[0]->GetNbinsY();
  TH1D* nominal[binsizex];
  TH1D* upvar[binsizex];
  TH1D* dwvar[binsizex];
  TH1D *temp, *temp1, *temp2;
  double var_nom, var_up, var_dw;
  double syst_up, syst_dw;
  double test1, test2, test3;
  
/*
  for(int i=1;i<=binsizex;i++){
    for(int j=0;j<9;j++) nominal[j] = bkg2D[j]->GetBinContent(i,1);
  }
*/
  for(int j=0;j<9;j++){
    nominal[j] = bkg2D[j]->ProjectionX(Form("bin%d",j),1,1);
    if(j==0) temp = (TH1D*)nominal[0]->Clone();
    else    temp->Add(nominal[j]);
  }
  var_nom = temp->Integral();

  for(int i=2;i<=binsizey;i++){
    for(int j=0;j<9;j++){
      upvar[j] = bkg2D[j]->ProjectionX(Form("binup%d",j),i,i);
      dwvar[j] = bkg2D[j]->ProjectionX(Form("bindw%d",j),i+1,i+1);
      if(j==0) {
        temp1 = (TH1D*)upvar[j]->Clone();
        temp2 = (TH1D*)dwvar[j]->Clone();
      } else  {
        temp1->Add(upvar[j]);
        temp2->Add(dwvar[j]);
      }

/*
      test1 = nominal[j]->Integral();
      test2 = upvar[j]->Integral();
      test3 = dwvar[j]->Integral();
      test2 = test2-test1;
      test3 = test3-test1;
      cout << histname[j] << endl;
      if(test2 > 0) {
        if(test3 > 0) cout << "Total number of events with "<< i/2 << " th systematic source : " << test1 << " +" << abs(test2) << ", +" << abs(test3) << endl;
        else            cout << "Total number of events with "<< i/2 << " th systematic source : " << test1 << " +" << abs(test2) << ", -" << abs(test3) << endl;
      } else {
        if(test3 > 0) cout << "Total number of events with "<< i/2 << " th systematic source : " << test1 << " +" << abs(test3) << ", -" << abs(test2) << endl;
        else            cout << "Total number of events with "<< i/2 << " th systematic source : " << test1 << " -" << abs(test2) << ", -" << abs(test3) << endl;
      }
*/

      //cout << "total " << j << ", " << temp1->Integral() << ", " << temp2->Integral() << endl;
      //for(int ii=1;ii<=binsizex;ii++) cout << ii << ", " << upvar[j]->GetBinContent(ii) << endl;
      //for(int ii=1;ii<=binsizex;ii++) cout << ii << ", " << dwvar[j]->GetBinContent(ii) << endl;
      //if(j==2 || j==3){
      //cout << "up variation" << endl;
      //for(int ii=1;ii<=binsizex;ii++){if(isnan(upvar[j]->GetBinContent(ii))) cout << ii << endl;}
      //cout << "down variation" << endl;
      //for(int ii=1;ii<=binsizex;ii++){if(isnan(dwvar[j]->GetBinContent(ii))) cout << ii << endl;}
      //}
      //for(int ii=1;ii<=binsizex;ii++){if(isnan(upvar[j]->GetBinContent(ii))||isnan(dwvar[j]->GetBinContent(ii))) cout << j << endl;}
    }
    var_up = temp1->Integral();
    var_dw = temp2->Integral();
    
    syst_up = var_up-var_nom;
    syst_dw = var_dw-var_nom;
    if(syst_up > 0) {
      if(syst_dw > 0) cout << "Total number of events with "<< i/2 << " th systematic source : " << var_nom << " +" << abs(syst_up) << ", +" << abs(syst_dw) << endl;
      else            cout << "Total number of events with "<< i/2 << " th systematic source : " << var_nom << " +" << abs(syst_up) << ", -" << abs(syst_dw) << endl;
    } else {
      if(syst_dw > 0) cout << "Total number of events with "<< i/2 << " th systematic source : " << var_nom << " +" << abs(syst_dw) << ", -" << abs(syst_up) << endl;
      else            cout << "Total number of events with "<< i/2 << " th systematic source : " << var_nom << " -" << abs(syst_up) << ", -" << abs(syst_dw) << endl;
    }
    //cout << "Total number of events : " << var_nom << " +- " << var_up-var_nom << ", " << var_dw-var_nom << endl;
    i++;
  }
}
void SystComparison(){
  CalcComparison("2016", "Val_QCD");
  CalcComparison("2017", "Val_QCD");
  CalcComparison("2018", "Val_QCD");
  //CalcComparison("2016", "Pre");
  //CalcComparison("2016", "CR_Top16_V");
  //CalcComparison("2017", "CR_NonIso_0b_RMTdPhi_1Boost");
  //CalcComparison("2017", "CR_NonIso_b_RMTdPhi_1Boost");
  //CalcComparison("2018", "CR_NonIso_RMT");
  //CalcComparison("BlindData_2016", "SR_Had_V_b_45j");
  //CalcComparison("BlindData_2017", "SR_Had_V_b_45j");
  //CalcComparison("BlindData_2018", "SR_Had_V_b_45j");
  //CalcComparison("2017", "CR_NonIso_RdPhi");
  //CalcComparison("2016", "CR_QCD16_V");
  //CalcComparison("2016", "CR_Top16_Top");

  //CalcComparison("2016", "");
  //CalcComparison("2017", "");
  //CalcComparison("2018", "");
}
