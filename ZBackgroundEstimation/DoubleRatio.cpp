void CalcComparison(TString period="2016"){
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TString locate = "/Users/huhchanggi/temp/220512/run_2022_05_16";
  TFile* file1 = TFile::Open(locate+".root");

  TString histname[2][10];
  TString path = "/Counts_vs_MRR2Bin/Syst_vs_MRR2Bin/";
  //TString path = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/";
  histname[0][0] = path+"DYToLL_"+period+"_CR_1PhoInv";
  histname[0][1] = path+"GJets_"+period+"_CR_1PhoInv";
  histname[0][2] = path+"Higgs_"+period+"_CR_1PhoInv";
  histname[0][3] = path+"Multiboson_"+period+"_CR_1PhoInv";
  histname[0][4] = path+"Multijet_"+period+"_CR_1PhoInv";
  histname[0][5] = path+"Top_"+period+"_CR_1PhoInv";
  histname[0][6] = path+"TT_powheg_pythia8_"+period+"_CR_1PhoInv";
  histname[0][7] = path+"WToLNu_"+period+"_CR_1PhoInv";
  histname[0][8] = path+"ZToNuNu_"+period+"_CR_1PhoInv";
  histname[0][9] = path+"Data_"+period+"_CR_1PhoInv";
  histname[1][0] = path+"DYToLL_"+period+"_CR_2LepInv";
  histname[1][1] = path+"GJets_"+period+"_CR_2LepInv";
  histname[1][2] = path+"Higgs_"+period+"_CR_2LepInv";
  histname[1][3] = path+"Multiboson_"+period+"_CR_2LepInv";
  histname[1][4] = path+"Multijet_"+period+"_CR_2LepInv";
  histname[1][5] = path+"Top_"+period+"_CR_2LepInv";
  histname[1][6] = path+"TT_powheg_pythia8_"+period+"_CR_2LepInv";
  histname[1][7] = path+"WToLNu_"+period+"_CR_2LepInv";
  histname[1][8] = path+"ZToNuNu_"+period+"_CR_2LepInv";
  histname[1][9] = path+"Data_"+period+"_CR_2LepInv";

  TH2D* bkg2D[2][10];
  TH1D* htemp[2][10];
  for(int j=0;j<2;j++){
    for(int i=0; i<10;i++) {
      //cout << histname[j][i] << endl;
      bkg2D[j][i] = (TH2D*)file1->Get(histname[j][i]);
      if(bkg2D[j][i] == NULL) {
        htemp[j][i] = new TH1D(Form("bin%d%d",i,j), "", 1,0,1000);
        htemp[j][i]->SetBinContent(1,0);
        htemp[j][i]->SetBinError(1,0);
        continue;
      }
      htemp[j][i] = (TH1D*)bkg2D[j][i]->ProjectionX(Form("bin%d%d",i,j),1,1);
      htemp[j][i]->Rebin(6);
      if(htemp[j][i]->GetNbinsX() != 1) cout << "# of bin isn't 6" << endl;
    }
  }
  htemp[0][9]->Add(htemp[0][0],-1);
  htemp[0][9]->Add(htemp[0][2],-1);
  htemp[0][9]->Add(htemp[0][3],-1);
  htemp[0][9]->Add(htemp[0][4],-1);
  htemp[0][9]->Add(htemp[0][5],-1);
  htemp[0][9]->Add(htemp[0][6],-1);
  htemp[0][9]->Add(htemp[0][7],-1);
  htemp[0][9]->Add(htemp[0][8],-1);
  htemp[0][9]->Divide(htemp[0][1]);

  htemp[1][9]->Add(htemp[1][1],-1);
  htemp[1][9]->Add(htemp[1][2],-1);
  htemp[1][9]->Add(htemp[1][3],-1);
  htemp[1][9]->Add(htemp[1][4],-1);
  htemp[1][9]->Add(htemp[1][5],-1);
  htemp[1][9]->Add(htemp[1][6],-1);
  htemp[1][9]->Add(htemp[1][7],-1);
  htemp[1][9]->Add(htemp[1][8],-1);
  htemp[1][9]->Divide(htemp[1][0]);

  TH1D* hRatio = (TH1D*)htemp[0][9]->Clone();
  hRatio->Divide(htemp[1][9]);

  cout << period << "Double ratio : " << hRatio->GetBinContent(1) << ", Error : " << hRatio->GetBinError(1) << endl;

}
void DoubleRatio(){
  CalcComparison("2016");
  CalcComparison("2017");
  CalcComparison("2018");
}
