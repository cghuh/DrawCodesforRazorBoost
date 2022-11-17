void save(TString name="QCD"){
  TFile* file1 = TFile::Open("2016_"+name+".root");
  TFile* file2 = TFile::Open("2016APV_"+name+".root");
  TFile* file3 = TFile::Open("2017_"+name+".root");
  TFile* file4 = TFile::Open("2018_"+name+".root");


  TH2D* h2016[6];
  TH2D* h2016APV[6];
  TH2D* h2017[6];
  TH2D* h2018[6];

  h2016[0] = (TH2D*)file1->Get("btag_eff_b_loose");
  h2016[1] = (TH2D*)file1->Get("btag_eff_c_loose");
  h2016[2] = (TH2D*)file1->Get("btag_eff_l_loose");
  h2016[3] = (TH2D*)file1->Get("btag_eff_b_medium");
  h2016[4] = (TH2D*)file1->Get("btag_eff_c_medium");
  h2016[5] = (TH2D*)file1->Get("btag_eff_l_medium");
  h2016APV[0] = (TH2D*)file1->Get("btag_eff_b_loose");
  h2016APV[1] = (TH2D*)file1->Get("btag_eff_c_loose");
  h2016APV[2] = (TH2D*)file1->Get("btag_eff_l_loose");
  h2016APV[3] = (TH2D*)file1->Get("btag_eff_b_medium");
  h2016APV[4] = (TH2D*)file1->Get("btag_eff_c_medium");
  h2016APV[5] = (TH2D*)file1->Get("btag_eff_l_medium");
  h2017[0] = (TH2D*)file2->Get("btag_eff_b_loose");
  h2017[1] = (TH2D*)file2->Get("btag_eff_c_loose");
  h2017[2] = (TH2D*)file2->Get("btag_eff_l_loose");
  h2017[3] = (TH2D*)file2->Get("btag_eff_b_medium");
  h2017[4] = (TH2D*)file2->Get("btag_eff_c_medium");
  h2017[5] = (TH2D*)file2->Get("btag_eff_l_medium");
  h2018[0] = (TH2D*)file3->Get("btag_eff_b_loose");
  h2018[1] = (TH2D*)file3->Get("btag_eff_c_loose");
  h2018[2] = (TH2D*)file3->Get("btag_eff_l_loose");
  h2018[3] = (TH2D*)file3->Get("btag_eff_b_medium");
  h2018[4] = (TH2D*)file3->Get("btag_eff_c_medium");
  h2018[5] = (TH2D*)file3->Get("btag_eff_l_medium");

  h2016[0]->SetName("btag_eff_b_loose_2016");
  h2016[1]->SetName("btag_eff_c_loose_2016");
  h2016[2]->SetName("btag_eff_l_loose_2016");
  h2016[3]->SetName("btag_eff_b_medium_2016");
  h2016[4]->SetName("btag_eff_c_medium_2016");
  h2016[5]->SetName("btag_eff_l_medium_2016");
  h2016APV[0]->SetName("btag_eff_b_loose_2016APV");
  h2016APV[1]->SetName("btag_eff_c_loose_2016APV");
  h2016APV[2]->SetName("btag_eff_l_loose_2016APV");
  h2016APV[3]->SetName("btag_eff_b_medium_2016APV");
  h2016APV[4]->SetName("btag_eff_c_medium_2016APV");
  h2016APV[5]->SetName("btag_eff_l_medium_2016APV");
  h2017[0]->SetName("btag_eff_b_loose_2017");
  h2017[1]->SetName("btag_eff_c_loose_2017");
  h2017[2]->SetName("btag_eff_l_loose_2017");
  h2017[3]->SetName("btag_eff_b_medium_2017");
  h2017[4]->SetName("btag_eff_c_medium_2017");
  h2017[5]->SetName("btag_eff_l_medium_2017");
  h2018[0]->SetName("btag_eff_b_loose_2018");
  h2018[1]->SetName("btag_eff_c_loose_2018");
  h2018[2]->SetName("btag_eff_l_loose_2018");
  h2018[3]->SetName("btag_eff_b_medium_2018");
  h2018[4]->SetName("btag_eff_c_medium_2018");
  h2018[5]->SetName("btag_eff_l_medium_2018");

  TFile* output = new TFile(name+".root", "recreate");
  h2016[0]->Write();
  h2016[1]->Write();
  h2016[2]->Write();
  h2016[3]->Write();
  h2016[4]->Write();
  h2016[5]->Write();
  h2016APV[0]->Write();
  h2016APV[1]->Write();
  h2016APV[2]->Write();
  h2016APV[3]->Write();
  h2016APV[4]->Write();
  h2016APV[5]->Write();
  h2017[0]->Write();
  h2017[1]->Write();
  h2017[2]->Write();
  h2017[3]->Write();
  h2017[4]->Write();
  h2017[5]->Write();
  h2018[0]->Write();
  h2018[1]->Write();
  h2018[2]->Write();
  h2018[3]->Write();
  h2018[4]->Write();
  h2018[5]->Write();

}

void btag(){
   save("QCD");
   save("ST");
   save("TT");
   save("WJetsToLNu");
   //save("SMS-T5ttcc");
}
