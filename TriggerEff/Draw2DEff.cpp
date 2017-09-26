void canvas(TString hname, TCanvas* ctemp ){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  cout << hname << endl;
  TFile* f1 = TFile::Open("SingleMuon_TriggerEfficiency.root");
  TFile* f2 = TFile::Open("SingleElectron_TriggerEfficiency.root");
  TFile* f3 = TFile::Open("MET_TriggerEfficiency.root");
  
  TH2D* h_HT_pre_eff_SM = (TH2D*)f1->Get(hname);
  TH2D* h_HT_pre_eff_SE = (TH2D*)f2->Get(hname);
  TH2D* h_HT_pre_eff_SL = (TH2D*)f3->Get(hname);
  

  TH2D* htemp1 = (TH2D*)f1->Get(hname);
  TH2D* htemp2 = (TH2D*)f2->Get(hname);
  TH2D* htemp3 = (TH2D*)f3->Get(hname);
  htemp1->SetName("HLTeff_diff_ele_muon_"+hname);
  htemp2->SetName("HLTeff_diff_muon_MET"+hname);
  htemp3->SetName("HLTeff_diff_ele_MET_"+hname);
  TH2D* htemp4 = (TH2D*)f1->Get(hname);
  TH2D* htemp5 = (TH2D*)f2->Get(hname);
  TH2D* htemp6 = (TH2D*)f3->Get(hname);
  htemp4->SetName("HLTeff_diff_ratio_ele_muon_"+hname);
  htemp5->SetName("HLTeff_diff_ratio_muon_MET"+hname);
  htemp6->SetName("HLTeff_diff_ratio_ele_MET_"+hname);
  TH2D* htemp7 = (TH2D*)f1->Get(hname);
  TH2D* htemp8 = (TH2D*)f2->Get(hname);
  TH2D* htemp9 = (TH2D*)f3->Get(hname);
  htemp7->SetName("HLTeff_ratio_ele_muon_"+hname);
  htemp8->SetName("HLTeff_ratio_muon_MET"+hname);
  htemp9->SetName("HLTeff_ratio_ele_MET_"+hname);

  Double_t SM_value, SE_value, MET_value;
  Double_t SM_Error, SE_Error, MET_Error;

  cout << h_HT_pre_eff_SM->GetXaxis()->GetNbins() << ", " << h_HT_pre_eff_SM->GetYaxis()->GetNbins() << endl;

  for(int i=1;i<=h_HT_pre_eff_SM->GetXaxis()->GetNbins();i++){
    cout << i << endl;
    for(int j=1;j<=h_HT_pre_eff_SM->GetYaxis()->GetNbins();j++){
      SM_value = h_HT_pre_eff_SM->GetBinContent(i,j);
      SE_value = h_HT_pre_eff_SE->GetBinContent(i,j);
      MET_value = h_HT_pre_eff_SL->GetBinContent(i,j);
      SM_Error = h_HT_pre_eff_SM->GetBinError(i,j);
      SE_Error = h_HT_pre_eff_SE->GetBinError(i,j);
      MET_Error = h_HT_pre_eff_SL->GetBinError(i,j);

      if(SM_value == 0){
        htemp1->SetBinContent(i,j,0);
        htemp2->SetBinContent(i,j,0);
        htemp4->SetBinContent(i,j,0);
        htemp5->SetBinContent(i,j,0);
        htemp6->SetBinContent(i,j,0);
        htemp7->SetBinContent(i,j,0);
        htemp3->SetBinContent(i,j,0);
        htemp8->SetBinContent(i,j,0);
        htemp9->SetBinContent(i,j,0);
      }

      cout << " SE  = " << SE_value;
      cout << ", SM  = " << SM_value;
      cout << ", MET = " << MET_value << endl;

      cout << " SE  = " << SE_Error;
      cout << ", SM  = " << SM_Error;
      cout << ", MET = " << MET_Error << endl;

      cout << " SE - MET = " << SE_value-MET_value;
      cout << ", SM - MET = " << SM_value-MET_value;
      cout << ", SE - SM  = " << SE_value-SM_value << endl;

      cout << " SE - MET = " <<  (SE_Error-MET_Error)/(SE_value-MET_value);
      cout << ", SM - MET = " << (SM_Error-MET_Error)/(SM_value-MET_value);
      cout << ", SE - SM  = " << (SE_Error-SM_Error)/(SE_value-SM_value) << endl;

      htemp3->SetBinContent(i,j,SE_Error-MET_Error);
      htemp8->SetBinContent(i,j,(SE_Error-MET_Error)/(SE_value*MET_value));
      htemp9->SetBinContent(i,j,(abs(SE_value-MET_value))/(SE_Error+MET_Error)*2);
      if(SM_value == 0) continue;
      htemp2->SetBinContent(i,j,SM_Error-MET_Error);
      htemp6->SetBinContent(i,j,(SM_Error-MET_Error)/(SM_value*MET_value));
      htemp7->SetBinContent(i,j,(abs(SM_value-MET_value))/(SM_Error+MET_Error)*2);
      htemp1->SetBinContent(i,j,SE_Error-SM_Error);
      htemp4->SetBinContent(i,j,(SE_Error-SM_Error)/(SE_value*SM_value));
      htemp5->SetBinContent(i,j,(abs(SE_value-SM_value))/(SE_Error+SM_Error)*2);
    }
  }   
  //float xmin = htemp1->GetXaxis()->GetBinLowEdge(htemp1->GetXaxis()->GetFirst());
  string text = "#scale[0.7]{Error_{#epsilon_{e}}-Error_{#epsilon_{#mu}}}";
  //string text = "#scale[0.7]{Error_{#epsilon_{e}}-Error_{#epsilon_{MET}}}";
  //string text = "#scale[0.7]{Error_{#epsilon_{#mu}}-Error_{#epsilon_{MET}}}";
  TLatex* era_lat = new TLatex(700,1050, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.06);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);

  TCanvas* ctemp1 = new TCanvas("ctemp1","",750,750);
  ctemp1->SetGridx();
  ctemp1->SetGridy();
  htemp1->Draw("colz");
  era_lat->Draw();
  ctemp1->SaveAs("fig/HLTeff_diff_ele_muon_"+hname+".png");
  htemp4->Draw("colz");
  text = "#scale[0.7]{#frac{Error_{#epsilon_{e}}-Error_{#epsilon_{#mu}}}{#epsilon_{e}#epsilon_{#mu}}}";
  era_lat = new TLatex(700,1050, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.06);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();
  ctemp1->SaveAs("fig/HLTeff_ratio_ele_muon_"+hname+".png");
  htemp5->Draw("colz");
  text = "#scale[0.7]{#frac{2|#epsilon_{e}-#epsilon_{#mu}|}{Error_{#epsilon_{e}}+Error_{#epsilon_{#mu}}}}";
  era_lat = new TLatex(700,1050, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.06);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();
  ctemp1->SaveAs("fig/HLTeff_diff_ratio_ele_muon_"+hname+".png");
  htemp2->Draw("colz");
  text = "#scale[0.7]{Error_{#epsilon_{#mu}}-Error_{#epsilon_{MET}}}";
  era_lat = new TLatex(700,1050, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.06);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();
  ctemp1->SaveAs("fig/HLTeff_diff_muon_MET_"+hname+".png");
  htemp6->Draw("colz");
  text = "#scale[0.7]{#frac{Error_{#epsilon_{#mu}}-Error_{#epsilon_{MET}}}{#epsilon_{#mu}#epsilon_{MET}}}";
  era_lat = new TLatex(700,1050, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.06);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();
  ctemp1->SaveAs("fig/HLTeff_ratio_muon_MET_"+hname+".png");
  htemp7->Draw("colz");
  text = "#scale[0.7]{#frac{2|#epsilon_{#mu}-#epsilon_{MET}|}{Error_{#epsilon_{#mu}}+Error_{#epsilon_{MET}}}}";
  era_lat = new TLatex(700,1050, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.06);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();
  ctemp1->SaveAs("fig/HLTeff_diff_ratio_muon_MET_"+hname+".png");
  htemp3->Draw("colz");
  text = "#scale[0.7]{Error_{#epsilon_{e}}-Error_{#epsilon_{MET}}}";
  era_lat = new TLatex(700,1050, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.06);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();
  ctemp1->SaveAs("fig/HLTeff_diff_ele_MET_"+hname+".png");
  htemp8->Draw("colz");
  text = "#scale[0.7]{#frac{Error_{#epsilon_{e}}-Error_{#epsilon_{MET}}}{#epsilon_{e}#epsilon_{MET}}}";
  era_lat = new TLatex(700,1050, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.06);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();
  ctemp1->SaveAs("fig/HLTeff_ratio_ele_MET_"+hname+".png");
  htemp9->Draw("colz");
  text = "#scale[0.7]{#frac{2|#epsilon_{e}-#epsilon_{MET}|}{Error_{#epsilon_{e}}+Error_{#epsilon_{MET}}}}";
  era_lat = new TLatex(700,1050, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.06);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();
  ctemp1->SaveAs("fig/HLTeff_diff_ratio_ele_MET_"+hname+".png");

}

void Draw2DEff(){

  TCanvas* c1;
  TString hname;

  hname = "h_HT_j1pt_pre_eff";
  canvas(hname, c1);
}
