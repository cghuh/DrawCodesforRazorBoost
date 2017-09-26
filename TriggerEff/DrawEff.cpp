void canvas(TString hname, TCanvas* ctemp ){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  cout << hname << endl;
  TFile* f1 = TFile::Open("SingleMuon_TriggerEfficiency.root");
  TFile* f2 = TFile::Open("SingleElectron_TriggerEfficiency.root");
  TFile* f3 = TFile::Open("MET_TriggerEfficiency.root");
  
  TGraphAsymmErrors* h_HT_pre_eff_SM = (TGraphAsymmErrors*)f1->Get(hname);
  TGraphAsymmErrors* h_HT_pre_eff_SE = (TGraphAsymmErrors*)f2->Get(hname);
  TGraphAsymmErrors* h_HT_pre_eff_SL = (TGraphAsymmErrors*)f3->Get(hname);
  
  ctemp = new TCanvas("ctemp","",750,750);
  ctemp->SetGridx();
  ctemp->SetGridy();
  h_HT_pre_eff_SE->GetYaxis()->SetRangeUser(0,1);

  if(hname.Contains("h_HT"))h_HT_pre_eff_SE->GetXaxis()->SetTitle("H_{T} [GeV]");

  h_HT_pre_eff_SE->GetYaxis()->SetTitle("#epsilon_{HLT_PFHT800 or HLT_PFHT900 or HLT_AK8PFJet450}");
  //if(hname.Contains("pre_eff_HT")) h_HT_pre_eff_SE->GetXaxis()->SetRangeUser(0,2000);

  h_HT_pre_eff_SM->SetLineColor(kBlack);
  h_HT_pre_eff_SM->SetMarkerColor(kBlack);
  h_HT_pre_eff_SM->SetMarkerStyle(20);
  h_HT_pre_eff_SE->SetLineColor(kRed);
  h_HT_pre_eff_SE->SetMarkerColor(kRed);
  h_HT_pre_eff_SE->GetXaxis()->SetRangeUser(0,70);
  h_HT_pre_eff_SE->SetMarkerStyle(20);
  h_HT_pre_eff_SL->SetLineColor(kBlue);
  h_HT_pre_eff_SL->SetMarkerColor(kBlue);
  h_HT_pre_eff_SL->SetMarkerStyle(20);

  h_HT_pre_eff_SE->Draw("AP");
  h_HT_pre_eff_SM->Draw("Psame");
  h_HT_pre_eff_SL->Draw("Psame");

  leg = new TLegend(0.50,0.15,0.9,0.3);
  leg->SetHeader("Pass trigger pre-selection");
  leg->AddEntry(h_HT_pre_eff_SM, "SingleMuon (IsoMu24, 1 TightMuon)", "lp");
  leg->AddEntry(h_HT_pre_eff_SE, "SingleElectron (Ele23+Ele27, 1 TightEle)", "lp");
  leg->AddEntry(h_HT_pre_eff_SL, "MET (PFMET120, no Lepton)", "lp");
  leg->Draw();

  TText *t = new TText();
  if(hname == "h_HT_j1pt_pre_eff_HT"){
    //const char* label[100] = {"[400, 500]", "[500, 600]","[600, 700]","[700, 750]","[750, 800]","[800, 850]","[850, 900]","[900, 950]","[950, 1000]","[1000, 1500]","[400, 500]", "[500, 600]","[600, 700]","[700, 750]","[750, 800]","[800, 850]","[850, 900]","[900, 950]","[950, 1000]","[1000, 1500]","[400, 500]", "[500, 600]","[600, 700]","[700, 750]","[750, 800]","[800, 850]","[850, 900]","[900, 950]","[950, 1000]","[1000, 1500]","[400, 500]", "[500, 600]","[600, 700]","[700, 750]","[750, 800]","[800, 850]","[850, 900]","[900, 950]","[950, 1000]","[1000, 1500]","[400, 500]", "[500, 600]","[600, 700]","[700, 750]","[750, 800]","[800, 850]","[850, 900]","[900, 950]","[950, 1000]","[1000, 1500]","[400, 500]", "[500, 600]","[600, 700]","[700, 750]","[750, 800]","[800, 850]","[850, 900]","[900, 950]","[950, 1000]","[1000, 1500]","[400, 500]", "[500, 600]","[600, 700]","[700, 750]","[750, 800]","[800, 850]","[850, 900]","[900, 950]","[950, 1000]","[1000, 1500]"};
    //cout << h_HT_pre_eff_SE->GetXaxis()->GetNbins() << endl;
    //for(int i=1;i<=h_HT_pre_eff_SE->GetXaxis()->GetNbins();i++){
    //  h_HT_pre_eff_SE->GetXaxis()->SetBinLabel(i, label[i-1]);
    //}
  //h_HT_pre_eff_SE->GetXaxis()->SetLabelSize(0.03);
  //h_HT_pre_eff_SE->GetXaxis()->SetNdivisions(502);
  //h_HT_pre_eff_SE->GetXaxis()->LabelsOption("v");
  //t->SetTextAlign(33);
  //t->SetTextSize(0.015);
  //t->SetTextFont(42);
  //for (Int_t i=0;i<h_HT_pre_eff_SE->GetXaxis()->GetNbins();i++) {
  //    t->DrawText(h_HT_pre_eff_SE->GetXaxis()->GetBinCenter(i),0,label[i]);
  // }
    float xmax = h_HT_pre_eff_SE->GetXaxis()->GetBinUpEdge(h_HT_pre_eff_SE->GetXaxis()->GetLast());
    //string text = "#scale[0.41]{[200,300]  [300, 400]  [400, 450]  [450, 500]  [500, 550]  [550, 600]  [600, 1000] }";
    string text = "#scale[0.41]{    [200,300]  [300, 400]  [400, 450]  [450, 500]  [500, 550]  [550, 600]  [600, 1000]}";
    TLatex* era_lat = new TLatex(xmax,0.5, text.c_str());
    era_lat->SetTextAlign(32);
    era_lat->SetTextSize(0.06);
    era_lat->SetTextFont(42);
    era_lat->SetLineWidth(2);
    era_lat->Draw();
  }
  ctemp->SaveAs("fig/HLTeff_"+hname+".png");
  TH1F* h1 = new TH1F("h1", ";H_{T} [GeV]",70,0,70);
  for(int i=1;i<=70;i++){
   h1->GetBinContent(i, h_HT_pre_eff_SM->GetErrorY(i)-h_HT_pre_eff_SE->GetErrorY(i));
  }
  h1->Draw();
  ctemp->SaveAs("fig/HLTeff_diff_"+hname+".png");

}

void DrawEff(){

  TCanvas* c1;
  TString hname;

  hname = "h_HT_pre_eff_HT";
  canvas(hname, c1);
  hname = "h_HT_j1pt_pre_eff_HT";
  canvas(hname, c1);
}
