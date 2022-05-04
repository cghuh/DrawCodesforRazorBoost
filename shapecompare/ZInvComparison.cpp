void DrawComparison(TString period="2016", TString region="Signal"){
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TString locate = "/Users/huhchanggi/temp/210604/run_2021_06_01";
  TFile* file1 = TFile::Open(locate+".root");
  locate = "/Users/huhchanggi/temp/210901/run_2021_08_30";
  TFile* file2 = TFile::Open(locate+".root");

  TString histname[3];
  TString path = "/Counts_vs_MRR2Bins/Syst_vs_MRR2Bins/";
  histname[0] = path+"ZToNuNu_"+period+"_Val_"+region;
  histname[1] = path+"ZToNuNu_"+period+"_Val_"+region;
  histname[2] = path+"ZToNuNu_"+period+"_Val_"+region+"_L";

  TH2D* bkg2D[3];
  bkg2D[0] = (TH2D*)file1->Get(histname[0]);
  for(int j=1; j<=2; j++) bkg2D[j] = (TH2D*)file2->Get(histname[j]);

  int binsize = bkg2D[0]->GetNbinsX();
  double bin[binsize+1];
  string name;

  for(int i=1;i<=bkg2D[0]->GetNbinsX()+1;i++) bin[i-1] = bkg2D[0]->GetXaxis()->GetBinLowEdge(i);

  TH1D* bkg[3];
  for(int i=0; i<3; i++) {
    name = "bkg_"+to_string(i);
    bkg[i] = new TH1D(name.c_str(), "", binsize, bin);
    for(int k=1;k<=bkg2D[0]->GetNbinsX();k++) {
      if(bkg2D[i] == NULL) {
        bkg[i]->SetBinContent(k, 0);
        bkg[i]->SetBinError(k, 0);
      }
      else {
        bkg[i]->SetBinContent(k, bkg2D[i]->GetBinContent(k,1));
        bkg[i]->SetBinError(k,   bkg2D[i]->GetBinError(k,1));
        //cout << i << ", bin : " << k << ", Error : " << bkg2D[i]->GetBinError(k,1) << endl;
      }
    }
  }

  TCanvas* c1 = new TCanvas("c1", "", 1000, 1000);
  TPad* pad1 = new TPad("pad1", "The pad 1", 0.,0.2,1.,1.0,0.);
  TPad* pad2 = new TPad("pad2", "The pad 2", 0.,0.0,1.,0.278,10);
  pad1->Draw();
  pad2->Draw();

  pad1->cd();
  //pad1->SetGridx();
  //pad1->SetGridy();
  pad1->SetLogy();
  bkg[0]->SetMarkerColor(kRed);
  bkg[0]->SetLineColor(kRed);
  bkg[1]->SetMarkerColor(kBlack);
  bkg[1]->SetLineColor(kBlack);
  bkg[1]->SetMarkerStyle(20);
  bkg[2]->SetMarkerColor(kGreen);
  bkg[2]->SetLineColor(kGreen);
  bkg[2]->SetMarkerStyle(20);
  bkg[0]->GetYaxis()->SetRangeUser(1,1e4);

  float xmin = bkg[0]->GetXaxis()->GetBinLowEdge(bkg[0]->GetXaxis()->GetFirst());
  float xmax = bkg[0]->GetXaxis()->GetBinUpEdge(bkg[0]->GetXaxis()->GetLast());

  bkg[0]->Draw("HIST");

  string text = "CMS #scale[0.7]{#font[52]{Work in progress}}";
  TLatex* cms_lat = new TLatex(xmin, 1.5e4, text.c_str());
  cms_lat->SetLineWidth(2);
  cms_lat->Draw();
  text = "#scale[0.7]{137 fb^{-1} (13 TeV)}";
  TLatex* era_lat = new TLatex(xmax*0.97,1.5e4, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.04);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();
  bkg[1]->Draw("EPsame");
  bkg[2]->Draw("EPsame");

  //pad1->Update();
  char legentry[100];
   
  auto leg = new TLegend(0.45,0.65,0.90,0.90);
  //leg->SetBorderSize(0);
  leg->SetHeader(period+", "+region+" Validation region");
  leg->SetNColumns(2);
  sprintf(legentry, "%.1f",bkg[0]->Integral());
  leg->AddEntry(bkg[0],"Z(vv) MC","L");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",bkg[1]->Integral());
  leg->AddEntry(bkg[1],"#gamma estimate","L P E");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",bkg[2]->Integral());
  leg->AddEntry(bkg[2],"W(lv) estimate","L P E");
  leg->AddEntry((TObject*)0,legentry,"");
  leg->Draw();
 
  TH1D* compare1 = (TH1D*)bkg[1]->Clone(); 
  TH1D* compare2 = (TH1D*)bkg[2]->Clone(); 
  compare1->Divide(bkg[0]);
  compare2->Divide(bkg[0]);
  //for(int i=1;i<=compare1->GetNbinsX();i++) cout << i << ", " << compare1->GetBinContent(i) << endl;
  pad2->cd();
  compare1->SetStats(0);
  compare2->SetStats(0);
  //pad2->SetGridx();
  //pad2->SetGridy();
  //compare1->GetXaxis()->SetTitle("M_{R} X R^{2}");
  //compare1->GetXaxis()->SetTitleOffset(0.5);
  //compare1->GetXaxis()->SetTitleSize(0.08);
  compare1->GetXaxis()->SetLabelSize(0.06);
  compare1->GetYaxis()->SetLabelSize(0.09);
  compare1->GetYaxis()->SetTitle("Estimate/MC");
  compare1->GetYaxis()->SetTitleSize(0.1);
  compare1->GetYaxis()->SetTitleOffset(0.5);
  compare1->GetYaxis()->SetRangeUser(0.5,1.5);
  compare1->Draw("EP");
  compare2->Draw("EPsame");

  c1->SaveAs("fig/"+period+"_compare_ZinvEst_"+region+"Val.png");
}
void ZInvComparison(){
  DrawComparison("2016", "Signal");
  DrawComparison("2016", "QCD");
  DrawComparison("2017", "Signal");
  DrawComparison("2017", "QCD");
  DrawComparison("2018", "Signal");
  DrawComparison("2018", "QCD");
}
