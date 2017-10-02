void DrawComparison(TString NtupleName="TTX"){
  //gStyle->SetOptStat(0);
  gStyle->SetOptStat("nemr");
  gStyle->SetOptTitle(0);
  TString locate = "/gatbawi/palgongsan/chuh/susy/susy170830/";
  TFile* file1 = TFile::Open(locate+NtupleName+".root");
  //TFile* file2 = TFile::Open("ssb_calEle.root");

  TString name;
  if(NtupleName=="TT") name="T";
  if(NtupleName=="WJet")name="W";
  if(NtupleName=="Multijet") name="Q";

  TH1F* varR2S = (TH1F*)file1->Get("R2_S");
  varR2S->SetName("S Region");
  //varR2S->Rebin(20);
  varR2S->Scale(1./varR2S->Integral());
  TH1F* varR2C;
  if(NtupleName=="TT")  varR2C = (TH1F*)file1->Get("R2_T");
  if(NtupleName=="WJet") varR2C = (TH1F*)file1->Get("R2_W");
  if(NtupleName=="Multijet")  varR2C = (TH1F*)file1->Get("R2_Q");
  varR2C->SetName(name+" Region");
  //varR2C->Rebin(20);
  varR2C->Scale(1./varR2C->Integral());

  TCanvas* c1 = new TCanvas("c1", "", 700, 700);
  TPad* pad1 = new TPad("pad1", "The pad 1", 0.,0.2,1.,1.0,0.);
  TPad* pad2 = new TPad("pad2", "The pad 2", 0.,0.0,1.,0.278,10);
  pad1->Draw();
  pad2->Draw();

  pad1->cd();
  //pad1->SetGridx();
  //pad1->SetGridy();
  //pad1->SetLogy();
  varR2S->SetMarkerColor(kBlack);
  varR2S->GetYaxis()->SetRangeUser(0,0.8);
  varR2S->SetLineColor(kBlack);
  varR2S->SetLineWidth(2);
  varR2S->GetXaxis()->SetLabelSize(0);
  varR2S->GetYaxis()->SetTitle("A.U.");
  varR2S->SetMarkerStyle(20);
  varR2S->Draw("HIST");
  varR2C->SetLineColor(kRed);
  varR2C->SetLineWidth(2);
  varR2C->Draw("HISTsames");

  pad1->Update();
  TPaveStats *st = (TPaveStats*)varR2C->FindObject("stats");
  st->SetY1NDC(0.7);
  st->SetY2NDC(0.54);
   
  leg = new TLegend(0.55,0.78,0.77,0.88);
  leg->SetBorderSize(0);
  leg->AddEntry(varR2S, "S region", "l");
  leg->AddEntry(varR2C, name+" region", "l");
  leg->Draw();
 
  TH1F* comparePt = (TH1F*)varR2S->Clone(); 
  comparePt->Divide(varR2C);
  pad2->cd();
  comparePt->SetStats(0);
  pad2->SetGridx();
  pad2->SetGridy();
  comparePt->SetXTitle("R2");
  comparePt->GetYaxis()->SetRangeUser(0,2);
  comparePt->GetXaxis()->SetTitleOffset(0.5);
  comparePt->GetXaxis()->SetTitleSize(0.08);
  comparePt->GetXaxis()->SetLabelSize(0.06);
  comparePt->GetYaxis()->SetLabelSize(0.09);
  comparePt->GetYaxis()->SetTitleSize(0.);
  comparePt->Draw("EP");

  TH1F* varMRS = (TH1F*)file1->Get("MR_S");
  varMRS->SetName("S Region");
  //varMRS->Rebin(10);
  varMRS->Scale(1./varMRS->Integral());
  TH1F* varMRC;
  if(NtupleName=="TT")  varMRC = (TH1F*)file1->Get("MR_T");
  if(NtupleName=="WJet") varMRC = (TH1F*)file1->Get("MR_W");
  if(NtupleName=="Multijet")  varMRC = (TH1F*)file1->Get("MR_Q");
  varMRC->SetName(name+" Region");
  //varMRC->Rebin(10);
  varMRC->Scale(1./varMRC->Integral());

  TCanvas* c2 = new TCanvas("c2", "", 700, 700);
  pad1 = new TPad("pad1", "The pad 1", 0.,0.2,1.,1.0,0.);
  pad2 = new TPad("pad2", "The pad 2", 0.,0.0,1.,0.278,10);
  pad1->Draw();
  pad2->Draw();

  pad1->cd();
  //pad1->SetGridx();
  //pad1->SetGridy();
  //pad1->SetLogy();
  varMRS->SetMarkerColor(kBlack);
  varMRS->GetYaxis()->SetRangeUser(0,0.8);
  varMRS->SetLineColor(kBlack);
  varMRS->SetLineWidth(2);
  varMRS->GetXaxis()->SetLabelSize(0);
  varMRS->GetYaxis()->SetTitle("A.U.");
  varMRS->SetMarkerStyle(20);
  varMRS->Draw("HIST");
  varMRC->SetLineColor(kRed);
  varMRC->SetLineWidth(2);
  varMRC->Draw("HISTsames");

  pad1->Update();
  st = (TPaveStats*)varMRC->FindObject("stats");
  st->SetY1NDC(0.7);
  st->SetY2NDC(0.54);
   
  leg = new TLegend(0.55,0.78,0.77,0.88);
  leg->SetBorderSize(0);
  leg->AddEntry(varMRS, "S region", "l");
  leg->AddEntry(varMRC, name+" region", "l");
  leg->Draw();
 
  comparePt = (TH1F*)varMRS->Clone(); 
  comparePt->Divide(varMRC);
  pad2->cd();
  comparePt->SetStats(0);
  pad2->SetGridx();
  pad2->SetGridy();
  comparePt->SetXTitle("MR");
  comparePt->GetYaxis()->SetRangeUser(0,2);
  comparePt->GetXaxis()->SetTitleOffset(0.5);
  comparePt->GetXaxis()->SetTitleSize(0.08);
  comparePt->GetXaxis()->SetLabelSize(0.06);
  comparePt->GetYaxis()->SetLabelSize(0.09);
  comparePt->GetYaxis()->SetTitleSize(0.);
  comparePt->Draw("EP");

  c1->SaveAs("fig/R2Shape_Svs"+name+".png");
  c2->SaveAs("fig/MRShape_Svs"+name+".png");
}
void ShapeComparison(){
  DrawComparison("TT");
  DrawComparison("WJet");
  DrawComparison("Multijet");
}
