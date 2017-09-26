#include <iostream>
#include <string>
#include <sstream>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TColor.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"

using namespace std;

Double_t* getVariableBinEdges(int num_entries, Double_t* tmp_array)
{
  Double_t* my_array = new Double_t[num_entries];
  for (int i = 0; i != num_entries; ++i) {
    my_array[i] = tmp_array[i];
 
    //cout << "bin edge " << i << " : " << my_array[i] << endl;
  }
   return my_array;
}

void fillefferr2D(TH2D* h, TH2D* hp, TH2D* heu, TH2D* hel)
{
  double total, passed, cen, eup, elow;
  for (int i=1; i<h->GetXaxis()->GetNbins()+1; i++) {
    for (int j=1; j<h->GetXaxis()->GetNbins()+1; j++) {
      passed = hp->GetBinContent(i, j);
      total = h->GetBinContent(i, j);
      if (total == 0) continue;
      cen = (double)passed / (double)total;
      eup = TEfficiency::ClopperPearson(total, passed, 0.687, true);
      elow = TEfficiency::ClopperPearson(total, passed, 0.687, false);
      heu->SetBinContent(i, j, eup-cen);
      hel->SetBinContent(i, j, cen-elow);
    }
  }
}

TGraphAsymmErrors* filleff(TH1D *h, TH1D *hp) 
{
  TGraphAsymmErrors *he = new TGraphAsymmErrors(hp, h, "cl=0.683 b(1,1) mode");
  int neffpoints = he->GetN();
  //int neffpoints = he->GetXaxis()->GetNbins();
  //cout << neffpoints << endl;
  //cout << i << "  np: " << neffpointsMR << endl;
  for (int p=0; p<neffpoints; p++) {
    he->SetPointEXhigh(p, 0);
    he->SetPointEXlow(p, 0);
  }
  return he;
}

TH1D* ChangeBinning(TH2D* h){
  TH1::SetDefaultSumw2();
  TH1D* htemp = new TH1D("htemp","",70, 0, 70);
  const char* label[71] = {"[400, 500]", "[500, 600]","[600, 700]","[700, 750]","[750, 800]","[800, 850]","[850, 900]","[900, 950]","[950, 1000]","[1000, 1500]","[400, 500]", "[500, 600]","[600, 700]","[700, 750]","[750, 800]","[800, 850]","[850, 900]","[900, 950]","[950, 1000]","[1000, 1500]","[400, 500]", "[500, 600]","[600, 700]","[700, 750]","[750, 800]","[800, 850]","[850, 900]","[900, 950]","[950, 1000]","[1000, 1500]","[400, 500]", "[500, 600]","[600, 700]","[700, 750]","[750, 800]","[800, 850]","[850, 900]","[900, 950]","[950, 1000]","[1000, 1500]","[400, 500]", "[500, 600]","[600, 700]","[700, 750]","[750, 800]","[800, 850]","[850, 900]","[900, 950]","[950, 1000]","[1000, 1500]","[400, 500]", "[500, 600]","[600, 700]","[700, 750]","[750, 800]","[800, 850]","[850, 900]","[900, 950]","[950, 1000]","[1000, 1500]","[400, 500]", "[500, 600]","[600, 700]","[700, 750]","[750, 800]","[800, 850]","[850, 900]","[900, 950]","[950, 1000]","[1000, 1500]"}; 
  
  int count=1;
  //cout << h->GetNbinsX() << ", " << h->GetNbinsY() << endl;
  for(int j=1;j<=h->GetNbinsY();j++){
    for(int i=1;i<=h->GetNbinsX();i++){
      htemp->SetBinContent(count,h->GetBinContent(i,j));
      htemp->GetXaxis()->SetBinLabel( count, label[count-1]);
      count++;
    }
  }
  return htemp;
}

void calceff(TString ntuplePtl = "SingleElectron")
{

  TH1::SetDefaultSumw2();
  TFile* file1 = TFile::Open("/gatbawi/palgongsan/chuh/susy/susy170801/"+ntuplePtl+".root");

  TH1D *h_HT_pre  = (TH1D*)file1->Get("h_HT_pre");
  TH1D *h_HT_pre_pass = (TH1D*)file1->Get("h_HT_pre_pass");

  TH2D *h_HT_j1pt_pre = (TH2D*)file1->Get("h_HT_j1pt_pre");
  //TH2D *h_HT_j1pt_pre_eff_up = (TH2D*)file1->Get("h_HT_j1pt_pre");
  //TH2D *h_HT_j1pt_pre_eff_low = (TH2D*)file1->Get("h_HT_j1pt_pre");
  TH2D *h_HT_j1pt_pre_pass = (TH2D*)file1->Get("h_HT_j1pt_pre_pass");
  TH1D *h_HT_j1pt_pre_1D = ChangeBinning(h_HT_j1pt_pre);
  TH1D *h_HT_j1pt_pre_pass_1D = ChangeBinning(h_HT_j1pt_pre_pass);
  //cout << h_HT_j1pt_pre_1D->GetXaxis()->GetNbins() << endl;

  TGraphAsymmErrors* h_HT_pre_eff_HT = filleff(h_HT_pre, h_HT_pre_pass);
  h_HT_pre_eff_HT->SetName("h_HT_pre_eff_HT");

  TGraphAsymmErrors* h_HT_j1pt_pre_eff_HT = filleff(h_HT_j1pt_pre_1D, h_HT_j1pt_pre_pass_1D);
  h_HT_j1pt_pre_eff_HT->SetName("h_HT_j1pt_pre_eff_HT");
  //cout << h_HT_j1pt_pre_eff_HT->GetXaxis()->GetNbins() << endl;

  // HT-j1pt
  TH2D* h_HT_j1pt_pre_eff = (TH2D*)h_HT_j1pt_pre_pass->Clone();
  h_HT_j1pt_pre_eff->Sumw2();
  h_HT_j1pt_pre_eff->SetName("h_HT_j1pt_pre_eff");
  h_HT_j1pt_pre_eff->Divide(h_HT_j1pt_pre);
  //fillefferr2D(h_HT_j1pt_pre, h_HT_j1pt_pre_pass, h_HT_j1pt_pre_eff_up, h_HT_j1pt_pre_eff_low);

  TFile* output = new TFile(ntuplePtl+"_TriggerEfficiency.root", "recreate");
  h_HT_pre_eff_HT->Write();
  h_HT_j1pt_pre_eff_HT->Write();
  h_HT_j1pt_pre_eff->Write();

  TCanvas* c1 = new TCanvas("ctemp","",750,750);
  h_HT_j1pt_pre_eff->Draw("colz");
  c1->SaveAs("fig/HLTeff_"+ntuplePtl+".png");

}

void eff(){
  calceff("SingleElectron");
  calceff("SingleMuon");
  calceff("MET");
}
