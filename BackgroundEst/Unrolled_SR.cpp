#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TCanvas.h"
#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include "TString.h"
#include "TROOT.h"
#include "TStyle.h"
#include <TChain.h>
#include <math.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>
#include <THStack.h>
using namespace std;

void UnrolledPlots(){

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();

  delete gROOT->FindObject("c1");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();
  TString dir = "/Users/huhchanggi/temp/220512/run_2022_05_16.root";
  TFile* file0 = TFile::Open(dir);

  TString path = "/Counts_vs_MRR2/Syst_vs_MRR2/";

  TString histname[29];
	histname[0]  = "_SR_Had_1htop";
	histname[1]  = "_SR_Had_2htop";
	histname[2]  = "_SR_Had_V_b_45j";
	histname[3]  = "_SR_Had_V_b_6j";
	histname[4]  = "_SR_Had_1V_0b_34j";
	histname[5]  = "_SR_Had_1V_0b_5j";
	histname[6]  = "_SR_Had_2V_0b_24j";
	histname[7]  = "_SR_Had_2V_0b_5j";
	histname[8]  = "_SR_Had_H_b_45j";
	histname[9]  = "_SR_Had_H_b_6j";
	histname[10] = "_SR_Had_2H_b_6j";
	histname[11] = "_SR_Had_HV_b_6j";
	histname[12] = "_SR_Had_1H_0b_34j";
	histname[13] = "_SR_Had_1H_0b_5j";
	histname[14] = "_SR_Had_2H_0b_34j";
	histname[15] = "_SR_Had_2H_0b_5j";
	histname[16] = "_SR_Had_HV_0b_24j";
	histname[17] = "_SR_Had_HV_0b_5j";
	histname[18] = "_SR_Lep_1htop";
	histname[19] = "_SR_Lep_V_b";
	histname[20] = "_SR_Lep_V_0b";
	histname[21] = "_SR_Lep_H_b";
	histname[22] = "_SR_Lep_H_0b";
	histname[23] = "_SR_Leptop_0htop";
	histname[24] = "_SR_Leptop_1htop";
	histname[25] = "_SR_Lepjet_0V_24j";
	histname[26] = "_SR_Lepjet_0V_5j";
	histname[27] = "_SR_Lepjet_1V_24j";
	histname[28] = "_SR_Lepjet_1V_5j";

  TString period[3];
  period[0] = "_BlindData_2016";
  period[1] = "_BlindData_2017";
  period[2] = "_BlindData_2018";

  TH2D* hempty = new TH2D("empty", "", 6,0,6,83,0,83);

  TH2D* hist[10];
  int NBins_had = 64;
  int NBins_H   = 96;
  int NBins_lep = 72;
  int NBins = NBins_had + NBins_H + NBins_lep;
  NBins = 174;
  hist[0] = new TH2D("multijet", ";Bin;variation", 174, 0, 174, 83, 0, 83);
  for(int i=1;i<10;i++) hist[i] = (TH2D*)hist[0]->Clone(); 
  TH2D* hhtemp[12][3];
  TH2D* htemp[10];

  for(int i=0;i<29;i++){
    for(int j=0;j<3;j++){
      hhtemp[0][j] = (TH2D*)file0->Get(path+"Multijet"+period[j]+histname[i]);
      hhtemp[1][j] = (TH2D*)file0->Get(path+"TT_powheg_pythia8"+period[j]+histname[i]);
      hhtemp[2][j] = (TH2D*)file0->Get(path+"WToLNu"+period[j]+histname[i]);
      hhtemp[3][j] = (TH2D*)file0->Get(path+"ZToNuNu"+period[j]+histname[i]);
      hhtemp[4][j] = (TH2D*)file0->Get(path+"Multiboson"+period[j]+histname[i]);
      hhtemp[5][j] = (TH2D*)file0->Get(path+"DYToLL"+period[j]+histname[i]);
      hhtemp[6][j] = (TH2D*)file0->Get(path+"GJets"+period[j]+histname[i]);
      hhtemp[7][j] = (TH2D*)file0->Get(path+"Top"+period[j]+histname[i]);
      hhtemp[8][j] = (TH2D*)file0->Get(path+"Higgs"+period[j]+histname[i]);
      hhtemp[9][j] = (TH2D*)file0->Get(path+"T5ttcc"+period[j]+histname[i]);
      hhtemp[10][j] = (TH2D*)file0->Get(path+"T5qqqqVV"+period[j]+histname[i]);
      hhtemp[11][j] = (TH2D*)file0->Get(path+"T2bW"+period[j]+histname[i]);
    }
    for(int j=0;j<7;j++){
      if(hhtemp[j][0] != NULL) {
        htemp[j] = (TH2D*)hhtemp[j][0]->Clone();
        if(hhtemp[j][1] != NULL) htemp[j]->Add(hhtemp[j][1]);
        if(hhtemp[j][2] != NULL) htemp[j]->Add(hhtemp[j][2]);
      }
      else if(hhtemp[j][1] != NULL) {
        htemp[j] = (TH2D*)hhtemp[j][1]->Clone();
        if(hhtemp[j][2] != NULL) htemp[j]->Add(hhtemp[j][2]);
      }
      else if(hhtemp[j][2] != NULL) htemp[j] = (TH2D*)hhtemp[j][2]->Clone();
      else htemp[j] = (TH2D*)hempty->Clone();
    }
    for(int j=9;j<12;j++){
      if(hhtemp[j][0] != NULL) {
        htemp[j-2] = (TH2D*)hhtemp[j][0]->Clone();
        if(hhtemp[j][1] != NULL) htemp[j-2]->Add(hhtemp[j][1]);
        if(hhtemp[j][2] != NULL) htemp[j-2]->Add(hhtemp[j][2]);
      }
      else if(hhtemp[j][1] != NULL) {
        htemp[j-2] = (TH2D*)hhtemp[j][1]->Clone();
        if(hhtemp[j][2] != NULL) htemp[j-2]->Add(hhtemp[j][2]);
      }
      else if(hhtemp[j][2] != NULL) htemp[j-2] = (TH2D*)hhtemp[j][2]->Clone();
      else htemp[j-2] = (TH2D*)hempty->Clone();
    }
    htemp[1]->Add(hhtemp[7][0]);
    htemp[1]->Add(hhtemp[7][1]);
    htemp[1]->Add(hhtemp[7][2]);
    htemp[4]->Add(hhtemp[8][0]);
    htemp[4]->Add(hhtemp[8][1]);
    htemp[4]->Add(hhtemp[8][2]);
    
    for(int x=1;x<=6;x++){
      for(int y=1;y<=83;y++){
        for(int z=0;z<10;z++){
          hist[z]->SetBinContent(i*6+x, y, htemp[z]->GetBinContent(x,y));
          hist[z]->SetBinError(i*6+x, y, htemp[z]->GetBinError(x,y));
        }
      }
    }
  }

  TH1D* histo[10];
  histo[0] = new TH1D("hmultijet", ";Bin;Counts", 174, 0, 174);
  for(int i=1;i<10;i++) histo[i]=(TH1D*)histo[0]->Clone();

  double error;

  for(int i=0;i<10;i++){
    for(int x=1;x<=174;x++){
      error = 0;
      histo[i]->SetBinContent(x, hist[i]->GetBinContent(x,1));
      error += (hist[i]->GetBinError(x,1)*hist[i]->GetBinError(x,1));
      for(int y=2;y<=83;y++) {
        if(hist[i]->GetBinContent(x,1)-hist[i]->GetBinContent(x,y) > 0) error += (hist[i]->GetBinContent(x,1)-hist[i]->GetBinContent(x,y))*(hist[i]->GetBinContent(x,1)-hist[i]->GetBinContent(x,y));
        //if(hist[i]->GetBinContent(x,1)-hist[i]->GetBinContent(x,y) < 0) error += (hist[i]->GetBinContent(x,1)-hist[i]->GetBinContent(x,y))*(hist[i]->GetBinContent(x,1)-hist[i]->GetBinContent(x,y));
      }
      histo[i]->SetBinError(x, sqrt(error));
    }
  }

 	histo[0]->SetFillColor(607);
 	histo[1]->SetFillColor(418);
 	histo[2]->SetFillColor(633);
 	histo[3]->SetFillColor(433);
 	histo[4]->SetFillColor(922);
 	histo[5]->SetFillColor(593);
 	histo[6]->SetFillColor(800);

 	histo[7]->SetLineColor(633);
 	histo[7]->SetLineStyle(2);
 	histo[7]->SetLineWidth(5);
 	histo[8]->SetLineColor(619);
 	histo[8]->SetLineStyle(2);
 	histo[8]->SetLineWidth(5);
 	histo[9]->SetLineColor(803);
 	histo[9]->SetLineStyle(2);
 	histo[9]->SetLineWidth(5);

  TH1D* total;
  total = (TH1D*)histo[0]->Clone();
  total->Add(histo[1]);
  total->Add(histo[2]);
  total->Add(histo[3]);
  total->Add(histo[4]);
  total->Add(histo[5]);
  total->Add(histo[6]);
  total->SetFillColor(13);
  total->SetFillStyle(3235);

  THStack *hs1 = new THStack("hs1","");
  hs1->Add(histo[6]);
  hs1->Add(histo[5]);
  hs1->Add(histo[4]);
  hs1->Add(histo[3]);
  hs1->Add(histo[2]);
  hs1->Add(histo[1]);
  hs1->Add(histo[0]);


  TCanvas* c1 = new TCanvas("c1", "", 1920, 1080);
  c1->Divide(1,2);
  TPad *canvas_up = (TPad*)c1->GetListOfPrimitives()->FindObject("c1_1");
  TPad *canvas_dw = (TPad*)c1->GetListOfPrimitives()->FindObject("c1_2");
  double up_height     = 0.83; // please tune so that the upper figures size will meet your requirement
  double dw_correction = 1.290;//1.390 // please tune so that the smaller canvas size will work in your environment

  double dw_height    = (1. - up_height) * dw_correction;

  // set pad size
  canvas_up->SetPad(0., 1 - up_height, 1., 1.);
  canvas_dw->SetPad(0., 0., 1., dw_height);

  canvas_up->SetLogy();
  canvas_up->SetFrameFillColor(0);
  canvas_up->SetFillColor(0);
  canvas_dw->SetFillColor(0);
  canvas_dw->SetFrameFillColor(0);
  // set top margin 0 for bottom figure
  canvas_dw->SetTopMargin(0.03);//0.03
  canvas_dw->SetLeftMargin(0.1);
  canvas_dw->SetBottomMargin(0.3);
  /**************************************************************************************************/    
  canvas_up->cd();
  hs1->SetMaximum(1e5);
  hs1->SetMinimum(1e-1);
  hs1->Draw("hist");
  histo[7]->Draw("sameHIST");
  histo[8]->Draw("sameHIST");
  histo[9]->Draw("sameHIST");
  total->Draw("SAME E2");

  TLine* line;
  string text;
  TLatex* lat;
  for(int i=0;i<29;i++){
    line = new TLine(i*6, 1e-1, i*6, 2.2e5);
    line->SetLineStyle(2);
    line->Draw("same");
    text = (string)histname[i];
    text = text.substr(text.find("SR"));
    lat = new TLatex(6*i+5, 2.5e2, text.c_str());
    lat->SetTextSize(0.015);
    lat->SetTextAngle(90);
    lat->Draw();
  }

  hs1->GetXaxis()->SetTitle("Bin Number");
  hs1->GetYaxis()->SetTitle("Counts");

  text = "CMS #scale[0.7]{#font[52]{Work in progress Run2}}";
  TLatex* cms_lat = new TLatex(1, 3e5, text.c_str());
  cms_lat->SetLineWidth(2);
  cms_lat->Draw();
  text = "#scale[0.7]{137 fb^{-1} (13 TeV)}";
  TLatex* era_lat = new TLatex(174*0.95,3e5, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.04);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();

  TLegend* leg = new TLegend(0.49,0.71,0.89,0.89);

  leg->SetNColumns(4);
  leg->SetTextSize(0.03);

  leg->AddEntry(histo[0],  "Multijet", "f");
  leg->AddEntry(histo[1],  "t#bar{t} + single top", "f");
  leg->AddEntry(histo[2],  "W(#rightarrowl#nu)", "f");
  leg->AddEntry(histo[3],  "Z(#rightarrow#nu#nu)", "f");
  leg->AddEntry(histo[4],  "VV(V) + t#bar{t}(X) + Higgs", "f");
  leg->AddEntry(histo[5],  "Z#rightarrowll", "f");
  leg->AddEntry(histo[6],  "#gamma+jets", "f");
  leg->AddEntry((TObject*)0,"","");
  leg->AddEntry(histo[7], "T5ttcc", "l");
  leg->AddEntry(histo[8], "T5qqqqVV", "l");
  leg->AddEntry(histo[9], "T2bW", "l");

  leg->SetFillColor(0);
  //leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->Draw();

  canvas_dw->cd();                                                                                                                                                                                                                          
  gPad->SetGridy();
  TH1D* dividend1 = (TH1D*)total->Clone();
  for(int i=1;i<=total->GetNbinsX();i++){
    dividend1->SetBinContent(i,1);
    if(total->GetBinContent(i)==0) dividend1->SetBinError(i,0);
    else dividend1->SetBinError(i,total->GetBinError(i)/total->GetBinContent(i));
  }
  dividend1->SetFillColor(12);
  dividend1->SetFillStyle(3004);
  dividend1->SetMarkerStyle(1);
  dividend1->GetYaxis()->SetTitle("");
  dividend1->GetXaxis()->SetTitle("BinNumber");
  dividend1->GetXaxis()->SetTitleOffset(0.3);
  dividend1->GetXaxis()->SetTitleSize(0.15);
  dividend1->GetXaxis()->SetLabelSize(0.);
  dividend1->GetYaxis()->SetLabelSize(0.06);
  dividend1->SetMaximum(3);
  dividend1->SetMinimum(-1);
  dividend1->GetYaxis()->SetNdivisions(505);

  dividend1->Draw("E2");

  for(int i=0;i<29;i++){
    line = new TLine(i*6, -1, i*6, 3);
    line->SetLineStyle(2);
    line->Draw("same");
  }

  c1->SaveAs("MRR2_AllsignalRegions.png"); 

}

void Unrolled_SR(){
  UnrolledPlots();
}
