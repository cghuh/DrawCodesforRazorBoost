TH1D* BinMerging(TH2D* target){
	TH1::SetDefaultSumw2();
	TH1D* hist = new TH1D("hist", "", 22,0,22);
	TH1D* htemp1 = new TH1D("htemp1", "", 2,0,2);
	TH1D* htemp2 = new TH1D("htemp2", "", 3,0,3);
	const char *label[22] = {"[0.08, 0.12]", "[0.12, 0.16]", "[0.16, 0.24]", "[0.24, 0.50]", "[0.50, 1.50]", "[0.08, 0.12]", "[0.12, 0.16]", "[0.16, 0.24]", "[0.24, 0.50]", "[0.50, 1.50]", "[0.08, 0.12]", "[0.12, 0.16]", "[0.16, 0.24]", "[0.24, 0.50]", "[0.50, 1.50]", "[0.08, 0.12]", "[0.12, 0.16]", "[0.16, 0.24]", "[0.24, 1.50]", "[0.08, 0.12]", "[0.12, 0.16]", "[0.16, 1.50]"};
	int count=1,k=1;
	for(int i=1;i<=target->GetNbinsX();i++){
		for(int j=1;j<=target->GetNbinsY();j++){
			if(k == 19) {
				htemp1->SetBinContent(1, target->GetBinContent(i,j));
				htemp1->SetBinError(1,target->GetBinError(i,j));
			}
			else if(k == 20) {
				htemp1->SetBinContent(2, target->GetBinContent(i,j));
				htemp1->SetBinError(2,target->GetBinError(i,j));
				htemp1->Rebin(2); 
				hist->SetBinContent(count, htemp1->GetBinContent(1));
				hist->SetBinError(count, htemp1->GetBinError(1));
				hist->GetXaxis()->SetBinLabel( count, label[count-1]);
				count++;
			}
			else if(k == 23) {
				htemp2->SetBinContent(1, target->GetBinContent(i,j));
				htemp2->SetBinError(1,target->GetBinError(i,j));
			}
			else if(k == 24) {
				htemp2->SetBinContent(2, target->GetBinContent(i,j));
				htemp2->SetBinError(2,target->GetBinError(i,j));
			}
			else if(k == 25) {
				htemp2->SetBinContent(3, target->GetBinContent(i,j));
				htemp2->SetBinError(3,target->GetBinError(i,j));
				htemp2->Rebin(3); 
				hist->SetBinContent(count, htemp2->GetBinContent(1));
				hist->SetBinError(count, htemp2->GetBinError(1));
				hist->GetXaxis()->SetBinLabel( count, label[count-1]);
				count++;
			}
			else {
				hist->SetBinContent(count,target->GetBinContent(i,j));
				hist->SetBinError(count,target->GetBinError(i,j));
				hist->GetXaxis()->SetBinLabel( count, label[count-1]);
				count++;
			}
			k++;
		}
	}
  delete gROOT->FindObject("htemp1");
  delete gROOT->FindObject("htemp2");
	return hist;
}
void DrawR2MR_Blind(TH1D* hMRR2_Multijet,TH1D* hMRR2_TTST,TH1D* hMRR2_Wjet,TH1D* hMRR2_Zjet,TH1D* hMRR2_DYToLL,TH1D* hMRR2_Gjet,TH1D* hMRR2_Multiboson,TH1D* sig0,TH1D* sig1,TH1D* sig2,TH1D* sig3,TString histname,TString njet,TString channel,TString target){

  TString name;
  if(histname.Contains("_S")) name = "S";
  if(histname.Contains("_s")) name = "s";
  if(histname.Contains("_Q")) name = "Q";
  if(histname.Contains("_q")) name = "q";
  if(histname.Contains("_T")) name = "T";
  if(histname.Contains("_W")) name = "W";
  if(histname.Contains("_Z")) name = "Z";
  if(histname.Contains("_G")) name = "G";
  if(histname.Contains("_P")) name = "P";
  if(target == "W") name="S_Wjet";
  if(target == "Top") name="S_Topjet";

  hMRR2_Multijet->SetFillColor(619);
  hMRR2_TTST->SetFillColor(633);
  hMRR2_Wjet->SetFillColor(418);
  hMRR2_Zjet->SetFillColor(401);
  hMRR2_DYToLL->SetFillColor(40);
  hMRR2_Gjet->SetFillColor(803);
  hMRR2_Multiboson->SetFillColor(601);

  sig0->SetLineColor(862);
  sig1->SetLineColor(841);
  sig2->SetLineColor(403);
  sig3->SetLineColor(603);

  sig0->SetLineWidth(2);
  sig1->SetLineWidth(2);
  sig2->SetLineWidth(2);
  sig3->SetLineWidth(2);

  THStack *hs1 = new THStack("hs1","");
  hs1->Add(hMRR2_Multiboson);
  hs1->Add(hMRR2_Wjet);
  hs1->Add(hMRR2_Zjet);
  hs1->Add(hMRR2_DYToLL);
  hs1->Add(hMRR2_Gjet);
  hs1->Add(hMRR2_Multijet);
  hs1->Add(hMRR2_TTST);

  TCanvas* canvas_1 = new TCanvas("canvas_1", "", 800, 800);
  canvas_1->SetLogy();
  canvas_1->SetGridx();
  canvas_1->SetGridy();
  hs1->SetMaximum(1e7);
  hs1->SetMinimum(1e-3);
  //hs1->GetXaxis()->SetTitle("MR_R2");
  hs1->Draw("HIST");
  sig0->Draw("histsame");
  sig1->Draw("histsame");
  sig2->Draw("histsame");
  sig3->Draw("histsame");
  hs1->GetXaxis()->SetLabelOffset(999);

  float xmin = hMRR2_TTST->GetXaxis()->GetBinLowEdge(hMRR2_TTST->GetXaxis()->GetFirst());
  float xmax = hMRR2_TTST->GetXaxis()->GetBinUpEdge(hMRR2_TTST->GetXaxis()->GetLast());

  //string text = "CMS Simulation #scale[0.7]{#font[52]{Work in progress 2017}}";
  string text = "CMS Simulation #scale[0.7]{#font[52]{Work in progress Run2}}";
  TLatex* cms_lat = new TLatex(xmin, 3.1e7, text.c_str());
  cms_lat->SetTextSize(0.04);
  cms_lat->SetLineWidth(2);
  cms_lat->Draw();
  if(channel == "TOP_") text = "#scale[0.7]{41.5 fb^{-1} (13 TeV)}";
  //else                  text = "#scale[0.7]{41.5 fb^{-1} (13 TeV)}";
  else                  text = "#scale[0.7]{137 fb^{-1} (13 TeV)}";
  TLatex* era_lat = new TLatex(xmax,4.0e7, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.03);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();

  text = "#scale[0.56]{   [800, 1000]   [1000, 1200]   [1200, 1600]  [1600, 2000]   [2000, 4000]}";
  era_lat = new TLatex(xmax,3.e4, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.052);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();

  auto leg = new TLegend(0.20,0.71,0.9,0.9);

  char legentry[10];

  leg->SetNColumns(4);
  leg->SetTextSize(0.03);
  //if(channel == "TOP_") leg->SetHeader("TOP ana, "+name+" region, "+njet);
  //else                  leg->SetHeader("W ana, "+name+" region, "+njet);
  sprintf(legentry, "%.2f",sig1->Integral());
  leg->AddEntry(sig1,  "TChiWZ_600_200", "l");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.2f",sig0->Integral());
  leg->AddEntry(sig0,  "TChiWZ_700_200", "l");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.2f",sig2->Integral());
  leg->AddEntry(sig2,  "TChiWZ_800_300", "l");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.2f",sig3->Integral());
  leg->AddEntry(sig3,  "TChiWZ_900_200", "l");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.2f",hMRR2_Multijet->Integral());
  leg->AddEntry(hMRR2_Multijet,  "Multijet", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.2f",hMRR2_TTST->Integral());
  leg->AddEntry(hMRR2_TTST,  "t#bar{t} + single top", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.2f",hMRR2_Wjet->Integral());
  leg->AddEntry(hMRR2_Wjet,  "W(#rightarrowl#nu)", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.2f",hMRR2_Zjet->Integral());
  leg->AddEntry(hMRR2_Zjet,  "Z(#rightarrow#nu#nu)", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.2f",hMRR2_DYToLL->Integral());
  leg->AddEntry(hMRR2_DYToLL,  "DYToLL", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.2f",hMRR2_Gjet->Integral());
  leg->AddEntry(hMRR2_Gjet,  "#gamma + jet", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.2f",hMRR2_Multiboson->Integral());
  leg->AddEntry(hMRR2_Multiboson,  "VV(V) + t#bar{t}(X)", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  leg->Draw();

  canvas_1->SaveAs("./plot/"+channel+"MRR2_"+name+njet+".png");

}
void DrawR2MR(TH1D* hMRR2_Multijet,TH1D* hMRR2_TTST,TH1D* hMRR2_Wjet,TH1D* hMRR2_Zjet,TH1D* hMRR2_DYToLL,TH1D* hMRR2_Gjet,TH1D* hMRR2_Multiboson,TH1D* hMRR2_data,TString histname,TString njet,TString channel){

  TString name;
  if(histname.Contains("_S")) name = "S";
  if(histname.Contains("_s")) name = "s";
  if(histname.Contains("_Q")) name = "Q";
  if(histname.Contains("_q")) name = "q";
  if(histname.Contains("_T")) name = "T";
  if(histname.Contains("_W")) name = "W";
  if(histname.Contains("_Z")) name = "Z";
  if(histname.Contains("_G")) name = "G";
  if(histname.Contains("_P")) name = "P";

  hMRR2_Multijet->SetFillColor(619);
  hMRR2_TTST->SetFillColor(633);
  hMRR2_Wjet->SetFillColor(418);
  hMRR2_Zjet->SetFillColor(401);
  hMRR2_DYToLL->SetFillColor(40);
  hMRR2_Gjet->SetFillColor(803);
  hMRR2_Multiboson->SetFillColor(601);

  THStack *hs1 = new THStack("hs1","");
  hs1->Add(hMRR2_Multiboson);
  hs1->Add(hMRR2_Wjet);
  hs1->Add(hMRR2_Zjet);
  hs1->Add(hMRR2_DYToLL);
  hs1->Add(hMRR2_Gjet);
  hs1->Add(hMRR2_Multijet);
  hs1->Add(hMRR2_TTST);

  TCanvas* canvas_1 = new TCanvas("canvas_1", "", 800, 800);
  canvas_1->Divide(1,2);
  TPad *canvas_up = (TPad*)canvas_1->GetListOfPrimitives()->FindObject("canvas_1_1");
  TPad *canvas_dw = (TPad*)canvas_1->GetListOfPrimitives()->FindObject("canvas_1_2");

  double up_height     = 0.7; // please tune so that the upper figures size will meet your requirement
  double dw_correction = 1.210;//1.390 // please tune so that the smaller canvas size will work in your environment
  double dw_height    = (1. - up_height) * dw_correction;

  canvas_up->SetPad(0., 1.-up_height, 1., 1.);
  canvas_dw->SetPad(0., 0., 1., dw_height);
  canvas_up->SetLogy();
  canvas_up->SetFrameFillColor(0);
  canvas_up->SetFillColor(0);
  canvas_dw->SetFillColor(0);
  canvas_dw->SetFrameFillColor(0);
  canvas_up->SetLeftMargin(0.1);
  canvas_dw->SetTopMargin(0.03);//0.03
  canvas_dw->SetLeftMargin(0.1);
  canvas_dw->SetBottomMargin(0.3);
  canvas_up->SetGridx();
  canvas_up->SetGridy();

  canvas_up->cd();
  hMRR2_data->SetMarkerColor(kBlack);
  hMRR2_data->SetLineColor(1);
  hMRR2_data->SetMarkerStyle(20);
  //hMRR2_data->SetMarkerSize(0.);
  hMRR2_data->SetMaximum(1e7);
  hMRR2_data->SetMinimum(1e-1);
  hs1->SetMaximum(1e7);
  hs1->SetMinimum(1e-1);
  //hs1->GetXaxis()->SetTitle("MR_R2");
  hMRR2_data->GetXaxis()->SetTitle("MR_R2");
  hMRR2_data->GetXaxis()->SetLabelOffset(999);
  hs1->Draw("HIST");
  hMRR2_data->Draw("EPsame");
  //hs1->GetXaxis()->SetLabelOffset(999);

  float xmin = hMRR2_data->GetXaxis()->GetBinLowEdge(hMRR2_data->GetXaxis()->GetFirst());
  float xmax = hMRR2_data->GetXaxis()->GetBinUpEdge(hMRR2_data->GetXaxis()->GetLast());

  string text = "CMS Simulation #scale[0.7]{#font[52]{Work in progress 2017}}";
  TLatex* cms_lat = new TLatex(xmin, 3.0e7, text.c_str());
  cms_lat->SetTextSize(0.04);
  cms_lat->SetLineWidth(2);
  cms_lat->Draw();
  if(channel == "TOP_") text = "#scale[0.7]{TOP ana, 41.5 fb^{-1} (13 TeV)}";
  else                  text = "#scale[0.7]{W ana, 41.5 fb^{-1} (13 TeV)}";
  TLatex* era_lat = new TLatex(xmax,3.4e7, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.03);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();

  text = "#scale[0.7]{   [800, 1000]   [1000, 1200]   [1200, 1600]  [1600, 2000]   [2000, 4000]}";
  era_lat = new TLatex(xmax,3.e4, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.052);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();

  auto leg = new TLegend(0.40,0.71,0.9,0.9);

  char legentry[10];

  leg->SetNColumns(4);
  leg->SetTextSize(0.03);
  if(channel == "TOP_") leg->SetHeader("TOP ana, "+name+" region, "+njet);
  else                  leg->SetHeader("W ana, "+name+" region, "+njet);
  sprintf(legentry, "%.1f",hMRR2_data->Integral());
  leg->AddEntry(hMRR2_data, "data", "lep");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_Multijet->Integral());
  leg->AddEntry(hMRR2_Multijet,  "Multijet", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_TTST->Integral());
  leg->AddEntry(hMRR2_TTST,  "t#bar{t} + single top", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_Wjet->Integral());
  leg->AddEntry(hMRR2_Wjet,  "W(#rightarrowl#nu)", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_Zjet->Integral());
  leg->AddEntry(hMRR2_Zjet,  "Z(#rightarrow#nu#nu)", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_DYToLL->Integral());
  leg->AddEntry(hMRR2_DYToLL,  "DYToLL", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_Gjet->Integral());
  leg->AddEntry(hMRR2_Gjet,  "#gamma + jet", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_Multiboson->Integral());
  leg->AddEntry(hMRR2_Multiboson,  "VV(V) + t#bar{t}(X)", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  leg->Draw();

	TH1D* dividend  = (TH1D*)hMRR2_data->Clone();
	TH1D* dividend1 = (TH1D*)hMRR2_data->Clone();
  TH1D* Tot_MC    = (TH1D*)hMRR2_Multijet->Clone();
  Tot_MC->Add(hMRR2_TTST);
  Tot_MC->Add(hMRR2_Wjet);
  Tot_MC->Add(hMRR2_Zjet);
  Tot_MC->Add(hMRR2_DYToLL);
  Tot_MC->Add(hMRR2_Gjet);
  Tot_MC->Add(hMRR2_Multiboson);
  //Tot_MC->Add();

	canvas_dw->cd();
  gPad->SetGridy();

	dividend->Divide(Tot_MC);

	dividend->SetLineColor(kBlack);
	dividend->SetMarkerColor(kBlack);
	dividend->SetMarkerStyle(20);
	dividend->SetMarkerSize(1.3);

	dividend->GetYaxis()->SetRangeUser(0.0,2.0);

	dividend->SetTitle(0);
	dividend->SetStats(0);

	dividend->GetYaxis()->SetLabelSize(0.09);
	dividend->GetXaxis()->SetLabelSize(0.09);

	dividend->GetYaxis()->SetNdivisions(502);

	dividend->GetXaxis()->SetTickLength(0.07);
	dividend->GetYaxis()->SetTickLength(0.03);

	dividend->GetYaxis()->SetTitle("Data / MC");

	dividend->GetXaxis()->SetTitleOffset(0.81);
	dividend->GetXaxis()->SetTitleSize(0.12);

	dividend->GetYaxis()->SetTitleOffset(0.35);
	dividend->GetYaxis()->SetTitleSize(0.12);//dividend->SetTitleSize(0.12,"Y");
  for(int i=1;i<=Tot_MC->GetNbinsX();i++){
    dividend1->SetBinContent(i,1);
    if(Tot_MC->GetBinContent(i)==0) dividend1->SetBinError(i,0);
    else dividend1->SetBinError(i,Tot_MC->GetBinError(i)/Tot_MC->GetBinContent(i));
  }
	dividend1->SetFillColor(13);
	//dividend1->SetMarkerColor(kWhite);
	dividend1->SetFillStyle(3235);
  dividend1->SetMarkerStyle(1);
	
	dividend->Draw("PE1");
	dividend1->Draw("sameE2");
  

  canvas_1->SaveAs("./plot/"+channel+"MRR2_"+name+njet+".png");

}
void ClosureTest(TString njet, TString channel){
  cout << njet + channel << endl;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();
  TString dir = "/eos/cms/store/user/chuh/RazorBoost/191104_"+channel+"PFHT1050_PFHTXXX_PFMETXXX_PFMHTXXX_IDTight/added/";
  TFile* file1 = TFile::Open(dir+"Multijet.root");
  TFile* file2 = TFile::Open(dir+"TTST.root");
  TFile* file3 = TFile::Open(dir+"WJet.root");
  TFile* file4 = TFile::Open(dir+"ZJet.root");
  TFile* file5 = TFile::Open(dir+"DYToLL.root");
  TFile* file6 = TFile::Open(dir+"Multiboson+TTX.root");
  TFile* file7 = TFile::Open(dir+"GJet.root");
  TFile* file0 = TFile::Open(dir+"data.root");
  TFile* files0 = TFile::Open(dir+"TChiWZ.root");
  TFile* files1 = TFile::Open(dir+"TChiWZ_600_200.root");
  TFile* files2 = TFile::Open(dir+"TChiWZ_800_300.root");
  TFile* files3 = TFile::Open(dir+"TChiWZ_900_200.root");

  //TFile* files0 = TFile::Open(dir+"T1tttt.root");
  //TFile* files1 = TFile::Open(dir+"T2tt.root");
  //TFile* files2 = TFile::Open(dir+"T5ttcc.root");

  TString histname_Wjet_S = "R2_MR_Wjet_S"+njet;
  TString histname_Topjet_S = "R2_MR_Topjet_S"+njet;
  if(TString(njet).Contains("nj")) histname_Topjet_S = "R2_MR_Topjet_Topjet_S"+njet;
  //TString histname_Wjet_S = "R2_MR_S"+njet;
  //TString histname_Topjet_S = "R2_MR_S"+njet;
  cout << histname_Topjet_S << endl;
  TString histname_Q = "R2_MR_Q"+njet;
  TString histname_s = "R2_MR_s"+njet;
  TString histname_q = "R2_MR_q"+njet;
  TString histname_T = "R2_MR_T"+njet;
  TString histname_W = "R2_MR_W"+njet;
  TString histname_Z = "R2_MR_Z"+njet;
  TString histname_G = "R2_MR_G"+njet;
  TString histname_P = "R2_MR_P";

  histname_Wjet_S = "R2_MR_WZc0b_S"+njet;
  histname_Wjet_S = "R2_MR_WZ0b_S_mDPhi"+njet;
  TH2D*   Multijet_Wjet_S = (TH2D*)file1->Get(histname_Wjet_S);
  TH2D*       TTST_Wjet_S = (TH2D*)file2->Get(histname_Wjet_S);
  TH2D*       Wjet_Wjet_S = (TH2D*)file3->Get(histname_Wjet_S);
  TH2D*       Zjet_Wjet_S = (TH2D*)file4->Get(histname_Wjet_S);
  TH2D*     DYToLL_Wjet_S = (TH2D*)file5->Get(histname_Wjet_S);
  TH2D* Multiboson_Wjet_S = (TH2D*)file6->Get(histname_Wjet_S);
  TH2D*       Gjet_Wjet_S = (TH2D*)file7->Get(histname_Wjet_S);
  TH2D*       data_Wjet_S = (TH2D*)file0->Get(histname_Wjet_S);
  TH2D*   Multijet_Topjet_S = (TH2D*)file1->Get(histname_Topjet_S);
  TH2D*       TTST_Topjet_S = (TH2D*)file2->Get(histname_Topjet_S);
  TH2D*       Wjet_Topjet_S = (TH2D*)file3->Get(histname_Topjet_S);
  TH2D*       Zjet_Topjet_S = (TH2D*)file4->Get(histname_Topjet_S);
  TH2D*     DYToLL_Topjet_S = (TH2D*)file5->Get(histname_Topjet_S);
  TH2D* Multiboson_Topjet_S = (TH2D*)file6->Get(histname_Topjet_S);
  TH2D*       Gjet_Topjet_S = (TH2D*)file7->Get(histname_Topjet_S);
  TH2D*       data_Topjet_S = (TH2D*)file0->Get(histname_Topjet_S);
  TH2D*       sig0_Wjet_S = (TH2D*)files0->Get(histname_Wjet_S);
  TH2D*       sig1_Wjet_S = (TH2D*)files1->Get(histname_Wjet_S);
  TH2D*       sig2_Wjet_S = (TH2D*)files2->Get(histname_Wjet_S);
  TH2D*       sig3_Wjet_S = (TH2D*)files3->Get(histname_Wjet_S);
  TH2D*       sig0_Topjet_S = (TH2D*)files0->Get(histname_Topjet_S);
  TH2D*       sig1_Topjet_S = (TH2D*)files1->Get(histname_Topjet_S);
  TH2D*       sig2_Topjet_S = (TH2D*)files2->Get(histname_Topjet_S);
  TH2D*       sig3_Topjet_S = (TH2D*)files3->Get(histname_Topjet_S);
  TH2D*   Multijet_s = (TH2D*)file1->Get(histname_s);
  TH2D*       TTST_s = (TH2D*)file2->Get(histname_s);
  TH2D*       Wjet_s = (TH2D*)file3->Get(histname_s);
  TH2D*       Zjet_s = (TH2D*)file4->Get(histname_s);
  TH2D*     DYToLL_s = (TH2D*)file5->Get(histname_s);
  TH2D* Multiboson_s = (TH2D*)file6->Get(histname_s);
  TH2D*       Gjet_s = (TH2D*)file7->Get(histname_s);
  TH2D*       data_s = (TH2D*)file0->Get(histname_s);
  TH2D*   Multijet_Q = (TH2D*)file1->Get(histname_Q);
  TH2D*       TTST_Q = (TH2D*)file2->Get(histname_Q);
  TH2D*       Wjet_Q = (TH2D*)file3->Get(histname_Q);
  TH2D*       Zjet_Q = (TH2D*)file4->Get(histname_Q);
  TH2D*     DYToLL_Q = (TH2D*)file5->Get(histname_Q);
  TH2D* Multiboson_Q = (TH2D*)file6->Get(histname_Q);
  TH2D*       Gjet_Q = (TH2D*)file7->Get(histname_Q);
  TH2D*       data_Q = (TH2D*)file0->Get(histname_Q);
  TH2D*   Multijet_q = (TH2D*)file1->Get(histname_q);
  TH2D*       TTST_q = (TH2D*)file2->Get(histname_q);
  TH2D*       Wjet_q = (TH2D*)file3->Get(histname_q);
  TH2D*       Zjet_q = (TH2D*)file4->Get(histname_q);
  TH2D*     DYToLL_q = (TH2D*)file5->Get(histname_q);
  TH2D* Multiboson_q = (TH2D*)file6->Get(histname_q);
  TH2D*       Gjet_q = (TH2D*)file7->Get(histname_q);
  TH2D*       data_q = (TH2D*)file0->Get(histname_q);
  TH2D*   Multijet_T = (TH2D*)file1->Get(histname_T);
  TH2D*       TTST_T = (TH2D*)file2->Get(histname_T);
  TH2D*       Wjet_T = (TH2D*)file3->Get(histname_T);
  TH2D*       Zjet_T = (TH2D*)file4->Get(histname_T);
  TH2D*     DYToLL_T = (TH2D*)file5->Get(histname_T);
  TH2D* Multiboson_T = (TH2D*)file6->Get(histname_T);
  TH2D*       Gjet_T = (TH2D*)file7->Get(histname_T);
  TH2D*       data_T = (TH2D*)file0->Get(histname_T);
  TH2D*   Multijet_W = (TH2D*)file1->Get(histname_W);
  TH2D*       TTST_W = (TH2D*)file2->Get(histname_W);
  TH2D*       Wjet_W = (TH2D*)file3->Get(histname_W);
  TH2D*       Zjet_W = (TH2D*)file4->Get(histname_W);
  TH2D*     DYToLL_W = (TH2D*)file5->Get(histname_W);
  TH2D* Multiboson_W = (TH2D*)file6->Get(histname_W);
  TH2D*       Gjet_W = (TH2D*)file7->Get(histname_W);
  TH2D*       data_W = (TH2D*)file0->Get(histname_W);
  TH2D*   Multijet_Z = (TH2D*)file1->Get(histname_Z);
  TH2D*       TTST_Z = (TH2D*)file2->Get(histname_Z);
  TH2D*       Wjet_Z = (TH2D*)file3->Get(histname_Z);
  TH2D*       Zjet_Z = (TH2D*)file4->Get(histname_Z);
  TH2D*     DYToLL_Z = (TH2D*)file5->Get(histname_Z);
  TH2D* Multiboson_Z = (TH2D*)file6->Get(histname_Z);
  TH2D*       Gjet_Z = (TH2D*)file7->Get(histname_Z);
  TH2D*       data_Z = (TH2D*)file0->Get(histname_Z);
  TH2D*   Multijet_G = (TH2D*)file1->Get(histname_G);
  TH2D*       TTST_G = (TH2D*)file2->Get(histname_G);
  TH2D*       Wjet_G = (TH2D*)file3->Get(histname_G);
  TH2D*       Zjet_G = (TH2D*)file4->Get(histname_G);
  TH2D*     DYToLL_G = (TH2D*)file5->Get(histname_G);
  TH2D* Multiboson_G = (TH2D*)file6->Get(histname_G);
  TH2D*       Gjet_G = (TH2D*)file7->Get(histname_G);
  TH2D*       data_G = (TH2D*)file0->Get(histname_G);
  TH2D*   Multijet_P = (TH2D*)file1->Get(histname_P);
  TH2D*       TTST_P = (TH2D*)file2->Get(histname_P);
  TH2D*       Wjet_P = (TH2D*)file3->Get(histname_P);
  TH2D*       Zjet_P = (TH2D*)file4->Get(histname_P);
  TH2D*     DYToLL_P = (TH2D*)file5->Get(histname_P);
  TH2D* Multiboson_P = (TH2D*)file6->Get(histname_P);
  TH2D*       Gjet_P = (TH2D*)file7->Get(histname_P);
  TH2D*       data_P = (TH2D*)file0->Get(histname_P);

  TH1D* hMRR2_Wjet_S_Multijet = new TH1D("MRR2_Wjet_S_Multijet", "", 22,0,22);
  TH1D* hMRR2_Wjet_S_TTST = new TH1D("MRR2_Wjet_S_TTST", "", 22,0,22);
  TH1D* hMRR2_Wjet_S_Wjet = new TH1D("MRR2_Wjet_S_Wjet", "", 22,0,22);
  TH1D* hMRR2_Wjet_S_Zjet = new TH1D("MRR2_Wjet_S_Zjet", "", 22,0,22);
  TH1D* hMRR2_Wjet_S_DYToLL = new TH1D("MRR2_Wjet_S_DYToLL", "", 22,0,22);
  TH1D* hMRR2_Wjet_S_Multiboson = new TH1D("MRR2_Wjet_S_Multiboson", "", 22,0,22);
  TH1D* hMRR2_Wjet_S_Gjet = new TH1D("MRR2_Wjet_S_Gjet", "", 22,0,22);
  TH1D* hMRR2_Wjet_S_data = new TH1D("MRR2_Wjet_S_data", "", 22,0,22);

  TH1D* hMRR2_Topjet_S_Multijet = new TH1D("MRR2_Topjet_S_Multijet", "", 22,0,22);
  TH1D* hMRR2_Topjet_S_TTST = new TH1D("MRR2_Topjet_S_TTST", "", 22,0,22);
  TH1D* hMRR2_Topjet_S_Wjet = new TH1D("MRR2_Topjet_S_Wjet", "", 22,0,22);
  TH1D* hMRR2_Topjet_S_Zjet = new TH1D("MRR2_Topjet_S_Zjet", "", 22,0,22);
  TH1D* hMRR2_Topjet_S_DYToLL = new TH1D("MRR2_Topjet_S_DYToLL", "", 22,0,22);
  TH1D* hMRR2_Topjet_S_Multiboson = new TH1D("MRR2_Topjet_S_Multiboson", "", 22,0,22);
  TH1D* hMRR2_Topjet_S_Gjet = new TH1D("MRR2_Topjet_S_Gjet", "", 22,0,22);
  TH1D* hMRR2_Topjet_S_data = new TH1D("MRR2_Topjet_S_data", "", 22,0,22);

  TH1D* hMRR2_Wjet_S_sig0 = new TH1D("MRR2_Wjet_S_sig0", "", 22,0,22);
  TH1D* hMRR2_Wjet_S_sig1 = new TH1D("MRR2_Wjet_S_sig1", "", 22,0,22);
  TH1D* hMRR2_Wjet_S_sig2 = new TH1D("MRR2_Wjet_S_sig2", "", 22,0,22);
  TH1D* hMRR2_Wjet_S_sig3 = new TH1D("MRR2_Wjet_S_sig3", "", 22,0,22);
  TH1D* hMRR2_Topjet_S_sig0 = new TH1D("MRR2_Topjet_S_sig0", "", 22,0,22);
  TH1D* hMRR2_Topjet_S_sig1 = new TH1D("MRR2_Topjet_S_sig1", "", 22,0,22);
  TH1D* hMRR2_Topjet_S_sig2 = new TH1D("MRR2_Topjet_S_sig2", "", 22,0,22);
  TH1D* hMRR2_Topjet_S_sig3 = new TH1D("MRR2_Topjet_S_sig3", "", 22,0,22);

  TH1D* hMRR2_s_Multijet = new TH1D("MRR2_s_Multijet", "", 22,0,22);
  TH1D* hMRR2_s_TTST = new TH1D("MRR2_s_TTST", "", 22,0,22);
  TH1D* hMRR2_s_Wjet = new TH1D("MRR2_s_Wjet", "", 22,0,22);
  TH1D* hMRR2_s_Zjet = new TH1D("MRR2_s_Zjet", "", 22,0,22);
  TH1D* hMRR2_s_DYToLL = new TH1D("MRR2_s_DYToLL", "", 22,0,22);
  TH1D* hMRR2_s_Multiboson = new TH1D("MRR2_s_Multiboson", "", 22,0,22);
  TH1D* hMRR2_s_Gjet = new TH1D("MRR2_s_Gjet", "", 22,0,22);
  TH1D* hMRR2_s_data = new TH1D("MRR2_s_data", "", 22,0,22);

  TH1D* hMRR2_Q_Multijet = new TH1D("MRR2_Q_Multijet", "", 22,0,22);
  TH1D* hMRR2_Q_TTST = new TH1D("MRR2_Q_TTST", "", 22,0,22);
  TH1D* hMRR2_Q_Wjet = new TH1D("MRR2_Q_Wjet", "", 22,0,22);
  TH1D* hMRR2_Q_Zjet = new TH1D("MRR2_Q_Zjet", "", 22,0,22);
  TH1D* hMRR2_Q_DYToLL = new TH1D("MRR2_Q_DYToLL", "", 22,0,22);
  TH1D* hMRR2_Q_Multiboson = new TH1D("MRR2_Q_Multiboson", "", 22,0,22);
  TH1D* hMRR2_Q_Gjet = new TH1D("MRR2_Q_Gjet", "", 22,0,22);
  TH1D* hMRR2_Q_data = new TH1D("MRR2_Q_data", "", 22,0,22);

  TH1D* hMRR2_q_Multijet = new TH1D("MRR2_q_Multijet", "", 22,0,22);
  TH1D* hMRR2_q_TTST = new TH1D("MRR2_q_TTST", "", 22,0,22);
  TH1D* hMRR2_q_Wjet = new TH1D("MRR2_q_Wjet", "", 22,0,22);
  TH1D* hMRR2_q_Zjet = new TH1D("MRR2_q_Zjet", "", 22,0,22);
  TH1D* hMRR2_q_DYToLL = new TH1D("MRR2_q_DYToLL", "", 22,0,22);
  TH1D* hMRR2_q_Multiboson = new TH1D("MRR2_q_Multiboson", "", 22,0,22);
  TH1D* hMRR2_q_Gjet = new TH1D("MRR2_q_Gjet", "", 22,0,22);
  TH1D* hMRR2_q_data = new TH1D("MRR2_q_data", "", 22,0,22);

  TH1D* hMRR2_T_Multijet = new TH1D("MRR2_T_Multijet", "", 22,0,22);
  TH1D* hMRR2_T_TTST = new TH1D("MRR2_T_TTST", "", 22,0,22);
  TH1D* hMRR2_T_Wjet = new TH1D("MRR2_T_Wjet", "", 22,0,22);
  TH1D* hMRR2_T_Zjet = new TH1D("MRR2_T_Zjet", "", 22,0,22);
  TH1D* hMRR2_T_DYToLL = new TH1D("MRR2_T_DYToLL", "", 22,0,22);
  TH1D* hMRR2_T_Multiboson = new TH1D("MRR2_T_Multiboson", "", 22,0,22);
  TH1D* hMRR2_T_Gjet = new TH1D("MRR2_T_Gjet", "", 22,0,22);
  TH1D* hMRR2_T_data = new TH1D("MRR2_T_data", "", 22,0,22);

  TH1D* hMRR2_W_Multijet = new TH1D("MRR2_W_Multijet", "", 22,0,22);
  TH1D* hMRR2_W_TTST = new TH1D("MRR2_W_TTST", "", 22,0,22);
  TH1D* hMRR2_W_Wjet = new TH1D("MRR2_W_Wjet", "", 22,0,22);
  TH1D* hMRR2_W_Zjet = new TH1D("MRR2_W_Zjet", "", 22,0,22);
  TH1D* hMRR2_W_DYToLL = new TH1D("MRR2_W_DYToLL", "", 22,0,22);
  TH1D* hMRR2_W_Multiboson = new TH1D("MRR2_W_Multiboson", "", 22,0,22);
  TH1D* hMRR2_W_Gjet = new TH1D("MRR2_W_Gjet", "", 22,0,22);
  TH1D* hMRR2_W_data = new TH1D("MRR2_W_data", "", 22,0,22);

  TH1D* hMRR2_Z_Multijet = new TH1D("MRR2_Z_Multijet", "", 22,0,22);
  TH1D* hMRR2_Z_TTST = new TH1D("MRR2_Z_TTST", "", 22,0,22);
  TH1D* hMRR2_Z_Wjet = new TH1D("MRR2_Z_Wjet", "", 22,0,22);
  TH1D* hMRR2_Z_Zjet = new TH1D("MRR2_Z_Zjet", "", 22,0,22);
  TH1D* hMRR2_Z_DYToLL = new TH1D("MRR2_Z_DYToLL", "", 22,0,22);
  TH1D* hMRR2_Z_Multiboson = new TH1D("MRR2_Z_Multiboson", "", 22,0,22);
  TH1D* hMRR2_Z_Gjet = new TH1D("MRR2_Z_Gjet", "", 22,0,22);
  TH1D* hMRR2_Z_data = new TH1D("MRR2_Z_data", "", 22,0,22);

  TH1D* hMRR2_G_Multijet = new TH1D("MRR2_G_Multijet", "", 22,0,22);
  TH1D* hMRR2_G_TTST = new TH1D("MRR2_G_TTST", "", 22,0,22);
  TH1D* hMRR2_G_Wjet = new TH1D("MRR2_G_Wjet", "", 22,0,22);
  TH1D* hMRR2_G_Zjet = new TH1D("MRR2_G_Zjet", "", 22,0,22);
  TH1D* hMRR2_G_DYToLL = new TH1D("MRR2_G_DYToLL", "", 22,0,22);
  TH1D* hMRR2_G_Multiboson = new TH1D("MRR2_G_Multiboson", "", 22,0,22);
  TH1D* hMRR2_G_Gjet = new TH1D("MRR2_G_Gjet", "", 22,0,22);
  TH1D* hMRR2_G_data = new TH1D("MRR2_G_data", "", 22,0,22);

  TH1D* hMRR2_P_Multijet = new TH1D("MRR2_P_Multijet", "", 22,0,22);
  TH1D* hMRR2_P_TTST = new TH1D("MRR2_P_TTST", "", 22,0,22);
  TH1D* hMRR2_P_Wjet = new TH1D("MRR2_P_Wjet", "", 22,0,22);
  TH1D* hMRR2_P_Zjet = new TH1D("MRR2_P_Zjet", "", 22,0,22);
  TH1D* hMRR2_P_DYToLL = new TH1D("MRR2_P_DYToLL", "", 22,0,22);
  TH1D* hMRR2_P_Multiboson = new TH1D("MRR2_P_Multiboson", "", 22,0,22);
  TH1D* hMRR2_P_Gjet = new TH1D("MRR2_P_Gjet", "", 22,0,22);
  TH1D* hMRR2_P_data = new TH1D("MRR2_P_data", "", 22,0,22);

  hMRR2_Wjet_S_Multijet = BinMerging(Multijet_Wjet_S);
  hMRR2_Wjet_S_TTST = BinMerging(TTST_Wjet_S);
  hMRR2_Wjet_S_Wjet = BinMerging(Wjet_Wjet_S);
  hMRR2_Wjet_S_Zjet = BinMerging(Zjet_Wjet_S);
  hMRR2_Wjet_S_DYToLL = BinMerging(DYToLL_Wjet_S);
  hMRR2_Wjet_S_Multiboson = BinMerging(Multiboson_Wjet_S);
  hMRR2_Wjet_S_Gjet = BinMerging(Gjet_Wjet_S);
  hMRR2_Wjet_S_data = BinMerging(data_Wjet_S);
  hMRR2_Topjet_S_Multijet = BinMerging(Multijet_Topjet_S);
  hMRR2_Topjet_S_TTST = BinMerging(TTST_Topjet_S);
  hMRR2_Topjet_S_Wjet = BinMerging(Wjet_Topjet_S);
  hMRR2_Topjet_S_Zjet = BinMerging(Zjet_Topjet_S);
  hMRR2_Topjet_S_DYToLL = BinMerging(DYToLL_Topjet_S);
  hMRR2_Topjet_S_Multiboson = BinMerging(Multiboson_Topjet_S);
  hMRR2_Topjet_S_Gjet = BinMerging(Gjet_Topjet_S);
  hMRR2_Topjet_S_data = BinMerging(data_Topjet_S);

  hMRR2_Wjet_S_sig0 = BinMerging(sig0_Wjet_S);
  hMRR2_Wjet_S_sig1 = BinMerging(sig1_Wjet_S);
  hMRR2_Wjet_S_sig2 = BinMerging(sig2_Wjet_S);
  hMRR2_Wjet_S_sig3 = BinMerging(sig3_Wjet_S);
  hMRR2_Topjet_S_sig0 = BinMerging(sig0_Topjet_S);
  hMRR2_Topjet_S_sig1 = BinMerging(sig1_Topjet_S);
  hMRR2_Topjet_S_sig2 = BinMerging(sig2_Topjet_S);
  hMRR2_Topjet_S_sig3 = BinMerging(sig3_Topjet_S);

  hMRR2_s_Multijet = BinMerging(Multijet_s);
  hMRR2_s_TTST = BinMerging(TTST_s);
  hMRR2_s_Wjet = BinMerging(Wjet_s);
  hMRR2_s_Zjet = BinMerging(Zjet_s);
  hMRR2_s_DYToLL = BinMerging(DYToLL_s);
  hMRR2_s_Multiboson = BinMerging(Multiboson_s);
  hMRR2_s_Gjet = BinMerging(Gjet_s);
  hMRR2_s_data = BinMerging(data_s);
  hMRR2_Q_Multijet = BinMerging(Multijet_Q);
  hMRR2_Q_TTST = BinMerging(TTST_Q);
  hMRR2_Q_Wjet = BinMerging(Wjet_Q);
  hMRR2_Q_Zjet = BinMerging(Zjet_Q);
  hMRR2_Q_DYToLL = BinMerging(DYToLL_Q);
  hMRR2_Q_Multiboson = BinMerging(Multiboson_Q);
  hMRR2_Q_Gjet = BinMerging(Gjet_Q);
  hMRR2_Q_data = BinMerging(data_Q);
  hMRR2_q_Multijet = BinMerging(Multijet_q);
  hMRR2_q_TTST = BinMerging(TTST_q);
  hMRR2_q_Wjet = BinMerging(Wjet_q);
  hMRR2_q_Zjet = BinMerging(Zjet_q);
  hMRR2_q_DYToLL = BinMerging(DYToLL_q);
  hMRR2_q_Multiboson = BinMerging(Multiboson_q);
  hMRR2_q_Gjet = BinMerging(Gjet_q);
  hMRR2_q_data = BinMerging(data_q);
  hMRR2_T_Multijet = BinMerging(Multijet_T);
  hMRR2_T_TTST = BinMerging(TTST_T);
  hMRR2_T_Wjet = BinMerging(Wjet_T);
  hMRR2_T_Zjet = BinMerging(Zjet_T);
  hMRR2_T_DYToLL = BinMerging(DYToLL_T);
  hMRR2_T_Multiboson = BinMerging(Multiboson_T);
  hMRR2_T_Gjet = BinMerging(Gjet_T);
  hMRR2_T_data = BinMerging(data_T);
  hMRR2_W_Multijet = BinMerging(Multijet_W);
  hMRR2_W_TTST = BinMerging(TTST_W);
  hMRR2_W_Wjet = BinMerging(Wjet_W);
  hMRR2_W_Zjet = BinMerging(Zjet_W);
  hMRR2_W_DYToLL = BinMerging(DYToLL_W);
  hMRR2_W_Multiboson = BinMerging(Multiboson_W);
  hMRR2_W_Gjet = BinMerging(Gjet_W);
  hMRR2_W_data = BinMerging(data_W);
  hMRR2_Z_Multijet = BinMerging(Multijet_Z);
  hMRR2_Z_TTST = BinMerging(TTST_Z);
  hMRR2_Z_Wjet = BinMerging(Wjet_Z);
  hMRR2_Z_Zjet = BinMerging(Zjet_Z);
  hMRR2_Z_DYToLL = BinMerging(DYToLL_Z);
  hMRR2_Z_Multiboson = BinMerging(Multiboson_Z);
  hMRR2_Z_Gjet = BinMerging(Gjet_Z);
  hMRR2_Z_data = BinMerging(data_Z);
  hMRR2_G_Multijet = BinMerging(Multijet_G);
  hMRR2_G_TTST = BinMerging(TTST_G);
  hMRR2_G_Wjet = BinMerging(Wjet_G);
  hMRR2_G_Zjet = BinMerging(Zjet_G);
  hMRR2_G_DYToLL = BinMerging(DYToLL_G);
  hMRR2_G_Multiboson = BinMerging(Multiboson_G);
  hMRR2_G_Gjet = BinMerging(Gjet_G);
  hMRR2_G_data = BinMerging(data_G);
  hMRR2_P_Multijet = BinMerging(Multijet_P);
  hMRR2_P_TTST = BinMerging(TTST_P);
  hMRR2_P_Wjet = BinMerging(Wjet_P);
  hMRR2_P_Zjet = BinMerging(Zjet_P);
  hMRR2_P_DYToLL = BinMerging(DYToLL_P);
  hMRR2_P_Multiboson = BinMerging(Multiboson_P);
  hMRR2_P_Gjet = BinMerging(Gjet_P);
  hMRR2_P_data = BinMerging(data_P);

  DrawR2MR_Blind(hMRR2_Wjet_S_Multijet,hMRR2_Wjet_S_TTST,hMRR2_Wjet_S_Wjet,hMRR2_Wjet_S_Zjet,hMRR2_Wjet_S_DYToLL,hMRR2_Wjet_S_Gjet,hMRR2_Wjet_S_Multiboson,hMRR2_Wjet_S_sig0,hMRR2_Wjet_S_sig1,hMRR2_Wjet_S_sig2,hMRR2_Wjet_S_sig3,histname_Wjet_S,njet,channel, "W");
  //DrawR2MR_Blind(hMRR2_Topjet_S_Multijet,hMRR2_Topjet_S_TTST,hMRR2_Topjet_S_Wjet,hMRR2_Topjet_S_Zjet,hMRR2_Topjet_S_DYToLL,hMRR2_Topjet_S_Gjet,hMRR2_Topjet_S_Multiboson,hMRR2_Topjet_S_sig0,hMRR2_Topjet_S_sig1,hMRR2_Topjet_S_sig2,histname_Topjet_S,njet,channel, "TOP");
  //DrawR2MR_Blind(hMRR2_Wjet_S_Multijet,hMRR2_Wjet_S_TTST,hMRR2_Wjet_S_Wjet,hMRR2_Wjet_S_Zjet,hMRR2_Wjet_S_DYToLL,hMRR2_Wjet_S_Gjet,hMRR2_Wjet_S_Multiboson,histname_Wjet_S,njet,channel);
  DrawR2MR(hMRR2_s_Multijet,hMRR2_s_TTST,hMRR2_s_Wjet,hMRR2_s_Zjet,hMRR2_s_DYToLL,hMRR2_s_Gjet,hMRR2_s_Multiboson,hMRR2_s_data,histname_s,njet,channel);
  DrawR2MR(hMRR2_Q_Multijet,hMRR2_Q_TTST,hMRR2_Q_Wjet,hMRR2_Q_Zjet,hMRR2_Q_DYToLL,hMRR2_Q_Gjet,hMRR2_Q_Multiboson,hMRR2_Q_data,histname_Q,njet,channel);
  DrawR2MR(hMRR2_q_Multijet,hMRR2_q_TTST,hMRR2_q_Wjet,hMRR2_q_Zjet,hMRR2_q_DYToLL,hMRR2_q_Gjet,hMRR2_q_Multiboson,hMRR2_q_data,histname_q,njet,channel);
  DrawR2MR(hMRR2_T_Multijet,hMRR2_T_TTST,hMRR2_T_Wjet,hMRR2_T_Zjet,hMRR2_T_DYToLL,hMRR2_T_Gjet,hMRR2_T_Multiboson,hMRR2_T_data,histname_T,njet,channel);
  DrawR2MR(hMRR2_W_Multijet,hMRR2_W_TTST,hMRR2_W_Wjet,hMRR2_W_Zjet,hMRR2_W_DYToLL,hMRR2_W_Gjet,hMRR2_W_Multiboson,hMRR2_W_data,histname_W,njet,channel);
  DrawR2MR(hMRR2_Z_Multijet,hMRR2_Z_TTST,hMRR2_Z_Wjet,hMRR2_Z_Zjet,hMRR2_Z_DYToLL,hMRR2_Z_Gjet,hMRR2_Z_Multiboson,hMRR2_Z_data,histname_Z,njet,channel);
  DrawR2MR(hMRR2_G_Multijet,hMRR2_G_TTST,hMRR2_G_Wjet,hMRR2_G_Zjet,hMRR2_G_DYToLL,hMRR2_G_Gjet,hMRR2_G_Multiboson,hMRR2_G_data,histname_G,njet,channel);
  DrawR2MR(hMRR2_P_Multijet,hMRR2_P_TTST,hMRR2_P_Wjet,hMRR2_P_Zjet,hMRR2_P_DYToLL,hMRR2_P_Gjet,hMRR2_P_Multiboson,hMRR2_P_data,histname_P,njet,channel);
}

void unrolled(){
  ClosureTest("","");
  ClosureTest("_nj45","");
  ClosureTest("_nj6","");
  //ClosureTest("","TOP_");
}
