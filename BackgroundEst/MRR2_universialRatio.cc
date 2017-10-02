void ClosureTest(TString njet, TString channel){
  cout << njet + channel << endl;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();
  TString dir = "/uscms/homes/c/chuh/work/susy170922/";
  TFile* file1 = TFile::Open(dir+"Multijet.root");
  TFile* file2 = TFile::Open(dir+"TTST.root");
  TFile* file3 = TFile::Open(dir+"WJet.root");
  TFile* file4 = TFile::Open(dir+"ZJDY.root");
  TFile* file5 = TFile::Open(dir+"Multiboson.root");
  TFile* file6 = TFile::Open(dir+"GJet.root");
  TFile* file7 = TFile::Open(dir+"data.root");
  TFile* sfile1 = TFile::Open(dir+"T2tt_mStop_850_mLSP_100.root");
  TFile* sfile2 = TFile::Open(dir+"T5tttt_mGluino_1400_mLSP_300.root");
  TFile* sfile3 = TFile::Open(dir+"T5ttcc_mGluino_1400_mLSP_300.root");
  TFile* sfile4 = TFile::Open(dir+"T1tttt_mGluino_1400_mLSP_300.root");

  TString histname_S = "R2_MR"+channel+"_S"+njet;
  TString histname_Q = "R2_MR"+channel+"_Q"+njet;
  TString histname_s = "R2_MR"+channel+"_s"+njet;
  TString histname_q = "R2_MR"+channel+"_q"+njet;
  TString histname_T = "R2_MR"+channel+"_T"+njet;
  TString histname_W = "R2_MR"+channel+"_W"+njet;
  TString histname_Z = "R2_MR"+channel+"_Z"+njet;
  TString histname_G = "R2_MR"+channel+"_G"+njet;
  TString histname_PS = "R2_MR";
  if(channel=="_topana") histname_G = histname_S;

  TH2D*   Multijet_S = (TH2D*)file1->Get(histname_S);
  TH2D*       TTST_S = (TH2D*)file2->Get(histname_S);
  TH2D*       Wjet_S = (TH2D*)file3->Get(histname_S);
  TH2D*       Zjet_S = (TH2D*)file4->Get(histname_S);
  TH2D* Multiboson_S = (TH2D*)file5->Get(histname_S);
  TH2D*       Gjet_S = (TH2D*)file6->Get(histname_S);
  TH2D*       data_S = (TH2D*)file7->Get(histname_S);
  TH2D*   		T2tt_S = (TH2D*)sfile1->Get(histname_S);
  TH2D*   	T5tttt_S = (TH2D*)sfile2->Get(histname_S);
  TH2D*   	T5ttcc_S = (TH2D*)sfile3->Get(histname_S);
  TH2D*   	T1tttt_S = (TH2D*)sfile4->Get(histname_S);
  TH2D*   Multijet_s = (TH2D*)file1->Get(histname_s);
  TH2D*       TTST_s = (TH2D*)file2->Get(histname_s);
  TH2D*       Wjet_s = (TH2D*)file3->Get(histname_s);
  TH2D*       Zjet_s = (TH2D*)file4->Get(histname_s);
  TH2D* Multiboson_s = (TH2D*)file5->Get(histname_s);
  TH2D*       Gjet_s = (TH2D*)file6->Get(histname_s);
  TH2D*       data_s = (TH2D*)file7->Get(histname_s);
  TH2D*   		T2tt_s = (TH2D*)sfile1->Get(histname_s);
  TH2D*   	T5tttt_s = (TH2D*)sfile2->Get(histname_s);
  TH2D*   	T5ttcc_s = (TH2D*)sfile3->Get(histname_s);
  TH2D*   	T1tttt_s = (TH2D*)sfile4->Get(histname_s);
  TH2D*   Multijet_Q = (TH2D*)file1->Get(histname_Q);
  TH2D*       TTST_Q = (TH2D*)file2->Get(histname_Q);
  TH2D*       Wjet_Q = (TH2D*)file3->Get(histname_Q);
  TH2D*       Zjet_Q = (TH2D*)file4->Get(histname_Q);
  TH2D* Multiboson_Q = (TH2D*)file5->Get(histname_Q);
  TH2D*       Gjet_Q = (TH2D*)file6->Get(histname_Q);
  TH2D*       data_Q = (TH2D*)file7->Get(histname_Q);
  TH2D*   		T2tt_Q = (TH2D*)sfile1->Get(histname_Q);
  TH2D*   	T5tttt_Q = (TH2D*)sfile2->Get(histname_Q);
  TH2D*   	T5ttcc_Q = (TH2D*)sfile3->Get(histname_Q);
  TH2D*   	T1tttt_Q = (TH2D*)sfile4->Get(histname_Q);
  TH2D*   Multijet_q = (TH2D*)file1->Get(histname_q);
  TH2D*       TTST_q = (TH2D*)file2->Get(histname_q);
  TH2D*       Wjet_q = (TH2D*)file3->Get(histname_q);
  TH2D*       Zjet_q = (TH2D*)file4->Get(histname_q);
  TH2D* Multiboson_q = (TH2D*)file5->Get(histname_q);
  TH2D*       Gjet_q = (TH2D*)file6->Get(histname_q);
  TH2D*       data_q = (TH2D*)file7->Get(histname_q);
  TH2D*   		T2tt_q = (TH2D*)sfile1->Get(histname_q);
  TH2D*   	T5tttt_q = (TH2D*)sfile2->Get(histname_q);
  TH2D*   	T5ttcc_q = (TH2D*)sfile3->Get(histname_q);
  TH2D*   	T1tttt_q = (TH2D*)sfile4->Get(histname_q);
  TH2D*   Multijet_T = (TH2D*)file1->Get(histname_T);
  TH2D*       TTST_T = (TH2D*)file2->Get(histname_T);
  TH2D*       Wjet_T = (TH2D*)file3->Get(histname_T);
  TH2D*       Zjet_T = (TH2D*)file4->Get(histname_T);
  TH2D* Multiboson_T = (TH2D*)file5->Get(histname_T);
  TH2D*       Gjet_T = (TH2D*)file6->Get(histname_T);
  TH2D*       data_T = (TH2D*)file7->Get(histname_T);
  TH2D*   		T2tt_T = (TH2D*)sfile1->Get(histname_T);
  TH2D*   	T5tttt_T = (TH2D*)sfile2->Get(histname_T);
  TH2D*   	T5ttcc_T = (TH2D*)sfile3->Get(histname_T);
  TH2D*   	T1tttt_T = (TH2D*)sfile4->Get(histname_T);
  TH2D*   Multijet_W = (TH2D*)file1->Get(histname_W);
  TH2D*       TTST_W = (TH2D*)file2->Get(histname_W);
  TH2D*       Wjet_W = (TH2D*)file3->Get(histname_W);
  TH2D*       Zjet_W = (TH2D*)file4->Get(histname_W);
  TH2D* Multiboson_W = (TH2D*)file5->Get(histname_W);
  TH2D*				Gjet_W = (TH2D*)file6->Get(histname_W);
  TH2D*       data_W = (TH2D*)file7->Get(histname_W);
  TH2D*   		T2tt_W = (TH2D*)sfile1->Get(histname_W);
  TH2D*   	T5tttt_W = (TH2D*)sfile2->Get(histname_W);
  TH2D*   	T5ttcc_W = (TH2D*)sfile3->Get(histname_W);
  TH2D*   	T1tttt_W = (TH2D*)sfile4->Get(histname_W);
  TH2D*   Multijet_Z = (TH2D*)file1->Get(histname_Z);
  TH2D*       TTST_Z = (TH2D*)file2->Get(histname_Z);
  TH2D*       Wjet_Z = (TH2D*)file3->Get(histname_Z);
  TH2D*       Zjet_Z = (TH2D*)file4->Get(histname_Z);
  TH2D* Multiboson_Z = (TH2D*)file5->Get(histname_Z);
  TH2D*				Gjet_Z = (TH2D*)file6->Get(histname_Z);
  TH2D*       data_Z = (TH2D*)file7->Get(histname_Z);
  TH2D*   		T2tt_Z = (TH2D*)sfile1->Get(histname_Z);
  TH2D*   	T5tttt_Z = (TH2D*)sfile2->Get(histname_Z);
  TH2D*   	T5ttcc_Z = (TH2D*)sfile3->Get(histname_Z);
  TH2D*   	T1tttt_Z = (TH2D*)sfile4->Get(histname_Z);
  TH2D*   Multijet_G = (TH2D*)file1->Get(histname_G);
  TH2D*       TTST_G = (TH2D*)file2->Get(histname_G);
  TH2D*       Wjet_G = (TH2D*)file3->Get(histname_G);
  TH2D*       Zjet_G = (TH2D*)file4->Get(histname_G);
  TH2D* Multiboson_G = (TH2D*)file5->Get(histname_G);
  TH2D*				Gjet_G = (TH2D*)file6->Get(histname_G);
  TH2D*       data_G = (TH2D*)file7->Get(histname_G);
  TH2D*   		T2tt_G = (TH2D*)sfile1->Get(histname_G);
  TH2D*   	T5tttt_G = (TH2D*)sfile2->Get(histname_G);
  TH2D*   	T5ttcc_G = (TH2D*)sfile3->Get(histname_G);
  TH2D*   	T1tttt_G = (TH2D*)sfile4->Get(histname_G);
  TH2D*   Multijet_PS = (TH2D*)file1->Get(histname_PS);
  TH2D*       TTST_PS = (TH2D*)file2->Get(histname_PS);
  TH2D*       Wjet_PS = (TH2D*)file3->Get(histname_PS);
  TH2D*       Zjet_PS = (TH2D*)file4->Get(histname_PS);
  TH2D* Multiboson_PS = (TH2D*)file5->Get(histname_PS);
  TH2D*				Gjet_PS = (TH2D*)file6->Get(histname_PS);
  TH2D*       data_PS = (TH2D*)file7->Get(histname_PS);
  TH2D*   		T2tt_PS = (TH2D*)sfile1->Get(histname_PS);
  TH2D*   	T5tttt_PS = (TH2D*)sfile2->Get(histname_PS);
  TH2D*   	T5ttcc_PS = (TH2D*)sfile3->Get(histname_PS);
  TH2D*   	T1tttt_PS = (TH2D*)sfile4->Get(histname_PS);

  TH1F* hMRR2_S_Multijet = new TH1F("MRR2_S_Multijet", "", 25,0,25);
  TH1F* hMRR2_S_TTST = new TH1F("MRR2_S_TTST", "", 25,0,25);
  TH1F* hMRR2_S_Wjet = new TH1F("MRR2_S_Wjet", "", 25,0,25);
  TH1F* hMRR2_S_Zjet = new TH1F("MRR2_S_Zjet", "", 25,0,25);
  TH1F* hMRR2_S_Multiboson = new TH1F("MRR2_S_Multiboson", "", 25,0,25);
  TH1F* hMRR2_S_Gjet = new TH1F("MRR2_S_Gjet", "", 25,0,25);
  TH1F* hMRR2_S_data = new TH1F("MRR2_S_data", "", 25,0,25);
  TH1F* hMRR2_S_T2tt = new TH1F("MRR2_S_T2tt", "", 25,0,25);
  TH1F* hMRR2_S_T5tttt = new TH1F("MRR2_S_T5tttt", "", 25,0,25);
  TH1F* hMRR2_S_T5ttcc = new TH1F("MRR2_S_T5ttcc", "", 25,0,25);
  TH1F* hMRR2_S_T1tttt = new TH1F("MRR2_S_T1tttt", "", 25,0,25);

  TH1F* hMRR2_s_Multijet = new TH1F("MRR2_s_Multijet", "", 25,0,25);
  TH1F* hMRR2_s_TTST = new TH1F("MRR2_s_TTST", "", 25,0,25);
  TH1F* hMRR2_s_Wjet = new TH1F("MRR2_s_Wjet", "", 25,0,25);
  TH1F* hMRR2_s_Zjet = new TH1F("MRR2_s_Zjet", "", 25,0,25);
  TH1F* hMRR2_s_Multiboson = new TH1F("MRR2_s_Multiboson", "", 25,0,25);
  TH1F* hMRR2_s_Gjet = new TH1F("MRR2_s_Gjet", "", 25,0,25);
  TH1F* hMRR2_s_data = new TH1F("MRR2_s_data", "", 25,0,25);
  TH1F* hMRR2_s_T2tt = new TH1F("MRR2_s_T2tt", "", 25,0,25);
  TH1F* hMRR2_s_T5tttt = new TH1F("MRR2_s_T5tttt", "", 25,0,25);
  TH1F* hMRR2_s_T5ttcc = new TH1F("MRR2_s_T5ttcc", "", 25,0,25);
  TH1F* hMRR2_s_T1tttt = new TH1F("MRR2_s_T1tttt", "", 25,0,25);

  TH1F* hMRR2_Q_Multijet = new TH1F("MRR2_Q_Multijet", "", 25,0,25);
  TH1F* hMRR2_Q_TTST = new TH1F("MRR2_Q_TTST", "", 25,0,25);
  TH1F* hMRR2_Q_Wjet = new TH1F("MRR2_Q_Wjet", "", 25,0,25);
  TH1F* hMRR2_Q_Zjet = new TH1F("MRR2_Q_Zjet", "", 25,0,25);
  TH1F* hMRR2_Q_Multiboson = new TH1F("MRR2_Q_Multiboson", "", 25,0,25);
  TH1F* hMRR2_Q_Gjet = new TH1F("MRR2_Q_Gjet", "", 25,0,25);
  TH1F* hMRR2_Q_data = new TH1F("MRR2_Q_data", "", 25,0,25);
  TH1F* hMRR2_Q_T2tt = new TH1F("MRR2_Q_T2tt", "", 25,0,25);
  TH1F* hMRR2_Q_T5tttt = new TH1F("MRR2_Q_T5tttt", "", 25,0,25);
  TH1F* hMRR2_Q_T5ttcc = new TH1F("MRR2_Q_T5ttcc", "", 25,0,25);
  TH1F* hMRR2_Q_T1tttt = new TH1F("MRR2_Q_T1tttt", "", 25,0,25);

  TH1F* hMRR2_q_Multijet = new TH1F("MRR2_q_Multijet", "", 25,0,25);
  TH1F* hMRR2_q_TTST = new TH1F("MRR2_q_TTST", "", 25,0,25);
  TH1F* hMRR2_q_Wjet = new TH1F("MRR2_q_Wjet", "", 25,0,25);
  TH1F* hMRR2_q_Zjet = new TH1F("MRR2_q_Zjet", "", 25,0,25);
  TH1F* hMRR2_q_Multiboson = new TH1F("MRR2_q_Multiboson", "", 25,0,25);
  TH1F* hMRR2_q_Gjet = new TH1F("MRR2_q_Gjet", "", 25,0,25);
  TH1F* hMRR2_q_data = new TH1F("MRR2_q_data", "", 25,0,25);
  TH1F* hMRR2_q_T2tt = new TH1F("MRR2_q_T2tt", "", 25,0,25);
  TH1F* hMRR2_q_T5tttt = new TH1F("MRR2_q_T5tttt", "", 25,0,25);
  TH1F* hMRR2_q_T5ttcc = new TH1F("MRR2_q_T5ttcc", "", 25,0,25);
  TH1F* hMRR2_q_T1tttt = new TH1F("MRR2_q_T1tttt", "", 25,0,25);

  TH1F* hMRR2_T_Multijet = new TH1F("MRR2_T_Multijet", "", 25,0,25);
  TH1F* hMRR2_T_TTST = new TH1F("MRR2_T_TTST", "", 25,0,25);
  TH1F* hMRR2_T_Wjet = new TH1F("MRR2_T_Wjet", "", 25,0,25);
  TH1F* hMRR2_T_Zjet = new TH1F("MRR2_T_Zjet", "", 25,0,25);
  TH1F* hMRR2_T_Multiboson = new TH1F("MRR2_T_Multiboson", "", 25,0,25);
  TH1F* hMRR2_T_Gjet = new TH1F("MRR2_T_Gjet", "", 25,0,25);
  TH1F* hMRR2_T_data = new TH1F("MRR2_T_data", "", 25,0,25);
  TH1F* hMRR2_T_T2tt = new TH1F("MRR2_T_T2tt", "", 25,0,25);
  TH1F* hMRR2_T_T5tttt = new TH1F("MRR2_T_T5tttt", "", 25,0,25);
  TH1F* hMRR2_T_T5ttcc = new TH1F("MRR2_T_T5ttcc", "", 25,0,25);
  TH1F* hMRR2_T_T1tttt = new TH1F("MRR2_T_T1tttt", "", 25,0,25);

  TH1F* hMRR2_W_Multijet = new TH1F("MRR2_W_Multijet", "", 25,0,25);
  TH1F* hMRR2_W_TTST = new TH1F("MRR2_W_TTST", "", 25,0,25);
  TH1F* hMRR2_W_Wjet = new TH1F("MRR2_W_Wjet", "", 25,0,25);
  TH1F* hMRR2_W_Zjet = new TH1F("MRR2_W_Zjet", "", 25,0,25);
  TH1F* hMRR2_W_Multiboson = new TH1F("MRR2_W_Multiboson", "", 25,0,25);
  TH1F* hMRR2_W_Gjet = new TH1F("MRR2_W_Gjet", "", 25,0,25);
  TH1F* hMRR2_W_data = new TH1F("MRR2_W_data", "", 25,0,25);
  TH1F* hMRR2_W_T2tt = new TH1F("MRR2_W_T2tt", "", 25,0,25);
  TH1F* hMRR2_W_T5tttt = new TH1F("MRR2_W_T5tttt", "", 25,0,25);
  TH1F* hMRR2_W_T5ttcc = new TH1F("MRR2_W_T5ttcc", "", 25,0,25);
  TH1F* hMRR2_W_T1tttt = new TH1F("MRR2_W_T1tttt", "", 25,0,25);

  TH1F* hMRR2_Z_Multijet = new TH1F("MRR2_Z_Multijet", "", 25,0,25);
  TH1F* hMRR2_Z_TTST = new TH1F("MRR2_Z_TTST", "", 25,0,25);
  TH1F* hMRR2_Z_Wjet = new TH1F("MRR2_Z_Wjet", "", 25,0,25);
  TH1F* hMRR2_Z_Zjet = new TH1F("MRR2_Z_Zjet", "", 25,0,25);
  TH1F* hMRR2_Z_Multiboson = new TH1F("MRR2_Z_Multiboson", "", 25,0,25);
  TH1F* hMRR2_Z_Gjet = new TH1F("MRR2_Z_Gjet", "", 25,0,25);
  TH1F* hMRR2_Z_data = new TH1F("MRR2_Z_data", "", 25,0,25);
  TH1F* hMRR2_Z_T2tt = new TH1F("MRR2_Z_T2tt", "", 25,0,25);
  TH1F* hMRR2_Z_T5tttt = new TH1F("MRR2_Z_T5tttt", "", 25,0,25);
  TH1F* hMRR2_Z_T5ttcc = new TH1F("MRR2_Z_T5ttcc", "", 25,0,25);
  TH1F* hMRR2_Z_T1tttt = new TH1F("MRR2_Z_T1tttt", "", 25,0,25);

  TH1F* hMRR2_G_Multijet = new TH1F("MRR2_G_Multijet", "", 25,0,25);
  TH1F* hMRR2_G_TTST = new TH1F("MRR2_G_TTST", "", 25,0,25);
  TH1F* hMRR2_G_Wjet = new TH1F("MRR2_G_Wjet", "", 25,0,25);
  TH1F* hMRR2_G_Zjet = new TH1F("MRR2_G_Zjet", "", 25,0,25);
  TH1F* hMRR2_G_Multiboson = new TH1F("MRR2_G_Multiboson", "", 25,0,25);
  TH1F* hMRR2_G_Gjet = new TH1F("MRR2_G_Gjet", "", 25,0,25);
  TH1F* hMRR2_G_data = new TH1F("MRR2_G_data", "", 25,0,25);
  TH1F* hMRR2_G_T2tt = new TH1F("MRR2_G_T2tt", "", 25,0,25);
  TH1F* hMRR2_G_T5tttt = new TH1F("MRR2_G_T5tttt", "", 25,0,25);
  TH1F* hMRR2_G_T5ttcc = new TH1F("MRR2_G_T5ttcc", "", 25,0,25);
  TH1F* hMRR2_G_T1tttt = new TH1F("MRR2_G_T1tttt", "", 25,0,25);

  TH1F* hMRR2_PS_Multijet = new TH1F("MRR2_PS_Multijet", "", 25,0,25);
  TH1F* hMRR2_PS_TTST = new TH1F("MRR2_PS_TTST", "", 25,0,25);
  TH1F* hMRR2_PS_Wjet = new TH1F("MRR2_PS_Wjet", "", 25,0,25);
  TH1F* hMRR2_PS_Zjet = new TH1F("MRR2_PS_Zjet", "", 25,0,25);
  TH1F* hMRR2_PS_Multiboson = new TH1F("MRR2_PS_Multiboson", "", 25,0,25);
  TH1F* hMRR2_PS_Gjet = new TH1F("MRR2_PS_Gjet", "", 25,0,25);
  TH1F* hMRR2_PS_data = new TH1F("MRR2_PS_data", "", 25,0,25);
  TH1F* hMRR2_PS_T2tt = new TH1F("MRR2_PS_T2tt", "", 25,0,25);
  TH1F* hMRR2_PS_T5tttt = new TH1F("MRR2_PS_T5tttt", "", 25,0,25);
  TH1F* hMRR2_PS_T5ttcc = new TH1F("MRR2_PS_T5ttcc", "", 25,0,25);
  TH1F* hMRR2_PS_T1tttt = new TH1F("MRR2_PS_T1tttt", "", 25,0,25);

  //hMRR2_S_data->GetXaxis()->SetLabelSize(0.01);
  //hMRR2_S_data->GetXaxis()->LabelsOption("v");

  const char *label[26] = {"[0.08, 0.12]", "[0.12, 0.16]", "[0.16, 0.24]", "[0.24, 0.50]", "[0.50, 1.00]", "[0.08, 0.12]", "[0.12, 0.16]", "[0.16, 0.24]", "[0.24, 0.50]", "[0.50, 1.00]", "[0.08, 0.12]", "[0.12, 0.16]", "[0.16, 0.24]", "[0.24, 0.50]", "[0.50, 1.00]", "[0.08, 0.12]", "[0.12, 0.16]", "[0.16, 0.24]", "[0.24, 0.50]", "[0.50, 1.00]", "[0.08, 0.12]", "[0.12, 0.16]", "[0.16, 0.24]", "[0.24, 0.50]", "[0.50, 1.00]"};

  int count=1;
  for(int i=3;i<=data_S->GetNbinsX();i++){
    for(int j=3;j<=data_S->GetNbinsY();j++){
      hMRR2_S_Multijet->SetBinContent(count,Multijet_S->GetBinContent(i,j));
      hMRR2_S_Multijet->SetBinError(count,Multijet_S->GetBinError(i,j));
      hMRR2_S_TTST->SetBinContent(count,TTST_S->GetBinContent(i,j));
      hMRR2_S_TTST->SetBinError(count,TTST_S->GetBinError(i,j));
      hMRR2_S_Wjet->SetBinContent(count,Wjet_S->GetBinContent(i,j));
      hMRR2_S_Wjet->SetBinError(count,Wjet_S->GetBinError(i,j));
      hMRR2_S_Zjet->SetBinContent(count,Zjet_S->GetBinContent(i,j));
      hMRR2_S_Zjet->SetBinError(count,Zjet_S->GetBinError(i,j));
      hMRR2_S_Multiboson->SetBinContent(count,Multiboson_S->GetBinContent(i,j));
      hMRR2_S_Multiboson->SetBinError(count,Multiboson_S->GetBinError(i,j));
      hMRR2_S_Gjet->SetBinContent(count,Gjet_S->GetBinContent(i,j));
      hMRR2_S_Gjet->SetBinError(count,Gjet_S->GetBinError(i,j));
      hMRR2_S_data->SetBinContent(count,data_S->GetBinContent(i,j));
      hMRR2_S_data->SetBinError(count,data_S->GetBinError(i,j));
      hMRR2_S_T2tt->SetBinContent(count,T2tt_S->GetBinContent(i,j));
      hMRR2_S_T2tt->SetBinError(count,T2tt_S->GetBinError(i,j));
      hMRR2_S_T5tttt->SetBinContent(count,T5tttt_S->GetBinContent(i,j));
      hMRR2_S_T5tttt->SetBinError(count,T5tttt_S->GetBinError(i,j));
      hMRR2_S_T5ttcc->SetBinContent(count,T5ttcc_S->GetBinContent(i,j));
      hMRR2_S_T5ttcc->SetBinError(count,T5ttcc_S->GetBinError(i,j));
      hMRR2_S_T1tttt->SetBinContent(count,T1tttt_S->GetBinContent(i,j));
      hMRR2_S_T1tttt->SetBinError(count,T1tttt_S->GetBinError(i,j));
      hMRR2_s_Multijet->SetBinContent(count,Multijet_s->GetBinContent(i,j));
      hMRR2_s_Multijet->SetBinError(count,Multijet_s->GetBinError(i,j));
      hMRR2_s_TTST->SetBinContent(count,TTST_s->GetBinContent(i,j));
      hMRR2_s_TTST->SetBinError(count,TTST_s->GetBinError(i,j));
      hMRR2_s_Wjet->SetBinContent(count,Wjet_s->GetBinContent(i,j));
      hMRR2_s_Wjet->SetBinError(count,Wjet_s->GetBinError(i,j));
      hMRR2_s_Zjet->SetBinContent(count,Zjet_s->GetBinContent(i,j));
      hMRR2_s_Zjet->SetBinError(count,Zjet_s->GetBinError(i,j));
      hMRR2_s_Multiboson->SetBinContent(count,Multiboson_s->GetBinContent(i,j));
      hMRR2_s_Multiboson->SetBinError(count,Multiboson_s->GetBinError(i,j));
      hMRR2_s_Gjet->SetBinContent(count,Gjet_s->GetBinContent(i,j));
      hMRR2_s_Gjet->SetBinError(count,Gjet_s->GetBinError(i,j));
      hMRR2_s_data->SetBinContent(count,data_s->GetBinContent(i,j));
      hMRR2_s_data->SetBinError(count,data_s->GetBinError(i,j));
      hMRR2_s_T2tt->SetBinContent(count,T2tt_s->GetBinContent(i,j));
      hMRR2_s_T2tt->SetBinError(count,T2tt_s->GetBinError(i,j));
      hMRR2_s_T5tttt->SetBinContent(count,T5tttt_s->GetBinContent(i,j));
      hMRR2_s_T5tttt->SetBinError(count,T5tttt_s->GetBinError(i,j));
      hMRR2_s_T5ttcc->SetBinContent(count,T5ttcc_s->GetBinContent(i,j));
      hMRR2_s_T5ttcc->SetBinError(count,T5ttcc_s->GetBinError(i,j));
      hMRR2_s_T1tttt->SetBinContent(count,T1tttt_s->GetBinContent(i,j));
      hMRR2_s_T1tttt->SetBinError(count,T1tttt_s->GetBinError(i,j));
      hMRR2_Q_Multijet->SetBinContent(count,Multijet_Q->GetBinContent(i,j));
      hMRR2_Q_Multijet->SetBinError(count,Multijet_Q->GetBinError(i,j));
      hMRR2_Q_TTST->SetBinContent(count,TTST_Q->GetBinContent(i,j));
      hMRR2_Q_TTST->SetBinError(count,TTST_Q->GetBinError(i,j));
      hMRR2_Q_Wjet->SetBinContent(count,Wjet_Q->GetBinContent(i,j));
      hMRR2_Q_Wjet->SetBinError(count,Wjet_Q->GetBinError(i,j));
      hMRR2_Q_Zjet->SetBinContent(count,Zjet_Q->GetBinContent(i,j));
      hMRR2_Q_Zjet->SetBinError(count,Zjet_Q->GetBinError(i,j));
      hMRR2_Q_Multiboson->SetBinContent(count,Multiboson_Q->GetBinContent(i,j));
      hMRR2_Q_Multiboson->SetBinError(count,Multiboson_Q->GetBinError(i,j));
      hMRR2_Q_Gjet->SetBinContent(count,Gjet_Q->GetBinContent(i,j));
      hMRR2_Q_Gjet->SetBinError(count,Gjet_Q->GetBinError(i,j));
      hMRR2_Q_data->SetBinContent(count,data_Q->GetBinContent(i,j));
      hMRR2_Q_data->SetBinError(count,data_Q->GetBinError(i,j));
      hMRR2_Q_T2tt->SetBinContent(count,T2tt_Q->GetBinContent(i,j));
      hMRR2_Q_T2tt->SetBinError(count,T2tt_Q->GetBinError(i,j));
      hMRR2_Q_T5tttt->SetBinContent(count,T5tttt_Q->GetBinContent(i,j));
      hMRR2_Q_T5tttt->SetBinError(count,T5tttt_Q->GetBinError(i,j));
      hMRR2_Q_T5ttcc->SetBinContent(count,T5ttcc_Q->GetBinContent(i,j));
      hMRR2_Q_T5ttcc->SetBinError(count,T5ttcc_Q->GetBinError(i,j));
      hMRR2_Q_T1tttt->SetBinContent(count,T1tttt_Q->GetBinContent(i,j));
      hMRR2_Q_T1tttt->SetBinError(count,T1tttt_Q->GetBinError(i,j));
      hMRR2_q_Multijet->SetBinContent(count,Multijet_q->GetBinContent(i,j));
      hMRR2_q_Multijet->SetBinError(count,Multijet_q->GetBinError(i,j));
      hMRR2_q_TTST->SetBinContent(count,TTST_q->GetBinContent(i,j));
      hMRR2_q_TTST->SetBinError(count,TTST_q->GetBinError(i,j));
      hMRR2_q_Wjet->SetBinContent(count,Wjet_q->GetBinContent(i,j));
      hMRR2_q_Wjet->SetBinError(count,Wjet_q->GetBinError(i,j));
      hMRR2_q_Zjet->SetBinContent(count,Zjet_q->GetBinContent(i,j));
      hMRR2_q_Zjet->SetBinError(count,Zjet_q->GetBinError(i,j));
      hMRR2_q_Multiboson->SetBinContent(count,Multiboson_q->GetBinContent(i,j));
      hMRR2_q_Multiboson->SetBinError(count,Multiboson_q->GetBinError(i,j));
      hMRR2_q_Gjet->SetBinContent(count,Gjet_q->GetBinContent(i,j));
      hMRR2_q_Gjet->SetBinError(count,Gjet_q->GetBinError(i,j));
      hMRR2_q_data->SetBinContent(count,data_q->GetBinContent(i,j));
      hMRR2_q_data->SetBinError(count,data_q->GetBinError(i,j));
      hMRR2_q_T2tt->SetBinContent(count,T2tt_q->GetBinContent(i,j));
      hMRR2_q_T2tt->SetBinError(count,T2tt_q->GetBinError(i,j));
      hMRR2_q_T5tttt->SetBinContent(count,T5tttt_q->GetBinContent(i,j));
      hMRR2_q_T5tttt->SetBinError(count,T5tttt_q->GetBinError(i,j));
      hMRR2_q_T5ttcc->SetBinContent(count,T5ttcc_q->GetBinContent(i,j));
      hMRR2_q_T5ttcc->SetBinError(count,T5ttcc_q->GetBinError(i,j));
      hMRR2_q_T1tttt->SetBinContent(count,T1tttt_q->GetBinContent(i,j));
      hMRR2_q_T1tttt->SetBinError(count,T1tttt_q->GetBinError(i,j));
      hMRR2_T_Multijet->SetBinContent(count,Multijet_T->GetBinContent(i,j));
      hMRR2_T_Multijet->SetBinError(count,Multijet_T->GetBinError(i,j));
      hMRR2_T_TTST->SetBinContent(count,TTST_T->GetBinContent(i,j));
      hMRR2_T_TTST->SetBinError(count,TTST_T->GetBinError(i,j));
      hMRR2_T_Wjet->SetBinContent(count,Wjet_T->GetBinContent(i,j));
      hMRR2_T_Wjet->SetBinError(count,Wjet_T->GetBinError(i,j));
      hMRR2_T_Zjet->SetBinContent(count,Zjet_T->GetBinContent(i,j));
      hMRR2_T_Zjet->SetBinError(count,Zjet_T->GetBinError(i,j));
      hMRR2_T_Multiboson->SetBinContent(count,Multiboson_T->GetBinContent(i,j));
      hMRR2_T_Multiboson->SetBinError(count,Multiboson_T->GetBinError(i,j));
      hMRR2_T_Gjet->SetBinContent(count,Gjet_T->GetBinContent(i,j));
      hMRR2_T_Gjet->SetBinError(count,Gjet_T->GetBinError(i,j));
      hMRR2_T_data->SetBinContent(count,data_T->GetBinContent(i,j));
      hMRR2_T_data->SetBinError(count,data_T->GetBinError(i,j));
      hMRR2_T_T2tt->SetBinContent(count,T2tt_T->GetBinContent(i,j));
      hMRR2_T_T2tt->SetBinError(count,T2tt_T->GetBinError(i,j));
      hMRR2_T_T5tttt->SetBinContent(count,T5tttt_T->GetBinContent(i,j));
      hMRR2_T_T5tttt->SetBinError(count,T5tttt_T->GetBinError(i,j));
      hMRR2_T_T5ttcc->SetBinContent(count,T5ttcc_T->GetBinContent(i,j));
      hMRR2_T_T5ttcc->SetBinError(count,T5ttcc_T->GetBinError(i,j));
      hMRR2_T_T1tttt->SetBinContent(count,T1tttt_T->GetBinContent(i,j));
      hMRR2_T_T1tttt->SetBinError(count,T1tttt_T->GetBinError(i,j));
      hMRR2_W_Multijet->SetBinContent(count,Multijet_W->GetBinContent(i,j));
      hMRR2_W_Multijet->SetBinError(count,Multijet_W->GetBinError(i,j));
      hMRR2_W_TTST->SetBinContent(count,TTST_W->GetBinContent(i,j));
      hMRR2_W_TTST->SetBinError(count,TTST_W->GetBinError(i,j));
      hMRR2_W_Wjet->SetBinContent(count,Wjet_W->GetBinContent(i,j));
      hMRR2_W_Wjet->SetBinError(count,Wjet_W->GetBinError(i,j));
      hMRR2_W_Zjet->SetBinContent(count,Zjet_W->GetBinContent(i,j));
      hMRR2_W_Zjet->SetBinError(count,Zjet_W->GetBinError(i,j));
      hMRR2_W_Multiboson->SetBinContent(count,Multiboson_W->GetBinContent(i,j));
      hMRR2_W_Multiboson->SetBinError(count,Multiboson_W->GetBinError(i,j));
      hMRR2_W_Gjet->SetBinContent(count,Gjet_W->GetBinContent(i,j));
      hMRR2_W_Gjet->SetBinError(count,Gjet_W->GetBinError(i,j));
      hMRR2_W_data->SetBinContent(count,data_W->GetBinContent(i,j));
      hMRR2_W_data->SetBinError(count,data_W->GetBinError(i,j));
      hMRR2_W_T2tt->SetBinContent(count,T2tt_W->GetBinContent(i,j));
      hMRR2_W_T2tt->SetBinError(count,T2tt_W->GetBinError(i,j));
      hMRR2_W_T5tttt->SetBinContent(count,T5tttt_W->GetBinContent(i,j));
      hMRR2_W_T5tttt->SetBinError(count,T5tttt_W->GetBinError(i,j));
      hMRR2_W_T5ttcc->SetBinContent(count,T5ttcc_W->GetBinContent(i,j));
      hMRR2_W_T5ttcc->SetBinError(count,T5ttcc_W->GetBinError(i,j));
      hMRR2_W_T1tttt->SetBinContent(count,T1tttt_W->GetBinContent(i,j));
      hMRR2_W_T1tttt->SetBinError(count,T1tttt_W->GetBinError(i,j));
      hMRR2_Z_Multijet->SetBinContent(count,Multijet_Z->GetBinContent(i,j));
      hMRR2_Z_Multijet->SetBinError(count,Multijet_Z->GetBinError(i,j));
      hMRR2_Z_TTST->SetBinContent(count,TTST_Z->GetBinContent(i,j));
      hMRR2_Z_TTST->SetBinError(count,TTST_Z->GetBinError(i,j));
      hMRR2_Z_Wjet->SetBinContent(count,Wjet_Z->GetBinContent(i,j));
      hMRR2_Z_Wjet->SetBinError(count,Wjet_Z->GetBinError(i,j));
      hMRR2_Z_Zjet->SetBinContent(count,Zjet_Z->GetBinContent(i,j));
      hMRR2_Z_Zjet->SetBinError(count,Zjet_Z->GetBinError(i,j));
      hMRR2_Z_Multiboson->SetBinContent(count,Multiboson_Z->GetBinContent(i,j));
      hMRR2_Z_Multiboson->SetBinError(count,Multiboson_Z->GetBinError(i,j));
      hMRR2_Z_Gjet->SetBinContent(count,Gjet_Z->GetBinContent(i,j));
      hMRR2_Z_Gjet->SetBinError(count,Gjet_Z->GetBinError(i,j));
      hMRR2_Z_data->SetBinContent(count,data_Z->GetBinContent(i,j));
      hMRR2_Z_data->SetBinError(count,data_Z->GetBinError(i,j));
      hMRR2_Z_T2tt->SetBinContent(count,T2tt_Z->GetBinContent(i,j));
      hMRR2_Z_T2tt->SetBinError(count,T2tt_Z->GetBinError(i,j));
      hMRR2_Z_T5tttt->SetBinContent(count,T5tttt_Z->GetBinContent(i,j));
      hMRR2_Z_T5tttt->SetBinError(count,T5tttt_Z->GetBinError(i,j));
      hMRR2_Z_T5ttcc->SetBinContent(count,T5ttcc_Z->GetBinContent(i,j));
      hMRR2_Z_T5ttcc->SetBinError(count,T5ttcc_Z->GetBinError(i,j));
      hMRR2_Z_T1tttt->SetBinContent(count,T1tttt_Z->GetBinContent(i,j));
      hMRR2_Z_T1tttt->SetBinError(count,T1tttt_Z->GetBinError(i,j));
      hMRR2_G_Multijet->SetBinContent(count,Multijet_G->GetBinContent(i,j));
      hMRR2_G_Multijet->SetBinError(count,Multijet_G->GetBinError(i,j));
      hMRR2_G_TTST->SetBinContent(count,TTST_G->GetBinContent(i,j));
      hMRR2_G_TTST->SetBinError(count,TTST_G->GetBinError(i,j));
      hMRR2_G_Wjet->SetBinContent(count,Wjet_G->GetBinContent(i,j));
      hMRR2_G_Wjet->SetBinError(count,Wjet_G->GetBinError(i,j));
      hMRR2_G_Zjet->SetBinContent(count,Zjet_G->GetBinContent(i,j));
      hMRR2_G_Zjet->SetBinError(count,Zjet_G->GetBinError(i,j));
      hMRR2_G_Multiboson->SetBinContent(count,Multiboson_G->GetBinContent(i,j));
      hMRR2_G_Multiboson->SetBinError(count,Multiboson_G->GetBinError(i,j));
      hMRR2_G_Gjet->SetBinContent(count,Gjet_G->GetBinContent(i,j));
      hMRR2_G_Gjet->SetBinError(count,Gjet_G->GetBinError(i,j));
      hMRR2_G_data->SetBinContent(count,data_G->GetBinContent(i,j));
      hMRR2_G_data->SetBinError(count,data_G->GetBinError(i,j));
      hMRR2_G_T2tt->SetBinContent(count,T2tt_G->GetBinContent(i,j));
      hMRR2_G_T2tt->SetBinError(count,T2tt_G->GetBinError(i,j));
      hMRR2_G_T5tttt->SetBinContent(count,T5tttt_G->GetBinContent(i,j));
      hMRR2_G_T5tttt->SetBinError(count,T5tttt_G->GetBinError(i,j));
      hMRR2_G_T5ttcc->SetBinContent(count,T5ttcc_G->GetBinContent(i,j));
      hMRR2_G_T5ttcc->SetBinError(count,T5ttcc_G->GetBinError(i,j));
      hMRR2_G_T1tttt->SetBinContent(count,T1tttt_G->GetBinContent(i,j));
      hMRR2_G_T1tttt->SetBinError(count,T1tttt_G->GetBinError(i,j));
      hMRR2_PS_Multijet->SetBinContent(count,Multijet_PS->GetBinContent(i,j));
      hMRR2_PS_Multijet->SetBinError(count,Multijet_PS->GetBinError(i,j));
      hMRR2_PS_TTST->SetBinContent(count,TTST_PS->GetBinContent(i,j));
      hMRR2_PS_TTST->SetBinError(count,TTST_PS->GetBinError(i,j));
      hMRR2_PS_Wjet->SetBinContent(count,Wjet_PS->GetBinContent(i,j));
      hMRR2_PS_Wjet->SetBinError(count,Wjet_PS->GetBinError(i,j));
      hMRR2_PS_Zjet->SetBinContent(count,Zjet_PS->GetBinContent(i,j));
      hMRR2_PS_Zjet->SetBinError(count,Zjet_PS->GetBinError(i,j));
      hMRR2_PS_Multiboson->SetBinContent(count,Multiboson_PS->GetBinContent(i,j));
      hMRR2_PS_Multiboson->SetBinError(count,Multiboson_PS->GetBinError(i,j));
      hMRR2_PS_Gjet->SetBinContent(count,Gjet_PS->GetBinContent(i,j));
      hMRR2_PS_Gjet->SetBinError(count,Gjet_PS->GetBinError(i,j));
      hMRR2_PS_data->SetBinContent(count,data_PS->GetBinContent(i,j));
      hMRR2_PS_data->SetBinError(count,data_PS->GetBinError(i,j));
      hMRR2_PS_T2tt->SetBinContent(count,T2tt_PS->GetBinContent(i,j));
      hMRR2_PS_T2tt->SetBinError(count,T2tt_PS->GetBinError(i,j));
      hMRR2_PS_T5tttt->SetBinContent(count,T5tttt_PS->GetBinContent(i,j));
      hMRR2_PS_T5tttt->SetBinError(count,T5tttt_PS->GetBinError(i,j));
      hMRR2_PS_T5ttcc->SetBinContent(count,T5ttcc_PS->GetBinContent(i,j));
      hMRR2_PS_T5ttcc->SetBinError(count,T5ttcc_PS->GetBinError(i,j));
      hMRR2_PS_T1tttt->SetBinContent(count,T1tttt_PS->GetBinContent(i,j));
      hMRR2_PS_T1tttt->SetBinError(count,T1tttt_PS->GetBinError(i,j));
      hMRR2_S_data->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_s_data->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_Q_data->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_q_data->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_T_data->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_W_data->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_Z_data->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_G_data->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_PS_data->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_S_Multijet->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_s_Multijet->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_Q_Multijet->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_q_Multijet->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_T_Multijet->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_W_Multijet->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_Z_Multijet->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_G_Multijet->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_PS_Multijet->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_S_TTST->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_s_TTST->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_Q_TTST->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_q_TTST->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_T_TTST->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_W_TTST->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_Z_TTST->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_G_TTST->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_PS_TTST->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_S_Wjet->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_s_Wjet->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_Q_Wjet->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_q_Wjet->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_T_Wjet->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_W_Wjet->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_Z_Wjet->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_G_Wjet->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_PS_Wjet->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_S_Zjet->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_s_Zjet->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_Q_Zjet->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_q_Zjet->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_T_Zjet->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_W_Zjet->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_Z_Zjet->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_G_Zjet->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_PS_Zjet->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_S_Multiboson->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_s_Multiboson->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_Q_Multiboson->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_q_Multiboson->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_T_Multiboson->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_W_Multiboson->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_Z_Multiboson->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_G_Multiboson->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_PS_Multiboson->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_S_T2tt->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_s_T2tt->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_Q_T2tt->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_q_T2tt->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_T_T2tt->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_W_T2tt->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_Z_T2tt->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_G_T2tt->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_PS_T2tt->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_S_T5tttt->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_s_T5tttt->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_Q_T5tttt->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_q_T5tttt->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_T_T5tttt->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_W_T5tttt->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_Z_T5tttt->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_G_T5tttt->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_PS_T5tttt->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_S_T5ttcc->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_s_T5ttcc->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_Q_T5ttcc->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_q_T5ttcc->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_T_T5ttcc->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_W_T5ttcc->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_Z_T5ttcc->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_G_T5ttcc->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_PS_T5ttcc->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_S_T1tttt->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_s_T1tttt->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_Q_T1tttt->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_q_T1tttt->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_T_T1tttt->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_W_T1tttt->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_Z_T1tttt->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_G_T1tttt->GetXaxis()->SetBinLabel( count, label[count-1]);
      hMRR2_PS_T1tttt->GetXaxis()->SetBinLabel( count, label[count-1]);
      count++;
    }
  }

  TH1F* NT_QCD  = (TH1F*)hMRR2_S_data->Clone();
  TH1F* NSp_QCD = (TH1F*)hMRR2_S_data->Clone();
  TH1F* NSp_T   = (TH1F*)hMRR2_S_data->Clone();
  TH1F* NSp_W   = (TH1F*)hMRR2_S_data->Clone();
  TH1F* NSp_Z   = (TH1F*)hMRR2_S_data->Clone();
  TH1F* NSp_G   = (TH1F*)hMRR2_S_data->Clone();
  TH1F* NSp     = (TH1F*)hMRR2_S_data->Clone();
  TH1F* NQp     = (TH1F*)hMRR2_S_data->Clone();
  TH1F* NQp1    = (TH1F*)hMRR2_S_data->Clone();
  TH1F* NQp2    = (TH1F*)hMRR2_S_data->Clone();

  NT_QCD = (TH1F*)hMRR2_Q_data->Clone();
  NT_QCD->Add(hMRR2_Q_Wjet, -1);
  NT_QCD->Add(hMRR2_Q_Zjet, -1);
  NT_QCD->Add(hMRR2_Q_Gjet, -1);
  NT_QCD->Add(hMRR2_Q_TTST, -1);
  NT_QCD->Add(hMRR2_Q_Multiboson, -1);
  //NT_QCD->Multiply(hMRR2_T_Multijet);
  //NT_QCD->Divide(hMRR2_Q_Multijet);
  NT_QCD->Scale(hMRR2_T_Multijet->Integral()/hMRR2_Q_Multijet->Integral());
  cout << NT_QCD->Integral() << endl;
  
  NSp_QCD = (TH1F*)hMRR2_Q_data->Clone();
  NSp_QCD->Add(hMRR2_Q_Zjet, -1);
  NSp_QCD->Add(hMRR2_Q_Gjet, -1);
  NSp_QCD->Add(hMRR2_Q_Wjet, -1);
  NSp_QCD->Add(hMRR2_Q_TTST, -1);
  NSp_QCD->Add(hMRR2_Q_Multiboson, -1);
  cout << "data in Q region - other bkg in Q region" << endl;
  for(int i=1;i<=NSp_QCD->GetNbinsX();i++) cout << NSp_QCD->GetBinContent(i) << ", ";
  cout << endl;
  //NSp_QCD->Multiply(hMRR2_s_Multijet);
  //NSp_QCD->Divide(hMRR2_Q_Multijet);
  NSp_QCD->Scale(hMRR2_s_Multijet->Integral()/hMRR2_Q_Multijet->Integral());
  cout << "ClosureTest in S' region from Q region" << endl;
  for(int i=1;i<=NSp_QCD->GetNbinsX();i++) cout << NSp_QCD->GetBinContent(i) << ", ";
  cout << endl;
  //cout << NSp_QCD->Integral() << endl;

  NSp_W = (TH1F*)hMRR2_W_data->Clone();
  NSp_W->Add(hMRR2_W_Zjet, -1);
  NSp_W->Add(hMRR2_W_Gjet, -1);
  NSp_W->Add(hMRR2_W_TTST, -1);
  NSp_W->Add(hMRR2_W_Multijet, -1);
  NSp_W->Add(hMRR2_W_Multiboson, -1);
  cout << "data in W region - other bkg in W region" << endl;
  for(int i=1;i<=NSp_W->GetNbinsX();i++) cout << NSp_W->GetBinContent(i) << ", ";
  cout << endl;
  //NSp_W->Multiply(hMRR2_s_Wjet);
  //NSp_W->Divide(hMRR2_W_Wjet);
  NSp_W->Scale(hMRR2_s_Wjet->Integral()/hMRR2_W_Wjet->Integral());
  cout << "ClosureTest in S' region from W region" << endl;
  for(int i=1;i<=NSp_W->GetNbinsX();i++) cout << NSp_W->GetBinContent(i) << ", ";
  cout << endl;
  //cout << NSp_W->Integral() << endl;

  NSp_T = (TH1F*)hMRR2_T_data->Clone();
  NSp_T->Add(hMRR2_T_Zjet, -1);
  NSp_T->Add(hMRR2_T_Gjet, -1);
  NSp_T->Add(hMRR2_T_Wjet, -1);
  //NSp_T->Add(hMRR2_T_Multijet, -1);
  NSp_T->Add(hMRR2_T_Multiboson, -1);
  NSp_T->Add(NT_QCD,-1);
  cout << "data in T region - other bkg in T region" << endl;
  for(int i=1;i<=NSp_T->GetNbinsX();i++) cout << NSp_T->GetBinContent(i) << ", ";
  cout << endl;
  //NSp_T->Multiply(hMRR2_s_TTST);
  //NSp_T->Divide(hMRR2_T_TTST);
  NSp_T->Scale(hMRR2_s_TTST->Integral()/hMRR2_T_TTST->Integral());
  cout << "ClosureTest in S' region from T region" << endl;
  for(int i=1;i<=NSp_T->GetNbinsX();i++) cout << NSp_T->GetBinContent(i) << ", ";
  cout << endl;
  //cout << NSp_T->Integral() << endl;

  NSp_Z = (TH1F*)hMRR2_Z_data->Clone();
  NSp_Z->Add(hMRR2_Z_TTST, -1);
  NSp_Z->Add(hMRR2_Z_Wjet, -1);
  NSp_Z->Add(hMRR2_Z_Gjet, -1);
  NSp_Z->Add(hMRR2_Z_Multijet, -1);
  NSp_Z->Add(hMRR2_Z_Multiboson, -1);

  cout << "data in Z region - other bkg in Z region" << endl;
  for(int i=1;i<=NSp_Z->GetNbinsX();i++) cout << NSp_Z->GetBinContent(i) << ", ";
  cout << endl;
  //cout << NSp_Z->Integral() << endl;
  //NSp_Z->Multiply(hMRR2_s_Zjet);
  //for(int i=1;i<=NSp_Z->GetNbinsX();i++) cout << NSp_Z->GetBinContent(i) << ", ";
  //cout << endl;
  //cout << NSp_Z->Integral() << endl;
  //NSp_Z->Divide(hMRR2_Z_Zjet);
  NSp_Z->Scale(hMRR2_s_Zjet->Integral()/hMRR2_Z_Zjet->Integral());
  cout << "ClosureTest in S' region from Z region" << endl;
  for(int i=1;i<=NSp_Z->GetNbinsX();i++) cout << NSp_Z->GetBinContent(i) << ", ";
  cout << endl;
  //cout << hMRR2_Z_Zjet->Integral() << "div" << endl;
  //for(int i=1;i<=NSp_Z->GetNbinsX();i++) cout << hMRR2_Z_Zjet->GetBinContent(i) << ", ";
  //cout << endl;
  //cout << NSp_Z->Integral() << endl;
  for(int i=1;i<=NSp_Z->GetNbinsX();i++) cout << NSp_Z->GetBinContent(i) << ", ";
  cout << endl;

  NSp_G = (TH1F*)hMRR2_G_data->Clone();
  NSp_G->Add(hMRR2_G_TTST, -1);
  NSp_G->Add(hMRR2_G_Wjet, -1);
  NSp_G->Add(hMRR2_G_Zjet, -1);
  NSp_G->Add(hMRR2_G_Multijet, -1);
  NSp_G->Add(hMRR2_G_Multiboson, -1);
  //NSp_G->Multiply(hMRR2_s_Gjet);
  //NSp_G->Divide(hMRR2_G_Gjet);
  NSp_G->Scale(hMRR2_s_Gjet->Integral()/hMRR2_G_Gjet->Integral());

  NSp = (TH1F*)NSp_QCD->Clone();
  NSp->Add(NSp_W);
  NSp->Add(NSp_T);
  NSp->Add(NSp_Z);
  //NSp->Add(NSp_G);
  NSp->Add(hMRR2_s_Gjet);
  NSp->Add(hMRR2_s_Multiboson);
  //cout << NSp->Integral() << endl;

  NQp  = (TH1F*)hMRR2_Q_data->Clone();
  NQp1 = (TH1F*)hMRR2_q_Multijet->Clone();
  NQp1->Add(hMRR2_q_TTST);
  NQp1->Add(hMRR2_q_Wjet);
  NQp1->Add(hMRR2_q_Zjet);
  NQp1->Add(hMRR2_q_Gjet);
  NQp1->Add(hMRR2_q_Multiboson);
  NQp2 = (TH1F*)hMRR2_Q_Multijet->Clone();
  NQp2->Add(hMRR2_Q_TTST);
  NQp2->Add(hMRR2_Q_Wjet);
  NQp2->Add(hMRR2_Q_Zjet);
  NQp2->Add(hMRR2_Q_Gjet);
  NQp2->Add(hMRR2_Q_Multiboson);
  //NQp->Multiply(NQp1);
  //NQp->Divide(NQp2);
  NQp->Scale(NQp1->Integral()/NQp2->Integral());

  hMRR2_S_Multijet->SetFillColor(619);
  hMRR2_S_TTST->SetFillColor(633);
  hMRR2_S_Wjet->SetFillColor(418);
  hMRR2_S_Zjet->SetFillColor(401);
  hMRR2_S_Gjet->SetFillColor(803);
  hMRR2_S_Multiboson->SetFillColor(601);

  TH1D* hs = (TH1D*)hMRR2_S_Multijet->Clone();
  hs->Add(hMRR2_S_TTST);
  hs->Add(hMRR2_S_Wjet);
  hs->Add(hMRR2_S_Zjet);
  hs->Add(hMRR2_S_Gjet);
  hs->Add(hMRR2_S_Multiboson);

  THStack *hs1 = new THStack("hs1","");
  hs1->Add(hMRR2_S_Multiboson);
  hs1->Add(hMRR2_S_Wjet);
  hs1->Add(hMRR2_S_Zjet);
  hs1->Add(hMRR2_S_Gjet);
  hs1->Add(hMRR2_S_Multijet);
  hs1->Add(hMRR2_S_TTST);

  TCanvas* canvas_1 = new TCanvas("canvas_1", "", 800, 800);
  canvas_1->Divide(1,2);
  TPad *canvas_up = (TPad*)canvas_1->GetListOfPrimitives()->FindObject("canvas_1_1");
  TPad *canvas_dw = (TPad*)canvas_1->GetListOfPrimitives()->FindObject("canvas_1_2");

  double up_height     = 0.7; // please tune so that the upper figures size will meet your requirement
  double dw_correction = 1.210;//1.390 // please tune so that the smaller canvas size will work in your environment
  double dw_height    = (1. - up_height) * dw_correction;

  canvas_up->SetPad(0., 0., 1., 1.);
  canvas_dw->SetPad(0., 0., 1., 0.);
  canvas_up->SetLogy();
  canvas_up->SetFrameFillColor(0);
  canvas_up->SetFillColor(0);
  canvas_dw->SetFillColor(0);
  canvas_dw->SetFrameFillColor(0);
  canvas_up->SetLeftMargin(0.1);
  canvas_dw->SetTopMargin(0.03);//0.03
  canvas_dw->SetLeftMargin(0.1);
  canvas_dw->SetBottomMargin(0.3);

  canvas_up->cd();
  hMRR2_S_data->SetMarkerColor(kBlack);
  hMRR2_S_data->SetLineColor(kBlack);
  hMRR2_S_data->SetMarkerStyle(20);
  hMRR2_S_data->SetMarkerSize(1.3);
  hs1->SetMaximum(1e3);
  hs1->SetMinimum(1e-1);
  hs1->Draw("HIST");
  hs1->GetXaxis()->SetLabelOffset(999);
  //hs1->GetYaxis()->SetTitle("Events/bin");
  //hMRR2_S_data->Draw("EPsame");
  hMRR2_S_T2tt->SetLineColor(403);
  hMRR2_S_T5tttt->SetLineColor(862);
  hMRR2_S_T5ttcc->SetLineColor(12);
  hMRR2_S_T1tttt->SetLineColor(841);
  hMRR2_S_T2tt->Draw("HISTSAME");
  hMRR2_S_T5tttt->Draw("HISTSAME");
  hMRR2_S_T5ttcc->Draw("HISTSAME");
  hMRR2_S_T1tttt->Draw("HISTSAME");

  float xmin = hMRR2_S_data->GetXaxis()->GetBinLowEdge(hMRR2_S_data->GetXaxis()->GetFirst());
  float xmax = hMRR2_S_data->GetXaxis()->GetBinUpEdge(hMRR2_S_data->GetXaxis()->GetLast());

  string text = "CMS Simulation #scale[0.7]{#font[52]{Work in progress 2016}}";
  TLatex* cms_lat = new TLatex(xmin, 2.0e3, text.c_str());
  cms_lat->SetTextSize(0.04);
  cms_lat->SetLineWidth(2);
  cms_lat->Draw();
  text = "#scale[0.7]{W ana, 35.9 fb^{-1} (13 TeV)}";
  if(channel=="_topana")text = "#scale[0.7]{Top ana, 35.9 fb^{-1} (13 TeV)}";
  TLatex* era_lat = new TLatex(xmax,2.4e3, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.03);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();
  text = "#scale[0.7]{   [800, 1000]   [1000, 1200]   [1200, 1600]  [1600, 2000]   [2000, 4000]}";
  era_lat = new TLatex(xmax,3.e1, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.042);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();

  leg = new TLegend(0.45,0.71,0.9,0.9);
  //if(TString(njet).Contains("nj35")) leg->SetHeader("S' Closure test, 3 #leq jet < 5");
  //else if(TString(njet).Contains("nj35")) leg->SetHeader("S' Closure test, 6 #leq jet");
  //else leg->SetHeader("S' Closure test");

  char legentry[10];

  leg->SetNColumns(4);
  leg->SetTextSize(0.02);
  //sprintf(legentry, "%.1f",hMRR2_S_data->Integral());
  //leg->AddEntry(hMRR2_S_data, "data", "lep");
  //leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_S_Multijet->Integral());
  leg->AddEntry(hMRR2_S_Multijet,  "Multijet", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_S_TTST->Integral());
  leg->AddEntry(hMRR2_S_TTST,  "t#bar{t} + single top", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_S_Wjet->Integral());
  leg->AddEntry(hMRR2_S_Wjet,  "W(#rightarrowl#nu)", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_S_Zjet->Integral());
  leg->AddEntry(hMRR2_S_Zjet,  "Z(#rightarrow#nu#nu)", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_S_Gjet->Integral());
  leg->AddEntry(hMRR2_S_Gjet,  "#gamma + jet", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_S_Multiboson->Integral());
  leg->AddEntry(hMRR2_S_Multiboson,  "VV(V) + t#bar{t}(X)", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_S_T2tt->Integral());
  leg->AddEntry(hMRR2_S_T2tt, "T2tt", "l");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_S_T5tttt->Integral());
  leg->AddEntry(hMRR2_S_T5tttt, "T5tttt", "l");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_S_T5ttcc->Integral());
  leg->AddEntry(hMRR2_S_T5ttcc, "T5ttcc", "l");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_S_T1tttt->Integral());
  leg->AddEntry(hMRR2_S_T1tttt, "T1tttt", "l");
  leg->AddEntry((TObject*)0,legentry,"");
  leg->Draw();
  
  hMRR2_Q_Multijet->SetFillColor(619);
  hMRR2_Q_TTST->SetFillColor(633);
  hMRR2_Q_Wjet->SetFillColor(418);
  hMRR2_Q_Zjet->SetFillColor(401);
  hMRR2_Q_Gjet->SetFillColor(803);
  hMRR2_Q_Multiboson->SetFillColor(601);

  hs = (TH1D*)hMRR2_Q_Multijet->Clone();
  hs->Add(hMRR2_Q_TTST);
  hs->Add(hMRR2_Q_Wjet);
  hs->Add(hMRR2_Q_Zjet);
  hs->Add(hMRR2_Q_Gjet);
  hs->Add(hMRR2_Q_Multiboson);

  hs1 = new THStack("hs1","");
  hs1->Add(hMRR2_Q_Multiboson);
  hs1->Add(hMRR2_Q_Wjet);
  hs1->Add(hMRR2_Q_Zjet);
  hs1->Add(hMRR2_Q_Gjet);
  hs1->Add(hMRR2_Q_Multijet);
  hs1->Add(hMRR2_Q_TTST);

  TCanvas* canvas_2 = new TCanvas("canvas_2", "", 800, 800);
  canvas_2->Divide(1,2);
  canvas_up = (TPad*)canvas_2->GetListOfPrimitives()->FindObject("canvas_2_1");
  canvas_dw = (TPad*)canvas_2->GetListOfPrimitives()->FindObject("canvas_2_2");

  up_height     = 0.7; // please tune so that the upper figures size will meet your requirement
  dw_correction = 1.210;//1.390 // please tune so that the smaller canvas size will work in your environment
  dw_height    = (1. - up_height) * dw_correction;

  canvas_up->SetPad(0., 1 - up_height, 1., 1.);
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

  canvas_up->cd();
  hMRR2_Q_data->SetMarkerColor(kBlack);
  hMRR2_Q_data->SetLineColor(kBlack);
  hMRR2_Q_data->SetMarkerStyle(20);
  hMRR2_Q_data->SetMarkerSize(1.3);
  hs1->SetMaximum(1e4);
  hs1->SetMinimum(1e-1);
  hs1->Draw("HIST");
  //hs1->GetYaxis()->SetTitle("Events/bin");
  hMRR2_Q_data->Draw("EPsame");
  hMRR2_Q_T2tt->SetLineColor(403);
  hMRR2_Q_T5tttt->SetLineColor(862);
  hMRR2_Q_T5ttcc->SetLineColor(12);
  hMRR2_Q_T1tttt->SetLineColor(841);
  hMRR2_Q_T2tt->Draw("HISTSAME");
  hMRR2_Q_T5tttt->Draw("HISTSAME");
  hMRR2_Q_T5ttcc->Draw("HISTSAME");
  hMRR2_Q_T1tttt->Draw("HISTSAME");

  xmin = hMRR2_Q_data->GetXaxis()->GetBinLowEdge(hMRR2_Q_data->GetXaxis()->GetFirst());
  xmax = hMRR2_Q_data->GetXaxis()->GetBinUpEdge(hMRR2_Q_data->GetXaxis()->GetLast());

  text = "CMS #scale[0.7]{#font[52]{Work in progress 2016}}";
  cms_lat = new TLatex(xmin, 2.2e4, text.c_str());
  cms_lat->SetTextSize(0.07);
  cms_lat->SetLineWidth(2);
  cms_lat->Draw();
  text = "#scale[0.7]{W ana, 35.9 fb^{-1} (13 TeV)}";
  if(channel=="_topana")text = "#scale[0.7]{Top ana, 35.9 fb^{-1} (13 TeV)}";
  era_lat = new TLatex(xmax,2.6e4, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.06);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();
  text = "#scale[0.7]{   [800, 1000]   [1000, 1200]   [1200, 1600]  [1600, 2000]   [2000, 4000]}";
  era_lat = new TLatex(xmax,1.e2, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.06);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();

  leg = new TLegend(0.45,0.62,0.9,0.9);
  //if(TString(njet).Contains("nj35")) leg->SetHeader("S' Closure test, 3 #leq jet < 5");
  //else if(TString(njet).Contains("nj35")) leg->SetHeader("S' Closure test, 6 #leq jet");
  //else leg->SetHeader("S' Closure test");

  //char legentry[10];

  leg->SetNColumns(4);
  leg->SetTextSize(0.03);
  sprintf(legentry, "%.1f",hMRR2_Q_data->Integral());
  leg->AddEntry(hMRR2_Q_data, "data", "lep");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_Q_Multijet->Integral());
  leg->AddEntry(hMRR2_Q_Multijet,  "Multijet", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_Q_TTST->Integral());
  leg->AddEntry(hMRR2_Q_TTST,  "t#bar{t} + single top", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_Q_Wjet->Integral());
  leg->AddEntry(hMRR2_Q_Wjet,  "W(#rightarrowl#nu)", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_Q_Zjet->Integral());
  leg->AddEntry(hMRR2_Q_Zjet,  "Z(#rightarrow#nu#nu)", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_Q_Gjet->Integral());
  leg->AddEntry(hMRR2_Q_Gjet,  "#gamma + jet", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_Q_Multiboson->Integral());
  leg->AddEntry(hMRR2_Q_Multiboson,  "VV(V) + t#bar{t}(X)", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_Q_T2tt->Integral());
  leg->AddEntry(hMRR2_Q_T2tt, "T2tt", "l");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_Q_T5tttt->Integral());
  leg->AddEntry(hMRR2_Q_T5tttt, "T5tttt", "l");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_Q_T5ttcc->Integral());
  leg->AddEntry(hMRR2_Q_T5ttcc, "T5ttcc", "l");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_Q_T1tttt->Integral());
  leg->AddEntry(hMRR2_Q_T1tttt, "T1tttt", "l");
  leg->AddEntry((TObject*)0,legentry,"");
  leg->Draw();
  
  compareSp = (TH1F*)hMRR2_Q_data->Clone(); 
  compareSp1 = (TH1F*)hMRR2_Q_data->Clone(); 
  compareSp->Divide(hs);

  canvas_dw->cd();
  gPad->SetGridy();
  compareSp->SetLineColor(kBlack);
  compareSp->SetMarkerColor(kBlack);
  compareSp->SetMarkerStyle(20);
  compareSp->SetMarkerSize(1.3);
  compareSp->GetYaxis()->SetRangeUser(0.0,2.0);
  compareSp->SetTitle(0);
  compareSp->SetStats(0);
  compareSp->GetYaxis()->SetLabelSize(0.09);
  compareSp->GetXaxis()->SetLabelSize(0.09);
  compareSp->GetXaxis()->LabelsOption("v");
  compareSp->GetYaxis()->SetNdivisions(502);
  compareSp->GetXaxis()->SetTickLength(0.07);
  compareSp->GetYaxis()->SetTickLength(0.03);
  compareSp->GetYaxis()->SetTitle("#frac{Data}{Estimate}");
  compareSp->GetXaxis()->SetTitleOffset(0.81);
  compareSp->GetXaxis()->SetTitleSize(0.12);
  compareSp->GetYaxis()->SetTitleOffset(0.35);
  compareSp->GetYaxis()->SetTitleSize(0.10);//dividend->SetTitleSize(0.12,"Y");
  for(int i=1;i<=hMRR2_Q_data->GetNbinsX();i++){
    compareSp1->SetBinContent(i,1);
    if(hMRR2_Q_data->GetBinContent(i)==0) compareSp1->SetBinError(i,0);
    else compareSp1->SetBinError(i,hs->GetBinError(i)/hs->GetBinContent(i));
  }
  compareSp1->SetFillColor(13);
  compareSp1->SetFillStyle(3235);
  compareSp1->SetMarkerStyle(1);
  compareSp->Draw("EP");
  compareSp1->Draw("E2same");


  hMRR2_T_Multijet->SetFillColor(619);
  hMRR2_T_TTST->SetFillColor(633);
  hMRR2_T_Wjet->SetFillColor(418);
  hMRR2_T_Zjet->SetFillColor(401);
  hMRR2_T_Gjet->SetFillColor(803);
  hMRR2_T_Multiboson->SetFillColor(601);

  hs = (TH1D*)hMRR2_T_Multijet->Clone();
  hs->Add(hMRR2_T_TTST);
  hs->Add(hMRR2_T_Wjet);
  hs->Add(hMRR2_T_Zjet);
  hs->Add(hMRR2_T_Gjet);
  hs->Add(hMRR2_T_Multiboson);

  hs1 = new THStack("hs1","");
  hs1->Add(hMRR2_T_Multiboson);
  hs1->Add(hMRR2_T_Wjet);
  hs1->Add(hMRR2_T_Gjet);
  hs1->Add(hMRR2_T_Multijet);
  hs1->Add(hMRR2_T_TTST);

  TCanvas* canvas_3 = new TCanvas("canvas_3", "", 800, 800);
  canvas_3->Divide(1,2);
  canvas_up = (TPad*)canvas_3->GetListOfPrimitives()->FindObject("canvas_3_1");
  canvas_dw = (TPad*)canvas_3->GetListOfPrimitives()->FindObject("canvas_3_2");

  up_height     = 0.7; // please tune so that the upper figures size will meet your requirement
  dw_correction = 1.210;//1.390 // please tune so that the smaller canvas size will work in your environment
  dw_height    = (1. - up_height) * dw_correction;

  canvas_up->SetPad(0., 1 - up_height, 1., 1.);
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

  canvas_up->cd();
  hMRR2_T_data->SetMarkerColor(kBlack);
  hMRR2_T_data->SetLineColor(kBlack);
  hMRR2_T_data->SetMarkerStyle(20);
  hMRR2_T_data->SetMarkerSize(1.3);
  hs1->SetMaximum(1e3);
  hs1->SetMinimum(1e-1);
  hs1->Draw("HIST");
  //hs1->GetYaxis()->SetTitle("Events/bin");
  hMRR2_T_data->Draw("EPsame");
  hMRR2_T_T2tt->SetLineColor(403);
  hMRR2_T_T5tttt->SetLineColor(862);
  hMRR2_T_T5ttcc->SetLineColor(12);
  hMRR2_T_T1tttt->SetLineColor(841);
  hMRR2_T_T2tt->Draw("HISTSAME");
  hMRR2_T_T5tttt->Draw("HISTSAME");
  hMRR2_T_T5ttcc->Draw("HISTSAME");
  hMRR2_T_T1tttt->Draw("HISTSAME");

  xmin = hMRR2_T_data->GetXaxis()->GetBinLowEdge(hMRR2_T_data->GetXaxis()->GetFirst());
  xmax = hMRR2_T_data->GetXaxis()->GetBinUpEdge(hMRR2_T_data->GetXaxis()->GetLast());

  text = "CMS #scale[0.7]{#font[52]{Work in progress 2016}}";
  cms_lat = new TLatex(xmin, 1.9e3, text.c_str());
  cms_lat->SetTextSize(0.07);
  cms_lat->SetLineWidth(2);
  cms_lat->Draw();
  text = "#scale[0.7]{W ana, 35.9 fb^{-1} (13 TeV)}";
  if(channel=="_topana")text = "#scale[0.7]{Top ana, 35.9 fb^{-1} (13 TeV)}";
  era_lat = new TLatex(xmax,2.3e3, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.06);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();
  text = "#scale[0.7]{   [800, 1000]   [1000, 1200]   [1200, 1600]  [1600, 2000]   [2000, 4000]}";
  era_lat = new TLatex(xmax,3.e1, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.06);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();

  leg = new TLegend(0.45,0.62,0.9,0.9);
  //if(TString(njet).Contains("nj35")) leg->SetHeader("S' Closure test, 3 #leq jet < 5");
  //else if(TString(njet).Contains("nj35")) leg->SetHeader("S' Closure test, 6 #leq jet");
  //else leg->SetHeader("S' Closure test");

  //char legentry[10];

  leg->SetNColumns(4);
  leg->SetTextSize(0.03);
  sprintf(legentry, "%.1f",hMRR2_T_data->Integral());
  leg->AddEntry(hMRR2_T_data, "data", "lep");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_T_Multijet->Integral());
  leg->AddEntry(hMRR2_T_Multijet,  "Multijet", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_T_TTST->Integral());
  leg->AddEntry(hMRR2_T_TTST,  "t#bar{t} + single top", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_T_Wjet->Integral());
  leg->AddEntry(hMRR2_T_Wjet,  "W(#rightarrowl#nu)", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_T_Zjet->Integral());
  leg->AddEntry(hMRR2_T_Zjet,  "Z(#rightarrow#nu#nu)", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_T_Gjet->Integral());
  leg->AddEntry(hMRR2_T_Gjet,  "#gamma + jet", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_T_Multiboson->Integral());
  leg->AddEntry(hMRR2_T_Multiboson,  "VV(V) + t#bar{t}(X)", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_T_T2tt->Integral());
  leg->AddEntry(hMRR2_T_T2tt, "T2tt", "l");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_T_T5tttt->Integral());
  leg->AddEntry(hMRR2_T_T5tttt, "T5tttt", "l");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_T_T5ttcc->Integral());
  leg->AddEntry(hMRR2_T_T5ttcc, "T5ttcc", "l");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_T_T1tttt->Integral());
  leg->AddEntry(hMRR2_T_T1tttt, "T1tttt", "l");
  leg->AddEntry((TObject*)0,legentry,"");
  leg->Draw();
  
  compareSp = (TH1F*)hMRR2_T_data->Clone(); 
  compareSp1 = (TH1F*)hMRR2_T_data->Clone(); 
  compareSp->Divide(hs);

  canvas_dw->cd();
  gPad->SetGridy();
  compareSp->SetLineColor(kBlack);
  compareSp->SetMarkerColor(kBlack);
  compareSp->SetMarkerStyle(20);
  compareSp->SetMarkerSize(1.3);
  compareSp->GetYaxis()->SetRangeUser(0.0,2.0);
  compareSp->SetTitle(0);
  compareSp->SetStats(0);
  compareSp->GetYaxis()->SetLabelSize(0.09);
  compareSp->GetXaxis()->SetLabelSize(0.09);
  compareSp->GetXaxis()->LabelsOption("v");
  compareSp->GetYaxis()->SetNdivisions(502);
  compareSp->GetXaxis()->SetTickLength(0.07);
  compareSp->GetYaxis()->SetTickLength(0.03);
  compareSp->GetYaxis()->SetTitle("#frac{Data}{Estimate}");
  compareSp->GetXaxis()->SetTitleOffset(0.81);
  compareSp->GetXaxis()->SetTitleSize(0.12);
  compareSp->GetYaxis()->SetTitleOffset(0.35);
  compareSp->GetYaxis()->SetTitleSize(0.10);//dividend->SetTitleSize(0.12,"Y");
  for(int i=1;i<=hMRR2_T_data->GetNbinsX();i++){
    compareSp1->SetBinContent(i,1);
    if(hMRR2_T_data->GetBinContent(i)==0) compareSp1->SetBinError(i,0);
    else compareSp1->SetBinError(i,hs->GetBinError(i)/hs->GetBinContent(i));
  }
  compareSp1->SetFillColor(13);
  compareSp1->SetFillStyle(3235);
  compareSp1->SetMarkerStyle(1);
  compareSp->Draw("EP");
  compareSp1->Draw("E2same");


  hMRR2_W_Multijet->SetFillColor(619);
  hMRR2_W_TTST->SetFillColor(633);
  hMRR2_W_Wjet->SetFillColor(418);
  hMRR2_W_Zjet->SetFillColor(401);
  hMRR2_W_Gjet->SetFillColor(803);
  hMRR2_W_Multiboson->SetFillColor(601);

  hs = (TH1D*)hMRR2_W_Multijet->Clone();
  hs->Add(hMRR2_W_TTST);
  hs->Add(hMRR2_W_Wjet);
  hs->Add(hMRR2_W_Zjet);
  hs->Add(hMRR2_W_Gjet);
  hs->Add(hMRR2_W_Multiboson);

  hs1 = new THStack("hs1","");
  hs1->Add(hMRR2_W_Multiboson);
  hs1->Add(hMRR2_W_Wjet);
  hs1->Add(hMRR2_W_Zjet);
  hs1->Add(hMRR2_W_Gjet);
  hs1->Add(hMRR2_W_Multijet);
  hs1->Add(hMRR2_W_TTST);

  TCanvas* canvas_4 = new TCanvas("canvas_4", "", 800, 800);
  canvas_4->Divide(1,2);
  canvas_up = (TPad*)canvas_4->GetListOfPrimitives()->FindObject("canvas_4_1");
  canvas_dw = (TPad*)canvas_4->GetListOfPrimitives()->FindObject("canvas_4_2");

  up_height     = 0.7; // please tune so that the upper figures size will meet your requirement
  dw_correction = 1.210;//1.390 // please tune so that the smaller canvas size will work in your environment
  dw_height    = (1. - up_height) * dw_correction;

  canvas_up->SetPad(0., 1 - up_height, 1., 1.);
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

  canvas_up->cd();
  hMRR2_W_data->SetMarkerColor(kBlack);
  hMRR2_W_data->SetLineColor(kBlack);
  hMRR2_W_data->SetMarkerStyle(20);
  hMRR2_W_data->SetMarkerSize(1.3);
  hs1->SetMaximum(1e3);
  hs1->SetMinimum(1e-1);
  hs1->Draw("HIST");
  //hs1->GetYaxis()->SetTitle("Events/bin");
  hMRR2_W_data->Draw("EPsame");
  hMRR2_W_T2tt->SetLineColor(403);
  hMRR2_W_T5tttt->SetLineColor(862);
  hMRR2_W_T5ttcc->SetLineColor(12);
  hMRR2_W_T1tttt->SetLineColor(841);
  hMRR2_W_T2tt->Draw("HISTSAME");
  hMRR2_W_T5tttt->Draw("HISTSAME");
  hMRR2_W_T5ttcc->Draw("HISTSAME");
  hMRR2_W_T1tttt->Draw("HISTSAME");

  xmin = hMRR2_W_data->GetXaxis()->GetBinLowEdge(hMRR2_W_data->GetXaxis()->GetFirst());
  xmax = hMRR2_W_data->GetXaxis()->GetBinUpEdge(hMRR2_W_data->GetXaxis()->GetLast());

  text = "CMS #scale[0.7]{#font[52]{Work in progress 2016}}";
  cms_lat = new TLatex(xmin, 1.9e3, text.c_str());
  cms_lat->SetTextSize(0.07);
  cms_lat->SetLineWidth(2);
  cms_lat->Draw();
  text = "#scale[0.7]{W ana, 35.9 fb^{-1} (13 TeV)}";
  if(channel=="_topana")text = "#scale[0.7]{Top ana, 35.9 fb^{-1} (13 TeV)}";
  era_lat = new TLatex(xmax,2.3e3, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.06);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();
  text = "#scale[0.7]{   [800, 1000]   [1000, 1200]   [1200, 1600]  [1600, 2000]   [2000, 4000]}";
  era_lat = new TLatex(xmax,3.e1, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.06);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();

  leg = new TLegend(0.45,0.62,0.9,0.9);
  //if(TString(njet).Contains("nj35")) leg->SetHeader("S' Closure test, 3 #leq jet < 5");
  //else if(TString(njet).Contains("nj35")) leg->SetHeader("S' Closure test, 6 #leq jet");
  //else leg->SetHeader("S' Closure test");

  //char legentry[10];

  leg->SetNColumns(4);
  leg->SetTextSize(0.03);
  sprintf(legentry, "%.1f",hMRR2_W_data->Integral());
  leg->AddEntry(hMRR2_W_data, "data", "lep");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_W_Multijet->Integral());
  leg->AddEntry(hMRR2_W_Multijet,  "Multijet", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_W_TTST->Integral());
  leg->AddEntry(hMRR2_W_TTST,  "t#bar{t} + single top", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_W_Wjet->Integral());
  leg->AddEntry(hMRR2_W_Wjet,  "W(#rightarrowl#nu)", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_W_Zjet->Integral());
  leg->AddEntry(hMRR2_W_Zjet,  "Z(#rightarrow#nu#nu)", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_W_Gjet->Integral());
  leg->AddEntry(hMRR2_W_Gjet,  "#gamma + jet", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_W_Multiboson->Integral());
  leg->AddEntry(hMRR2_W_Multiboson,  "VV(V) + t#bar{t}(X)", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_W_T2tt->Integral());
  leg->AddEntry(hMRR2_W_T2tt, "T2tt", "l");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_W_T5tttt->Integral());
  leg->AddEntry(hMRR2_W_T5tttt, "T5tttt", "l");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_W_T5ttcc->Integral());
  leg->AddEntry(hMRR2_W_T5ttcc, "T5ttcc", "l");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_W_T1tttt->Integral());
  leg->AddEntry(hMRR2_W_T1tttt, "T1tttt", "l");
  leg->AddEntry((TObject*)0,legentry,"");
  leg->Draw();
  
  compareSp = (TH1F*)hMRR2_W_data->Clone(); 
  compareSp1 = (TH1F*)hMRR2_W_data->Clone(); 
  compareSp->Divide(hs);

  canvas_dw->cd();
  gPad->SetGridy();
  compareSp->SetLineColor(kBlack);
  compareSp->SetMarkerColor(kBlack);
  compareSp->SetMarkerStyle(20);
  compareSp->SetMarkerSize(1.3);
  compareSp->GetYaxis()->SetRangeUser(0.0,2.0);
  compareSp->SetTitle(0);
  compareSp->SetStats(0);
  compareSp->GetYaxis()->SetLabelSize(0.09);
  compareSp->GetXaxis()->SetLabelSize(0.09);
  compareSp->GetXaxis()->LabelsOption("v");
  compareSp->GetYaxis()->SetNdivisions(502);
  compareSp->GetXaxis()->SetTickLength(0.07);
  compareSp->GetYaxis()->SetTickLength(0.03);
  compareSp->GetYaxis()->SetTitle("#frac{Data}{Estimate}");
  compareSp->GetXaxis()->SetTitleOffset(0.81);
  compareSp->GetXaxis()->SetTitleSize(0.12);
  compareSp->GetYaxis()->SetTitleOffset(0.35);
  compareSp->GetYaxis()->SetTitleSize(0.10);//dividend->SetTitleSize(0.12,"Y");
  for(int i=1;i<=hMRR2_W_data->GetNbinsX();i++){
    compareSp1->SetBinContent(i,1);
    if(hMRR2_W_data->GetBinContent(i)==0) compareSp1->SetBinError(i,0);
    else compareSp1->SetBinError(i,hs->GetBinError(i)/hs->GetBinContent(i));
  }
  compareSp1->SetFillColor(13);
  compareSp1->SetFillStyle(3235);
  compareSp1->SetMarkerStyle(1);
  compareSp->Draw("EP");
  compareSp1->Draw("E2same");


  hMRR2_Z_Multijet->SetFillColor(619);
  hMRR2_Z_TTST->SetFillColor(633);
  hMRR2_Z_Wjet->SetFillColor(418);
  hMRR2_Z_Zjet->SetFillColor(401);
  hMRR2_Z_Gjet->SetFillColor(803);
  hMRR2_Z_Multiboson->SetFillColor(601);

  hs = (TH1D*)hMRR2_Z_Multijet->Clone();
  hs->Add(hMRR2_Z_TTST);
  hs->Add(hMRR2_Z_Wjet);
  hs->Add(hMRR2_Z_Zjet);
  hs->Add(hMRR2_Z_Gjet);
  hs->Add(hMRR2_Z_Multiboson);

  hs1 = new THStack("hs1","");
  hs1->Add(hMRR2_Z_Multiboson);
  hs1->Add(hMRR2_Z_Wjet);
  hs1->Add(hMRR2_Z_Zjet);
  hs1->Add(hMRR2_Z_Gjet);
  hs1->Add(hMRR2_Z_Multijet);
  hs1->Add(hMRR2_Z_TTST);

  TCanvas* canvas_5 = new TCanvas("canvas_5", "", 800, 800);
  canvas_5->Divide(1,2);
  canvas_up = (TPad*)canvas_5->GetListOfPrimitives()->FindObject("canvas_5_1");
  canvas_dw = (TPad*)canvas_5->GetListOfPrimitives()->FindObject("canvas_5_2");

  up_height     = 0.7; // please tune so that the upper figures size will meet your requirement
  dw_correction = 1.210;//1.390 // please tune so that the smaller canvas size will work in your environment
  dw_height    = (1. - up_height) * dw_correction;

  canvas_up->SetPad(0., 1 - up_height, 1., 1.);
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

  canvas_up->cd();
  hMRR2_Z_data->SetMarkerColor(kBlack);
  hMRR2_Z_data->SetLineColor(kBlack);
  hMRR2_Z_data->SetMarkerStyle(20);
  hMRR2_Z_data->SetMarkerSize(1.3);
  hs1->SetMaximum(1e3);
  hs1->SetMinimum(1e-1);
  hs1->Draw("HIST");
  //hs1->GetYaxis()->SetTitle("Events/bin");
  hMRR2_Z_data->Draw("EPsame");
  hMRR2_Z_T2tt->SetLineColor(403);
  hMRR2_Z_T5tttt->SetLineColor(862);
  hMRR2_Z_T5ttcc->SetLineColor(12);
  hMRR2_Z_T1tttt->SetLineColor(841);
  hMRR2_Z_T2tt->Draw("HISTSAME");
  hMRR2_Z_T5tttt->Draw("HISTSAME");
  hMRR2_Z_T5ttcc->Draw("HISTSAME");
  hMRR2_Z_T1tttt->Draw("HISTSAME");

  xmin = hMRR2_Z_data->GetXaxis()->GetBinLowEdge(hMRR2_Z_data->GetXaxis()->GetFirst());
  xmax = hMRR2_Z_data->GetXaxis()->GetBinUpEdge(hMRR2_Z_data->GetXaxis()->GetLast());

  text = "CMS #scale[0.7]{#font[52]{Work in progress 2016}}";
  cms_lat = new TLatex(xmin, 1.9e3, text.c_str());
  cms_lat->SetTextSize(0.07);
  cms_lat->SetLineWidth(2);
  cms_lat->Draw();
  text = "#scale[0.7]{W ana, 35.9 fb^{-1} (13 TeV)}";
  if(channel=="_topana")text = "#scale[0.7]{Top ana, 35.9 fb^{-1} (13 TeV)}";
  era_lat = new TLatex(xmax,2.3e3, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.06);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();
  text = "#scale[0.7]{   [800, 1000]   [1000, 1200]   [1200, 1600]  [1600, 2000]   [2000, 4000]}";
  era_lat = new TLatex(xmax,3.e1, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.06);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();

  leg = new TLegend(0.45,0.62,0.9,0.9);
  //if(TString(njet).Contains("nj35")) leg->SetHeader("S' Closure test, 3 #leq jet < 5");
  //else if(TString(njet).Contains("nj35")) leg->SetHeader("S' Closure test, 6 #leq jet");
  //else leg->SetHeader("S' Closure test");

  //char legentry[10];

  leg->SetNColumns(4);
  leg->SetTextSize(0.03);
  sprintf(legentry, "%.1f",hMRR2_Z_data->Integral());
  leg->AddEntry(hMRR2_Z_data, "data", "lep");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_Z_Multijet->Integral());
  leg->AddEntry(hMRR2_Z_Multijet,  "Multijet", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_Z_TTST->Integral());
  leg->AddEntry(hMRR2_Z_TTST,  "t#bar{t} + single top", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_Z_Wjet->Integral());
  leg->AddEntry(hMRR2_Z_Wjet,  "W(#rightarrowl#nu)", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_Z_Zjet->Integral());
  leg->AddEntry(hMRR2_Z_Zjet,  "Z(#rightarrow#nu#nu)", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_Z_Gjet->Integral());
  leg->AddEntry(hMRR2_Z_Gjet,  "#gamma + jet", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_Z_Multiboson->Integral());
  leg->AddEntry(hMRR2_Z_Multiboson,  "VV(V) + t#bar{t}(X)", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_Z_T2tt->Integral());
  leg->AddEntry(hMRR2_Z_T2tt, "T2tt", "l");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_Z_T5tttt->Integral());
  leg->AddEntry(hMRR2_Z_T5tttt, "T5tttt", "l");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_Z_T5ttcc->Integral());
  leg->AddEntry(hMRR2_Z_T5ttcc, "T5ttcc", "l");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_Z_T1tttt->Integral());
  leg->AddEntry(hMRR2_Z_T1tttt, "T1tttt", "l");
  leg->AddEntry((TObject*)0,legentry,"");
  leg->Draw();
  
  compareSp = (TH1F*)hMRR2_Z_data->Clone(); 
  compareSp1 = (TH1F*)hMRR2_Z_data->Clone(); 
  compareSp->Divide(hs);

  canvas_dw->cd();
  gPad->SetGridy();
  compareSp->SetLineColor(kBlack);
  compareSp->SetMarkerColor(kBlack);
  compareSp->SetMarkerStyle(20);
  compareSp->SetMarkerSize(1.3);
  compareSp->GetYaxis()->SetRangeUser(0.0,2.0);
  compareSp->SetTitle(0);
  compareSp->SetStats(0);
  compareSp->GetYaxis()->SetLabelSize(0.09);
  compareSp->GetXaxis()->SetLabelSize(0.09);
  compareSp->GetXaxis()->LabelsOption("v");
  compareSp->GetYaxis()->SetNdivisions(502);
  compareSp->GetXaxis()->SetTickLength(0.07);
  compareSp->GetYaxis()->SetTickLength(0.03);
  compareSp->GetYaxis()->SetTitle("#frac{Data}{Estimate}");
  compareSp->GetXaxis()->SetTitleOffset(0.81);
  compareSp->GetXaxis()->SetTitleSize(0.12);
  compareSp->GetYaxis()->SetTitleOffset(0.35);
  compareSp->GetYaxis()->SetTitleSize(0.10);//dividend->SetTitleSize(0.12,"Y");
  for(int i=1;i<=hMRR2_Z_data->GetNbinsX();i++){
    compareSp1->SetBinContent(i,1);
    if(hMRR2_Z_data->GetBinContent(i)==0) compareSp1->SetBinError(i,0);
    else compareSp1->SetBinError(i,hs->GetBinError(i)/hs->GetBinContent(i));
  }
  compareSp1->SetFillColor(13);
  compareSp1->SetFillStyle(3235);
  compareSp1->SetMarkerStyle(1);
  compareSp->Draw("EP");
  compareSp1->Draw("E2same");


  hMRR2_G_Multijet->SetFillColor(619);
  hMRR2_G_TTST->SetFillColor(633);
  hMRR2_G_Wjet->SetFillColor(418);
  hMRR2_G_Zjet->SetFillColor(401);
  hMRR2_G_Gjet->SetFillColor(803);
  hMRR2_G_Multiboson->SetFillColor(601);

  hs = (TH1D*)hMRR2_G_Multijet->Clone();
  hs->Add(hMRR2_G_TTST);
  hs->Add(hMRR2_G_Wjet);
  hs->Add(hMRR2_G_Zjet);
  hs->Add(hMRR2_G_Gjet);
  hs->Add(hMRR2_G_Multiboson);

  hs1 = new THStack("hs1","");
  hs1->Add(hMRR2_G_Multiboson);
  hs1->Add(hMRR2_G_Wjet);
  hs1->Add(hMRR2_G_Zjet);
  hs1->Add(hMRR2_G_Gjet);
  hs1->Add(hMRR2_G_Multijet);
  hs1->Add(hMRR2_G_TTST);

  TCanvas* canvas_6 = new TCanvas("canvas_6", "", 800, 800);
  canvas_6->Divide(1,2);
  canvas_up = (TPad*)canvas_6->GetListOfPrimitives()->FindObject("canvas_6_1");
  canvas_dw = (TPad*)canvas_6->GetListOfPrimitives()->FindObject("canvas_6_2");

  up_height     = 0.7; // please tune so that the upper figures size will meet your requirement
  dw_correction = 1.210;//1.390 // please tune so that the smaller canvas size will work in your environment
  dw_height    = (1. - up_height) * dw_correction;

  canvas_up->SetPad(0., 1 - up_height, 1., 1.);
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

  canvas_up->cd();
  hMRR2_G_data->SetMarkerColor(kBlack);
  hMRR2_G_data->SetLineColor(kBlack);
  hMRR2_G_data->SetMarkerStyle(20);
  hMRR2_G_data->SetMarkerSize(1.3);
  hs1->SetMaximum(1e3);
  hs1->SetMinimum(1e-1);
  hs1->Draw("HIST");
  //hs1->GetYaxis()->SetTitle("Events/bin");
  hMRR2_G_data->Draw("EPsame");

  xmin = hMRR2_G_data->GetXaxis()->GetBinLowEdge(hMRR2_G_data->GetXaxis()->GetFirst());
  xmax = hMRR2_G_data->GetXaxis()->GetBinUpEdge(hMRR2_G_data->GetXaxis()->GetLast());

  text = "CMS #scale[0.7]{#font[52]{Work in progress 2016}}";
  cms_lat = new TLatex(xmin, 1.9e3, text.c_str());
  cms_lat->SetTextSize(0.07);
  cms_lat->SetLineWidth(2);
  cms_lat->Draw();
  text = "#scale[0.7]{W ana, 35.9 fb^{-1} (13 TeV)}";
  if(channel=="_topana")text = "#scale[0.7]{Top ana, 35.9 fb^{-1} (13 TeV)}";
  era_lat = new TLatex(xmax,2.3e3, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.06);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();
  text = "#scale[0.7]{   [800, 1000]   [1000, 1200]   [1200, 1600]  [1600, 2000]   [2000, 4000]}";
  era_lat = new TLatex(xmax,3.e1, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.06);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();

  leg = new TLegend(0.45,0.62,0.9,0.9);
  //if(TString(njet).Contains("nj35")) leg->SetHeader("S' Closure test, 3 #leq jet < 5");
  //else if(TString(njet).Contains("nj35")) leg->SetHeader("S' Closure test, 6 #leq jet");
  //else leg->SetHeader("S' Closure test");

  //char legentry[10];

  leg->SetNColumns(2);
  leg->SetTextSize(0.03);
  leg->AddEntry(hMRR2_G_data, "data", "lep");
  leg->AddEntry(hMRR2_G_Multijet,  "Multijet", "f");
  leg->AddEntry(hMRR2_G_TTST,  "t#bar{t} + single top", "f");
  leg->AddEntry(hMRR2_G_Wjet,  "W(#rightarrowl#nu)", "f");
  leg->AddEntry(hMRR2_G_Gjet,  "G(#rightarrow#nu#nu)", "f");
  leg->AddEntry(hMRR2_G_Gjet,  "#gamma + jet", "f");
  leg->AddEntry(hMRR2_G_Multiboson,  "VV(V) + t#bar{t}(X)", "f");
  leg->AddEntry(hMRR2_G_T2tt, "T2tt", "l");
  leg->AddEntry(hMRR2_G_T5tttt, "T5tttt", "l");
  leg->AddEntry(hMRR2_G_T5ttcc, "T5ttcc", "l");
  leg->AddEntry(hMRR2_G_T1tttt, "T1tttt", "l");
  leg->Draw();
  
  compareSp = (TH1F*)hMRR2_G_data->Clone(); 
  compareSp1 = (TH1F*)hMRR2_G_data->Clone(); 
  compareSp->Divide(hs);

  canvas_dw->cd();
  gPad->SetGridy();
  compareSp->SetLineColor(kBlack);
  compareSp->SetMarkerColor(kBlack);
  compareSp->SetMarkerStyle(20);
  compareSp->SetMarkerSize(1.3);
  compareSp->GetYaxis()->SetRangeUser(0.0,2.0);
  compareSp->SetTitle(0);
  compareSp->SetStats(0);
  compareSp->GetYaxis()->SetLabelSize(0.09);
  compareSp->GetXaxis()->SetLabelSize(0.09);
  compareSp->GetXaxis()->LabelsOption("v");
  compareSp->GetYaxis()->SetNdivisions(502);
  compareSp->GetXaxis()->SetTickLength(0.07);
  compareSp->GetYaxis()->SetTickLength(0.03);
  compareSp->GetYaxis()->SetTitle("#frac{Data}{Estimate}");
  compareSp->GetXaxis()->SetTitleOffset(0.81);
  compareSp->GetXaxis()->SetTitleSize(0.12);
  compareSp->GetYaxis()->SetTitleOffset(0.35);
  compareSp->GetYaxis()->SetTitleSize(0.10);//dividend->SetTitleSize(0.12,"Y");
  for(int i=1;i<=hMRR2_G_data->GetNbinsX();i++){
    compareSp1->SetBinContent(i,1);
    if(hMRR2_G_data->GetBinContent(i)==0) compareSp1->SetBinError(i,0);
    else compareSp1->SetBinError(i,hs->GetBinError(i)/hs->GetBinContent(i));
  }
  compareSp1->SetFillColor(13);
  compareSp1->SetFillStyle(3235);
  compareSp1->SetMarkerStyle(1);
  compareSp->Draw("EP");
  compareSp1->Draw("E2same");


  hMRR2_PS_Multijet->SetFillColor(619);
  hMRR2_PS_TTST->SetFillColor(633);
  hMRR2_PS_Wjet->SetFillColor(418);
  hMRR2_PS_Zjet->SetFillColor(401);
  hMRR2_PS_Gjet->SetFillColor(803);
  hMRR2_PS_Multiboson->SetFillColor(601);

  hs = (TH1D*)hMRR2_PS_Multijet->Clone();
  hs->Add(hMRR2_PS_TTST);
  hs->Add(hMRR2_PS_Wjet);
  hs->Add(hMRR2_PS_Zjet);
  hs->Add(hMRR2_PS_Gjet);
  hs->Add(hMRR2_PS_Multiboson);

  hs1 = new THStack("hs1","");
  hs1->Add(hMRR2_PS_Multiboson);
  hs1->Add(hMRR2_PS_Wjet);
  hs1->Add(hMRR2_PS_Zjet);
  hs1->Add(hMRR2_PS_Gjet);
  hs1->Add(hMRR2_PS_Multijet);
  hs1->Add(hMRR2_PS_TTST);

  TCanvas* canvas_7 = new TCanvas("canvas_7", "", 800, 800);
  canvas_7->Divide(1,2);
  canvas_up = (TPad*)canvas_7->GetListOfPrimitives()->FindObject("canvas_7_1");
  canvas_dw = (TPad*)canvas_7->GetListOfPrimitives()->FindObject("canvas_7_2");

  up_height     = 0.7; // please tune so that the upper figures size will meet your requirement
  dw_correction = 1.210;//1.390 // please tune so that the smaller canvas size will work in your environment
  dw_height    = (1. - up_height) * dw_correction;

  canvas_up->SetPad(0., 1 - up_height, 1., 1.);
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

  canvas_up->cd();
  hMRR2_PS_data->SetMarkerColor(kBlack);
  hMRR2_PS_data->SetLineColor(kBlack);
  hMRR2_PS_data->SetMarkerStyle(20);
  hMRR2_PS_data->SetMarkerSize(1.3);
  hs1->SetMaximum(1e6);
  hs1->SetMinimum(1);
  hs1->Draw("HIST");
  //hs1->GetYaxis()->SetTitle("Events/bin");
  hMRR2_PS_data->Draw("EPsame");

  xmin = hMRR2_PS_data->GetXaxis()->GetBinLowEdge(hMRR2_PS_data->GetXaxis()->GetFirst());
  xmax = hMRR2_PS_data->GetXaxis()->GetBinUpEdge(hMRR2_PS_data->GetXaxis()->GetLast());

  text = "CMS #scale[0.7]{#font[52]{Work in progress 2016}}";
  cms_lat = new TLatex(xmin, 2.3e6, text.c_str());
  cms_lat->SetTextSize(0.07);
  cms_lat->SetLineWidth(2);
  cms_lat->Draw();
  text = "#scale[0.7]{Preselection, 35.9 fb^{-1} (13 TeV)}";
  era_lat = new TLatex(xmax,2.9e6, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.06);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();
  text = "#scale[0.7]{   [800, 1000]   [1000, 1200]   [1200, 1600]  [1600, 2000]   [2000, 4000]}";
  era_lat = new TLatex(xmax,3.e3, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.06);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();

  leg = new TLegend(0.45,0.62,0.9,0.9);
  //if(TString(njet).Contains("nj35")) leg->SetHeader("S' Closure test, 3 #leq jet < 5");
  //else if(TString(njet).Contains("nj35")) leg->SetHeader("S' Closure test, 6 #leq jet");
  //else leg->SetHeader("S' Closure test");

  //char legentry[10];

  leg->SetNColumns(4);
  leg->SetTextSize(0.03);
  sprintf(legentry, "%.1f",hMRR2_PS_data->Integral());
  leg->AddEntry(hMRR2_PS_data, "data", "lep");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_PS_Multijet->Integral());
  leg->AddEntry(hMRR2_PS_Multijet,  "Multijet", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_PS_TTST->Integral());
  leg->AddEntry(hMRR2_PS_TTST,  "t#bar{t} + single top", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_PS_Wjet->Integral());
  leg->AddEntry(hMRR2_PS_Wjet,  "W(#rightarrowl#nu)", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_PS_Zjet->Integral());
  leg->AddEntry(hMRR2_PS_Zjet,  "Z(#rightarrow#nu#nu)", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_PS_Gjet->Integral());
  leg->AddEntry(hMRR2_PS_Gjet,  "#gamma + jet", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_PS_Multiboson->Integral());
  leg->AddEntry(hMRR2_PS_Multiboson,  "VV(V) + t#bar{t}(X)", "f");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_PS_T2tt->Integral());
  leg->AddEntry(hMRR2_PS_T2tt, "T2tt", "l");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_PS_T5tttt->Integral());
  leg->AddEntry(hMRR2_PS_T5tttt, "T5tttt", "l");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_PS_T5ttcc->Integral());
  leg->AddEntry(hMRR2_PS_T5ttcc, "T5ttcc", "l");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.1f",hMRR2_PS_T1tttt->Integral());
  leg->AddEntry(hMRR2_PS_T1tttt, "T1tttt", "l");
  leg->AddEntry((TObject*)0,legentry,"");
  leg->Draw();
  
  compareSp = (TH1F*)hMRR2_PS_data->Clone(); 
  compareSp1 = (TH1F*)hMRR2_PS_data->Clone(); 
  compareSp->Divide(hs);

  canvas_dw->cd();
  gPad->SetGridy();
  compareSp->SetLineColor(kBlack);
  compareSp->SetMarkerColor(kBlack);
  compareSp->SetMarkerStyle(20);
  compareSp->SetMarkerSize(1.3);
  compareSp->GetYaxis()->SetRangeUser(0.0,2.0);
  compareSp->SetTitle(0);
  compareSp->SetStats(0);
  compareSp->GetYaxis()->SetLabelSize(0.09);
  compareSp->GetXaxis()->SetLabelSize(0.09);
  compareSp->GetXaxis()->LabelsOption("v");
  compareSp->GetYaxis()->SetNdivisions(502);
  compareSp->GetXaxis()->SetTickLength(0.07);
  compareSp->GetYaxis()->SetTickLength(0.03);
  compareSp->GetYaxis()->SetTitle("#frac{Data}{Estimate}");
  compareSp->GetXaxis()->SetTitleOffset(0.81);
  compareSp->GetXaxis()->SetTitleSize(0.12);
  compareSp->GetYaxis()->SetTitleOffset(0.35);
  compareSp->GetYaxis()->SetTitleSize(0.10);//dividend->SetTitleSize(0.12,"Y");
  for(int i=1;i<=hMRR2_PS_data->GetNbinsX();i++){
    compareSp1->SetBinContent(i,1);
    if(hMRR2_PS_data->GetBinContent(i)==0) compareSp1->SetBinError(i,0);
    else compareSp1->SetBinError(i,hs->GetBinError(i)/hs->GetBinContent(i));
  }
  compareSp1->SetFillColor(13);
  compareSp1->SetFillStyle(3235);
  compareSp1->SetMarkerStyle(1);
  compareSp->Draw("EP");
  compareSp1->Draw("E2same");




  //Closure Test

  hs1 = new THStack("hs1","");
  NSp_QCD->SetLineColor(619);
  NSp_QCD->SetLineWidth(4);
  hs1->Add(NSp_QCD);
  NSp_W->SetLineColor(418);
  NSp_W->SetLineWidth(4);
  hs1->Add(NSp_W);
  NSp_Z->SetLineColor(401);
  NSp_Z->SetLineWidth(4);
  hs1->Add(NSp_Z);
  //NSp_G->SetLineColor(803);
  //NSp_G->SetLineWidth(4);
  //hs1->Add(NSp_G);
  NSp_T->SetLineColor(633);
  NSp_T->SetLineWidth(4);
  hs1->Add(NSp_T);

  TCanvas* comSp = new TCanvas("comSp", "", 800, 800);
  comSp->Divide(1,2);
  canvas_up = (TPad*)comSp->GetListOfPrimitives()->FindObject("comSp_1");
  canvas_dw = (TPad*)comSp->GetListOfPrimitives()->FindObject("comSp_2");

  up_height     = 0.7; // please tune so that the upper figures size will meet your requirement
  dw_correction = 1.210;//1.390 // please tune so that the smaller canvas size will work in your environment
  dw_height    = (1. - up_height) * dw_correction;

  canvas_up->SetPad(0., 1 - up_height, 1., 1.);
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

  canvas_up->cd();
  hMRR2_s_data->SetMarkerColor(kBlack);
  hMRR2_s_data->SetLineColor(kBlack);
  hMRR2_s_data->SetMarkerStyle(20);
  hMRR2_s_data->SetMarkerSize(1.3);
  //NSp->SetFillColor(kTeal);
  NSp->SetLineColor(kTeal);
  NSp->SetMarkerColor(kTeal);
  NSp->SetFillStyle(3003);
  NSp->SetMaximum(1e5);
  NSp->SetMinimum(1e-1);
  NSp->GetXaxis()->SetLabelSize(0.01);
  NSp->GetXaxis()->LabelsOption("v");
  //NSp->GetYaxis()->SetTitle("Events/bin");
  NSp->Draw("E2");
  hs1->Draw("HISTsame");
  hMRR2_s_data->Draw("EPsame");

  //xmin = hMRR2_S_data->GetXaxis()->GetBinLowEdge(hMRR2_S_data->GetXaxis()->GetFirst());
  //xmax = hMRR2_S_data->GetXaxis()->GetBinUpEdge(hMRR2_S_data->GetXaxis()->GetLast());

  text = "CMS #scale[0.7]{#font[52]{Work in progress 2016}}";
  cms_lat = new TLatex(xmin, 2.0e5, text.c_str());
  cms_lat->SetTextSize(0.07);
  cms_lat->SetLineWidth(2);
  cms_lat->Draw();
  text = "#scale[0.7]{W ana, 35.9 fb^{-1} (13 TeV)}";
  if(channel=="_topana")text = "#scale[0.7]{Top ana, 35.9 fb^{-1} (13 TeV)}";
  era_lat = new TLatex(xmax,2.2e5, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.06);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();
  text = "#scale[0.7]{   [800, 1000]   [1000, 1200]   [1200, 1600]  [1600, 2000]   [2000, 4000]}";
  era_lat = new TLatex(xmax,5.e2, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.06);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();

  leg = new TLegend(0.45,0.62,0.9,0.9);
  if(TString(njet).Contains("nj35")) leg->SetHeader("S' Closure test, 3 #leq jet < 5");
  else if(TString(njet).Contains("nj6")) leg->SetHeader("S' Closure test, 6 #leq jet");
  else leg->SetHeader("S' Closure test");

  //char legentry[10];

  leg->SetNColumns(2);
  leg->SetTextSize(0.03);
  sprintf(legentry, "%.2f",NSp_QCD->Integral());
  leg->AddEntry(NSp_QCD,  "Multijet estimate", "l");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.2f",NSp_T->Integral());
  leg->AddEntry(NSp_T,  "t#bar{t} + single top estimate", "l");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.2f",NSp_Z->Integral());
  leg->AddEntry(NSp_Z,  "Z(#rightarrow#nu#nu) estimate", "l");
  leg->AddEntry((TObject*)0,legentry,"");
  //sprintf(legentry, "%.2f",NSp_G->Integral());
  //leg->AddEntry(NSp_G,  "#gamma + jet estimate", "l");
  //leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.2f",NSp_W->Integral());
  leg->AddEntry(NSp_W,  "W(#rightarrowl#nu) estimate", "l");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.2f",hMRR2_s_data->Integral());
  leg->AddEntry(hMRR2_s_data, "data", "lp");
  leg->AddEntry((TObject*)0,legentry,"");
  leg->Draw();

  compareSp  = (TH1F*)hMRR2_s_data->Clone();
  compareSp1 = (TH1F*)hMRR2_s_data->Clone();
  compareSp->Divide(NSp);

  canvas_dw->cd();
  gPad->SetGridy();
  compareSp->SetLineColor(kBlack);
  compareSp->SetMarkerColor(kBlack);
  compareSp->SetMarkerStyle(20);
  compareSp->SetMarkerSize(1.3);
  compareSp->GetYaxis()->SetRangeUser(0.0,2.0);
  compareSp->SetTitle(0);
  compareSp->SetStats(0);
  compareSp->GetYaxis()->SetLabelSize(0.09);
  compareSp->GetXaxis()->SetLabelSize(0.09);
  compareSp->GetXaxis()->LabelsOption("v");
  compareSp->GetYaxis()->SetNdivisions(502);
  compareSp->GetXaxis()->SetTickLength(0.07);
  compareSp->GetYaxis()->SetTickLength(0.03);
  compareSp->GetYaxis()->SetTitle("#frac{Data}{Estimate}");
  compareSp->GetXaxis()->SetTitleOffset(0.81);
  compareSp->GetXaxis()->SetTitleSize(0.12);
  compareSp->GetYaxis()->SetTitleOffset(0.35);
  compareSp->GetYaxis()->SetTitleSize(0.10);//dividend->SetTitleSize(0.12,"Y");
  for(int i=1;i<=NSp->GetNbinsX();i++){
    compareSp1->SetBinContent(i,1);
    if(NSp->GetBinContent(i)==0) compareSp1->SetBinError(i,0);
    else compareSp1->SetBinError(i,NSp->GetBinError(i)/NSp->GetBinContent(i));
  }
  compareSp1->SetFillColor(13);
  compareSp1->SetFillStyle(3235);
  compareSp1->SetMarkerStyle(1);

  compareSp->Draw("EP");
  compareSp1->Draw("E2same");

  TCanvas* comQp = new TCanvas("comQp", "", 800, 800);
  comQp->Divide(1,2);
  TPad *canvsa_up = (TPad*)comQp->GetListOfPrimitives()->FindObject("comQp_1");
  TPad *canvsa_dw = (TPad*)comQp->GetListOfPrimitives()->FindObject("comQp_2");

  up_height     = 0.7; // please tune so that the upper figures size will meet your requirement
  dw_correction = 1.210;//1.390 // please tune so that the smaller canvsa size will work in your environment
  dw_height    = (1. - up_height) * dw_correction;

  canvsa_up->SetPad(0., 1 - up_height, 1., 1.);
  canvsa_dw->SetPad(0., 0., 1., dw_height);
  canvsa_up->SetLogy();
  canvsa_up->SetFrameFillColor(0);
  canvsa_up->SetFillColor(0);
  canvsa_dw->SetFillColor(0);
  canvsa_dw->SetFrameFillColor(0);
  canvsa_up->SetLeftMargin(0.1);
  canvsa_dw->SetTopMargin(0.03);//0.03
  canvsa_dw->SetLeftMargin(0.1);
  canvsa_dw->SetBottomMargin(0.3);
  canvsa_up->cd();
  hMRR2_q_data->SetMarkerColor(kBlack);
  hMRR2_q_data->SetLineColor(kBlack);
  hMRR2_q_data->SetMarkerStyle(20);
  hMRR2_q_data->SetMarkerSize(1.3);
  NQp->SetFillColor(kTeal);
  NQp->SetLineColor(kTeal);
  NQp->SetMarkerColor(kTeal);
  NQp->SetFillStyle(3003);
  NQp->GetXaxis()->SetLabelSize(0.01);
  NQp->GetXaxis()->LabelsOption("v");
  NQp->SetMaximum(1e4);
  NQp->SetMinimum(1e-1);
  //NQp->GetYaxis()->SetTitle("Events/bin");
  NQp->Draw("E2");
  hMRR2_q_data->Draw("EPsame");

  xmin = NQp->GetXaxis()->GetBinLowEdge(NQp->GetXaxis()->GetFirst());
  xmax = NQp->GetXaxis()->GetBinUpEdge(NQp->GetXaxis()->GetLast());

  text = "CMS #scale[0.7]{#font[52]{Work in progress 2016}}";
  cms_lat = new TLatex(xmin, 2.0e4, text.c_str());
  cms_lat->SetTextSize(0.07);
  cms_lat->SetLineWidth(2);
  cms_lat->Draw();
  text = "#scale[0.7]{W ana, 35.9 fb^{-1} (13 TeV)}";
  if(channel=="_topana")text = "#scale[0.7]{Top ana, 35.9 fb^{-1} (13 TeV)}";
  era_lat = new TLatex(xmax,2.2e4, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.06);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();
  text = "#scale[0.7]{   [800, 1000]   [1000, 1200]   [1200, 1600]  [1600, 2000]   [2000, 4000]}";
  era_lat = new TLatex(xmax,1.e2, text.c_str());
  era_lat->SetTextAlign(32);
  era_lat->SetTextSize(0.06);
  era_lat->SetTextFont(42);
  era_lat->SetLineWidth(2);
  era_lat->Draw();

  leg = new TLegend(0.45,0.70,0.9,0.9);
  if(TString(njet).Contains("nj35")) leg->SetHeader("Q' Closure test, 3 #leq jet < 5");
  else if(TString(njet).Contains("nj6")) leg->SetHeader("Q' Closure test, 6 #leq jet");
  else leg->SetHeader("Q' Closure test");

  leg->SetNColumns(2);
  leg->SetTextSize(0.03);
  sprintf(legentry, "%.2f",NQp->Integral());
  leg->AddEntry(NQp,  "Total background estimate", "LPE");
  leg->AddEntry((TObject*)0,legentry,"");
  sprintf(legentry, "%.2f",hMRR2_q_data->Integral());
  leg->AddEntry(hMRR2_q_data, "data", "lp");
  leg->AddEntry((TObject*)0,legentry,"");
  leg->Draw();

  TH1F* compareQp  = (TH1F*)hMRR2_q_data->Clone();
  TH1F* compareQp1 = (TH1F*)hMRR2_q_data->Clone();
  compareQp->Divide(NQp);
  canvsa_dw->cd();
  gPad->SetGridy();
  compareQp->SetLineColor(kBlack);
  compareQp->SetMarkerColor(kBlack);
  compareQp->SetMarkerStyle(20);
  compareQp->SetMarkerSize(1.3);
  compareQp->GetYaxis()->SetRangeUser(0.0,2.0);
  compareQp->SetTitle(0);
  compareQp->SetStats(0);
  compareQp->GetYaxis()->SetLabelSize(0.09);
  compareQp->GetXaxis()->SetLabelSize(0.09);
  compareQp->GetXaxis()->LabelsOption("v");
  compareQp->GetYaxis()->SetNdivisions(502);
  compareQp->GetXaxis()->SetTickLength(0.07);
  compareQp->GetYaxis()->SetTickLength(0.03);
  compareQp->GetYaxis()->SetTitle("#frac{Data}{Estimate}");
  compareQp->GetXaxis()->SetTitleOffset(0.81);
  compareQp->GetXaxis()->SetTitleSize(0.12);
  compareQp->GetYaxis()->SetTitleOffset(0.35);
  compareQp->GetYaxis()->SetTitleSize(0.10);//dividend->SetTitleSize(0.12,"Y");
  for(int i=1;i<=NQp->GetNbinsX();i++){
    compareQp1->SetBinContent(i,1);
    if(NQp->GetBinContent(i)==0) compareQp1->SetBinError(i,0);
    else compareQp1->SetBinError(i,NQp->GetBinError(i)/NQp->GetBinContent(i));
  }
  compareQp1->SetFillColor(13);
  compareQp1->SetFillStyle(3235);
  compareQp1->SetMarkerStyle(1);

  compareQp->Draw("EP");
  compareQp1->Draw("E2same");


  TString Dir = "./plot2/";
  comSp->SaveAs(Dir+"ClosureTest_MRR2_s"+channel+njet+".png");
  comQp->SaveAs(Dir+"ClosureTest_MRR2_q"+channel+njet+".png");

}

void MRR2_universialRatio(){
  ClosureTest("","");
  ClosureTest("_nj35","");
  ClosureTest("_nj6","");
  ClosureTest("","_topana");
}
