tuple<vector<int>, vector<double>> LastBinMerging(TH1D* bkg, TH1D* sig){
  TH1::SetDefaultSumw2();
  vector<int> bins;
  vector<double> val;
  bins.clear();
  val.clear();
  double pre=0, nbkg=0, nsig=0, value=0;
  for(double i=bkg->GetNbinsX();i>0;i--){
    nbkg += bkg->GetBinContent(i);
    nsig += sig->GetBinContent(i);
    value = nsig/(TMath::Sqrt(nsig+nbkg));
    //value = sqrt(2*log(1+nsig/nbkg)-2*nsig);
    if(value < pre && nbkg > 1) {
      bins.push_back(i+1);
      val.push_back(pre);
      nbkg = bkg->GetBinContent(i); nsig = sig->GetBinContent(i); value = nsig/(TMath::Sqrt(nsig+nbkg)); pre = value;
    } else pre = value;
    //cout << i << ", " << nbkg << ", " << nsig << ", " << value << endl;
  }
  return make_tuple(bins, val);
}

std::vector<int> comb(int N, int K) {
  std::string bitmask(K, 1);
  bitmask.resize(N, 0);
  std::vector<int> store;
  do {
    for (int i = 0; i < N; ++i) {
      if (bitmask[i]) store.push_back(i+1);
    }
  } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
  return store;
}

//tuple<vector<int>, vector<double>> BinGrouping(TH1D* bkg, TH1D* sig){
vector<int> BinGrouping(TH1D* bkg, TH1D* sig){
//void BinGrouping(TH1D* bkg, TH1D* sig){
  TH1::SetDefaultSumw2();
  vector<int> bins;
  vector<int> val;
  vector<int> tmpbins;
  bins.clear();
  tmpbins.clear();
  val.clear();
  double store=0;
  double pre=0, nbkg=0, nsig=0, value=0;
  bool reset = false;

  const int N = bkg->GetNbinsX();
  for(int i=1;i<N;i++){
    tmpbins = comb(N-1, i);
    reset = false;
    val.clear();
    value = 0;
    for(int j=0;j<tmpbins.size();j++) {
      if(j==0 || reset) {
        nbkg = bkg->Integral(1,tmpbins.at(j));
        nsig = sig->Integral(1,tmpbins.at(j));
      } else {
        if((j+1)%i == 0) {
          nbkg = bkg->Integral(tmpbins.at(j-1),N);
          nsig = sig->Integral(tmpbins.at(j-1),N);
        } else{
          nbkg = bkg->Integral(tmpbins.at(j-1),tmpbins.at(j));
          nsig = sig->Integral(tmpbins.at(j-1),tmpbins.at(j));
        }
      }
      if(nbkg < 1) value = -999999; 
      else         value += nsig/(TMath::Sqrt(nsig+nbkg));
      //else         value += sqrt(2*log(1+nsig/nbkg)-2*nsig);
      if(isnan(value)) cout << "warning, nan value" << endl;
      val.push_back(tmpbins.at(j));

      if((j+1)%i == 0) {
        if(store == value) {
          bins = bins.size() > val.size() ? val : bins;
        } else if(store < value) {
          store = value;
          bins = val;
          val.clear();
        }
        reset = true; value = 0;
      } else{
        reset = false;
      }
      //cout << " " << bins.at(j);
      //if((j+1)%i == 0) cout << endl;
    }
  }
  cout << store << endl;
  //return make_tuple(bins, val);
  return bins;
}


void CalcSign(TString region, TString sample){
  cout << endl << "region : " << region << endl;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();
  TString dir = "/Users/chuh/Dropbox/Analysis/razor/";
  TFile* file0 = TFile::Open(dir+sample);


  TString histname;
  TString path = "/Background";
  histname = path+"/MRR2_bkg_"+region;
  //histname = path+"/MRR2_bkg_"+region+"_new";

  TH2D* bkg2D;
  TH1D* bkg;

  bkg2D = (TH2D*)file0->Get(histname);
  int binsize = bkg2D->GetNbinsX();
  double bin[binsize+1];

  for(int i=1;i<=bkg2D->GetNbinsX()+1;i++) bin[i-1] = bkg2D->GetXaxis()->GetBinLowEdge(i);

  bkg = new TH1D("bkg", "", binsize, bin);
  for(int k=1;k<=bkg2D->GetNbinsX();k++) {
    if(bkg2D == NULL) {
      bkg->SetBinContent(k, 0);
      bkg->SetBinError(k, 0);
    }
    else {
      bkg->SetBinContent(k, bkg2D->GetBinContent(k,1));
      bkg->SetBinError(k,   bkg2D->GetBinError(k,1));
    }
  }

/* 
  //Signal name
  T5qqqqWH
  T5bbbbZH
  T5ttcc
  T6ttZH
  TChiWZ
  TChiWW
  GluinoRPV
  SbottomRPV
*/

  TH1D* sig;
  TH1D* tmp;
  file0->cd("Signal2");
  TIter keyList(gDirectory->GetListOfKeys());
  TKey *key;
  TH2D* h1;
  TString name;
  bool flag = true;
  int count=0;

  while ((key = (TKey*)keyList())) {
    delete gROOT->FindObject("sigtemp");
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH2D")) continue;
    h1 = (TH2D*)key->ReadObj();
    name = h1->GetName();
    if(!name.Contains(region)) continue;
    if(name.Contains("_new")) continue;
    //if(!name.Contains("_new")) continue;

		tmp = new TH1D("sigtemp", "", binsize, bin);
		for(int k=1;k<=h1->GetNbinsX();k++) {
		  tmp->SetBinContent(k, h1->GetBinContent(k,1));
		  tmp->SetBinError(k,   h1->GetBinError(k,1));
		}
    count++;
    if(flag) {
      sig = (TH1D*)tmp->Clone();
      sig->SetName("Signal");
      flag = false;
    } else {
      sig->Add(tmp);
    }
  }

  //cout << count << " number of mass bin" << endl;
  sig->Scale(1./count);

  tuple<vector<int>, vector<double>> Merge;
  vector<int> bins;
  vector<double> valu;
  bins.clear();valu.clear();
  if(bkg->Integral() == 0) {cout << "No background events in this region" << endl; return 0;}
  else cout << bkg->Integral()  << " events in this region" << endl;
  if(sig->Integral() == 0) {cout << "No signal events in this region" << endl; return 0;}
  else cout << sig->Integral()  << " events in this region" << endl;

/*
  cout << "Last Bin Merging" << endl;
  Merge = LastBinMerging(bkg, sig);
  bins = get<0>(Merge);
  valu = get<1>(Merge);

  cout << bins.size() << ", " << valu.size() << endl;
  for(int i=0;i<bins.size();++i) cout << bins.at(i) << ", ";
  cout << endl;
  for(int i=0;i<valu.size();++i) cout << valu.at(i) << ", ";
  cout << endl;
*/
  double Value[30] = {0., 100., 200., 300., 400., 500., 600., 700., 800., 900., 1000., 1100., 1200., 1300., 1400., 3000.};
  //double Value[30] = {0, 30, 40, 50, 75, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 3000.};

  cout << "Bin Grouping" << endl;
  bins.clear();
  bins = BinGrouping(bkg, sig);
 
  if(bins.size() == 0) cout << "use big bin" << endl;
  cout << bins.size()+1 << endl;
  for(int i=0;i<bins.size();++i) cout << bins.at(i) << ", ";
  cout << endl << "0, ";
  for(int i=0;i<bins.size();++i) cout << Value[bins.at(i)] << ", ";
  cout << "3000" << endl;

  file0->Close();
}

void bin_optimize(){
  //TString SR[30] = {"SR_Had_1H_0b_34j", "SR_Had_1H_0b_5j", "SR_Had_1V_0b_34j", "SR_Had_1V_0b_5j", "SR_Had_1htop", "SR_Had_2H_0b_34j", "SR_Had_2H_0b_5j", "SR_Had_2H_b_6j", "SR_Had_2V_0b_24j", "SR_Had_2V_0b_5j", "SR_Had_2htop", "SR_Had_HV_0b_24j", "SR_Had_HV_0b_5j", "SR_Had_HV_b_6j", "SR_Had_H_b_45j", "SR_Had_H_b_6j", "SR_Had_V_b_45j", "SR_Had_V_b_6j", "SR_Lep_1htop", "SR_Lep_H_0b", "SR_Lep_H_b", "SR_Lep_V_0b", "SR_Lep_V_b", "SR_Lepjet_0V_24j", "SR_Lepjet_0V_5j", "SR_Lepjet_1V_24j", "SR_Lepjet_1V_5j", "SR_Leptop_0htop", "SR_Leptop_1htop"};
  TString SR[30] = {"SR_Had_H_0b_34j", "SR_Had_H_0b_5j", "SR_Had_1V_0b_34j", "SR_Had_1V_0b_5j", "SR_Had_1htop", "SR_Had_2V_0b_24j", "SR_Had_2V_0b_5j", "SR_Had_2htop", "SR_Had_HV_0b_24j", "SR_Had_HV_0b_5j", "SR_Had_HV_b_6j", "SR_Had_H_b_45j", "SR_Had_H_b_6j", "SR_Had_V_b_45j", "SR_Had_V_b_6j", "SR_Lep_1htop", "SR_Lep_H_0b", "SR_Lep_H_b", "SR_Lep_V_0b", "SR_Lep_V_b", "SR_Lepjet_0V_24j", "SR_Lepjet_0V_5j", "SR_Lepjet_1V_24j", "SR_Lepjet_1V_5j", "SR_Leptop_0htop", "SR_Leptop_1htop"};
  TString sample = "230131/run_2023_03_16.root";
  for(int i=0;i<26;i++){
  //for(int i=0;i<5;i++){
    CalcSign(SR[i], sample);
  }
  //CalcSign(SR[3], sample);
  //CalcSign(SR[4], sample);
}
