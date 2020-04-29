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
TH1D* BigBin(TH1D* target){
	TH1::SetDefaultSumw2();
	//it is only for MC in control region
	//hist scanning 
	TH1D* hist = (TH1D*)target->Clone();
	TH1D* htemp = new TH1D("htemp","",3,0,3);
	std::vector<int> emptybin;
	int bin,j=1;
	double temp;
	for(int i=1;i<=target->GetNbinsX();i++){
		if(target->GetBinContent(i)<0.1) emptybin.push_back(i);
	}
  bool trig = false;
	if(emptybin.size() != 0) {
		//for(unsigned i = emptybin.size() - 1; emptybin.size() > i; --i){
		for(unsigned i = 0; emptybin.size() > i; ++i){
      if(trig) {trig = false; continue;}
      if(emptybin.size()==1){
			  htemp = new TH1D("htemp","",2,0,2);
        htemp->SetBinContent(2, target->GetBinContent(emptybin[i]));
        htemp->SetBinContent(1, target->GetBinContent(emptybin[i]-1));
        htemp->SetBinError(2, target->GetBinError(emptybin[i]));
        htemp->SetBinError(1, target->GetBinError(emptybin[i]-1));
			  htemp->Rebin(2);
			  hist->SetBinContent(emptybin[i], htemp->GetBinContent(1));
		  	hist->SetBinError(emptybin[i], htemp->GetBinError(1));
			  hist->SetBinContent(emptybin[i]-1, htemp->GetBinContent(1));
		  	hist->SetBinError(emptybin[i]-1, htemp->GetBinError(1));
        cout<< emptybin[i] << ", " << htemp->GetBinContent(1) << ", " << htemp->GetBinError(1) << endl;
      }
      else{
        if(emptybin[i+1]-emptybin[i]==1){
        trig = true;
			  htemp = new TH1D("htemp","",3,0,3);
        htemp->SetBinContent(3, target->GetBinContent(emptybin[i]));
        htemp->SetBinContent(2, target->GetBinContent(emptybin[i]-1));
        htemp->SetBinContent(1, target->GetBinContent(emptybin[i]-2));
        htemp->SetBinError(2, target->GetBinError(emptybin[i]));
        htemp->SetBinError(1, target->GetBinError(emptybin[i]-1));
        htemp->SetBinError(1, target->GetBinError(emptybin[i]-1));
			  htemp->Rebin(2);
			  hist->SetBinContent(emptybin[i], htemp->GetBinContent(1));
		  	hist->SetBinError(emptybin[i], htemp->GetBinError(1));
			  hist->SetBinContent(emptybin[i]-1, htemp->GetBinContent(1));
		  	hist->SetBinError(emptybin[i]-1, htemp->GetBinError(1));
			  hist->SetBinContent(emptybin[i]-2, htemp->GetBinContent(1));
		  	hist->SetBinError(emptybin[i]-2, htemp->GetBinError(1));
        cout<< emptybin[i] << ", " << htemp->GetBinContent(1) << ", " << htemp->GetBinError(1) << endl;
      }
      else{
			  htemp = new TH1D("htemp","",2,0,2);
        htemp->SetBinContent(2, target->GetBinContent(emptybin[i]));
        htemp->SetBinContent(1, target->GetBinContent(emptybin[i]-1));
        htemp->SetBinError(2, target->GetBinError(emptybin[i]));
        htemp->SetBinError(1, target->GetBinError(emptybin[i]-1));
			  htemp->Rebin(2);
			  hist->SetBinContent(emptybin[i], htemp->GetBinContent(1));
		  	hist->SetBinError(emptybin[i], htemp->GetBinError(1));
			  hist->SetBinContent(emptybin[i]-1, htemp->GetBinContent(1));
		  	hist->SetBinError(emptybin[i]-1, htemp->GetBinError(1));
        cout<< emptybin[i] << ", " << htemp->GetBinContent(1) << ", " << htemp->GetBinError(1) << endl;
      }
      }
		}
	}
  delete gROOT->FindObject("htemp");
	return hist;
}
TH1D* Extrapolation(TH1D* target, TH1D* loose, TString channel, TString njet, TString sample){
	TH1::SetDefaultSumw2();
	TH1D* htemp = (TH1D*)target->Clone();
  //htemp->Clone(target);
	htemp->Divide(loose);
	TH1D* hResult = (TH1D*)target->Clone();
  //hResult->Clone(target);

	std::vector<TF1*> Gaus;
	std::vector<double> mean, Error;

	for(int i=1;i<=htemp->GetNbinsX();i++){
    Gaus.push_back(new TF1((std::string("Gaus_"+to_string(i))).c_str(), "gaus", -5,5));
  	mean.push_back(htemp->GetBinContent(i));
	  Error.push_back(htemp->GetBinError(i));
    //if(htemp->GetBinContent(i) == 0) noeve = true;
	}
	for(int i=0;i<Gaus.size();++i) Gaus[i]->SetParameters(1,mean[i],Error[i]);
	TF1* f1 = new TF1("f1", "gaus(0)+gaus(3)+gaus(6)+gaus(9)+gaus(12)+gaus(15)+gaus(18)+gaus(21)+gaus(24)+gaus(27)+gaus(30)+gaus(33)+gaus(36)+gaus(39)+gaus(42)+gaus(45)+gaus(48)+gaus(51)+gaus(54)+gaus(57)+gaus(60)+gaus(63)", -5,5);
  int temp;
  //cout << "fill parameter" << endl;
	/*for(int i=0;i<f1->GetNpar();i++){
    temp = (int)i/3;
		     if(i%3==0) f1->SetParameter(i,1);
		else if(i%3==1) f1->SetParameter(i,Gaus[temp+1]->GetParameter(1));
		else if(i%3==2) f1->SetParameter(i,Gaus[temp+1]->GetParameter(2));
	}*/
	//f1->SetParameter(0,1);
	for(int i=0;i<Gaus.size();++i){
    f1->SetParameter(3*i+0,Gaus[i]->GetParameter(0));
    f1->SetParameter(3*i+1,Gaus[i]->GetParameter(1));
    f1->SetParameter(3*i+2,Gaus[i]->GetParameter(2));
	}
	const int nq = 100000;
	Double_t xq[nq], yq[nq];
	double nomin=0, error=0, error1=0, error2=0;
	for(int i=0;i<nq+1;i++) xq[i] = (Double_t)(i-(Double_t)nq/2)/nq*6;
	f1->GetQuantiles(nq-35000,yq,xq);
	for(int i=0;i<nq;i++){
		if(yq[i] > 0.1745 && yq[i] < 0.1755) error1 = xq[i]; 
		if(yq[i] > 0.499 && yq[i] < 0.501) nomin  = xq[i]; 
		if(yq[i] > 0.8245 && yq[i] < 0.8255) error2 = xq[i]; 
	}
	if(error1 && error2) error = (error1+error2)/2;
	else                 error = error1+error2;
  cout << nomin << ", " << error1 << ", " << error2 << ", " << error << endl;
	/*
	TF1Convolution* conv;
	TF1* temp;

	conv = new TF1Convolution(gaus[0], gaus[1],-10,10,true);
	conv->SetNofPointsFFT(1000);

	for(int i=2; i<Gaus.size();++i){
		temp = new TF1("temp", *conv, -5,5,conv->GetNpar());
		conv = new TF1Convolution(temp, Gaus[i],-10,10,true);
		conv->SetNofPointsFFT(1000);
	}

	TF1* f1 = new TF1("f1", *conv, -5,5,conv->GetNpar());
	for(int i=0;i<conv->GetNpar();i++){
		if(i==0) f1->SetParameter(i,1);
		else if(i%2==1) f1->SetParameter(i,Gaus[i/2]->GetParameter(1));
		else if(i%2==0) f1->SetParameter(i,Gaus[i/2-1]->GetParameter(2));
		//cout << f1->GetParameter(i) << endl;
	}
	f1->SetNpx(1000);
	//this is temporary I'll use gaussian fit results(mean, sigma)
	const int nq = 100000;
	Double_t xq[nq], yq[nq];
	double nomin, error, error1=0, error2=0;
	for(int i=0;i<nq+1;i++) xq[i] = (Double_t)(i-(Double_t)nq/2.)/nq*6;
	f1->GetQuantiles(nq-40000,yq,xq);
	for(int i=0;i<nq;i++){
		if(yq[i] > 0.1749 && yq[i] < 0.1751)  nomin = xq[i]; 
		if(yq[i] > 0.4999 && yq[i] < 0.5001) error1 = xq[i]; 
		if(yq[i] > 0.8249 && yq[i] < 0.8251) error2 = xq[i]; 
	}
	if(error1 && error2) error = (error1+error2)/2.;
	else                 error = error1+error2;
	//This is gaussian fit results from convolution's
	TH1D* h1 = new TH1D("h1", "", 100,-3,3);
	h1->FillRandom("f1", 1e6);
	TF1* Gaus = new TF1("Gaus","gaus(0)",-3,3);
	Gaus->SetParameters(1,0.5,0.5);
	h1->Fit(Gaus);
	nomin = Gaus->GetParameter(1);
	error = Gaus->GetParameter(2);
	*/
  int count=0;
	for(int i=1;i<=htemp->GetNbinsX();i++){
		if(htemp->GetBinContent(i) != 0){
			hResult->SetBinContent(i, htemp->GetBinContent(i));
			hResult->SetBinError(i, htemp->GetBinError(i));
		}
    else {
  //cout << nomin << ", " << error1 << ", " << error2 << ", " << error << endl;
			hResult->SetBinContent(i, nomin);
			hResult->SetBinError(i, error);
      count++;
		}
	}
	hResult->Multiply(loose);
	//for(int i=1;i<=htemp->GetNbinsX();i++) cout << hResult->GetBinContent(i) << ", " << hResult->GetBinError(i) << endl;
  TCanvas* c1 = new TCanvas("c1","",700,700);
  f1->SetNpx(10000);
  f1->Draw();
  for(int i=0;i<Gaus.size();++i) Gaus[i]->Draw("same");
  //TFile* f1 = new TFile("./plot/"+target->GetTitle()+".root","recreate");
  if(count>0)c1->SaveAs("./plot/Extrapolation"+channel+njet+sample+".png","recreate");
	return hResult;
}
pair<double, double> PurityFit(TH1D* chisoFake, TH1D* chisoPrompt, TH1D* chiso){
	TH1::SetDefaultSumw2();
  cout << chiso->GetNbinsX() << endl;
	TObjArray *mc;
	TFractionFitter* fit;
	TH1F* result;
	double purity=0, error=0;
	mc = new TObjArray(3);
	mc->Add(chisoFake);
	mc->Add(chisoPrompt);
	fit = new TFractionFitter(chiso,mc);
	fit->Constrain(0,0.,0.5);
	fit->Constrain(1,0.5,1.);
	double prompt_value=0., prompt_error=0., fake_value=0., fake_error=0.;
	int status = fit->Fit();
	if(status == 0){
		result = (TH1F*)fit->GetPlot();
		fit->GetResult(0, fake_value, fake_error);
		fit->GetResult(1, prompt_value, prompt_error);
		//cout << "prmpt : " << prompt_value << " fake_value : " << fake_value << endl;
		purity = (prompt_value/(prompt_value+fake_value));
		double temp = TMath::Power(fake_value,2)/(TMath::Power(fake_value,2) + TMath::Power(prompt_value,2));
		error = TMath::Sqrt(temp*(TMath::Power(prompt_error,2)+TMath::Power(fake_error,2)));
	}
  TFile* f1 = new TFile("./plot/test.root", "update");
  result->Write();
  f1->Close();
  cout << purity << ", " << error << endl;
	return make_pair(purity,error);
}

tuple<TH1D*, TH1D*, TH1D*, TH1D*> Purity2D(TFile* f0, TFile* f1, TString njet){
	TH1::SetDefaultSumw2();
	//template - fake
	TH2D* htemp0_EB = (TH2D*)f0->Get("ChIsoTemplate_Fake_g_R2_EB");
	TH2D* htemp0_EE = (TH2D*)f0->Get("ChIsoTemplate_Fake_g_R2_EE");
	TH2D* htemp1_EB = (TH2D*)f1->Get("ChIsoTemplate_Fake_g_R2_EB_MC");
	TH2D* htemp1_EE = (TH2D*)f1->Get("ChIsoTemplate_Fake_g_R2_EE_MC");
	//template - prompt
	TH2D* htemp2_EB = (TH2D*)f1->Get("ChIsoTemplate_Prompt_g_R2_EB");
	TH2D* htemp2_EE = (TH2D*)f1->Get("ChIsoTemplate_Prompt_g_R2_EE");
	//fitting target
	TH2D* htemp3_EB = (TH2D*)f0->Get("R2_ChIso_gNoIso_EB"+njet);
	TH2D* htemp3_EE = (TH2D*)f0->Get("R2_ChIso_gNoIso_EE"+njet);

	std::vector<TH1D*> chisoFake_R2_EB;
	std::vector<TH1D*> chisoFakeMC_R2_EB;
	std::vector<TH1D*> chisoPrompt_R2_EB;
	std::vector<TH1D*> chiso_R2_EB;
	std::vector<TH1D*> chisoFake_R2_EE;
	std::vector<TH1D*> chisoFakeMC_R2_EE;
	std::vector<TH1D*> chisoPrompt_R2_EE;
	std::vector<TH1D*> chiso_R2_EE;

	for(int i=1;i<=htemp0_EB->GetNbinsX();i++){
    chisoFake_R2_EB.push_back(new TH1D((std::string("chisoFake_R2_EB_"+to_string(i))).c_str(), "CHIso", 20,0,1));
    chisoFakeMC_R2_EB.push_back(new TH1D((std::string("chisoFakeMC_R2_EB_"+to_string(i))).c_str(), "CHIso", 20,0,1));
    chisoPrompt_R2_EB.push_back(new TH1D((std::string("chisoPrompt_R2_EB_"+to_string(i))).c_str(), "CHIso", 20,0,1));
    chiso_R2_EB.push_back(new TH1D((std::string("chiso_R2_EB_"+to_string(i))).c_str(), "CHIso", 20,0,1));
    chisoFake_R2_EE.push_back(new TH1D((std::string("chisoFake_R2_EE_"+to_string(i))).c_str(), "CHIso", 20,0,1));
    chisoFakeMC_R2_EE.push_back(new TH1D((std::string("chisoFakeMC_R2_EE_"+to_string(i))).c_str(), "CHIso", 20,0,1));
    chisoPrompt_R2_EE.push_back(new TH1D((std::string("chisoPrompt_R2_EE_"+to_string(i))).c_str(), "CHIso", 20,0,1));
    chiso_R2_EE.push_back(new TH1D((std::string("chiso_R2_EE_"+to_string(i))).c_str(), "CHIso", 20,0,1));
	}
	for(int i=1;i<=htemp0_EB->GetNbinsX();i++){
	  for(int j=1;j<=htemp0_EB->GetNbinsY();j++){
      chisoFake_R2_EB[i-1]->SetBinContent(j,htemp0_EB->GetBinContent(i,j));
      chisoFake_R2_EB[i-1]->SetBinError(j,htemp0_EB->GetBinError(i,j));
      chisoFakeMC_R2_EB[i-1]->SetBinContent(j,htemp1_EB->GetBinContent(i,j));
      chisoFakeMC_R2_EB[i-1]->SetBinError(j,htemp1_EB->GetBinError(i,j));
      chisoPrompt_R2_EB[i-1]->SetBinContent(j,htemp2_EB->GetBinContent(i,j));
      chisoPrompt_R2_EB[i-1]->SetBinError(j,htemp2_EB->GetBinError(i,j));
      chiso_R2_EB[i-1]->SetBinContent(j,htemp3_EB->GetBinContent(i,j));
      chiso_R2_EB[i-1]->SetBinError(j,htemp3_EB->GetBinError(i,j));
      chisoFake_R2_EE[i-1]->SetBinContent(j,htemp0_EE->GetBinContent(i,j));
      chisoFake_R2_EE[i-1]->SetBinError(j,htemp0_EE->GetBinError(i,j));
      chisoFakeMC_R2_EE[i-1]->SetBinContent(j,htemp1_EE->GetBinContent(i,j));
      chisoFakeMC_R2_EE[i-1]->SetBinError(j,htemp1_EE->GetBinError(i,j));
      chisoPrompt_R2_EE[i-1]->SetBinContent(j,htemp2_EE->GetBinContent(i,j));
      chisoPrompt_R2_EE[i-1]->SetBinError(j,htemp2_EE->GetBinError(i,j));
      chiso_R2_EE[i-1]->SetBinContent(j,htemp3_EE->GetBinContent(i,j));
      chiso_R2_EE[i-1]->SetBinError(j,htemp3_EE->GetBinError(i,j));
    }
  }
	std::vector<double> purity_R2_EB;
	std::vector<double> purity_R2_EE;
	std::vector<double> purity_R2_EB_MC;
	std::vector<double> purity_R2_EE_MC;
	std::vector<double> error_R2_EB;
	std::vector<double> error_R2_EE;
	std::vector<double> error_R2_EB_MC;
	std::vector<double> error_R2_EE_MC;
	pair<double, double> result;

  cout << "Chiso Fit" << endl;
	for(int i=0;i<chisoFake_R2_EB.size();++i){
		result = PurityFit(chisoFake_R2_EB[i], chisoPrompt_R2_EB[i], chiso_R2_EB[i]);
		purity_R2_EB.push_back(result.first); error_R2_EB.push_back(result.second);
		result = PurityFit(chisoFakeMC_R2_EB[i], chisoPrompt_R2_EB[i], chiso_R2_EB[i]);
		purity_R2_EB_MC.push_back(result.first); error_R2_EB_MC.push_back(result.second);
		result = PurityFit(chisoFake_R2_EE[i], chisoPrompt_R2_EE[i], chiso_R2_EE[i]);
		purity_R2_EE.push_back(result.first); error_R2_EE.push_back(result.second);
		result = PurityFit(chisoFakeMC_R2_EE[i], chisoPrompt_R2_EE[i], chiso_R2_EE[i]);
		purity_R2_EE_MC.push_back(result.first); error_R2_EE_MC.push_back(result.second);
	}
  cout << "end Fit" << endl;

	TH1D* hist_R2_EB = new TH1D("hist_R2_EB", "", 5, 0, 5);
	TH1D* hist_R2_EE = new TH1D("hist_R2_EE", "", 5, 0, 5);
	TH1D* hist_R2_EB_MC = new TH1D("hist_R2_EB_MC", "", 5, 0, 5);
	TH1D* hist_R2_EE_MC = new TH1D("hist_R2_EE_MC", "", 5, 0, 5);

	for(int i=1;i<=5;i++){
		hist_R2_EB->SetBinContent(i,purity_R2_EB[i-1]);
		hist_R2_EB->SetBinError(i,error_R2_EB[i-1]);
		hist_R2_EB_MC->SetBinContent(i,purity_R2_EB_MC[i-1]);
		hist_R2_EB_MC->SetBinError(i,error_R2_EB_MC[i-1]);
		hist_R2_EE->SetBinContent(i,purity_R2_EE[i-1]);
		hist_R2_EE->SetBinError(i,error_R2_EE[i-1]);
		hist_R2_EE_MC->SetBinContent(i,purity_R2_EE_MC[i-1]);
		hist_R2_EE_MC->SetBinError(i,error_R2_EE_MC[i-1]);
	}

	return make_tuple(hist_R2_EB, hist_R2_EE, hist_R2_EB_MC, hist_R2_EE_MC);
}

tuple<TH1D*, TH1D*, TH1D*, TH1D*> Purity(TFile* f0, TFile* f1, TString njet){
	TH1::SetDefaultSumw2();
	//template - fake
	TH3D* htemp0_EB = (TH3D*)f0->Get("ChIsoTemplate_Fake_g_EB");
	TH3D* htemp0_EE = (TH3D*)f0->Get("ChIsoTemplate_Fake_g_EE");
	TH3D* htemp1_EB = (TH3D*)f1->Get("ChIsoTemplate_Fake_g_EB_MC");
	TH3D* htemp1_EE = (TH3D*)f1->Get("ChIsoTemplate_Fake_g_EE_MC");
	//template - prompt
	TH3D* htemp2_EB = (TH3D*)f1->Get("ChIsoTemplate_Prompt_g_EB");
	TH3D* htemp2_EE = (TH3D*)f1->Get("ChIsoTemplate_Prompt_g_EE");
	//fitting target
	TH3D* htemp3_EB = (TH3D*)f0->Get("R2_MR_ChIso_GNoIso_EB"+njet);
	TH3D* htemp3_EE = (TH3D*)f0->Get("R2_MR_ChIso_GNoIso_EE"+njet);
  TH3D* htemp;

  std::vector<TH3D*> R2_EB;
  std::vector<TH3D*> R2_EE;
  std::vector<TH3D*> MR_EB;
  std::vector<TH3D*> MR_EE;
	const char *name[8] = {"chisoFake_EB", "chisoFake_EE", "chisoFakeMC_EB", "chisoFakeMC_EE", "chisoPrompt_EB", "chisoPrompt_EE", "chiso_EB", "chiso_EE"};
  htemp = (TH3D*)htemp0_EB->Clone();
  htemp->RebinX(5,"htemp");
  R2_EB.push_back((TH3D*)htemp->Clone());
  htemp = (TH3D*)htemp1_EB->Clone();
  htemp->RebinX(5,"htemp");
  R2_EB.push_back((TH3D*)htemp->Clone());
  htemp = (TH3D*)htemp2_EB->Clone();
  htemp->RebinX(5,"htemp");
  R2_EB.push_back((TH3D*)htemp->Clone());
  htemp = (TH3D*)htemp3_EB->Clone();
  htemp->RebinX(5,"htemp");
  R2_EB.push_back((TH3D*)htemp->Clone());
  htemp = (TH3D*)htemp0_EE->Clone();
  htemp->RebinX(5,"htemp");
  R2_EE.push_back((TH3D*)htemp->Clone());
  htemp = (TH3D*)htemp1_EE->Clone();
  htemp->RebinX(5,"htemp");
  R2_EE.push_back((TH3D*)htemp->Clone());
  htemp = (TH3D*)htemp2_EE->Clone();
  htemp->RebinX(5,"htemp");
  R2_EE.push_back((TH3D*)htemp->Clone());
  htemp = (TH3D*)htemp3_EE->Clone();
  htemp->RebinX(5,"htemp");
  R2_EE.push_back((TH3D*)htemp->Clone());
  htemp = (TH3D*)htemp0_EB->Clone();
  htemp->RebinX(5,"htemp");
  MR_EB.push_back((TH3D*)htemp->Clone());
  htemp = (TH3D*)htemp1_EB->Clone();
  htemp->RebinY(5,"htemp");
  MR_EB.push_back((TH3D*)htemp->Clone());
  htemp = (TH3D*)htemp2_EB->Clone();
  htemp->RebinY(5,"htemp");
  MR_EB.push_back((TH3D*)htemp->Clone());
  htemp = (TH3D*)htemp3_EB->Clone();
  htemp->RebinY(5,"htemp");
  MR_EB.push_back((TH3D*)htemp->Clone());
  htemp = (TH3D*)htemp0_EE->Clone();
  htemp->RebinY(5,"htemp");
  MR_EE.push_back((TH3D*)htemp->Clone());
  htemp = (TH3D*)htemp1_EE->Clone();
  htemp->RebinY(5,"htemp");
  MR_EE.push_back((TH3D*)htemp->Clone());
  htemp = (TH3D*)htemp2_EE->Clone();
  htemp->RebinY(5,"htemp");
  MR_EE.push_back((TH3D*)htemp->Clone());
  htemp = (TH3D*)htemp3_EE->Clone();
  htemp->RebinY(5,"htemp");
  MR_EE.push_back((TH3D*)htemp->Clone());

	std::vector<TH1D*> chisoFake_R2_EB;
	std::vector<TH1D*> chisoFakeMC_R2_EB;
	std::vector<TH1D*> chisoPrompt_R2_EB;
	std::vector<TH1D*> chiso_R2_EB;
	std::vector<TH1D*> chisoFake_R2_EE;
	std::vector<TH1D*> chisoFakeMC_R2_EE;
	std::vector<TH1D*> chisoPrompt_R2_EE;
	std::vector<TH1D*> chiso_R2_EE;
	std::vector<TH1D*> chisoFake_MR_EB;
	std::vector<TH1D*> chisoFakeMC_MR_EB;
	std::vector<TH1D*> chisoPrompt_MR_EB;
	std::vector<TH1D*> chiso_MR_EB;
	std::vector<TH1D*> chisoFake_MR_EE;
	std::vector<TH1D*> chisoFakeMC_MR_EE;
	std::vector<TH1D*> chisoPrompt_MR_EE;
	std::vector<TH1D*> chiso_MR_EE;

	//Make 1D hist
	for(int i=1;i<=5;i++){
      chisoFake_R2_EB.push_back(R2_EB[0]->ProjectionZ(name[0],1,1,i,i,"e"));
      chisoFakeMC_R2_EB.push_back(R2_EB[1]->ProjectionZ(name[1],1,1,i,i,"e"));
      chisoPrompt_R2_EB.push_back(R2_EB[2]->ProjectionZ(name[2],1,1,i,i,"e"));
      chiso_R2_EB.push_back(R2_EB[3]->ProjectionZ(name[3],1,1,i,i,"e"));
      chisoFake_R2_EE.push_back(R2_EE[0]->ProjectionZ(name[0],1,1,i,i,"e"));
      chisoFakeMC_R2_EE.push_back(R2_EE[1]->ProjectionZ(name[0],1,1,i,i,"e"));
      chisoPrompt_R2_EE.push_back(R2_EE[2]->ProjectionZ(name[0],1,1,i,i,"e"));
      chiso_R2_EE.push_back(R2_EE[3]->ProjectionZ(name[0],1,1,i,i,"e"));
      chisoFake_MR_EB.push_back(MR_EB[0]->ProjectionZ(name[0],i,i,1,1,"e"));
      chisoFakeMC_MR_EB.push_back(MR_EB[1]->ProjectionZ(name[0],i,i,1,1,"e"));
      chisoPrompt_MR_EB.push_back(MR_EB[2]->ProjectionZ(name[0],i,i,1,1,"e"));
      chiso_MR_EB.push_back(MR_EB[3]->ProjectionZ(name[0],i,i,1,1,"e"));
      chisoFake_MR_EE.push_back(MR_EE[0]->ProjectionZ(name[0],i,i,1,1,"e"));
      chisoFakeMC_MR_EE.push_back(MR_EE[1]->ProjectionZ(name[0],i,i,1,1,"e"));
      chisoPrompt_MR_EE.push_back(MR_EE[2]->ProjectionZ(name[0],i,i,1,1,"e"));
      chiso_MR_EE.push_back(MR_EE[3]->ProjectionZ(name[0],i,i,1,1,"e"));
	}
	std::vector<double> purity_R2_EB;
	std::vector<double> purity_R2_EE;
	std::vector<double> purity_R2_EB_MC;
	std::vector<double> purity_R2_EE_MC;
	std::vector<double> error_R2_EB;
	std::vector<double> error_R2_EE;
	std::vector<double> error_R2_EB_MC;
	std::vector<double> error_R2_EE_MC;
	std::vector<double> purity_MR_EB;
	std::vector<double> purity_MR_EE;
	std::vector<double> purity_MR_EB_MC;
	std::vector<double> purity_MR_EE_MC;
	std::vector<double> error_MR_EB;
	std::vector<double> error_MR_EE;
	std::vector<double> error_MR_EB_MC;
	std::vector<double> error_MR_EE_MC;
	pair<double, double> result;

  cout << "Chiso Fit" << endl;
	for(int i=0;i<chisoFake_R2_EB.size();++i){
		result = PurityFit(chisoFake_R2_EB[i], chisoPrompt_R2_EB[i], chiso_R2_EB[i]);
		purity_R2_EB.push_back(result.first); error_R2_EB.push_back(result.second);
		result = PurityFit(chisoFakeMC_R2_EB[i], chisoPrompt_R2_EB[i], chiso_R2_EB[i]);
		purity_R2_EB_MC.push_back(result.first); error_R2_EB_MC.push_back(result.second);
		result = PurityFit(chisoFake_R2_EE[i], chisoPrompt_R2_EE[i], chiso_R2_EE[i]);
		purity_R2_EE.push_back(result.first); error_R2_EE.push_back(result.second);
		result = PurityFit(chisoFakeMC_R2_EE[i], chisoPrompt_R2_EE[i], chiso_R2_EE[i]);
		purity_R2_EE_MC.push_back(result.first); error_R2_EE_MC.push_back(result.second);
		result = PurityFit(chisoFake_MR_EB[i], chisoPrompt_MR_EB[i], chiso_MR_EB[i]);
		purity_MR_EB.push_back(result.first); error_MR_EB.push_back(result.second);
		result = PurityFit(chisoFakeMC_MR_EB[i], chisoPrompt_MR_EB[i], chiso_MR_EB[i]);
		purity_MR_EB_MC.push_back(result.first); error_MR_EB_MC.push_back(result.second);
		result = PurityFit(chisoFake_MR_EE[i], chisoPrompt_MR_EE[i], chiso_MR_EE[i]);
		purity_MR_EE.push_back(result.first); error_MR_EE.push_back(result.second);
		result = PurityFit(chisoFakeMC_MR_EE[i], chisoPrompt_MR_EE[i], chiso_MR_EE[i]);
		purity_MR_EE_MC.push_back(result.first); error_MR_EE_MC.push_back(result.second);
	}
  cout << "end Fit" << endl;

	TH1D* hist_R2_EB = new TH1D("hist_R2_EB", "", 5, 0, 5);
	TH1D* hist_R2_EE = new TH1D("hist_R2_EE", "", 5, 0, 5);
	TH1D* hist_R2_EB_MC = new TH1D("hist_R2_EB_MC", "", 5, 0, 5);
	TH1D* hist_R2_EE_MC = new TH1D("hist_R2_EE_MC", "", 5, 0, 5);
	TH1D* hist_MR_EB = new TH1D("hist_MR_EB", "", 5, 0, 5);
	TH1D* hist_MR_EE = new TH1D("hist_MR_EE", "", 5, 0, 5);
	TH1D* hist_MR_EB_MC = new TH1D("hist_MR_EB_MC", "", 5, 0, 5);
	TH1D* hist_MR_EE_MC = new TH1D("hist_MR_EE_MC", "", 5, 0, 5);

	for(int i=1;i<=5;i++){
		hist_R2_EB->SetBinContent(i,purity_R2_EB[i-1]);
		hist_R2_EB->SetBinError(i,error_R2_EB[i-1]);
		hist_R2_EB_MC->SetBinContent(i,purity_R2_EB_MC[i-1]);
		hist_R2_EB_MC->SetBinError(i,error_R2_EB_MC[i-1]);
		hist_R2_EE->SetBinContent(i,purity_R2_EE[i-1]);
		hist_R2_EE->SetBinError(i,error_R2_EE[i-1]);
		hist_R2_EE_MC->SetBinContent(i,purity_R2_EE_MC[i-1]);
		hist_R2_EE_MC->SetBinError(i,error_R2_EE_MC[i-1]);
		hist_MR_EB->SetBinContent(i,purity_MR_EB[i-1]);
		hist_MR_EB->SetBinError(i,error_MR_EB[i-1]);
		hist_MR_EB_MC->SetBinContent(i,purity_MR_EB_MC[i-1]);
		hist_MR_EB_MC->SetBinError(i,error_MR_EB_MC[i-1]);
		hist_MR_EE->SetBinContent(i,purity_MR_EE[i-1]);
		hist_MR_EE->SetBinError(i,error_MR_EE[i-1]);
		hist_MR_EE_MC->SetBinContent(i,purity_MR_EE_MC[i-1]);
		hist_MR_EE_MC->SetBinError(i,error_MR_EE_MC[i-1]);
	}

	return make_tuple(hist_R2_EB, hist_R2_EE, hist_MR_EB, hist_MR_EE);
}

pair<TH1D*, TH1D*> Fraction(TFile* f1){
	TH1::SetDefaultSumw2();
	TH3D* htemp_EB = (TH3D*)f1->Get("R2_MR_IsDirect_G_EB");
	TH3D* htemp_EE = (TH3D*)f1->Get("R2_MR_IsDirect_G_EE");
	TH1D* hist_EB = new TH1D("hist_EB", "", 25,0,25);
	TH1D* hist_EE = new TH1D("hist_EE", "", 25,0,25);
	double frag, prom;
	double frag_e, prom_e, bin=1;
	const char *name[2] = {"htemp_EB", "htemp_EE"};
	//Make 1D hist
	for(int i=1;i<=htemp_EB->GetNbinsX();i++){
		for(int j=1;j<=htemp_EB->GetNbinsY();j++){
			frag   = htemp_EB->GetBinContent(i,j,1);
			frag_e = htemp_EB->GetBinError(i,j,1);
			prom   = htemp_EB->GetBinContent(i,j,2);
			prom_e = htemp_EB->GetBinError(i,j,2);
      //hist_EB->SetBinContent(bin,prom/(prom+frag));
			//hist_EB->SetBinError(bin,frag*sqrt(TMath::Power(frag_e,2) + TMath::Power(prom_e,2)));
			frag = htemp_EE->GetBinContent(i,j,1);
			frag_e = htemp_EE->GetBinError(i,j,1);
			prom = htemp_EE->GetBinContent(i,j,2);
			prom_e = htemp_EE->GetBinError(i,j,2);
			//hist_EE->SetBinContent(bin,prom/(prom+frag));
			//hist_EE->SetBinError(bin,frag*sqrt(TMath::Power(frag_e,2) + TMath::Power(prom_e,2)));
      hist_EB->SetBinContent(bin,1);
			hist_EB->SetBinError(bin,0.01);
      hist_EE->SetBinContent(bin,1);
			hist_EE->SetBinError(bin,0.01);
      bin++;
      //cout << hist_EB->GetBinContent(bin) << ", " << hist_EE->GetBinContent(bin) << endl;
		}
	}
	return make_pair(hist_EB, hist_EE);
}

TH1D* CalcBackground(TH1D* data, TH1D* Bkg1, TH1D* Bkg2, TH1D* Bkg3, TH1D* Multiboson, TH1D* DY, TH1D* GJet, TH1D* Target_S, TH1D* Target_Target){
	TH1::SetDefaultSumw2();
  TH1D* htemp = (TH1D*)data->Clone();
  htemp->Add(Bkg1, -1);
  htemp->Add(Bkg2, -1);
  htemp->Add(Bkg3, -1);
  htemp->Add(Multiboson, -1);
  htemp->Add(DY, -1);
  htemp->Add(GJet, -1);
  htemp->Divide(Target_Target);

  TH1D* result = (TH1D*)Target_S->Clone();
  result->Multiply(htemp);
/*
  double value, Error;

  for(int i=1;i<=Target_S->GetNbinsX();i++){
    if(htemp->GetBinContent(i) == 0) {
      value = htemp->GetBinContent(i-1);
      Error = htemp->GetBinError(i-1);
      if(value == 0 && i < 20){
        value = htemp->GetBinContent(i-2);
        Error = htemp->GetBinError(i-2);
      }
      if(i==20){
        value = htemp->GetBinContent(i+1);
        Error = htemp->GetBinError(i+1);
      }
    }
    result->SetBinContent(i,Target_S->GetBinContent(i)*value);
    result->SetBinError(i,sqrt(TMath::Power(Target_S->GetBinError(i),2)+TMath::Power(Error,2)));
  }
*/
  return result;
}

TH1D* CalcZinv(TH1D* data_EB, TH1D* data_EE, TH1D* purity_EB, TH1D* purity_EE, TH1D* frac_EB, TH1D* frac_EE, TH1D* Zinv_data, TH1D* Zinv_MC, TH1D* GJet_data, TH1D* GJet_MC, TH1D* Zinv_S, TH1D* GJet_G){
	TH1::SetDefaultSumw2();
  TH1D* htemp_EB = (TH1D*)data_EB->Clone();
  for(int i=1;i<=25;i++) cout << htemp_EB->GetBinContent(i) << endl;
  htemp_EB->Multiply(purity_EB);
  for(int i=1;i<=25;i++) cout << htemp_EB->GetBinContent(i) << endl;
  htemp_EB->Multiply(frac_EB);
  for(int i=1;i<=25;i++) cout << htemp_EB->GetBinContent(i) << endl;
  htemp_EB->Scale((Zinv_data->Integral()/Zinv_MC->Integral())/(GJet_data->Integral()/GJet_MC->Integral()));
  for(int i=1;i<=25;i++) cout << htemp_EB->GetBinContent(i) << endl;
  TH1D* htemp_EE = (TH1D*)data_EE->Clone();
  htemp_EE->Multiply(purity_EE);
  htemp_EE->Multiply(frac_EE);
  htemp_EE->Scale((Zinv_data->Integral()/Zinv_MC->Integral())/(GJet_data->Integral()/GJet_MC->Integral()));

  TH1D* htemp = (TH1D*)Zinv_S->Clone();
  htemp->Divide(GJet_G);
  TH1D* result_EB = (TH1D*)Zinv_S->Clone();
  TH1D* result_EE = (TH1D*)Zinv_S->Clone();
  TH1D* result = (TH1D*)Zinv_S->Clone();
  double value1=0, Error1=0;
  double value2=0, Error2=0;

  for(int i=1;i<=Zinv_S->GetNbinsX();i++){
    if(i==19) {
      value1 = htemp_EB->GetBinContent(19);
      value1 += htemp_EB->GetBinContent(20);
      Error1 = sqrt(TMath::Power(htemp_EB->GetBinError(19),2) + TMath::Power(htemp_EB->GetBinError(20),2));
      value2 = htemp_EE->GetBinContent(19);
      value2 += htemp_EE->GetBinContent(20);
      Error2 = sqrt(TMath::Power(htemp_EE->GetBinError(19),2) + TMath::Power(htemp_EE->GetBinError(20),2));
    }
    else if(i==20) {
      value1 = htemp_EB->GetBinContent(21);
      Error1 = htemp_EB->GetBinError(21);
      value2 = htemp_EE->GetBinContent(21);
      Error2 = htemp_EE->GetBinError(21);
    }
    else if(i==21) {
      value1 = htemp_EB->GetBinContent(22);
      Error1 = htemp_EB->GetBinError(22);
      value2 = htemp_EE->GetBinContent(22);
      Error2 = htemp_EE->GetBinError(22);
    }
    else if(i==22) {
       value1 = htemp_EB->GetBinContent(23);
      value1 += htemp_EB->GetBinContent(24);
      value1 += htemp_EB->GetBinContent(25);
      Error1 = sqrt(TMath::Power(htemp_EB->GetBinError(23),2) + TMath::Power(htemp_EB->GetBinError(24),2) + TMath::Power(htemp_EB->GetBinError(25),2));
       value2 = htemp_EE->GetBinContent(23);
      value2 += htemp_EE->GetBinContent(24);
      value2 += htemp_EE->GetBinContent(25);
      Error2 = sqrt(TMath::Power(htemp_EE->GetBinError(23),2) + TMath::Power(htemp_EE->GetBinError(24),2) + TMath::Power(htemp_EE->GetBinError(25),2));
    }
    else {
      value1 = htemp_EB->GetBinContent(i);
      Error1 = htemp_EB->GetBinError(i);
      value2 = htemp_EE->GetBinContent(i);
      Error2 = htemp_EE->GetBinError(i);
    }
    result_EB->SetBinContent(i,htemp->GetBinContent(i)*value1);
    result_EB->SetBinError(i,sqrt(TMath::Power(htemp->GetBinError(i),2)+TMath::Power(Error1,2)));
    cout << "Z EB, " << htemp->GetBinContent(i) << ", " << value1 << endl;
    cout << "Z EE, " << htemp->GetBinContent(i) << ", " << value2 << endl;
/*
    if(data_EE->GetBinContent(i) == 0) {
      value2 = data_EE->GetBinContent(i-1);
      Error2 = data_EE->GetBinError(i-1);
      if(value2 == 0 && i < 20){
        value2 = data_EE->GetBinContent(i-2);
        Error2 = data_EE->GetBinError(i-2);
      }
      if(i==20 || i==16){
        value2 = data_EE->GetBinContent(i+1);
        Error2 = data_EE->GetBinError(i+1);
      }
    }*/
    result_EE->SetBinContent(i,Zinv_S->GetBinContent(i)*value2);
    result_EE->SetBinError(i,sqrt(TMath::Power(Zinv_S->GetBinError(i),2)+TMath::Power(Error2,2)));
  }
  //result->Add(result_EB,result_EE);
  return result_EB;
}

void Draw(TH1D* data, TH1D* Multijet, TH1D* TTST, TH1D* WJet, TH1D* ZJet, TH1D* Other, TString njet, TString channel, TString region, bool Zinv){
	TH1::SetDefaultSumw2();

	TH1D* total = (TH1D*)TTST->Clone();
	total->Add(Multijet);
	total->Add(WJet);
	total->Add(ZJet);
	total->Add(Other);

	auto hs1 = new THStack("hs1","");
	Other->SetLineColor(853);
	Other->SetLineWidth(4);
	hs1->Add(Other);
	WJet->SetLineColor(418);
	WJet->SetLineWidth(4);
	hs1->Add(WJet);
	ZJet->SetLineColor(429);
	ZJet->SetLineWidth(4);
	hs1->Add(ZJet);
	Multijet->SetLineColor(619);
	Multijet->SetLineWidth(4);
	hs1->Add(Multijet);
	TTST->SetLineColor(633);
	TTST->SetLineWidth(4);
	hs1->Add(TTST);

	TCanvas* comSp = new TCanvas("comSp", "", 700, 700);
	comSp->Divide(1,2);
	TPad* canvas_up = (TPad*)comSp->GetListOfPrimitives()->FindObject("comSp_1");
	TPad* canvas_dw = (TPad*)comSp->GetListOfPrimitives()->FindObject("comSp_2");

	double up_height     = 0.7; // please tune so that the upper figures size will meet your requirement
	double dw_correction = 1.210;//1.390 // please tune so that the smaller canvas size will work in your environment
	double dw_height    = (1. - up_height) * dw_correction;

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
	data->SetMarkerColor(kBlack);
	data->SetLineColor(kBlack);
	data->SetMarkerStyle(20);
	data->SetMarkerSize(1.3);
	//total->SetFillColor(kGray);
	total->SetLineColor(kGray);
	total->SetMarkerColor(kGray);
	total->SetFillStyle(3004);
	total->SetMaximum(1e5);
	total->SetMinimum(1e-1);
	total->GetXaxis()->SetLabelSize(0.1);
	total->GetXaxis()->LabelsOption("v");
	//total->GetYaxis()->SetTitle("Events/bin");
	total->Draw("E2");
	hs1->Draw("HISTsame");
	data->Draw("EPsame");

  float xmin = data->GetXaxis()->GetBinLowEdge(data->GetXaxis()->GetFirst());
  float xmax = data->GetXaxis()->GetBinUpEdge(data->GetXaxis()->GetLast());

	string text = "CMS #scale[0.7]{#font[52]{Work in progress 2017}}";
	auto cms_lat = new TLatex(xmin, 2.0e5, text.c_str());
	cms_lat->SetTextSize(0.07);
	cms_lat->SetLineWidth(2);
	cms_lat->Draw();
	text = "#scale[0.7]{W ana,  41.5 fb^{-1} (13 TeV)}";
	if(channel=="_TOP") text = "#scale[0.7]{Top ana, 41.5 fb^{-1} (13 TeV)}";
	auto era_lat = new TLatex(xmax,2.2e5, text.c_str());
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

  auto leg = new TLegend(0.45,0.71,0.9,0.9);;
	if(region == "s"){
		leg = new TLegend(0.45,0.62,0.9,0.9);
		if(TString(njet).Contains("nj45"))     leg->SetHeader("Signal-like validation region, Wn45 final state");
		else if(TString(njet).Contains("nj6")) leg->SetHeader("Signal-like validation region, Wn6 final state");
		else leg->SetHeader("Signal-like validation region, Top final state");
	}
	else if(region == "q"){
		leg = new TLegend(0.45,0.62,0.9,0.9);
		if(TString(njet).Contains("nj45"))     leg->SetHeader("Multijet validation region, Wn45 final state");
		else if(TString(njet).Contains("nj6")) leg->SetHeader("Multijet validation region, Wn6 final state");
		else leg->SetHeader("Multijet validation region, Top final state");
	}
	else{
		leg = new TLegend(0.45,0.62,0.9,0.9);
		if(TString(njet).Contains("nj45"))     leg->SetHeader("Wn45 final state");
		else if(TString(njet).Contains("nj6")) leg->SetHeader("Wn6 final state");
		else leg->SetHeader("Top final state");
	}

	char legentry[10];

	leg->SetNColumns(2);
	leg->SetTextSize(0.03);
	sprintf(legentry, "%.2f",Multijet->Integral());
	leg->AddEntry(Multijet,  "Multijet estimate", "l");
	leg->AddEntry((TObject*)0,legentry,"");
	sprintf(legentry, "%.2f",TTST->Integral());
	leg->AddEntry(TTST,  "t#bar{t} + single top estimate", "l");
	leg->AddEntry((TObject*)0,legentry,"");
	sprintf(legentry, "%.2f",WJet->Integral());
	leg->AddEntry(WJet,  "W(#rightarrowl#nu) estimate", "l");
	leg->AddEntry((TObject*)0,legentry,"");
	sprintf(legentry, "%.2f",ZJet->Integral());
	leg->AddEntry(ZJet,  "Z(#rightarrow#nu#nu) estimate", "l");
	leg->AddEntry((TObject*)0,legentry,"");
	sprintf(legentry, "%.2f",Other->Integral());
	leg->AddEntry(Other,  "Other", "l");
	leg->AddEntry((TObject*)0,legentry,"");
	sprintf(legentry, "%.2f",data->Integral());
	leg->AddEntry(data, "data", "lp");
	leg->AddEntry((TObject*)0,legentry,"");
	leg->Draw();

	TH1D* compare  = (TH1D*)data->Clone();
	TH1D* compare1 = (TH1D*)data->Clone();
	compare->Divide(total);

	canvas_dw->cd();
	gPad->SetGridy();
	compare->SetLineColor(kBlack);
	compare->SetMarkerColor(kBlack);
	compare->SetMarkerStyle(20);
	compare->SetMarkerSize(1.3);
	compare->GetYaxis()->SetRangeUser(0.0,2.0);
	compare->SetTitle(0);
	compare->SetStats(0);
	compare->GetYaxis()->SetLabelSize(0.09);
	compare->GetXaxis()->SetLabelSize(0.09);
	compare->GetXaxis()->LabelsOption("v");
	compare->GetYaxis()->SetNdivisions(502);
	compare->GetXaxis()->SetTickLength(0.07);
	compare->GetYaxis()->SetTickLength(0.03);
	compare->GetYaxis()->SetTitle("#frac{Data}{Estimate}");
	compare->GetXaxis()->SetTitleOffset(0.81);
	compare->GetXaxis()->SetTitleSize(0.12);
	compare->GetYaxis()->SetTitleOffset(0.35);
	compare->GetYaxis()->SetTitleSize(0.10);//dividend->SetTitleSize(0.12,"Y");
	for(int i=1;i<=total->GetNbinsX();i++){
		compare1->SetBinContent(i,1);
		if(total->GetBinContent(i)==0) compare1->SetBinError(i,0);
		else compare1->SetBinError(i,total->GetBinError(i)/total->GetBinContent(i));
	}
	compare1->SetFillColor(13);
	compare1->SetFillStyle(3235);
	compare1->SetMarkerStyle(1);

	compare->Draw("EP");
	compare1->Draw("E2same");

	TString Dir = "./plot/";
	if(Zinv) comSp->SaveAs(Dir+"BkgEst_"+region+channel+njet+".png");
  else     comSp->SaveAs(Dir+"BkgEst_L_"+region+channel+njet+".png");
}

void ClosureTest(TString njet, TString channel){
	cout << njet + channel << endl;
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	TString dir = "/eos/cms/store/user/chuh/RazorBoost/190823"+channel+"_PFHT1050_PFHTXXX_PFMETXXX_PFMHTXXX_IDTight_oldCR/added/";
	TFile* file[10];
	file[0] = TFile::Open(dir+"data.root");
	file[1] = TFile::Open(dir+"Multijet.root");
	file[2] = TFile::Open(dir+"TTST.root");
	file[3] = TFile::Open(dir+"WJet.root");
	file[4] = TFile::Open(dir+"ZJet.root");
	file[5] = TFile::Open(dir+"Multiboson+TTX.root");
	file[6] = TFile::Open(dir+"DYToLL.root");
	file[7] = TFile::Open(dir+"GJet.root");
	//file1 = TFile::Open(dir+"T2tt_mStop_850_mLSP_100.root");
	//file2 = TFile::Open(dir+"T5tttt_mGluino_1400_mLSP_300.root");
	//file3 = TFile::Open(dir+"T5ttcc_mGluino_1400_mLSP_300.root");
	//file4 = TFile::Open(dir+"T1tttt_mGluino_1400_mLSP_300.root");
	TString histname[20];
	histname[0] = "R2_MR_S"+njet;
	histname[1] = "R2_MR_s"+njet;
	histname[2] = "R2_MR_Q"+njet;
	histname[3] = "R2_MR_q"+njet;
	histname[4] = "R2_MR_T"+njet;
	histname[5] = "R2_MR_W"+njet;
	histname[6] = "R2_MR_Z"+njet;
	histname[7] = "R2_MR_G"+njet;
	histname[8] = "R2_MR_L"+njet;
	histname[9] = "R2_MR_looseS"+njet;
  if(channel=="_TOP") histname[9] = "R2_MR_loose3S";
	histname[10] = "R2_MR_looses"+njet;
  if(channel=="_TOP") histname[10] = "R2_MR_loose3s";
	histname[11] = "R2_MR_looseq"+njet;
  if(channel=="_TOP") histname[11] = "R2_MR_loose3q";

	std::vector<TH1D*> hMRR2_S;
	std::vector<TH1D*> hMRR2_s;
	std::vector<TH1D*> hMRR2_Q;
	std::vector<TH1D*> hMRR2_q;
	std::vector<TH1D*> hMRR2_T;
	std::vector<TH1D*> hMRR2_W;
	std::vector<TH1D*> hMRR2_Z;
	std::vector<TH1D*> hMRR2_G;
	std::vector<TH1D*> hMRR2_L;
	std::vector<TH1D*> hMRR2_looseS;
	std::vector<TH1D*> hMRR2_looses;
	std::vector<TH1D*> hMRR2_looseq;

	for(int i=0;i<8;++i){
		hMRR2_S.push_back(new TH1D((std::string("MRR2_S_"+to_string(i))).c_str(), ";R^{2};Events/Bins", 22,0,22));
		hMRR2_s.push_back(new TH1D((std::string("MRR2_s_"+to_string(i))).c_str(), ";R^{2};Events/Bins", 22,0,22));
		hMRR2_Q.push_back(new TH1D((std::string("MRR2_Q_"+to_string(i))).c_str(), ";R^{2};Events/Bins", 22,0,22));
		hMRR2_q.push_back(new TH1D((std::string("MRR2_q_"+to_string(i))).c_str(), ";R^{2};Events/Bins", 22,0,22));
		hMRR2_T.push_back(new TH1D((std::string("MRR2_T_"+to_string(i))).c_str(), ";R^{2};Events/Bins", 22,0,22));
		hMRR2_W.push_back(new TH1D((std::string("MRR2_W_"+to_string(i))).c_str(), ";R^{2};Events/Bins", 22,0,22));
		hMRR2_Z.push_back(new TH1D((std::string("MRR2_Z_"+to_string(i))).c_str(), ";R^{2};Events/Bins", 22,0,22));
		hMRR2_G.push_back(new TH1D((std::string("MRR2_G_"+to_string(i))).c_str(), ";R^{2};Events/Bins", 22,0,22));
		hMRR2_L.push_back(new TH1D((std::string("MRR2_L_"+to_string(i))).c_str(), ";R^{2};Events/Bins", 22,0,22));
		hMRR2_looseS.push_back(new TH1D((std::string("MRR2_looseS_"+to_string(i))).c_str(), ";R^{2};Events/Bins", 22,0,22));
		hMRR2_looses.push_back(new TH1D((std::string("MRR2_looses_"+to_string(i))).c_str(), ";R^{2};Events/Bins", 22,0,22));
		hMRR2_looseq.push_back(new TH1D((std::string("MRR2_looseq_"+to_string(i))).c_str(), ";R^{2};Events/Bins", 22,0,22));
	}

  cout << "Bin Merging Start" << endl;
  TFile* f1  = new TFile("./plot/results"+channel+njet+".root","recreate");
  f1->mkdir("mergedbin");
  f1->mkdir("Bigbin");
  f1->mkdir("extrapolation");
  f1->mkdir("purity");
  f1->mkdir("fraction");
  f1->mkdir("bkgest");
  f1->mkdir("zinvest");
  f1->cd("mergedbin");

	for(int i=0;i<8;++i){
		hMRR2_S[i] = BinMerging((TH2D*)file[i]->Get(histname[0]));
		hMRR2_s[i] = BinMerging((TH2D*)file[i]->Get(histname[1]));
		hMRR2_Q[i] = BinMerging((TH2D*)file[i]->Get(histname[2]));
		hMRR2_q[i] = BinMerging((TH2D*)file[i]->Get(histname[3]));
		hMRR2_T[i] = BinMerging((TH2D*)file[i]->Get(histname[4]));
		hMRR2_W[i] = BinMerging((TH2D*)file[i]->Get(histname[5]));
		hMRR2_Z[i] = BinMerging((TH2D*)file[i]->Get(histname[6]));
		hMRR2_G[i] = BinMerging((TH2D*)file[i]->Get(histname[7]));
		hMRR2_L[i] = BinMerging((TH2D*)file[i]->Get(histname[8]));
		hMRR2_looseS[i] = BinMerging((TH2D*)file[i]->Get(histname[9]));
		hMRR2_looses[i] = BinMerging((TH2D*)file[i]->Get(histname[10]));
		hMRR2_looseq[i] = BinMerging((TH2D*)file[i]->Get(histname[11]));

		hMRR2_S[i]->Write();
		hMRR2_s[i]->Write();
		hMRR2_Q[i]->Write();
		hMRR2_q[i]->Write();
		hMRR2_T[i]->Write();
		hMRR2_W[i]->Write();
		hMRR2_Z[i]->Write();
		hMRR2_G[i]->Write();
		hMRR2_L[i]->Write();
		hMRR2_looseS[i]->Write();
		hMRR2_looses[i]->Write();
		hMRR2_looseq[i]->Write();
	}
  //hMRR2_G_EB = BinMerging((TH2D*)file[0]->Get("R2_MR_G_EB"+njet));
  //hMRR2_G_EE = BinMerging((TH2D*)file[0]->Get("R2_MR_G_EE"+njet));
  TH1D* hMRR2_G_EB = new TH1D("MRR2_G_EB", ";R^{2};Events/Bins", 25,0,25);
  TH1D* hMRR2_G_EE = new TH1D("MRR2_G_EE", ";R^{2};Events/Bins", 25,0,25);
  TH2D* hMRR2_EB = (TH2D*)file[0]->Get("MR_R2_G_EB"+njet);
  TH2D* hMRR2_EE = (TH2D*)file[0]->Get("MR_R2_G_EE"+njet);
  int bin=0;
  for(int i=1;i<=5;i++){
     for(int j=1;j<=5;j++){
       bin++;
       hMRR2_G_EB->SetBinContent(bin,hMRR2_EB->GetBinContent(i,j));
       hMRR2_G_EB->SetBinError(bin,hMRR2_EB->GetBinError(i,j));
       hMRR2_G_EE->SetBinContent(bin,hMRR2_EE->GetBinContent(i,j));
       hMRR2_G_EE->SetBinError(bin,hMRR2_EE->GetBinError(i,j));
     }
  }
  hMRR2_G_EB->Write();
  hMRR2_G_EE->Write();

  cout << "Big bin for MC Start" << endl;
  f1->cd("Bigbin");

	TH1D* hMRR2_Q_M = BigBin((TH1D*)hMRR2_Q[1]);
	TH1D* hMRR2_T_T = BigBin((TH1D*)hMRR2_T[2]);
	TH1D* hMRR2_W_W = BigBin((TH1D*)hMRR2_W[3]);
	TH1D* hMRR2_G_G = BigBin((TH1D*)hMRR2_G[7]);
	TH1D* hMRR2_L_W = BigBin((TH1D*)hMRR2_L[3]);

  hMRR2_Q_M->Write();
  hMRR2_T_T->Write();
  hMRR2_W_W->Write();
  hMRR2_G_G->Write();
  hMRR2_L_W->Write();

  cout << "extrapolation for Signal region Start" << endl;
  f1->cd("extrapolation");

  TH1D* hMRR2_SO = (TH1D*)hMRR2_S[5]->Clone();
  hMRR2_SO->Add((TH1D*)hMRR2_S[6]);
  hMRR2_SO->Add((TH1D*)hMRR2_S[7]);
  TH1D* hMRR2_O_looseS = (TH1D*)hMRR2_looseS[5]->Clone();
  hMRR2_O_looseS->Add((TH1D*)hMRR2_looseS[6]);
  hMRR2_O_looseS->Add((TH1D*)hMRR2_looseS[7]);
	TH1D* hMRR2_S_M = Extrapolation((TH1D*)hMRR2_S[1], (TH1D*)hMRR2_looseS[1], channel, njet, "QCD");
	TH1D* hMRR2_S_T = Extrapolation((TH1D*)hMRR2_S[2], (TH1D*)hMRR2_looseS[2], channel, njet, "TT");
	TH1D* hMRR2_S_W = Extrapolation((TH1D*)hMRR2_S[3], (TH1D*)hMRR2_looseS[3], channel, njet, "WJet");
	TH1D* hMRR2_S_Z = Extrapolation((TH1D*)hMRR2_S[4], (TH1D*)hMRR2_looseS[4], channel, njet, "ZJet");
  TH1D* hMRR2_S_O = Extrapolation((TH1D*)hMRR2_SO,   (TH1D*)hMRR2_O_looseS, channel, njet, "Other");

  TH1D* hMRR2_sO = (TH1D*)hMRR2_s[5]->Clone();
  hMRR2_sO->Add((TH1D*)hMRR2_s[6]);
  hMRR2_sO->Add((TH1D*)hMRR2_s[7]);
  TH1D* hMRR2_O_looses = (TH1D*)hMRR2_looses[5]->Clone();
  hMRR2_O_looses->Add((TH1D*)hMRR2_looses[6]);
  hMRR2_O_looses->Add((TH1D*)hMRR2_looses[7]);
	TH1D* hMRR2_s_M = Extrapolation((TH1D*)hMRR2_s[1], (TH1D*)hMRR2_looses[1], channel, njet, "QCD");
	TH1D* hMRR2_s_T = Extrapolation((TH1D*)hMRR2_s[2], (TH1D*)hMRR2_looses[2], channel, njet, "TT");
	TH1D* hMRR2_s_W = Extrapolation((TH1D*)hMRR2_s[3], (TH1D*)hMRR2_looses[3], channel, njet, "WJet");
	TH1D* hMRR2_s_Z = Extrapolation((TH1D*)hMRR2_s[4], (TH1D*)hMRR2_looses[4], channel, njet, "ZJet");
  TH1D* hMRR2_s_O = Extrapolation((TH1D*)hMRR2_sO,   (TH1D*)hMRR2_O_looses, channel, njet, "Other");

  TH1D* hMRR2_qO = (TH1D*)hMRR2_q[5]->Clone();
  hMRR2_qO->Add((TH1D*)hMRR2_q[6]);
  hMRR2_qO->Add((TH1D*)hMRR2_q[7]);
  TH1D* hMRR2_O_looseq = (TH1D*)hMRR2_looseq[5]->Clone();
  hMRR2_O_looseq->Add((TH1D*)hMRR2_looseq[6]);
  hMRR2_O_looseq->Add((TH1D*)hMRR2_looseq[7]);
	TH1D* hMRR2_q_M = Extrapolation((TH1D*)hMRR2_q[1], (TH1D*)hMRR2_looseq[1], channel, njet, "QCD");
	TH1D* hMRR2_q_T = Extrapolation((TH1D*)hMRR2_q[2], (TH1D*)hMRR2_looseq[2], channel, njet, "TT");
	TH1D* hMRR2_q_W = Extrapolation((TH1D*)hMRR2_q[3], (TH1D*)hMRR2_looseq[3], channel, njet, "WJet");
	TH1D* hMRR2_q_Z = Extrapolation((TH1D*)hMRR2_q[4], (TH1D*)hMRR2_looseq[4], channel, njet, "ZJet");
  TH1D* hMRR2_q_O = Extrapolation((TH1D*)hMRR2_qO,   (TH1D*)hMRR2_O_looseq, channel, njet, "Other");

  hMRR2_S_M->Write();
  hMRR2_S_T->Write();
  hMRR2_S_W->Write();
  hMRR2_S_Z->Write();
  hMRR2_S_O->Write();
  hMRR2_s_M->Write();
  hMRR2_s_T->Write();
  hMRR2_s_W->Write();
  hMRR2_s_Z->Write();
  hMRR2_s_O->Write();
  hMRR2_q_M->Write();
  hMRR2_q_T->Write();
  hMRR2_q_W->Write();
  hMRR2_q_Z->Write();
  hMRR2_q_O->Write();

  cout << "purity calculation" << endl;
  f1->cd("purity");

	tuple<TH1D*, TH1D*, TH1D*, TH1D*> purity = Purity2D(file[0], file[7], njet);
	TH1D* purity_R2_EB = std::get<0>(purity);
	TH1D* purity_R2_EE = std::get<1>(purity);
	TH1D* purity_MR_EB = std::get<2>(purity);
	TH1D* purity_MR_EE = std::get<3>(purity);

  TH1D* purity_EB = new TH1D("purity_EB","", 25,0,25);
  TH1D* purity_EE = new TH1D("purity_EE","", 25,0,25);
  for(int i=1;i<=25;i++){
    purity_EB->SetBinContent(i,purity_R2_EB->GetBinContent((int)(i+4)/5));
    purity_EE->SetBinContent(i,purity_R2_EE->GetBinContent((int)(i+4)/5));
    purity_EB->SetBinError(i,purity_R2_EB->GetBinError((int)(i+4)/5));
    purity_EE->SetBinError(i,purity_R2_EE->GetBinError((int)(i+4)/5));
    //purity_EB->SetBinContent(i,0.7);
    //purity_EE->SetBinContent(i,0.6);
    //purity_EB->SetBinError(i,0.1);
    //purity_EE->SetBinError(i,0.1);
  }
  //for(int i=1;i<=25;i++) cout << purity_EB->GetBinContent(i) << ", " << purity_EE->GetBinContent(i) << endl;
  purity_EB->Write();
  purity_EE->Write();

  cout << "Fraction calculation" << endl;
  f1->cd("fraction");

	pair<TH1D*, TH1D*> frac = Fraction(file[7]);
	TH1D* frac_EB = frac.first;
	TH1D* frac_EE = frac.second;

  frac_EB->Write();
  frac_EE->Write();

  cout << "Zinv bkg calculation" << endl;
  f1->cd("zinvest");

  TH1D* BkgEst_S_Z = CalcZinv(hMRR2_G_EB, hMRR2_G_EE, purity_EB, purity_EE, frac_EB, frac_EE, hMRR2_Z[0], hMRR2_Z[6], hMRR2_G[0], hMRR2_G[7], hMRR2_S[4], hMRR2_G_G);
  TH1D* BkgEst_s_Z = CalcZinv(hMRR2_G_EB, hMRR2_G_EE, purity_EB, purity_EE, frac_EB, frac_EE, hMRR2_Z[0], hMRR2_Z[6], hMRR2_G[0], hMRR2_G[7], hMRR2_s[4], hMRR2_G_G);
  TH1D* BkgEst_q_Z = CalcZinv(hMRR2_G_EB, hMRR2_G_EE, purity_EB, purity_EE, frac_EB, frac_EE, hMRR2_Z[0], hMRR2_Z[6], hMRR2_G[0], hMRR2_G[7], hMRR2_q[4], hMRR2_G_G);

  BkgEst_S_Z->Write();
  BkgEst_s_Z->Write();
  BkgEst_q_Z->Write();

  cout << "Bkg calculation" << endl;
  f1->cd("bkgest");

  TH1D* BkgEst_S_M = CalcBackground(hMRR2_Q[0], hMRR2_Q[2], hMRR2_Q[3], hMRR2_Q[4], hMRR2_Q[5], hMRR2_Q[6], hMRR2_Q[7], hMRR2_S_M , hMRR2_Q_M);
  TH1D* BkgEst_s_M = CalcBackground(hMRR2_Q[0], hMRR2_Q[2], hMRR2_Q[3], hMRR2_Q[4], hMRR2_Q[5], hMRR2_Q[6], hMRR2_Q[7], hMRR2_s_M, hMRR2_Q_M);
  TH1D* BkgEst_q_M = CalcBackground(hMRR2_Q[0], hMRR2_Q[2], hMRR2_Q[3], hMRR2_Q[4], hMRR2_Q[5], hMRR2_Q[6], hMRR2_Q[7], hMRR2_q_M, hMRR2_Q_M);
  TH1D* BkgEst_S_T = CalcBackground(hMRR2_T[0], hMRR2_T[1], hMRR2_T[3], hMRR2_T[4], hMRR2_T[5], hMRR2_T[6], hMRR2_T[7], hMRR2_S_T , hMRR2_T_T);
  TH1D* BkgEst_s_T = CalcBackground(hMRR2_T[0], hMRR2_T[1], hMRR2_T[3], hMRR2_T[4], hMRR2_T[5], hMRR2_T[6], hMRR2_T[7], hMRR2_s_T, hMRR2_T_T);
  TH1D* BkgEst_q_T = CalcBackground(hMRR2_T[0], hMRR2_T[1], hMRR2_T[3], hMRR2_T[4], hMRR2_T[5], hMRR2_T[6], hMRR2_T[7], hMRR2_q_T, hMRR2_T_T);
  TH1D* BkgEst_S_W = CalcBackground(hMRR2_W[0], hMRR2_W[1], hMRR2_W[2], hMRR2_W[4], hMRR2_W[5], hMRR2_W[6], hMRR2_W[7], hMRR2_S_W , hMRR2_W_W);
  TH1D* BkgEst_s_W = CalcBackground(hMRR2_W[0], hMRR2_W[1], hMRR2_W[2], hMRR2_W[4], hMRR2_W[5], hMRR2_W[6], hMRR2_W[7], hMRR2_s_W, hMRR2_W_W);
  TH1D* BkgEst_q_W = CalcBackground(hMRR2_W[0], hMRR2_W[1], hMRR2_W[2], hMRR2_W[4], hMRR2_W[5], hMRR2_W[6], hMRR2_W[7], hMRR2_q_W, hMRR2_W_W);
  TH1D* BkgEst_S_L = CalcBackground(hMRR2_L[0], hMRR2_L[1], hMRR2_L[2], hMRR2_L[4], hMRR2_L[5], hMRR2_L[6], hMRR2_L[7], hMRR2_S_W , hMRR2_L_W);
  TH1D* BkgEst_s_L = CalcBackground(hMRR2_L[0], hMRR2_L[1], hMRR2_L[2], hMRR2_L[4], hMRR2_L[5], hMRR2_L[6], hMRR2_L[7], hMRR2_s_W, hMRR2_L_W);
  TH1D* BkgEst_q_L = CalcBackground(hMRR2_L[0], hMRR2_L[1], hMRR2_L[2], hMRR2_L[4], hMRR2_L[5], hMRR2_L[6], hMRR2_L[7], hMRR2_q_W, hMRR2_L_W);

  //TH1D* BkgEst_S_M = CalcBackground(hMRR2_Q[0], hMRR2_Q[2], hMRR2_Q[3], hMRR2_Q[4], hMRR2_Q[5], hMRR2_Q[6], hMRR2_Q[7], hMRR2_S[1], hMRR2_Q_M);
  //TH1D* BkgEst_S_T = CalcBackground(hMRR2_T[0], hMRR2_T[1], hMRR2_T[3], hMRR2_T[4], hMRR2_T[5], hMRR2_T[6], hMRR2_T[7], hMRR2_S[2], hMRR2_T_T);
  //TH1D* BkgEst_S_W = CalcBackground(hMRR2_W[0], hMRR2_W[1], hMRR2_W[2], hMRR2_W[4], hMRR2_W[5], hMRR2_W[6], hMRR2_W[7], hMRR2_S[3], hMRR2_W_W);
  //TH1D* BkgEst_S_L = CalcBackground(hMRR2_L[0], hMRR2_L[1], hMRR2_L[2], hMRR2_L[4], hMRR2_L[5], hMRR2_L[6], hMRR2_L[7], hMRR2_S[3], hMRR2_L_W);

  BkgEst_S_M->Write();
  BkgEst_s_M->Write();
  BkgEst_q_M->Write();
  BkgEst_S_T->Write();
  BkgEst_s_T->Write();
  BkgEst_q_T->Write();
  BkgEst_S_W->Write();
  BkgEst_s_W->Write();
  BkgEst_q_W->Write();
  BkgEst_S_L->Write();
  BkgEst_s_L->Write();
  BkgEst_q_L->Write();

  cout << "Draw" << endl;

	Draw(hMRR2_S[0], BkgEst_S_M, BkgEst_S_T, BkgEst_S_W, BkgEst_S_Z, hMRR2_S_O, njet, channel, "S", true);
	Draw(hMRR2_s[0], BkgEst_s_M, BkgEst_s_T, BkgEst_s_W, BkgEst_s_Z, hMRR2_s_O, njet, channel, "s", true);
	Draw(hMRR2_q[0], BkgEst_q_M, BkgEst_q_T, BkgEst_q_W, BkgEst_q_Z, hMRR2_q_O, njet, channel, "q", true);
	Draw(hMRR2_S[0], BkgEst_S_M, BkgEst_S_T, BkgEst_S_W, BkgEst_S_L, hMRR2_S_O, njet, channel, "S", false);
	Draw(hMRR2_s[0], BkgEst_s_M, BkgEst_s_T, BkgEst_s_W, BkgEst_s_L, hMRR2_s_O, njet, channel, "s", false);
	Draw(hMRR2_q[0], BkgEst_q_M, BkgEst_q_T, BkgEst_q_W, BkgEst_q_L, hMRR2_q_O, njet, channel, "q", false);
  f1->Close();
}

void BkgEst(){
	TH1::SetDefaultSumw2();
	ClosureTest("_nj45","");
	ClosureTest("_nj6","");
	ClosureTest("","_TOP");
}
