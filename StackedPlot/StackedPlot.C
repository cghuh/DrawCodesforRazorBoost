#include <algorithm>
#include <iostream>
#include <TEnv.h>
#include <TFile.h>
#include <string.h>
#include <TH1F.h>
#include <TH2D.h>
#include <map>
#include <TKey.h>
#include <Riostream.h>
#include <THStack.h>
#include <vector>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TPavesText.h>
#include <TPad.h>
#include <Rtypes.h>
#include <TStyle.h>
using namespace std;

/////////////////////////// global variable/////////////////////////////////////////////////////

int RebinVarNum;
int NumOfKeys;
int Num;/// Total Number of Root Files
char RootFilePath[200];//in RootFilepath.cfg

/////////////// define functions ///////////////////////////////////////////////////////////////

void MakeControlPlots (int ColorCode[][15], int ColorCodePlus[], TFile *File[], map<TString , TH1F *> HIST[]);

map<TString, TH1F *> MCSTACK ( char RootFileName[][200], map<TString , TH1F *> HIST[] );
//This function make ControlPlot eps file and return the tot_mc histogram.

void SaveYields( char RootFileName[][200], map<TString, TH1F *> HIST[], map<TString, TH1F *> MCtot );

std::string underscore2space(std::string text) {
	for(std::string::iterator it = text.begin(); it != text.end(); ++it) {
		if(*it == '_') {
			*it = ' ';
		}
	}   
	return text;
}


////////////// main code //////////////////////////////////////////////////////////////////////
int main(){


	/////////////////////////////////////load root file name,path and color code ////////////////  
	FILE *fp1 = fopen("RootFilePath.cfg","r");
	fscanf(fp1,"%s",RootFilePath);
	fclose(fp1);

	cout<<"******************************************************************************"<<endl;
	cout<<"RootFilePath = "<<RootFilePath<<endl;


	FILE *fp2 = fopen("Settings.cfg","r");

	char na1[1000][100];
	char na2[1000][100];
	int na3[1000];
	int ii=0;

	while(fscanf(fp2,"%s %s %d\n",na1[ii],na2[ii],&na3[ii])!=EOF) {
		ii++;
	}  ////////////////Read NumOfFile


	Num=ii;
	fclose(fp2);
	cout<<"Number of Files = "<<Num<<endl;



	FILE *fp3 = fopen("Settings.cfg","r");

	char RootFileName[Num][200];
	int ColorCode[Num][15];
	int ColorCodePlus[Num];


	for (int j=0; j<Num; j++) {
		//fscanf(fp3,"%s %s %d\n",RootFileName[j],ColorCode[j],&ColorCodePlus[j]);
		fscanf(fp3,"%s %d %d\n",RootFileName[j],&ColorCode[j],&ColorCodePlus[j]);

	}/////// read RootFileName, ColorCode


	fclose(fp3);

	char RootFileFullPath[Num][300];
	for(int i=0;i<Num;i++){
		strcpy(RootFileFullPath[i],RootFilePath);
	}
	for(int i=0;i<Num;i++){
		strcat(RootFileFullPath[i],RootFileName[i]);
		strcat(RootFileFullPath[i],".root");
	}
	for (int i=0;i<Num;i++){
		cout<<RootFileFullPath[i]<<endl;

		//cout<<ColorCode[i];
		//cout<<ColorCodePlus[i];
	}
	cout<<"******************************************************************************"<<endl;

	map<TString , TH1F *> HIST[Num];// data & each MC's histogram container [rootfilenum][keynum](array+map)
	map<TString , TH1F *> MCtot;//tot MC container
	MCtot.clear();
	
	/************ read root files ***************************/
	TFile *File[Num];

	for(int i=0;i<Num;i++){
		File[i] = new TFile(RootFileFullPath[i]);
	}   

	cout << " ////////////////////////////Successfully read Data & MC Files "<<endl;
	/********************************************************/	

	
	
	//////////// run  functions //////////////////////////////
	//														
	//
		MakeControlPlots(ColorCode,ColorCodePlus,File,HIST);

		MCtot = MCSTACK(RootFileName,HIST);
	
		//SaveYields(RootFileName,HIST,MCtot);
	//
	//	
	//////////////////////////////////////////////////////////

	//histogram container test in this main code
	//string kk;
	//cout<<MCtot["AnLepEnergy"]->Integral()<<endl;
	//cout<<HIST[0]["AnLepEnergy"]->Integral()<<endl;
	//cin>> kk;
}

void MakeControlPlots (int ColorCode[][15], int ColorCodePlus[], TFile *File[], map<TString , TH1F *> HIST[]){

	/***** Colorcode ***********************************/

	//Color_t color[Num];
	int color[Num] = {1, 12, 862, 841, 403, 803, 803, 601, 433, 418, 401,619,633};

    /*
	for(int i=0;i<Num;i++){
		if(strcmp(ColorCode[i],"kWhite")==0) color[i] =kWhite+ColorCodePlus[i];
		else if(strcmp(ColorCode[i],"kBlack")==0) color[i] =kBlack+ColorCodePlus[i];
		else if(strcmp(ColorCode[i],"kGray")==0) color[i] =kGray+ColorCodePlus[i];
		else if(strcmp(ColorCode[i],"kRed")==0) color[i] =kRed+ColorCodePlus[i];
		else if(strcmp(ColorCode[i],"kGreen")==0) color[i] =kGreen+ColorCodePlus[i];
		else if(strcmp(ColorCode[i],"kBlue")==0) color[i] =kBlue+ColorCodePlus[i];
		else if(strcmp(ColorCode[i],"kYellow")==0) color[i] =kYellow+ColorCodePlus[i];
		else if(strcmp(ColorCode[i],"kMagenta")==0) color[i] =kMagenta+ColorCodePlus[i];
		else if(strcmp(ColorCode[i],"kCyan")==0) color[i] =kCyan+ColorCodePlus[i];
		else if(strcmp(ColorCode[i],"kOrange")==0) color[i] =kOrange+ColorCodePlus[i];
		else if(strcmp(ColorCode[i],"kSpring")==0) color[i] =kSpring+ColorCodePlus[i];
		else if(strcmp(ColorCode[i],"kTeal")==0) color[i]=kTeal+ColorCodePlus[i];
		else if(strcmp(ColorCode[i],"kAzure")==0) color[i] =kAzure+ColorCodePlus[i];
		else if(strcmp(ColorCode[i],"kViolet")==0) color[i] =kViolet+ColorCodePlus[i];
		else if(strcmp(ColorCode[i],"kPink")==0) color[i] =kPink+ColorCodePlus[i];
		else
    color[i] = (Color_t)(ColorCode[i]+ColorCodePlus[i]);
	}
*/


	/***** make histogram **************************************************************************/

	for(int i=0;i<Num;i++){

		HIST[i].clear();
	}

	string histkey, histclonname, addname;
	TObject *obj;

	File[0]->cd();

	TDirectory *current_sourcedir = gDirectory;
	TIter nextkey( current_sourcedir->GetListOfKeys() );
	TKey *key;

	while ( (key= (TKey*)nextkey()))
	{

		obj = key->ReadObj();
   	histkey = key->GetName();

/*
    if(obj->InheritsFrom("TH2")){
      if(TString(histkey).Contains("HT"))continue;
      if(TString(histkey).Contains("btag"))continue;
      if(TString(histkey).Contains("read_speed"))continue;
      if(TString(histkey).Contains("run"))continue;
      if(TString(histkey).Contains("Gluino"))continue;
      if(TString(histkey).Contains("weight"))continue;
      if(TString(histkey).Contains("xsec"))continue;
      if(TString(histkey).Contains("trigger"))continue;
      if(TString(histkey).Contains("_nj"))continue;
			for(int i=0;i<Num;i++){
        File[i]->cd();
        gStyle->SetOptStat(0);
        TCanvas* ctemp = new TCanvas("ctemp","",1000,1000); 
        //gPad->SetLogz();
		  	TH2D *htemp = (TH2D*)gDirectory->Get(key->GetName())->Clone(histkey.c_str());
        htemp->Draw("colz");
	    	string Dir="./Histogram/2D/";
	    	string png=".png";
	FILE *fp3 = fopen("Settings.cfg","r");

	char RootFileName[Num][200];
	int ColorCode[Num][15];
	int ColorCodePlus[Num];


	for (int j=0; j<Num; j++) {
		//fscanf(fp3,"%s %s %d\n",RootFileName[j],ColorCode[j],&ColorCodePlus[j]);
		fscanf(fp3,"%s %d %d\n",RootFileName[j],&ColorCode[j],&ColorCodePlus[j]);

	}/////// read RootFileName, ColorCode

	    	string CanvasName=Dir+histkey+"_"+(string)RootFileName[i]+png;
	    	ctemp->SaveAs(CanvasName.c_str());
      }
      continue;
     }
  */  
		File[0]->cd();

    if(obj->InheritsFrom("TH2")) continue;
    if(histkey=="nvtx")continue;
    if(histkey=="nvtx_rw")continue;
    if(histkey=="pileup_weight")continue;
    if(histkey=="pileup_weight_up")continue;
    if(histkey=="pileup_weight_down")continue;
    if(histkey=="pileup_data")continue;
    if(histkey=="pileup_data_up")continue;
    if(histkey=="pileup_data_down")continue;
    if(histkey=="pileup_mc")continue;
    if(histkey=="totweight")continue;
    if(histkey=="nb")continue;
    if(histkey=="nw")continue;
    if(histkey=="ht_gen")continue;
    if(histkey=="h_njet_S_3Jet")continue;
    if(histkey=="h_njet_S_MRR2")continue;
    if(histkey=="jet1_AK4Puppi_pt_syst0")continue;
    if(histkey=="jet1_AK8Puppi_pt_syst0")continue;
    if(histkey=="GluinoMass")continue;
    if(histkey=="LSPMass")continue;
    if(histkey=="StopMass")continue;
    if(histkey=="counts")continue;
    if(TString(histkey).Contains("read_speed"))continue;
    if(TString(histkey).Contains("run"))continue;
    if(TString(histkey).Contains("trigger"))continue;
    if(TString(histkey).Contains("pass"))continue;
    if(TString(histkey).Contains("syst"))continue;
    if(TString(histkey).Contains("pre"))continue;
    if(TString(histkey).Contains("j1"))continue;
    if(TString(histkey).Contains("HT"))continue;
    if(TString(histkey).Contains("Gen"))continue;
    if(TString(histkey).Contains("fakerate"))continue;
    if(TString(histkey).Contains("Up"))continue;
    if(TString(histkey).Contains("Down"))continue;
    if(TString(histkey).Contains("MRR2"))continue;
    if(TString(histkey).Contains("ChargedIso"))continue;

		if( obj->InheritsFrom("TH1"))
		{


			File[0]->cd();
			TH1F *Hist = (TH1F*)gDirectory->Get(key->GetName())->Clone(histkey.c_str());
			HIST[0][histkey]=Hist;

			HIST[0][histkey]->SetLineColor(color[0]);
			HIST[0][histkey]->SetMarkerStyle(20);
			HIST[0][histkey]->SetMarkerSize(1.3);
			HIST[0][histkey]->SetMarkerColor(color[0]);

			File[1]->cd();
			Hist = (TH1F*)gDirectory->Get(key->GetName())->Clone(histkey.c_str());
			HIST[1][histkey]=Hist;
			HIST[1][histkey]->SetLineColor(color[1]);
			HIST[1][histkey]->SetLineStyle(2);
			HIST[1][histkey]->SetLineWidth(4);

			File[2]->cd();
			Hist = (TH1F*)gDirectory->Get(key->GetName())->Clone(histkey.c_str());
			HIST[2][histkey]=Hist;
			HIST[2][histkey]->SetLineColor(color[2]);
			HIST[2][histkey]->SetLineStyle(2);
			HIST[2][histkey]->SetLineWidth(4);

			File[3]->cd();
			Hist = (TH1F*)gDirectory->Get(key->GetName())->Clone(histkey.c_str());
			HIST[3][histkey]=Hist;
			HIST[3][histkey]->SetLineColor(color[3]);
			HIST[3][histkey]->SetLineStyle(2);
			HIST[3][histkey]->SetLineWidth(4);

			File[4]->cd();
			Hist = (TH1F*)gDirectory->Get(key->GetName())->Clone(histkey.c_str());
			HIST[4][histkey]=Hist;
			HIST[4][histkey]->SetLineColor(color[4]);
			HIST[4][histkey]->SetLineStyle(2);
			HIST[4][histkey]->SetLineWidth(4);

			for(int i=5;i<Num;i++)

			{
				File[i]->cd();

				TH1F *Hist = (TH1F*)gDirectory->Get(key->GetName())->Clone(histkey.c_str());
				HIST[i][histkey]=Hist;

				HIST[i][histkey]->SetLineColor(color[i]);
				HIST[i][histkey]->SetFillColor(color[i]);
				HIST[i][histkey]->SetFillStyle(1001);

			}

		}// if
                else continue;

	}// while


	cout << endl << "************* All Histograms loaded "<<endl;
}

map<TString, TH1F *> MCSTACK (char RootFileName[][200], map<TString , TH1F *> HIST[]){

	map<TString, TH1F *>totmc;
	
	char SigName[200];
	strcpy (SigName,RootFileName[Num-1]);
	strcat (SigName,"");

	TFile *STACK_File = new TFile("STACKHIST.root","RECREATE");

	FILE *NBins = fopen("NBins.txt","w");

	string histname;
	string histtitle;
	string stackname;
	int histnum=0;


	/************* read Rebin.cfg *******************************/    

	FILE *rebin = fopen("Rebin.cfg","r");

	char ka1[150][100];
	int ka2[150];
	int iii=0;

	while(fscanf(rebin,"%s %d\n",ka1[iii],&ka2[iii])!=EOF) {
		iii++;
	}  ////////////////Read RebinVarNum 


	RebinVarNum=iii;
	fclose(rebin);



	FILE *rrebin = fopen("Rebin.cfg","r");

	char RebinVar[RebinVarNum][30];
	int RebinSize[RebinVarNum];


	for (int j=0;j<RebinVarNum;j++)
	{
		fscanf(rrebin,"%s %d\n",RebinVar[j],&RebinSize[j]);
	}

	fclose(rrebin);



	for( std::map<TString, TH1F *>::iterator histoIter = HIST[0].begin(); histoIter != HIST[0].end(); ++histoIter)

	{

		histname = (*histoIter).first;
    //cout << histname << endl;

		/*********** Lets rebin  *********************************/
		for(int j=0;j<RebinVarNum;j++){
			if (strstr(histname.c_str(),RebinVar[j])!=NULL){

				for(int i=0;i<Num;i++){

					HIST[i][histname]->Rebin(RebinSize[j]);

					//cout<<histname<<" ** "<<RebinVar[j]<<" Rebinsize = "<<RebinSize[j]<<" ** Successful rebinning"<<endl;
				}
			}
		}
    if(TString(histname).Contains("tau21")){
     if(!TString(histname).Contains("notau21")||TString(histname).Contains("tau21_notau21")){
			for(int i=0;i<Num;i++){
        HIST[i][histname]->Rebin(5);
        }
       }
     }

		histtitle = HIST[0][histname]->GetXaxis()->GetTitle();

		THStack *MCStack = new THStack(histname.c_str(),histtitle.c_str());
		for(int k=5;k<Num;k++) {MCStack->Add(HIST[k][histname]);}//stack MCs

		/*********************   save the MC stack -> root files        **********/
		//STACK_File->cd();
		//MCStack->Write();////////////save the stacked histogram
		/**************************************************************************/

		histnum++;
		

    gStyle->SetOptTitle(kFALSE);
		TCanvas *c1 = new TCanvas("c1","",1000,1000);
    gPad->SetGridx();
    gPad->SetGridy();

		TH1F *Tot_MC = (TH1F*)HIST[3][histname]->Clone();
		for(int i=5;i<Num;i++) {Tot_MC->Add(HIST[i][histname]);}

		Tot_MC->SetFillColor(13);
		Tot_MC->SetFillStyle(3235);
		
		totmc[histname]=Tot_MC;

		double MCmax = Tot_MC->GetMaximum();
		double datamax = HIST[0][histname]->GetMaximum();
		double max = (MCmax>datamax)? MCmax : datamax;

		//MCStack->SetMaximum(max*1.2);
		//MCStack->SetMinimum(max*0.00005);
    //if(TString(histname).Contains("_nw")||TString(histname).Contains("_nb")||TString(histname).Contains("_nj")){
      MCStack->SetMaximum(1e8);
      MCStack->SetMinimum(1e-2);
    //}

		/***************** Canvas Setting  ***************************************/

    c1->Divide(1,2);

		TPad *canvas_up = (TPad*)c1->GetListOfPrimitives()->FindObject("c1_1");
		TPad *canvas_dw = (TPad*)c1->GetListOfPrimitives()->FindObject("c1_2");

		// define the size
		double up_height     = 0.7; // please tune so that the upper figures size will meet your requirement
		double dw_correction = 1.290;//1.390 // please tune so that the smaller canvas size will work in your environment
		//double font_size_dw  = 10; // please tune the font size parameter for bottom figure


		double dw_height    = (1. - up_height) * dw_correction;
    if(TString(histname.c_str()).Contains("_S")) {up_height = 1;dw_height=0;}

		// set pad size
		canvas_up->SetPad(0., 1 - up_height, 1., 1.);
		canvas_dw->SetPad(0., 0., 1., dw_height);
		//canvas_up->SetBottomMargin(0.3);  
		//canvas_dw->SetGridy();

    canvas_up->SetLogy();
		canvas_up->SetFrameFillColor(0);
		canvas_up->SetFillColor(0);
		canvas_dw->SetFillColor(0);
		canvas_dw->SetFrameFillColor(0);
    if(!(TString(histname.c_str()).Contains("_S"))) canvas_up->SetLeftMargin(0.1);
		// set top margin 0 for bottom figure
		canvas_dw->SetTopMargin(0.03);//0.03
		canvas_dw->SetLeftMargin(0.1);
		canvas_dw->SetBottomMargin(0.3);
		/**************************************************************************************************/    
		canvas_up->cd();
		
		MCStack->Draw("HIST");
	   MCStack->GetXaxis()->SetTitle(histtitle.c_str());
    if(!(TString(histname.c_str()).Contains("_S"))) MCStack->GetXaxis()->SetLabelOffset(999);
    if(TString(histname.c_str()).Contains("_S")){
	  	//MCStack->GetXaxis()->SetTitleOffset(0.6);
	  	//MCStack->GetXaxis()->SetTitleSize(0.70);
    }
    //cout << MCStack->GetXaxis()->GetTitle() << endl;
		MCStack->GetYaxis()->SetTitle("Events");
		MCStack->GetYaxis()->SetTitleSize(0.07);
		MCStack->GetYaxis()->SetTitleOffset(0.6);
		
		c1->Update();

    float xmin = Tot_MC->GetXaxis()->GetBinLowEdge(Tot_MC->GetXaxis()->GetFirst());
    float xmax = Tot_MC->GetXaxis()->GetBinUpEdge(Tot_MC->GetXaxis()->GetLast());
    //float ymin = Tot_MC->GetMinimum();
    //float ymax = Tot_MC->GetMaximum();

    string text = "CMS #scale[0.7]{#font[52]{Work in progress 2016}}";
    TLatex* cms_lat = new TLatex(xmin, 6e8, text.c_str());
    cms_lat->SetLineWidth(2);
    cms_lat->Draw();
    text = "#scale[0.7]{W ana, 35.9 fb^{-1} (13 TeV)}";
    if(TString(histname.c_str()).Contains("topana"))text = "#scale[0.7]{Top ana, 35.9 fb^{-1} (13 TeV)}";
    TLatex* era_lat = new TLatex(xmax*0.95,6e8, text.c_str());
    era_lat->SetTextAlign(32);
    era_lat->SetTextSize(0.04);
    era_lat->SetTextFont(42);
    era_lat->SetLineWidth(2);
    era_lat->Draw();


    if(!(TString(histname.c_str()).Contains("_S"))) HIST[0][histname]->Draw("samePE1");
		HIST[1][histname]->Draw("sameHIST");
		HIST[2][histname]->Draw("sameHIST");
		HIST[3][histname]->Draw("sameHIST");
		HIST[4][histname]->Draw("sameHIST");


		Tot_MC->Draw("SAME E2");
		/***************** legend settings *************************************************/


		TLegend *leg = new TLegend(0.44,0.51,0.89,0.89);
		leg->SetLineColor(kBlack);
		leg->SetFillColor(10);
		leg->SetLineWidth(0);

   	//TString histkey = key->GetName();

    if(TString(histname.c_str()).Contains("_S")){
      if(TString(histname.c_str()).Contains("nj6"))leg->SetHeader("S region, 6 #leq jet");
      else if(TString(histname.c_str()).Contains("nj35"))leg->SetHeader("S region, 3 #leq jet #leq 5");
      else leg->SetHeader("S region");
    }
    else if(TString(histname.c_str()).Contains("_s")){
      if(TString(histname.c_str()).Contains("nj6"))leg->SetHeader("S' region, 6 #leq jet");
      else if(TString(histname.c_str()).Contains("nj35"))leg->SetHeader("S' region, 3 #leq jet #leq 5");
      else leg->SetHeader("S' region");
    }
    else if(TString(histname.c_str()).Contains("_q")){
      if(TString(histname.c_str()).Contains("nj6"))leg->SetHeader("Q' region, 6 #leq jet");
      else if(TString(histname.c_str()).Contains("nj35"))leg->SetHeader("Q' region, 3 #leq jet #leq 5");
      else leg->SetHeader("Q' region");
    }
    else if(TString(histname.c_str()).Contains("_Q")){
      if(TString(histname.c_str()).Contains("nodPhi")){
        if(TString(histname.c_str()).Contains("nj6"))leg->SetHeader("Q region, nodPhi, 6 #leq jet");
        else if(TString(histname.c_str()).Contains("nj35"))leg->SetHeader("Q region, nodPhi, 3 #leq jet #leq 5");
        else leg->SetHeader("Q region, nodPhi");
        }
      else if(TString(histname.c_str()).Contains("nonb")){
        if(TString(histname.c_str()).Contains("nj6"))leg->SetHeader("Q region, nonb, 6 #leq jet");
        else if(TString(histname.c_str()).Contains("nj35"))leg->SetHeader("Q region, nonb, 3 #leq jet #leq 5");
        else leg->SetHeader("Q region, nonb");
        }
      else if(TString(histname.c_str()).Contains("notau21")){
        if(TString(histname.c_str()).Contains("nj6"))leg->SetHeader("Q region, notau21, 6 #leq jet");
        else if(TString(histname.c_str()).Contains("nj35"))leg->SetHeader("Q region, notau21, 3 #leq jet #leq 5");
        else leg->SetHeader("Q region, notau21");
        }
      else {
        if(TString(histname.c_str()).Contains("nj6"))leg->SetHeader("Q region, 6 #leq jet");
        else if(TString(histname.c_str()).Contains("nj35"))leg->SetHeader("Q region, 3 #leq jet #leq 5");
        else leg->SetHeader("Q region");
        }
    }
    else if(TString(histname.c_str()).Contains("_T")){
      if(TString(histname.c_str()).Contains("noMT")){
        if(TString(histname.c_str()).Contains("nj6"))leg->SetHeader("T region, noMT, 6 #leq jet");
        else if(TString(histname.c_str()).Contains("nj35"))leg->SetHeader("T region, noMT, 3 #leq jet #leq 5");
        else leg->SetHeader("T region, noMT");
        }
      else {
        if(TString(histname.c_str()).Contains("nj6"))leg->SetHeader("T region, 6 #leq jet");
        else if(TString(histname.c_str()).Contains("nj35"))leg->SetHeader("T region, 3 #leq jet #leq 5");
        else leg->SetHeader("T region");
        }
    }
    else if(TString(histname.c_str()).Contains("_W")){
      if(TString(histname.c_str()).Contains("noMT")){
        if(TString(histname.c_str()).Contains("nj6"))leg->SetHeader("W region, noMT, 6 #leq jet");
        else if(TString(histname.c_str()).Contains("nj35"))leg->SetHeader("W region, noMT, 3 #leq jet #leq 5");
        else leg->SetHeader("W region, noMT");
        }
      else if(TString(histname.c_str()).Contains("nonb")){
        if(TString(histname.c_str()).Contains("nj6"))leg->SetHeader("W region, nonb, 6 #leq jet");
        else if(TString(histname.c_str()).Contains("nj35"))leg->SetHeader("W region, nonb, 3 #leq jet #leq 5");
        else leg->SetHeader("W region, nonb");
        }
      else {
        if(TString(histname.c_str()).Contains("nj6"))leg->SetHeader("W region, 6 #leq jet");
        else if(TString(histname.c_str()).Contains("nj35"))leg->SetHeader("W region, 3 #leq jet #leq 5");
        else leg->SetHeader("W region");
        }
    }
    else if(TString(histname.c_str()).Contains("_Z")){
      if(TString(histname.c_str()).Contains("nj6"))leg->SetHeader("Z region, 6 #leq jet");
      else if(TString(histname.c_str()).Contains("nj35"))leg->SetHeader("Z region, 3 #leq jet #leq 5");
      else leg->SetHeader("Z region");
    }

		char legentry[100];

		//leg->SetHeader("The legend title");

		leg->SetNColumns(4);
    leg->SetTextSize(0.03);

    if(!(TString(histname.c_str()).Contains("_S"))){
		  sprintf(legentry, "%.2f",HIST[0][histname]->Integral());
		  leg->AddEntry(HIST[0][histname],RootFileName[0],"L P E");
		  leg->AddEntry((TObject*)0,legentry,"");
    }

		sprintf(legentry, "%.2f",HIST[1][histname]->Integral());
		//leg->AddEntry(HIST[1][histname],RootFileName[1],"L");
		leg->AddEntry(HIST[1][histname],"T5ttcc","L");
		leg->AddEntry((TObject*)0,legentry,"");

		sprintf(legentry, "%.2f",HIST[2][histname]->Integral());
		//leg->AddEntry(HIST[2][histname],RootFileName[2],"L");
		leg->AddEntry(HIST[2][histname],"T5tttt","L");
		leg->AddEntry((TObject*)0,legentry,"");

		sprintf(legentry, "%.2f",HIST[3][histname]->Integral());
		//leg->AddEntry(HIST[3][histname],RootFileName[3],"L");
		leg->AddEntry(HIST[3][histname],"T1tttt","L");
		leg->AddEntry((TObject*)0,legentry,"");

		sprintf(legentry, "%.2f",HIST[4][histname]->Integral());
		//leg->AddEntry(HIST[4][histname],RootFileName[4],"L");
		leg->AddEntry(HIST[4][histname],"T2tt","L");
		leg->AddEntry((TObject*)0,legentry,"");

		for(int i=5;i<Num;i++){
      //cout << RootFileName[i] << endl;
			sprintf(legentry, "%.2f",HIST[i][histname]->Integral());
			if((TString)RootFileName[i]=="TT") leg->AddEntry(HIST[i][histname],"t#bar{t}","f");
			else if((TString)RootFileName[i]=="ST") leg->AddEntry(HIST[i][histname],"Top","f");
			else if((TString)RootFileName[i]=="DYToLL") leg->AddEntry(HIST[i][histname],"Drell-Yan","f");
			else if((TString)RootFileName[i]=="GJet") leg->AddEntry(HIST[i][histname],"#gamma + jets","f");
			else if((TString)RootFileName[i]=="Multiboson") leg->AddEntry(HIST[i][histname],"VV(V) + t#bar{t}(X)","f");
			else if((TString)RootFileName[i]=="ZJet") leg->AddEntry(HIST[i][histname],"Z(#rightarrow#nu#nu)","f");
			else if((TString)RootFileName[i]=="WJet") leg->AddEntry(HIST[i][histname],"W(#rightarrowl#nu)","f");
		  else leg->AddEntry(HIST[i][histname],RootFileName[i],"f");
			leg->AddEntry((TObject*)0,legentry,"");
		}

		//sprintf(legentry, "%.2f",HIST[Num-1][histname]->Integral());
		//leg->AddEntry(HIST[Num-1][histname],SigName,"f");
		//leg->AddEntry((TObject*)0,legentry,"");

		sprintf(legentry, "%.2f",Tot_MC->Integral());
		leg->AddEntry(Tot_MC,"Tot_MC","");
		leg->AddEntry((TObject*)0,legentry,"");

		leg->Draw("SAME");

    //TText *t1 = new TText(0.2,0.8,"CMS Work in progress  35.9 fb^{-1} (13 TeV)");
    //t1->SetX(0);
    //t1->SetY(7e8);
    //t1->SetTextAlign(12);
    //t1->Draw();
		//HIST[0][histname]->Sumw2();
		//HIST[1][histname]->Sumw2();
		//HIST[2][histname]->Sumw2();
		//Tot_MC->Sumw2();



		TH1F* dividend = new TH1F(*HIST[0][histname]);
		TH1F* dividend1 = new TH1F(*HIST[0][histname]);
    
    if(!(TString(histname.c_str()).Contains("_S"))){
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
		

		/*****************************     Xaxis Title Settings ******************************************/
		//if (strstr(histname.c_str(),"Energy")!=NULL) dividend->GetXaxis()->SetTitle("Energy (GeV)");
		//else if (strstr(histname.c_str(),"pt")!=NULL) dividend->GetXaxis()->SetTitle("pT (GeV)");
		//else if (strstr(histname.c_str(),"eta")!=NULL) dividend->GetXaxis()->SetTitle("#eta");
		//else if (strstr(histname.c_str(),"phi")!=NULL) dividend->GetXaxis()->SetTitle("#phi");
		//else if (strstr(histname.c_str(),"Mass")!=NULL) dividend->GetXaxis()->SetTitle("Mass (GeV)");
		//else if (strstr(histname.c_str(),"mass")!=NULL) dividend->GetXaxis()->SetTitle("Mass (GeV)");

		dividend->Draw("PE1");
		dividend1->Draw("sameE2");
   }


		/*************  save canvases  ************************/
		string dir="./Histogram/";
		string eps=".png";
		string CanvasName=dir+histname+eps;
		c1->SaveAs(CanvasName.c_str());
		/******************************************************/

    c1->SetName(histname.c_str());
    c1->Write();

		/**** NBins.txt************************/
		char Title[30];
		sprintf(Title,"%s",HIST[0][histname]->GetName());
		char Nbin[15];
		sprintf(Nbin,"%1.0d",HIST[0][histname]->GetNbinsX());
		char BWidth[15];
		sprintf(BWidth,"%1lf",HIST[0][histname]->GetBinWidth(0));
		fprintf(NBins,"Key = %s\nNum.of bin = %s\nBin Width = %s\n\n",Title,Nbin,BWidth);
		/***************************************/
		//string rara;
		//cout<<HIST[0][histname]->GetBinError();//0->underflow
		//cin>>rara;
	}////// Iterator loop
	
	NumOfKeys=histnum;
	cout<<"///////////////////////////////////"<<histnum<<" histograms is stacked"<<endl;

	fclose(NBins);

return totmc;

}


void SaveYields( char RootFileName[][200], map<TString, TH1F *> HIST[], map<TString, TH1F *> MCtot ){
//GetTitle(),GetName(),GetBinError()//0->underflow
	FILE *table = fopen("Yields.txt","w");

	string histname;
	string RooName[Num];
	
	map <TString , double > Integrals[Num];
	map <TString , double > Integrals_TOT_MC;
	map <TString , int > Nbins;

	for(int i=0;i<Num;i++){
		RooName[i]=RootFileName[i];
	}

	for(int i=0;i<Num;i++){
		RooName[i]=underscore2space(RooName[i]);
	}
	
	for(int i=0;i<Num;i++){
		cout<<RooName[i]<<endl;
	}
	for(int i=0;i<Num;i++)
	{

		for( std::map<TString, TH1F *>::iterator histoIter = HIST[0].begin(); histoIter != HIST[0].end(); ++histoIter)
		{	
			histname = (*histoIter).first;
			Integrals[i][histname]=HIST[i][histname]->Integral();
		}
	}

	for( std::map<TString, TH1F *>::iterator histoIter = HIST[0].begin(); histoIter != HIST[0].end(); ++histoIter)
	{
		histname = (*histoIter).first;
		Integrals_TOT_MC[histname]=MCtot[histname]->Integral();
		Nbins[histname]=MCtot[histname]->GetNbinsX();
	}// save Entries, Nbins

	
	for( std::map<TString, TH1F *>::iterator histoIter = HIST[0].begin(); histoIter != HIST[0].end(); ++histoIter)
	{
		histname = (*histoIter).first;
		for(int i=0;i<Num;i++){
			HIST[i][histname]->Rebin(Nbins[histname]);
		}
		MCtot[histname]->Rebin(Nbins[histname]);

	}// NBins -> 1Bin // uncertainty calculation

	
	fprintf(table,"\\documentclass[10pt]{article}\n");
	fprintf(table,"\\usepackage{graphicx}\n");
	fprintf(table,"\\begin{document}\n");
	fprintf(table,"\\begin{center}\n");
	fprintf(table,"\\rotatebox{270}{\n");
	fprintf(table,"\\scalebox{0.7}{\n");
	fprintf(table,"\\begin{tabular}{");
	for(int i=0;i<9;i++){
		fprintf(table,"c");
	}
	fprintf(table,"}\n");
	fprintf(table,"\\hline\\hline\n");
	fprintf(table,"Selection & Di-lepton & Z Mass veto & Two jets & MET & 1 or more bTagging & 2 or more bTagging & exactly 2 bTagging & Top Reconstruction\n");
	fprintf(table,"\\\\\n");
	fprintf(table," & S1 & S2 & S3 & S4 & S5 & S6 & S7 & S8\n");
	fprintf(table,"\\\\\n");
	fprintf(table,"\\hline\n");
	
	for(int i=1;i<Num;i++){
		fprintf(table,"%s &%1.0lf$\\pm%1.1lf$ &%1.0lf$\\pm%1.1lf$ &%1.0lf$\\pm%1.1lf$ &%1.0lf$\\pm%1.1lf$ &%1.0lf$\\pm%1.1lf$ &%1.0lf$\\pm%1.1lf$ &%1.0lf$\\pm%1.1lf$ &%1.0lf$\\pm%1.1lf$\\\\\n",RooName[i].c_str(),Integrals[i]["_cf_Lep1pt_1_"],HIST[i]["_cf_Lep1pt_1_"]->GetBinError(1),Integrals[i]["_cf_Lep1pt_2_"],HIST[i]["_cf_Lep1pt_2_"]->GetBinError(1),Integrals[i]["_cf_Lep1pt_3_"],HIST[i]["_cf_Lep1pt_3_"]->GetBinError(1),Integrals[i]["_cf_Lep1pt_4_"],HIST[i]["_cf_Lep1pt_4_"]->GetBinError(1),Integrals[i]["_cf_Lep1pt_5_"],HIST[i]["_cf_Lep1pt_5_"]->GetBinError(1),Integrals[i]["_cf_Lep1pt_6_"],HIST[i]["_cf_Lep1pt_6_"]->GetBinError(1),Integrals[i]["_cf_Lep1pt_7_"],HIST[i]["_cf_Lep1pt_7_"]->GetBinError(1),Integrals[i]["TopMass"],HIST[i]["TopMass"]->GetBinError(1));
	if(i==Num-2) fprintf(table,"\\hline\n");
	}
	
	fprintf(table,"\\hline\\hline\n");
	
	
	fprintf(table,"MC total &%1.0lf$\\pm%1.1lf$ &%1.0lf$\\pm%1.1lf$ &%1.0lf$\\pm%1.1lf$ &%1.0lf$\\pm%1.1lf$ &%1.0lf$\\pm%1.1lf$ &%1.0lf$\\pm%1.1lf$ &%1.0lf$\\pm%1.1lf$ &%1.0lf$\\pm%1.1lf$\\\\\n",Integrals_TOT_MC["_cf_Lep1pt_1_"],MCtot["_cf_Lep1pt_1_"]->GetBinError(1),Integrals_TOT_MC["_cf_Lep1pt_2_"],MCtot["_cf_Lep1pt_2_"]->GetBinError(1),Integrals_TOT_MC["_cf_Lep1pt_3_"],MCtot["_cf_Lep1pt_3_"]->GetBinError(1),Integrals_TOT_MC["_cf_Lep1pt_4_"],MCtot["_cf_Lep1pt_4_"]->GetBinError(1),Integrals_TOT_MC["_cf_Lep1pt_5_"],MCtot["_cf_Lep1pt_5_"]->GetBinError(1),Integrals_TOT_MC["_cf_Lep1pt_6_"],MCtot["_cf_Lep1pt_6_"]->GetBinError(1),Integrals_TOT_MC["_cf_Lep1pt_7_"],MCtot["_cf_Lep1pt_7_"]->GetBinError(1),Integrals_TOT_MC["TopMass"],MCtot["TopMass"]->GetBinError(1));
	
	fprintf(table,"\\hline\n");
	
	fprintf(table,"Data &%1.0lf$\\pm%1.1lf$ &%1.0lf$\\pm%1.1lf$ &%1.0lf$\\pm%1.1lf$ &%1.0lf$\\pm%1.1lf$ &%1.0lf$\\pm%1.1lf$ &%1.0lf$\\pm%1.1lf$ &%1.0lf$\\pm%1.1lf$ &%1.0lf$\\pm%1.1lf$\\\\\n",Integrals[0]["_cf_Lep1pt_1_"],HIST[0]["_cf_Lep1pt_1_"]->GetBinError(1),Integrals[0]["_cf_Lep1pt_2_"],HIST[0]["_cf_Lep1pt_2_"]->GetBinError(1),Integrals[0]["_cf_Lep1pt_3_"],HIST[0]["_cf_Lep1pt_3_"]->GetBinError(1),Integrals[0]["_cf_Lep1pt_4_"],HIST[0]["_cf_Lep1pt_4_"]->GetBinError(1),Integrals[0]["_cf_Lep1pt_5_"],HIST[0]["_cf_Lep1pt_5_"]->GetBinError(1),Integrals[0]["_cf_Lep1pt_6_"],HIST[0]["_cf_Lep1pt_6_"]->GetBinError(1),Integrals[0]["_cf_Lep1pt_7_"],HIST[0]["_cf_Lep1pt_7_"]->GetBinError(1),Integrals[0]["TopMass"],HIST[0]["TopMass"]->GetBinError(1));

	fprintf(table,"\\hline\\hline\n");
	fprintf(table,"\\\\\n");
	fprintf(table,"\\end{tabular}\n");
	fprintf(table,"}\n");
	fprintf(table,"}\n");
	fprintf(table,"\\end{center}\n");
	fprintf(table,"\\end{document}\n");
	
	fclose(table);
/*
	cout<<Integrals[0]["AnLepEnergy"]<<endl<<Integrals[10]["AnLepEnergy"]<<endl<<Integrals_TOT_MC["AnLepEnergy"]<<endl;

	cout<<HIST[0]["AnLepEnergy"]->GetBinLowEdge(2)<<endl;

	cout<<Nbins["AnLepEnergy"]<<endl;
	cout<<HIST[0]["AnLepEnergy"]->GetNbinsX()<<endl<<MCtot["AnLepEnergy"]->GetNbinsX()<<endl;
	cout<<HIST[0]["Lep1pt"]->GetBinError(1)<<HIST[0]["_cf_Jet1pt_5_"]->GetBinError(1)<<endl;
	cout<<MCtot["Lep1pt_2b"]->Integral()<<endl<<MCtot["_cf_Lep1pt_7_"]->Integral(); */
}
