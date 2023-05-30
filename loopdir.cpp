void loopdir(TString dir="") {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gSystem->Exec("mkdir -p plot/"+dir+"/signal");
  gSystem->Exec("mkdir -p plot/"+dir+"/CR");
  TFile *f1 = TFile::Open("Plotter_out_2023_03_20.root");
  TString folder = "Counts_vs_MRR2";
  //TString folder = "Counts_vs_MRR2Bin";
  //TString folder = "BoostJetPt";
  //TString folder = "JetAK8Mass";
  //TString folder = "Counts_vs_newMRR2";
  //TString folder = "MRR2Bins";
  //TString folder = "NJetBins";
  f1->cd(folder);
  TIter keyList(gDirectory->GetListOfKeys());
  TKey *key;
  TCanvas* c1;
  TString name;
  while ((key = (TKey*)keyList())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TCanvas")) continue;
    c1 = (TCanvas*)key->ReadObj();
    name = c1->GetName();
    name.ReplaceAll(folder+"_","");
    if(!name.Contains("Ratio")) continue;
    if(name.Contains("Ele")) continue;
    if(name.Contains("Mu")) continue;
    if(name.Contains("MassTag")) continue;
    if(name.Contains("StackPlotSignal_"))       c1->SaveAs("plot/"+dir+"/signal/"+name+".png");
    else if(name.Contains("StackPlotHSignal_")) c1->SaveAs("plot/"+dir+"/signal/"+name+".png");
    else if(name.Contains("Boost"))             c1->SaveAs("plot/"+dir+"/CR/"+name+".png");
    else                                        c1->SaveAs("plot/"+dir+"/"+name+".png");
  }
}
