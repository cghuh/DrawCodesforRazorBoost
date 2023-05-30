void loopdir() {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gSystem->Exec("mkdir -p plot/error");
  TFile *f1 = TFile::Open("CF_Error.root");
  TIter keyList(gDirectory->GetListOfKeys());
  TKey *key;
  TCanvas* c1;
  TString name;
  while ((key = (TKey*)keyList())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TCanvas")) continue;
    c1 = (TCanvas*)key->ReadObj();
    name = c1->GetName();
    if(!TString(name).Contains("run2")) continue;
    //name.ReplaceAll(folder+"_","");
    c1->SaveAs("plot/error/"+name+".png");
  }
}
