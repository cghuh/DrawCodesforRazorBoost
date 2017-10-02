import sys
from string import *
from ROOT import TFile
from math import sqrt 

#def printCut(cut,info,signame="T1ttcc",option="counts"):
def printCut(cut,info,signame="T1tttt_mGluino_1400_mLSP_300",option="counts"):
    #ttj = 0
    wj = 0
    qcd = 0
    diboson = 0
    #top = 0
    tt = 0
    zjets = 0
    st = 0
    #zll = 0
    #triboson = 0
    ttX = 0
    #Wbb = 0
    #DYhad = 0
    DYJet = 0
    sig = 0
    data = 0
    for sample,t in info.iteritems():
        #print sample, t, cut
        #print "error"
        #print sample
        #print t
        #print cut
        #c = t[0][cut]
        c = t[cut]
        print c
        #w = t[1]
        effxs = c 
        #effxs = 0
        #if "TTJets" in sample:
        #    ttj = ttj + effxs
        if "WJet" in sample:
            wj = wj + effxs
        if "Multijet" in sample:
            qcd = qcd + effxs
        #if "WW" in sample:
        #    diboson = diboson + effxs
        #if "WZ" in sample:
        #    diboson = diboson + effxs
        #if "ZZ" in sample:
        #    diboson = diboson + effxs
        if "Multiboson" in sample:
            diboson = diboson + effxs
        if "TT.root" in sample:
            tt = tt + effxs
        if "ST" in sample:
            st = st + effxs
        if "ZJet" in sample:
            zjets = zjets + effxs
        if "DYToLL" in sample: 
            DYJet = DYJet + effxs
        #if "_WZZ_" in sample or "_WWZ" in sample or "_WWW_" in sample or "_WWGJets_" in sample or "_ZZZ" in sample:
        #    triboson = triboson + effxs
        #if "TTGJets_" in sample or "_TTbarW" in sample or "_ttbarZ_" in sample or "_TTWWJets_" in sample:
        if "GJet" in sample:
            ttX = ttX + effxs
        #if "_DYToCC_" in sample or "_DYToBB_" in sample:
        #    DYhad = DYhad + effxs
        #if "_Wbb_" in sample:
        #    Wbb = Wbb + effxs
        if signame in sample:
            sig = sig + effxs
        if "data" in sample:
            data = data + effxs
    total = tt + wj + qcd + diboson + st + zjets + ttX + DYJet
    if option == "counts":
        print "Cut %s \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d" % (cut,qcd,tt,wj,diboson,st,zjets,ttX,DYJet,sig,data)
        row = "%s & %.4g & %.4g & %.4g & %.4g & %.4g & %.4g & %.4g & %.4g & %d & %.4g & %.4g \\\\ \n" % (cut,qcd,tt,wj,diboson,st,zjets,ttX,DYJet,total,sig,data) 
        return row.replace("%","\%")
    
    elif option == "Percentage":
        if total == 0:
            total = 1
        print "Selection %s \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f" % (cut,qcd*100./total,tt*100./total,wj*100./total,diboson*100./total,st*100./total,zjets*100./total,ttX*100./total,DYJet*100./total)
        row = " & {0:.1%} & {1:.1%} & {2:.1%} & {3:.1%} & {4:.1%} & {5:.1%} & {6:.1%} & {7:.1%} \\\\ \n".format(qcd/total,tt/total,wj/total,diboson/total,st/total,zjets/total,ttX/total,DYJet/total )
        return row.replace("%","\%")

    
def getTTjetsComposition(cut,counts_allhad,counts_semilep,counts_dilep):
    allhad = counts_allhad[cut]
    semilep = counts_semilep[cut]
    dilep = counts_dilep[cut]
    total = allhad + semilep + dilep
    total = float(total)
    if total <= 0:
        print "Cut %s \t / \t / \t / " % cut
    else: 
        print "Cut %s \t %.2f \t %.2f \t %.2f" % (cut,allhad/total,semilep/total,dilep/total)


def makeCutflowTable(info,intlumi,cuts,signame):
    table = open("cutflow_TopAna.tex",'w')
    # get the table header stuff
    header = """\\begin{sidewaystable}[p]
\\centering
\\caption{Cutflow table, event counts are normalized to $%(lumi)s\\textrm{fb}^{-1}$. }
\\begin{tabular}{| l || c | c | c | c | c | c | c | c |}
\\hline
Selection & Multijet & TT & WJets & Multiboson & ST & ZJets & GJet & DYToLL  \\\\ \\hline
""" % {"lumi":intlumi}
    table.write(header)

    # put the body of the table
    for cut in cuts:
        #table.write( printCut(cut,info,signame,"counts") )
        table.write( printCut(cut,info,signame,"Percentage") )
        table.write("\\hline\n")
        
    # put the table closing stuff
    footer = """\\end{tabular}
\\label{tab:cutflow}
\\end{sidewaystable}
"""
    table.write(footer)
    
    table.close()

def allcuts():
    #cuts = ["NoCuts","Cleaning","Pileup","ISR","TopPt",
    #        "HCAL_noise","vertexg0","njetge3","HLT","jet1ptg200",
    #        "SIG","neleeq0","nmueq0","trackIso",
    #        "g1Mb0Ll","g1Mbg1W0Ll","1Mbg1W0Ll","g2Mbg1W0Ll","g1Mbg1W0Ll_mdPhiHat4","g1Mbg1W0Ll_mdPhiHatg4","g1Mb0Wg1uW0Ll",
    #        "0Lb0Ll","0Lbg1uW0Ll","0Lbg1uW0Ll_mdPhi0p3","0Lbg1uW0Ll_mdPhi0p5","0Lbg1uW0Ll_mdPhiHat4","0Lbg1uW0Ll_mdPhiHat5","0Lbg1W0Ll",
    #        "1Ll","g1Mb1Ll","g1Mbg1W1Ll",
    #        "g1Mbg1W1LlmT100","1Mbg1W1LlmT100","g2Mbg1W1LlmT100","g1Mbg1W1LlmT100_mdPhiHat4","g1Mbg1W1LlmT100_mdPhiHatg4","g1Mbg1W1LlmT",
    #        "0Lb1Ll","0Lbg1Y1Ll","0Lbg1Y1LlmT100",
    #        "0Lbg1Y1LlmT","0Lbg1Y1LlmT_mdPhiHat4","0Lbg1Y1LlmT_mdPhiHatg4",
    #        "2munoZmass","2mu","2mu0el","0Lb2mu0el","0Lbg1Y2mu0el","g1Mb2mu0el","g1Mbg1Y2mu0el",
    #        "2elnoZmass","2el","2el0mu","0Lb2el0mu","0Lbg1Y2el0mu","g1Mb2el0mu","g1Mbg1Y2el0mu",
    #        "2lnoZmass","2l","2l0ol","0Lb2l0ol","0Lbg1Y2l0ol","g1Mb2l0ol","g1Mbg1Y2l0ol",
    #        "g1Mbg1W0Ll_mdPhig0p3","g1Mbg1W0Ll_mdPhi0p3","g1Mbg1W1LlmT100_mdPhig0p3","g1Mbg1W1LlmT100_mdPhi0p3",
    #        "0Lbg1Y1LlmT_mdPhig0p3","0Lbg1Y1LlmT_mdPhi0p3",
    #        "g1Mbg1W0Ll_mdPhig0p5","g1Mbg1W0Ll_mdPhi0p5","g1Mbg1W1LlmT100_mdPhig0p5","g1Mbg1W1LlmT100_mdPhi0p5",
    #        "0Lbg1Y1LlmT_mdPhig0p5","0Lbg1Y1LlmT_mdPhi0p5",
    #        ]
    #cuts = ["NoCuts","GoodVertex","1 AK8Jet","nJet >= 3","MR > 800, R2 > 0.08", "Trigger",
    #        "No lepton", "nb >= 1", "nW >= 1", "mDPhi >= 0.4", "S",
    #        "No lepton", "nb = 0", "nW >= 1", "mDPhi < 0.25", "Q",
    #        "1 lepton", "nb >= 1", "nW >= 1", "mDPhi >= 0.4", "mT < 100", "T",
    #        "1 lepton", "nb = 0", "nW >= 1", "mDPhi >= 0.4", "30 <= mT < 100", "W",
    #        ]
    cuts = [
#"NoCuts", "Skim_1JetAK8",
# "Baseline_3Jet", "Baseline_MR_R2","Q_cut_HLT",
            #"X_cut_3Jet", "X_cut_MR_R2",
"X_cut_HLT",
						#"A_cut_0Ele", "A_cut_0Mu", "A_cut_1Top", 
"A_cut_mDPhi",
						#"X_cut_0Ele", "X_cut_0Mu", "X_cut_1aTop", 
"X_cut_mDPhi",
						#"V_cut_1Lep", "V_cut_1Top", "V_cut_mDPhi", 
"V_cut_MT",
						#"C_cut_1Lep", "C_cut_1mTop", "C_cut_mDPhi", 
"C_cut_MT",
						#"F_cut_2Lep", "F_cut_1mW", "F_cut_mDPhill", 
"F_cut_Mll",
"x_cut_mDPhi",
"a_cut_InvmDPhi",
            ]
    return cuts

if __name__ == '__main__':

    if len(sys.argv) < 2:
        print "Run as: python %s <samples file> " % (sys.argv[0])
        sys.exit()
    samplesfile = sys.argv[1]
    
    #datasets = open(samplesfile).readlines()
    datasets = open(samplesfile).read().splitlines()

    names = {}
    info = {}
    counts_allhad = {}
    counts_semilep = {}
    counts_dilep = {}
        
    # Integrated luminosity in fb-1s
    intlumi = 35.9

    # input information
    #inputdir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140522_noISR_btag_TopPt_newWtagger_eta2p4_wWtag_oldmass"
    #inputdir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140610_FullStatusReport"
    #inputdir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140729_preApp_comments"
    #inputdir = "/afs/cern.ch/work/n/nstrobbe/RazorBoost/GIT/Results/results_20140814"
    inputdir = "/gatbawi/palgongsan/chuh/susy/susy170830/"
    #analyzer = "rzrBoostMC"
    #signame = "T1ttcc_1000_325_300"
    signame = "T1tttt_mGluino_1400_mLSP_300"
    sigxs = 0.0243547
    
    # get signal info
    #fs = TFile.Open(inputdir + "/summary/" + analyzer + "_" + signame + ".root") #It should be fixed.
    fs = TFile.Open(inputdir + signame + ".root") #It should be fixed.
    histos = fs.Get("counts")
    countss = {}
    for bin in range(histos.GetNbinsX()):
        countss[histos.GetXaxis().GetBinLabel(bin+1)] = histos.GetBinContent(bin+1)
    #weights = sigxs/countss["NoCuts"]
    #info[signame] = (countss,weights)
    info[signame] = (countss)
    
    for d in datasets:
        #xsect = -1
        #totweight = -1
        #lumi = -1
        #d = strip(d)
        #d = split(d)
        #if len(d) >= 2:
        #    sampletype = d[1]
        #    if sampletype == "mc":
        #        lumi = intlumi
        #if len(d) >= 3:
        #    xsect = d[2]
        #if len(d) == 4:
        #    totweight = d[3]
        #names['xsect'] = xsect
        #names['totweight'] = totweight 
        #names['lumi'] = lumi
        #d = split(d[0], '%')
        # print d
        #user = d[0]
        #sample = d[1]
        sample = d
        #print sample
        samplew_ = sample.replace('/', '_')
        names['samplew_'] = samplew_

        #print samplew_
        # get the count:
        #f = TFile.Open(inputdir+d+".root")
        #f = TFile.Open(d)
        #f = TFile.Open("/home/cghuh3811/susy170201/added/data.root") 
        print d
        #d = split()
        #print d
        f = TFile.Open(d)
        histo = f.Get("counts")
        counts = {}
        for bin in range(histo.GetNbinsX()):
            counts[histo.GetXaxis().GetBinLabel(bin+1)] = histo.GetBinContent(bin+1)
        #print "sample, xsec, counts", samplew_, xsect, counts["NoCuts"]
        if counts["NoCuts"] == 0: continue
        #weight = float(xsect)/counts["NoCuts"]
        #info[samplew_] = (counts,weight)
        info[samplew_] = (counts)

        # now also get TTbar composition
        #if "_TTJets_" in samplew_:
        #    h_allhad = f.Get("counts_TTallhad")
        #    h_semilep = f.Get("counts_TTsemilep")
        #    h_dilep = f.Get("counts_TTdilep")
        #    for bin in range(h_allhad.GetNbinsX()):
        #        counts_allhad[h_allhad.GetXaxis().GetBinLabel(bin+1)] = h_allhad.GetBinContent(bin+1)
        #    for bin in range(h_semilep.GetNbinsX()):
        #        counts_semilep[h_semilep.GetXaxis().GetBinLabel(bin+1)] = h_semilep.GetBinContent(bin+1)
        #    for bin in range(h_dilep.GetNbinsX()):
        #        counts_dilep[h_dilep.GetXaxis().GetBinLabel(bin+1)] = h_dilep.GetBinContent(bin+1)

        f.Close()


    cuts = allcuts()
    print "Make cutflow table"

    print "Cut \t \t qcd \t tt \t wj \t diboson \t single top \t ZJet \t TTX \t DYJet \t signal \t data \n"
    #print "counts and percentage for %d fb-1 of data" % (intlumi)
    makeCutflowTable(info,intlumi,cuts,signame)
    print "No error in here"

    print "\n"
    
    #print "TTbar composition"
    #print "Cut \t allhad \t semilep \t dilep "
    #for cut in cuts:
    #    getTTjetsComposition(cut,counts_allhad,counts_semilep,counts_dilep)


