import ROOT
from array import array
import math
import os
def Fit(w, target='',input="",bin=''):
    #print(bin)
    w.factory(target)
    tmp='NLL'
    tmp+=bin
    like1='expr::'
    like2='("-log(@0)",likelihood'
    like3=')'
    w.factory(like1+tmp+like2+bin+like3)
    nll = w.function(tmp)
    minim = ROOT.RooMinimizer(nll)
    minim.setErrorLevel(0.5)
    minim.minimize("Minuit2","migrad")
    bestfitnll = nll.getVal()
    fitresult = minim.save() 
    c1 = ROOT.TCanvas("c1","",900,900)
    c1.SetRightMargin(0.13)
    hist=fitresult.correlationHist()
    hist.Draw("colz")

    bins="_corr_bin"
    format=".png"
    dir="corMatrix/"
    oput=dir+input+bins+bin+format
    c1.Print(oput)

    cfname=['cf_Q1','cf_Q2','cf_Q3','cf_Q4','cf_Q5','cf_Q6','cf_T1','cf_T2','cf_T3','cf_T4','cf_T5','cf_T6','cf_W1','cf_W2','cf_W3','cf_W4','cf_W5','cf_W6']
    format="_fitr.txt"
    dir="fitresult/"
    oput = dir+input+format
    with open(oput, 'a') as f:
        for i in range(0, 3):
            if((str.__contains__(input, 'noniso') or str.__contains__(input, 'lt')) and str.__contains__(str(cfname[6*i+int(bin)-1]), 'Q')):
                f.writelines(str(0))
                f.write('\n')
                f.writelines(str(0))
                f.write('\n')
            else:
                f.writelines(str(fitresult.floatParsFinal().find(str(cfname[6*i+int(bin)-1])).getValV()))
                f.write('\n')
                f.writelines(str(fitresult.floatParsFinal().find(str(cfname[6*i+int(bin)-1])).getError()))
                f.write('\n')

    del c1
    del fitresult
    del nll
    del minim
    del w
    del bestfitnll
    del hist


def Calculate(input='test'):
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    w = ROOT.RooWorkspace("w")
    ROOT.RooFit.SumW2Error(True)
    dir="histnum/"
    format=".txt"
    ifile=dir+input+format
    print(input)

    w.factory('expr::q1("cf_Q1*NQQ1*pow(SQQ1,uncQQ1)+cf_T1*NQT1*pow(SQT1,uncQT1)+cf_W1*NQW1*pow(SQW1,uncQW1)", NQQ1[1], NQT1[1], NQW1[1], SQQ1[1], SQT1[1], SQW1[1], cf_Q1[0.9,0.1,1.5], cf_T1[0.8,0,4], cf_W1[0.8,0,4], uncQQ1[0.1,-2,2], uncQT1[0.1,-2,2], uncQW1[0.1,-2,2])')
    w.factory('expr::q2("cf_Q2*NQQ2*pow(SQQ2,uncQQ2)+cf_T2*NQT2*pow(SQT2,uncQT2)+cf_W2*NQW2*pow(SQW2,uncQW2)", NQQ2[1], NQT2[1], NQW2[1], SQQ2[1], SQT2[1], SQW2[1], cf_Q2[0.9,0.1,1.5], cf_T2[0.8,0,4], cf_W2[0.8,0,4], uncQQ2[0.1,-2,2], uncQT2[0.1,-2,2], uncQW2[0.1,-2,2])')
    w.factory('expr::q3("cf_Q3*NQQ3*pow(SQQ3,uncQQ3)+cf_T3*NQT3*pow(SQT3,uncQT3)+cf_W3*NQW3*pow(SQW3,uncQW3)", NQQ3[1], NQT3[1], NQW3[1], SQQ3[1], SQT3[1], SQW3[1], cf_Q3[0.9,0.1,1.5], cf_T3[0.8,0,4], cf_W3[0.8,0,4], uncQQ3[0.1,-2,2], uncQT3[0.1,-2,2], uncQW3[0.1,-2,2])')
    w.factory('expr::q4("cf_Q4*NQQ4*pow(SQQ4,uncQQ4)+cf_T4*NQT4*pow(SQT4,uncQT4)+cf_W4*NQW4*pow(SQW4,uncQW4)", NQQ4[1], NQT4[1], NQW4[1], SQQ4[1], SQT4[1], SQW4[1], cf_Q4[0.9,0.1,1.5], cf_T4[0.8,0,4], cf_W4[0.8,0,4], uncQQ4[0.1,-2,2], uncQT4[0.1,-2,2], uncQW4[0.1,-2,2])')
    w.factory('expr::q5("cf_Q5*NQQ5*pow(SQQ5,uncQQ5)+cf_T5*NQT5*pow(SQT5,uncQT5)+cf_W5*NQW5*pow(SQW5,uncQW5)", NQQ5[1], NQT5[1], NQW5[1], SQQ5[1], SQT5[1], SQW5[1], cf_Q5[0.9,0.1,1.5], cf_T5[0.8,0,4], cf_W5[0.8,0,4], uncQQ5[0.1,-2,2], uncQT5[0.1,-2,2], uncQW5[0.1,-2,2])')
    w.factory('expr::q6("cf_Q6*NQQ6*pow(SQQ6,uncQQ6)+cf_T6*NQT6*pow(SQT6,uncQT6)+cf_W6*NQW6*pow(SQW6,uncQW6)", NQQ6[1], NQT6[1], NQW6[1], SQQ6[1], SQT6[1], SQW6[1], cf_Q6[0.9,0.1,1.5], cf_T6[0.8,0,4], cf_W6[0.8,0,4], uncQQ6[0.1,-2,2], uncQT6[0.1,-2,2], uncQW6[0.1,-2,2])')
    w.factory('expr::t1("cf_Q1*NTQ1*pow(STQ1,uncTQ1)+cf_T1*NTT1*pow(STT1,uncTT1)+cf_W1*NTW1*pow(STW1,uncTW1)", NTQ1[1], NTT1[1], NTW1[1], STQ1[1], STT1[1], STW1[1], cf_Q1[0.9,0.1,1.5], cf_T1[0.8,0,4], cf_W1[0.8,0,4], uncTQ1[0.1,-2,2], uncTT1[0.1,-2,2], uncTW1[0.1,-2,2])')
    w.factory('expr::t2("cf_Q2*NTQ2*pow(STQ2,uncTQ2)+cf_T2*NTT2*pow(STT2,uncTT2)+cf_W2*NTW2*pow(STW2,uncTW2)", NTQ2[1], NTT2[1], NTW2[1], STQ2[1], STT2[1], STW2[1], cf_Q2[0.9,0.1,1.5], cf_T2[0.8,0,4], cf_W2[0.8,0,4], uncTQ2[0.1,-2,2], uncTT2[0.1,-2,2], uncTW2[0.1,-2,2])')
    w.factory('expr::t3("cf_Q3*NTQ3*pow(STQ3,uncTQ3)+cf_T3*NTT3*pow(STT3,uncTT3)+cf_W3*NTW3*pow(STW3,uncTW3)", NTQ3[1], NTT3[1], NTW3[1], STQ3[1], STT3[1], STW3[1], cf_Q3[0.9,0.1,1.5], cf_T3[0.8,0,4], cf_W3[0.8,0,4], uncTQ3[0.1,-2,2], uncTT3[0.1,-2,2], uncTW3[0.1,-2,2])')
    w.factory('expr::t4("cf_Q4*NTQ4*pow(STQ4,uncTQ4)+cf_T4*NTT4*pow(STT4,uncTT4)+cf_W4*NTW4*pow(STW4,uncTW4)", NTQ4[1], NTT4[1], NTW4[1], STQ4[1], STT4[1], STW4[1], cf_Q4[0.9,0.1,1.5], cf_T4[0.8,0,4], cf_W4[0.8,0,4], uncTQ4[0.1,-2,2], uncTT4[0.1,-2,2], uncTW4[0.1,-2,2])')
    w.factory('expr::t5("cf_Q5*NTQ5*pow(STQ5,uncTQ5)+cf_T5*NTT5*pow(STT5,uncTT5)+cf_W5*NTW5*pow(STW5,uncTW5)", NTQ5[1], NTT5[1], NTW5[1], STQ5[1], STT5[1], STW5[1], cf_Q5[0.9,0.1,1.5], cf_T5[0.8,0,4], cf_W5[0.8,0,4], uncTQ5[0.1,-2,2], uncTT5[0.1,-2,2], uncTW5[0.1,-2,2])')
    w.factory('expr::t6("cf_Q6*NTQ6*pow(STQ6,uncTQ6)+cf_T6*NTT6*pow(STT6,uncTT6)+cf_W6*NTW6*pow(STW6,uncTW6)", NTQ6[1], NTT6[1], NTW6[1], STQ6[1], STT6[1], STW6[1], cf_Q6[0.9,0.1,1.5], cf_T6[0.8,0,4], cf_W6[0.8,0,4], uncTQ6[0.1,-2,2], uncTT6[0.1,-2,2], uncTW6[0.1,-2,2])')
    w.factory('expr::w1("cf_Q1*NWQ1*pow(SWQ1,uncWQ1)+cf_T1*NWT1*pow(SWT1,uncWT1)+cf_W1*NWW1*pow(SWW1,uncWW1)", NWQ1[1], NWT1[1], NWW1[1], SWQ1[1], SWT1[1], SWW1[1], cf_Q1[0.9,0.1,1.5], cf_T1[0.8,0,4], cf_W1[0.8,0,4], uncWQ1[0.1,-2,2], uncWT1[0.1,-2,2], uncWW1[0.1,-2,2])')
    w.factory('expr::w2("cf_Q2*NWQ2*pow(SWQ2,uncWQ2)+cf_T2*NWT2*pow(SWT2,uncWT2)+cf_W2*NWW2*pow(SWW2,uncWW2)", NWQ2[1], NWT2[1], NWW2[1], SWQ2[1], SWT2[1], SWW2[1], cf_Q2[0.9,0.1,1.5], cf_T2[0.8,0,4], cf_W2[0.8,0,4], uncWQ2[0.1,-2,2], uncWT2[0.1,-2,2], uncWW2[0.1,-2,2])')
    w.factory('expr::w3("cf_Q3*NWQ3*pow(SWQ3,uncWQ3)+cf_T3*NWT3*pow(SWT3,uncWT3)+cf_W3*NWW3*pow(SWW3,uncWW3)", NWQ3[1], NWT3[1], NWW3[1], SWQ3[1], SWT3[1], SWW3[1], cf_Q3[0.9,0.1,1.5], cf_T3[0.8,0,4], cf_W3[0.8,0,4], uncWQ3[0.1,-2,2], uncWT3[0.1,-2,2], uncWW3[0.1,-2,2])')
    w.factory('expr::w4("cf_Q4*NWQ4*pow(SWQ4,uncWQ4)+cf_T4*NWT4*pow(SWT4,uncWT4)+cf_W4*NWW4*pow(SWW4,uncWW4)", NWQ4[1], NWT4[1], NWW4[1], SWQ4[1], SWT4[1], SWW4[1], cf_Q4[0.9,0.1,1.5], cf_T4[0.8,0,4], cf_W4[0.8,0,4], uncWQ4[0.1,-2,2], uncWT4[0.1,-2,2], uncWW4[0.1,-2,2])')
    w.factory('expr::w5("cf_Q5*NWQ5*pow(SWQ5,uncWQ5)+cf_T5*NWT5*pow(SWT5,uncWT5)+cf_W5*NWW5*pow(SWW5,uncWW5)", NWQ5[1], NWT5[1], NWW5[1], SWQ5[1], SWT5[1], SWW5[1], cf_Q5[0.9,0.1,1.5], cf_T5[0.8,0,4], cf_W5[0.8,0,4], uncWQ5[0.1,-2,2], uncWT5[0.1,-2,2], uncWW5[0.1,-2,2])')
    w.factory('expr::w6("cf_Q6*NWQ6*pow(SWQ6,uncWQ6)+cf_T6*NWT6*pow(SWT6,uncWT6)+cf_W6*NWW6*pow(SWW6,uncWW6)", NWQ6[1], NWT6[1], NWW6[1], SWQ6[1], SWT6[1], SWW6[1], cf_Q6[0.9,0.1,1.5], cf_T6[0.8,0,4], cf_W6[0.8,0,4], uncWQ6[0.1,-2,2], uncWT6[0.1,-2,2], uncWW6[0.1,-2,2])')

    #18+54+54 (N_data, N_MC, Unc_MC)
    #Call numbers from histograms and put them in the PDFs
    N = list()
    with open(ifile) as f:
        while True:
            contents = f.readline()
            if not contents:
                break
            N.append(float(contents))

    varname=['NQ1','NQ2','NQ3','NQ4','NQ5','NQ6','NT1','NT2','NT3','NT4','NT5','NT6','NW1','NW2','NW3','NW4','NW5','NW6','NQQ1','NQQ2','NQQ3','NQQ4','NQQ5','NQQ6','NQT1','NQT2','NQT3','NQT4','NQT5','NQT6','NQW1','NQW2','NQW3','NQW4','NQW5','NQW6','NTQ1','NTQ2','NTQ3','NTQ4','NTQ5','NTQ6','NTT1','NTT2','NTT3','NTT4','NTT5','NTT6','NTW1','NTW2','NTW3','NTW4','NTW5','NTW6','NWQ1','NWQ2','NWQ3','NWQ4','NWQ5','NWQ6','NWT1','NWT2','NWT3','NWT4','NWT5','NWT6','NWW1','NWW2','NWW3','NWW4','NWW5','NWW6','SQQ1','SQQ2','SQQ3','SQQ4','SQQ5','SQQ6','SQT1','SQT2','SQT3','SQT4','SQT5','SQT6','SQW1','SQW2','SQW3','SQW4','SQW5','SQW6','STQ1','STQ2','STQ3','STQ4','STQ5','STQ6','STT1','STT2','STT3','STT4','STT5','STT6','STW1','STW2','STW3','STW4','STW5','STW6','SWQ1','SWQ2','SWQ3','SWQ4','SWQ5','SWQ6','SWT1','SWT2','SWT3','SWT4','SWT5','SWT6','SWW1','SWW2','SWW3','SWW4','SWW5','SWW6']
    uncname=['uncQQ1','uncQQ2','uncQQ3','uncQQ4','uncQQ5','uncQQ6','uncQT1','uncQT2','uncQT3','uncQT4','uncQT5','uncQT6','uncQW1','uncQW2','uncQW3','uncQW4','uncQW5','uncQW6','uncTQ1','uncTQ2','uncTQ3','uncTQ4','uncTQ5','uncTQ6','uncTT1','uncTT2','uncTT3','uncTT4','uncTT5','uncTT6','uncTW1','uncTW2','uncTW3','uncTW4','uncTW5','uncTW6','uncWQ1','uncWQ2','uncWQ3','uncWQ4','uncWQ5','uncWQ6','uncWT1','uncWT2','uncWT3','uncWT4','uncWT5','uncWT6','uncWW1','uncWW2','uncWW3','uncWW4','uncWW5','uncWW6']
    cfname=['cf_Q1','cf_Q2','cf_Q3','cf_Q4','cf_Q5','cf_Q6','cf_T1','cf_T2','cf_T3','cf_T4','cf_T5','cf_T6','cf_W1','cf_W2','cf_W3','cf_W4','cf_W5','cf_W6']


    for i in range(18, len(varname)):
        w.var(str(varname[i])).setVal(float(N[i]))

    for i in range(0, len(uncname)):
        if((str.__contains__(input, 'noniso') or str.__contains__(input, 'lt')) and str.__contains__(str(uncname[i]), 'Q')):
            w.var(str(uncname[i])).setVal(0)
            w.var(str(uncname[i])).setConstant(True)

    for i in range(0, 18):
        if((str.__contains__(input, 'noniso') or str.__contains__(input, 'lt')) and str.__contains__(str(cfname[i]), 'cf_Q')):
            w.var(str(cfname[i])).setVal(0)
            w.var(str(cfname[i])).setConstant(True)
            w.var(str(cfname[i])).removeError()
        elif((str.__contains__(input, 'qtw_20171BoostMRR2Bin_step0') and str.__contains__(str(cfname[i]), 'cf_Q1')) or (str.__contains__(input, 'qtw_20181BoostMRR2Bin') and str.__contains__(str(cfname[i]), 'cf_Q1')) or (str.__contains__(input, 'qtw_20172BoostMRR2Bin_step0') and str.__contains__(str(cfname[i]), 'cf_Q1')) or (str.__contains__(input, 'qtw_20182BoostMRR2Bin_step0') and str.__contains__(str(cfname[i]), 'cf_Q1'))):
            w.var(str(cfname[i])).setVal(0.5)
            w.var(str(cfname[i])).setMax(1.2)
            w.var(str(cfname[i])).setMin(0.0)
        elif((str.__contains__(input, 'qtw_20161BoostMRR2Bin_step0') and str.__contains__(str(cfname[i]), 'cf_Q1')) or (str.__contains__(input, 'qtw_20182BoostNJetBins_step0')) or (str.__contains__(input, 'qtw_20181BoostNJetBins')) or (str.__contains__(input, 'qtw_20161BoostNJetBins_step0') and str.__contains__(str(cfname[i]), 'cf_Q2')) or (str.__contains__(input, 'qtw_20171BoostNJetBins_step0') and (str.__contains__(str(cfname[i]), 'cf_Q2') or str.__contains__(str(cfname[i]), 'cf_Q3'))) or (str.__contains__(input, 'qtw_20172BoostMRR2Bin_step1') and str.__contains__(str(cfname[i]), 'cf_Q1'))):
            w.var(str(cfname[i])).setVal(0.7)
            w.var(str(cfname[i])).setMax(1.2)
            w.var(str(cfname[i])).setMin(0.0)
        elif(str.__contains__(input, 'qtw_20171BoostMRR2Bin') and str.__contains__(str(cfname[i]), 'cf_Q2')):
            w.var(str(cfname[i])).setVal(1.0)
            w.var(str(cfname[i])).setMax(2.0)
            w.var(str(cfname[i])).setMin(0.5)

    w.factory('Poisson::poisQ1(NQ1[19],q1)')
    w.factory('Poisson::poisQ2(NQ2[19],q2)')
    w.factory('Poisson::poisQ3(NQ3[19],q3)')
    w.factory('Poisson::poisQ4(NQ4[19],q4)')
    w.factory('Poisson::poisQ5(NQ5[19],q5)')
    w.factory('Poisson::poisQ6(NQ6[19],q6)')
    w.factory('Poisson::poisT1(NT1[19],t1)')
    w.factory('Poisson::poisT2(NT2[19],t2)')
    w.factory('Poisson::poisT3(NT3[19],t3)')
    w.factory('Poisson::poisT4(NT4[19],t4)')
    w.factory('Poisson::poisT5(NT5[19],t5)')
    w.factory('Poisson::poisT6(NT6[19],t6)')
    w.factory('Poisson::poisW1(NW1[19],w1)')
    w.factory('Poisson::poisW2(NW2[19],w2)')
    w.factory('Poisson::poisW3(NW3[19],w3)')
    w.factory('Poisson::poisW4(NW4[19],w4)')
    w.factory('Poisson::poisW5(NW5[19],w5)')
    w.factory('Poisson::poisW6(NW6[19],w6)')

    w.factory('Gaussian::UncQQ1(uncQQ1,0,1)')
    w.factory('Gaussian::UncQQ2(uncQQ2,0,1)')
    w.factory('Gaussian::UncQQ3(uncQQ3,0,1)')
    w.factory('Gaussian::UncQQ4(uncQQ4,0,1)')
    w.factory('Gaussian::UncQQ5(uncQQ5,0,1)')
    w.factory('Gaussian::UncQQ6(uncQQ6,0,1)')
    w.factory('Gaussian::UncQT1(uncQT1,0,1)')
    w.factory('Gaussian::UncQT2(uncQT2,0,1)')
    w.factory('Gaussian::UncQT3(uncQT3,0,1)')
    w.factory('Gaussian::UncQT4(uncQT4,0,1)')
    w.factory('Gaussian::UncQT5(uncQT5,0,1)')
    w.factory('Gaussian::UncQT6(uncQT6,0,1)')
    w.factory('Gaussian::UncQW1(uncQW1,0,1)')
    w.factory('Gaussian::UncQW2(uncQW2,0,1)')
    w.factory('Gaussian::UncQW3(uncQW3,0,1)')
    w.factory('Gaussian::UncQW4(uncQW4,0,1)')
    w.factory('Gaussian::UncQW5(uncQW5,0,1)')
    w.factory('Gaussian::UncQW6(uncQW6,0,1)')
    w.factory('Gaussian::UncTQ1(uncTQ1,0,1)')
    w.factory('Gaussian::UncTQ2(uncTQ2,0,1)')
    w.factory('Gaussian::UncTQ3(uncTQ3,0,1)')
    w.factory('Gaussian::UncTQ4(uncTQ4,0,1)')
    w.factory('Gaussian::UncTQ5(uncTQ5,0,1)')
    w.factory('Gaussian::UncTQ6(uncTQ6,0,1)')
    w.factory('Gaussian::UncTT1(uncTT1,0,1)')
    w.factory('Gaussian::UncTT2(uncTT2,0,1)')
    w.factory('Gaussian::UncTT3(uncTT3,0,1)')
    w.factory('Gaussian::UncTT4(uncTT4,0,1)')
    w.factory('Gaussian::UncTT5(uncTT5,0,1)')
    w.factory('Gaussian::UncTT6(uncTT6,0,1)')
    w.factory('Gaussian::UncTW1(uncTW1,0,1)')
    w.factory('Gaussian::UncTW2(uncTW2,0,1)')
    w.factory('Gaussian::UncTW3(uncTW3,0,1)')
    w.factory('Gaussian::UncTW4(uncTW4,0,1)')
    w.factory('Gaussian::UncTW5(uncTW5,0,1)')
    w.factory('Gaussian::UncTW6(uncTW6,0,1)')
    w.factory('Gaussian::UncWQ1(uncWQ1,0,1)')
    w.factory('Gaussian::UncWQ2(uncWQ2,0,1)')
    w.factory('Gaussian::UncWQ3(uncWQ3,0,1)')
    w.factory('Gaussian::UncWQ4(uncWQ4,0,1)')
    w.factory('Gaussian::UncWQ5(uncWQ5,0,1)')
    w.factory('Gaussian::UncWQ6(uncWQ6,0,1)')
    w.factory('Gaussian::UncWT1(uncWT1,0,1)')
    w.factory('Gaussian::UncWT2(uncWT2,0,1)')
    w.factory('Gaussian::UncWT3(uncWT3,0,1)')
    w.factory('Gaussian::UncWT4(uncWT4,0,1)')
    w.factory('Gaussian::UncWT5(uncWT5,0,1)')
    w.factory('Gaussian::UncWT6(uncWT6,0,1)')
    w.factory('Gaussian::UncWW1(uncWW1,0,1)')
    w.factory('Gaussian::UncWW2(uncWW2,0,1)')
    w.factory('Gaussian::UncWW3(uncWW3,0,1)')
    w.factory('Gaussian::UncWW4(uncWW4,0,1)')
    w.factory('Gaussian::UncWW5(uncWW5,0,1)')
    w.factory('Gaussian::UncWW6(uncWW6,0,1)')

    #for i in range(len(varname)):
    for i in range(0, 18):
        w.var(str(varname[i])).setVal(float(N[i]))

    dir="fitresult/"
    format="_fitr.txt"
    oput = dir+input+format
    if(os.path.isfile(oput)):
        os.remove(oput)
    Fit(w,'PROD::likelihood1(poisQ1,poisT1,poisW1,UncQQ1,UncQT1,UncQW1,UncTQ1,UncTT1,UncTW1,UncWQ1,UncWT1,UncWW1)',input,'1')
    Fit(w,'PROD::likelihood2(poisQ2,poisT2,poisW2,UncQQ2,UncQT2,UncQW2,UncTQ2,UncTT2,UncTW2,UncWQ2,UncWT2,UncWW2)',input,'2')
    Fit(w,'PROD::likelihood3(poisQ3,poisT3,poisW3,UncQQ3,UncQT3,UncQW3,UncTQ3,UncTT3,UncTW3,UncWQ3,UncWT3,UncWW3)',input,'3')
    Fit(w,'PROD::likelihood4(poisQ4,poisT4,poisW4,UncQQ4,UncQT4,UncQW4,UncTQ4,UncTT4,UncTW4,UncWQ4,UncWT4,UncWW4)',input,'4')
    Fit(w,'PROD::likelihood5(poisQ5,poisT5,poisW5,UncQQ5,UncQT5,UncQW5,UncTQ5,UncTT5,UncTW5,UncWQ5,UncWT5,UncWW5)',input,'5')
    if(str.__contains__(input, 'MRR2')):
        Fit(w,'PROD::likelihood6(poisQ6,poisT6,poisW6,UncQQ6,UncQT6,UncQW6,UncTQ6,UncTT6,UncTW6,UncWQ6,UncWT6,UncWW6)',input,'6')
    

def DrawPlot(input=''):
    #print(input)
    file0 = ROOT.TFile('/Users/chuh/Dropbox/Analysis/razor/230613/run_2023_06_09_step1.root')
    if(str.__contains__(input, 'MRR2')):
        hist = file0.Get('Counts_vs_MRR2Bin/Syst_vs_MRR2Bin/Data_2018_CR_L17_1Boost')
    else:
        hist = file0.Get('NJetBins/Data_2018_CR_L17_1Boost')

    binsize=hist.GetNbinsX()
    region=2
    if(str.__contains__(input, 'qtw')):
        region=3

    format="_fitr.txt"
    dir="fitresult/"
    oput = dir+input+format
    binx, binerr_x, cf_Q, cf_T, cf_W = array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' )  
    cf_Q_up, cf_T_up, cf_W_up = array( 'd' ), array( 'd' ), array( 'd' )

    for i in range(0,binsize):
        binx.append(hist.GetXaxis().GetBinCenter(i+1))
        binerr_x.append(hist.GetXaxis().GetBinWidth(i+1)/2)

    check=0
    with open(oput) as f:
        while True:
            contents = f.readline()
            if not contents:
                break
            if(check%6==0 or check%6==1): 
                if(check%2==0):
                    cf_Q.append(float(contents))
                else:
                    cf_Q_up.append(float(contents))
            elif(check%6==2 or check%6==3): 
                if(check%2==0):
                    cf_T.append(float(contents))
                else:
                    cf_T_up.append(float(contents))
            else: 
                if(check%2==0):
                    cf_W.append(float(contents))
                else:
                    cf_W_up.append(float(contents))
            check+=1

    if(str.__contains__(input, 'MRR2')):
        output = ROOT.TFile("CFs.root", "update")
    else:
        output = ROOT.TFile("NJet_CFs.root", "update")
    g1 = ROOT.TGraphAsymmErrors(binsize, binx, cf_Q, binerr_x, binerr_x, cf_Q_up, cf_Q_up)
    g2 = ROOT.TGraphAsymmErrors(binsize, binx, cf_T, binerr_x, binerr_x, cf_T_up, cf_T_up)
    g3 = ROOT.TGraphAsymmErrors(binsize, binx, cf_W, binerr_x, binerr_x, cf_W_up, cf_W_up)
    g1.SetMarkerStyle(21)
    g2.SetMarkerStyle(21)
    g3.SetMarkerStyle(21)
    g1.SetMarkerSize(1)
    g2.SetMarkerSize(1)
    g3.SetMarkerSize(1)
    g1.SetMarkerColor(1)
    g2.SetMarkerColor(2)
    g3.SetMarkerColor(4)
    g1.SetLineColor(1)
    g2.SetLineColor(2)
    g3.SetLineColor(4)
    format='_cf_Q'
    g1.SetName(input+format)
    format='_cf_T'
    g2.SetName(input+format)
    format='_cf_W'
    g3.SetName(input+format)
    g1.Write()
    g2.Write()
    g3.Write()
    output.Close()

    mg = ROOT.TMultiGraph()

    if(str.__contains__(input, 'qtw')):
        mg.Add(g1)
    mg.Add(g2)
    mg.Add(g3)

    c2 = ROOT.TCanvas("c2","",800,800)
    if(str.__contains__(input, 'MRR2')):
        mg.GetXaxis().SetTitle('M_{R}xR^{2}')
    else:
        mg.GetXaxis().SetTitle('N_{AK4jets}')
    mg.GetYaxis().SetTitle('Correction Factor')
    mg.SetMaximum(3)
    mg.SetMinimum(0)
    mg.Draw("AP")
    
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.03)
    latex.SetTextAlign(31)
    latex.DrawLatex(0.90, 0.93, "138 fb^{-1} (#sqrt{s} = 13 TeV)")
    latex.SetTextAlign(11)
    latex.DrawLatex(0.10,0.93,"CMS Work in Progress")

    if(str.__contains__(input, 'noniso')):
        leg = ROOT.TLegend(0.13,0.73,0.88,0.88)
    else:
        leg = ROOT.TLegend(0.13,0.68,0.88,0.88)
    leg.SetFillStyle(0)
    leg.SetLineColor(0)
    leg.SetTextSize(0.03)

    name =''
    if(str.__contains__(input, '2016')):
        name = '2016, '
    elif(str.__contains__(input, '2017')):
        name = '2017, '
    else:
        name = '2018, '

    if(str.__contains__(input, 'qtw')):
        name += 'Hadronic/Leptonic background estimation;'
    elif(str.__contains__(input, 'lt')):
        name += 'Z(#rightarrow#nu#nu)+jet estimation;'
    else:
        name += 'Non-Isolated lepton background estimation;'

    leg.SetHeader(name)

    if(str.__contains__(input, '1Boost')):
        leg.AddEntry(0, "1 boost jet final state", "")
    elif(str.__contains__(input, '2Boost')):
        leg.AddEntry(0, "#geq 2 boost jet final state", "")

    if(str.__contains__(input, 'qtw')):
        leg.AddEntry(g1,  "Multijet CF", "p")
    leg.AddEntry(g2,  "Top(TT+ST) CF", "p")
    leg.AddEntry(g3,  "W(#rightarrowl#nu)+jet CF", "p")
    leg.Draw()

    dir="fithist/"
    format=".png"
    oput = dir+input+format
    c2.Print(oput)
    file0.Close()

    del c2
    del mg

    res = list()
    res.append(g1)
    res.append(g2)
    res.append(g3)
    del g1
    del g2
    del g3
    return res

def compare(ctg='', var='', boost='', cf=''):
    if(str.__contains__(var, 'MRR2')):    
        file0 = ROOT.TFile('CFs.root')
        gname16 = ctg+'_2016'+boost+var+'_step0'+cf
        gname17 = ctg+'_2017'+boost+var+'_step0'+cf
        gname18 = ctg+'_2018'+boost+var+'_step0'+cf
    else:
        file0 = ROOT.TFile('NJet_CFs.root')
        gname16 = ctg+'_2016'+boost+var+'_step1'+cf
        gname17 = ctg+'_2017'+boost+var+'_step1'+cf
        gname18 = ctg+'_2018'+boost+var+'_step1'+cf

    g16 = file0.Get(gname16)
    g17 = file0.Get(gname17)
    g18 = file0.Get(gname18)

    g16.SetMarkerStyle(21)
    g17.SetMarkerStyle(21)
    g18.SetMarkerStyle(21)
    g16.SetMarkerSize(1)
    g17.SetMarkerSize(1)
    g18.SetMarkerSize(1)
    g16.SetMarkerColor(1)
    g17.SetMarkerColor(2)
    g18.SetMarkerColor(4)
    g16.SetLineColor(1)
    g17.SetLineColor(2)
    g18.SetLineColor(4)

    mg = ROOT.TMultiGraph()

    mg.Add(g16)
    mg.Add(g17)
    mg.Add(g18)

    c2 = ROOT.TCanvas("c2","",800,800)
    if(str.__contains__(var, 'MRR2')):
        mg.GetXaxis().SetTitle('M_{R}xR^{2}')
    else:
        mg.GetXaxis().SetTitle('N_{AK4jets}')
    mg.GetYaxis().SetTitle('Correction Factor')
    mg.SetMaximum(3)
    mg.SetMinimum(0)
    mg.Draw("AP")
    
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.03)
    latex.SetTextAlign(31)
    latex.DrawLatex(0.90, 0.93, "138 fb^{-1} (#sqrt{s} = 13 TeV)")
    latex.SetTextAlign(11)
    latex.DrawLatex(0.10,0.93,"CMS Work in Progress")

    if(str.__contains__(ctg, 'noniso')):
        leg = ROOT.TLegend(0.13,0.73,0.88,0.88)
    else:
        leg = ROOT.TLegend(0.13,0.68,0.88,0.88)
    leg.SetFillStyle(0)
    leg.SetLineColor(0)
    leg.SetTextSize(0.03)
    name =''
    if(str.__contains__(ctg, 'qtw')):
        name = 'Hadronic/Leptonic background estimation;'
    elif(str.__contains__(ctg, 'lt')):
        name = 'Z(#rightarrow#nu#nu)+jet estimation;'
    else:
        name = 'Non-Isolated lepton background estimation;'

    if(str.__contains__(cf, '_cf_Q')):
        name += ' Multijet CF'
    elif(str.__contains__(cf, '_cf_T')):
        name += ' Top(TT+ST) CF'
    else:
        name += ' W(#rightarrowl#nu)+jet CF'

    leg.SetHeader(name)
    if(str.__contains__(boost, '1Boost')):
        leg.AddEntry(0, "1 boost jet final state", "")
    elif(str.__contains__(boost, '2Boost')):
        leg.AddEntry(0, "#geq 2 boost jet final state", "")
    leg.AddEntry(g16,  "2016", "p")
    leg.AddEntry(g17,  "2017", "p")
    leg.AddEntry(g18,  "2018", "p")
    leg.Draw()

    dir="fithist/"
    format=".png"
    oput = dir+ctg+boost+var+cf+format
    c2.Print(oput)

"""
Calculate('qtw_20161BoostMRR2Bin_step0')
Calculate('qtw_20162BoostMRR2Bin_step0')
Calculate('qtw_20161BoostNJetBins_step0')
Calculate('qtw_20162BoostNJetBins_step0')
Calculate('lt_20161BoostMRR21vlBin_step0')
Calculate('lt_20162BoostMRR21vlBin_step0')
Calculate('lt_20161BoostNJetBins_step0')
Calculate('lt_20162BoostNJetBins_step0')
Calculate('noniso_2016MRR2Bin_step0')
Calculate('noniso_2016NJetBins_step0')
Calculate('qtw_20171BoostMRR2Bin_step0')
Calculate('qtw_20172BoostMRR2Bin_step0')
Calculate('qtw_20171BoostNJetBins_step0')
Calculate('qtw_20172BoostNJetBins_step0')
Calculate('lt_20171BoostMRR21vlBin_step0')
Calculate('lt_20172BoostMRR21vlBin_step0')
Calculate('lt_20171BoostNJetBins_step0')
Calculate('lt_20172BoostNJetBins_step0')
Calculate('noniso_2017MRR2Bin_step0')
Calculate('noniso_2017NJetBins_step0')
Calculate('qtw_20181BoostMRR2Bin_step0')
Calculate('qtw_20182BoostMRR2Bin_step0')
Calculate('qtw_20181BoostNJetBins_step0')
Calculate('qtw_20182BoostNJetBins_step0')
Calculate('lt_20181BoostMRR21vlBin_step0')
Calculate('lt_20182BoostMRR21vlBin_step0')
Calculate('lt_20181BoostNJetBins_step0')
Calculate('lt_20182BoostNJetBins_step0')
Calculate('noniso_2018MRR2Bin_step0')
Calculate('noniso_2018NJetBins_step0')

Calculate('qtw_20161BoostMRR2Bin_step1')
Calculate('qtw_20162BoostMRR2Bin_step1')
Calculate('qtw_20161BoostNJetBins_step1')
Calculate('qtw_20162BoostNJetBins_step1')
Calculate('lt_20161BoostMRR21vlBin_step1')
Calculate('lt_20162BoostMRR21vlBin_step1')
Calculate('lt_20161BoostNJetBins_step1')
Calculate('lt_20162BoostNJetBins_step1')
Calculate('noniso_2016MRR2Bin_step1')
Calculate('noniso_2016NJetBins_step1')
Calculate('qtw_20171BoostMRR2Bin_step1')
Calculate('qtw_20172BoostMRR2Bin_step1')
Calculate('qtw_20171BoostNJetBins_step1')
Calculate('qtw_20172BoostNJetBins_step1')
Calculate('lt_20171BoostMRR21vlBin_step1')
Calculate('lt_20172BoostMRR21vlBin_step1')
Calculate('lt_20171BoostNJetBins_step1')
Calculate('lt_20172BoostNJetBins_step1')
Calculate('noniso_2017MRR2Bin_step1')
Calculate('noniso_2017NJetBins_step1')
Calculate('qtw_20181BoostMRR2Bin_step1')
Calculate('qtw_20182BoostMRR2Bin_step1')
Calculate('qtw_20181BoostNJetBins_step1')
Calculate('qtw_20182BoostNJetBins_step1')
Calculate('lt_20181BoostMRR21vlBin_step1')
Calculate('lt_20182BoostMRR21vlBin_step1')
Calculate('lt_20181BoostNJetBins_step1')
Calculate('lt_20182BoostNJetBins_step1')
Calculate('noniso_2018MRR2Bin_step1')
Calculate('noniso_2018NJetBins_step1') 
"""

#"""
if(os.path.isfile('CFs.root')):
    os.remove('CFs.root')
MRR2 = list()
MRR2.append(DrawPlot('qtw_20161BoostMRR2Bin_step0'))
MRR2.append(DrawPlot('qtw_20162BoostMRR2Bin_step0'))
MRR2.append(DrawPlot('lt_20161BoostMRR21vlBin_step0'))
MRR2.append(DrawPlot('lt_20162BoostMRR21vlBin_step0'))
MRR2.append(DrawPlot('noniso_2016MRR2Bin_step0'))
MRR2.append(DrawPlot('qtw_20171BoostMRR2Bin_step0'))
MRR2.append(DrawPlot('qtw_20172BoostMRR2Bin_step0'))
MRR2.append(DrawPlot('lt_20171BoostMRR21vlBin_step0'))
MRR2.append(DrawPlot('lt_20172BoostMRR21vlBin_step0'))
MRR2.append(DrawPlot('noniso_2017MRR2Bin_step0'))
MRR2.append(DrawPlot('qtw_20181BoostMRR2Bin_step0'))
MRR2.append(DrawPlot('qtw_20182BoostMRR2Bin_step0'))
MRR2.append(DrawPlot('lt_20181BoostMRR21vlBin_step0'))
MRR2.append(DrawPlot('lt_20182BoostMRR21vlBin_step0'))
MRR2.append(DrawPlot('noniso_2018MRR2Bin_step0'))
MRR2.append(DrawPlot('qtw_20161BoostMRR2Bin_step1'))
MRR2.append(DrawPlot('qtw_20162BoostMRR2Bin_step1'))
MRR2.append(DrawPlot('lt_20161BoostMRR21vlBin_step1'))
MRR2.append(DrawPlot('lt_20162BoostMRR21vlBin_step1'))
MRR2.append(DrawPlot('noniso_2016MRR2Bin_step1'))
MRR2.append(DrawPlot('qtw_20171BoostMRR2Bin_step1'))
MRR2.append(DrawPlot('qtw_20172BoostMRR2Bin_step1'))
MRR2.append(DrawPlot('lt_20171BoostMRR21vlBin_step1'))
MRR2.append(DrawPlot('lt_20172BoostMRR21vlBin_step1'))
MRR2.append(DrawPlot('noniso_2017MRR2Bin_step1'))
MRR2.append(DrawPlot('qtw_20181BoostMRR2Bin_step1'))
MRR2.append(DrawPlot('qtw_20182BoostMRR2Bin_step1'))
MRR2.append(DrawPlot('lt_20181BoostMRR21vlBin_step1'))
MRR2.append(DrawPlot('lt_20182BoostMRR21vlBin_step1'))
MRR2.append(DrawPlot('noniso_2018MRR2Bin_step1'))

if(os.path.isfile('NJet_CFs.root')):
    os.remove('NJet_CFs.root')
NJet = list()
NJet.append(DrawPlot('qtw_20161BoostNJetBins_step0'))
NJet.append(DrawPlot('qtw_20162BoostNJetBins_step0'))
NJet.append(DrawPlot('lt_20161BoostNJetBins_step0'))
NJet.append(DrawPlot('lt_20162BoostNJetBins_step0'))
NJet.append(DrawPlot('noniso_2016NJetBins_step0'))
NJet.append(DrawPlot('qtw_20171BoostNJetBins_step0'))
NJet.append(DrawPlot('qtw_20172BoostNJetBins_step0'))
NJet.append(DrawPlot('lt_20171BoostNJetBins_step0'))
NJet.append(DrawPlot('lt_20172BoostNJetBins_step0'))
NJet.append(DrawPlot('noniso_2017NJetBins_step0'))
NJet.append(DrawPlot('qtw_20181BoostNJetBins_step0'))
NJet.append(DrawPlot('qtw_20182BoostNJetBins_step0'))
NJet.append(DrawPlot('lt_20181BoostNJetBins_step0'))
NJet.append(DrawPlot('lt_20182BoostNJetBins_step0'))
NJet.append(DrawPlot('noniso_2018NJetBins_step0'))
NJet.append(DrawPlot('qtw_20161BoostNJetBins_step1'))
NJet.append(DrawPlot('qtw_20162BoostNJetBins_step1'))
NJet.append(DrawPlot('lt_20161BoostNJetBins_step1'))
NJet.append(DrawPlot('lt_20162BoostNJetBins_step1'))
NJet.append(DrawPlot('noniso_2016NJetBins_step1'))
NJet.append(DrawPlot('qtw_20171BoostNJetBins_step1'))
NJet.append(DrawPlot('qtw_20172BoostNJetBins_step1'))
NJet.append(DrawPlot('lt_20171BoostNJetBins_step1'))
NJet.append(DrawPlot('lt_20172BoostNJetBins_step1'))
NJet.append(DrawPlot('noniso_2017NJetBins_step1'))
NJet.append(DrawPlot('qtw_20181BoostNJetBins_step1'))
NJet.append(DrawPlot('qtw_20182BoostNJetBins_step1'))
NJet.append(DrawPlot('lt_20181BoostNJetBins_step1'))
NJet.append(DrawPlot('lt_20182BoostNJetBins_step1'))
NJet.append(DrawPlot('noniso_2018NJetBins_step1'))

compare('qtw', 'MRR2Bin',  '1Boost', '_cf_Q')
compare('qtw', 'MRR2Bin',  '2Boost', '_cf_Q')
compare('qtw', 'MRR2Bin',  '1Boost', '_cf_T')
compare('qtw', 'MRR2Bin',  '2Boost', '_cf_T')
compare('qtw', 'MRR2Bin',  '1Boost', '_cf_W')
compare('qtw', 'MRR2Bin',  '2Boost', '_cf_W')
compare('qtw', 'NJetBins', '1Boost', '_cf_Q')
compare('qtw', 'NJetBins', '2Boost', '_cf_Q')
compare('qtw', 'NJetBins', '1Boost', '_cf_T')
compare('qtw', 'NJetBins', '2Boost', '_cf_T')
compare('qtw', 'NJetBins', '1Boost', '_cf_W')
compare('qtw', 'NJetBins', '2Boost', '_cf_W')

compare('lt', 'MRR21vlBin',  '1Boost', '_cf_T')
compare('lt', 'MRR21vlBin',  '2Boost', '_cf_T')
compare('lt', 'MRR21vlBin',  '1Boost', '_cf_W')
compare('lt', 'MRR21vlBin',  '2Boost', '_cf_W')
compare('lt', 'NJetBins', '1Boost', '_cf_T')
compare('lt', 'NJetBins', '2Boost', '_cf_T')
compare('lt', 'NJetBins', '1Boost', '_cf_W')
compare('lt', 'NJetBins', '2Boost', '_cf_W')

compare('noniso', 'MRR2Bin',  '', '_cf_T')
compare('noniso', 'NJetBins', '', '_cf_T')
compare('noniso', 'MRR2Bin',  '', '_cf_W')
compare('noniso', 'NJetBins', '', '_cf_W')
#"""