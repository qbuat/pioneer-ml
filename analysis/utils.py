import numpy as np
import ROOT
import array
from tabulate import tabulate

def get_scaling_factor():
        # scale factors due to cutoffs
    tau_muon = 2197. #ns
    
    f_50ps = 1 - np.exp(-0.05 / tau_muon)
    f_500ns = 1 - np.exp(-500 / tau_muon)
    print(f" 50ps Scale Factor: {f_50ps}\n500ns Scale Factor: {f_500ns}")
    return f_50ps

def print_cutflow(hist):
    cutflow_yields = []
    for i in range(hist.GetNbinsX()):
        cutflow_yields += [[
            i + 1,
            hist.GetXaxis().GetBinLabel(i+1),
            hist.GetBinContent(i+1)
        ]]
    table = tabulate(cutflow_yields, headers=['Step', 'Name', 'Entries'])
    print('Cutflow for {}'.format(hist.GetName()))
    print (table)


def cutflow_hist(cutflow_list, name):
    cutflow_bins = list(range(len(cutflow_list) + 1))
    hist = ROOT.TH1F(
        name + "_cutflow", "",
        len(cutflow_bins) - 1,
        array.array('f', cutflow_bins))
    for i_label, label in enumerate(cutflow_list):
        hist.GetXaxis().SetBinLabel(i_label + 1, label)
    return hist
