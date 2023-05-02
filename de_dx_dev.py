import uuid
import numpy as np
import ROOT
import os
import array
from core.reco import *
from core.variables import *
from core.io import load, write
import ROOT

# These should become class members
# Fiducial Volume and Time cuts:

theta_limits = [np.pi/3., 2 * np.pi / 3.]


timeWin = 200 #half width

atarWin = 1 # ATAR search time window
caloWin = 5 # Calo search time window

trMergeTime = 5
beamRejection = 5 # reject 5 ns after beam arrival

# scale factors due to cutoffs
tau_muon = 2197. #ns

f_50ps = 1 - np.exp(-0.05 / tau_muon)
f_500ns = 1 - np.exp(-500 / tau_muon)
print(f" 50ps Scale Factor: {f_50ps}\n500ns Scale Factor: {f_500ns}")

r_e_mu = 1.23e-4

# Scale factors due to lifetimes
t = np.linspace(0, 3000000, 3000001)
tau_pi = 26.03
tau_mu = tau_muon
N_pi = 30 * np.exp(- t / tau_pi )

N_mu = -tau_mu * tau_pi / (-tau_mu + tau_pi) * np.exp(-t / tau_mu) - tau_mu * tau_pi / (tau_mu - tau_pi) * np.exp(-t / tau_pi)

# from matplotlib import pyplot as plot
# plot.plot(t[:3000], N_mu[:3000])
# plot.plot(t[:3000], N_pi[:3000])
# plot.title("Number of Muons")
# plot.xlabel("time (ns)")
# plot.ylabel("n (a.u., not to scale)")

mu_frac = sum(N_mu[beamRejection:2 * timeWin]) / sum(N_mu)
pi_frac = np.exp(-beamRejection / tau_pi) - np.exp(-2 * timeWin / tau_pi)
print("Integration Window (ns): ", t[beamRejection], t[2 * timeWin])
print("Muon Lifetime Correction (numeric): ", mu_frac)
print("Pion Lifetime Correction: ", pi_frac)

def cutflow_hist(cutflow, name):
    cutflow_bins = list(range(len(cutflow) + 1))
    print (cutflow_bins)
    hist = ROOT.TH1F(
        name + "cutflow", "",
        len(cutflow_bins) - 1,
        array.array('f', cutflow_bins))
    for i_label, label in enumerate(cutflow):
        print(i_label, label)
        hist.GetXaxis().SetBinLabel(i_label + 1, label)
    return hist
    
                        
                        

def update_cutflow(hist, bin_number, weight=1):
    old = hist.GetBinContent(bin_number)
    hist.SetBinContent(bin_number, old + weight)
    
def InRange(val, rng):
    if val > rng[0] and val < rng[1]:
        return True
    else:
        return False
    
        

def type2key(evType):
    if (evType & 6156 == 6152):
        # Pi !DIF TAR, Mu DIF TAR
        return "mudif"
            
    elif (evType & 6156 == 4104):
        # PI !DIF, TAR, Mu TAR, !DIF
        return "mudar"
    elif (evType & 13 == 9):
        # Pi -> e nu !DIF TAR
        return "pienu"
    else:
        #Everything else
        return "other"
        
    
maxEvents = 1e5

keys = [
    "etmudif", "etmudar", "etpienu", "etother",
    
    "edEmudif", "edEmudar", "edEpienu", "edEother",
    "edE12mudif", "edE12mudar", "edE12pienu", "edE12other",
    "eRmudif", "eRmudar", "eRpienu", "eRother",
    
    "edE1mudif", "edE1mudar", "edE1pienu", "edE1other",
    "eR1mudif", "eR1mudar", "eR1pienu", "eR1other",
    "edE2mudif", "edE2mudar", "edE2pienu", "edE2other",
    "eR2mudif", "eR2mudar", "eR2pienu", "eR2other",
    
    "pidEdx1", "mudE12", "edE12",
    
    "mudEdx", "edEdx","mudEdx1", "edEdx1", "mudEdx2", "edEdx2",
    "mudtdx", "edtdx","mudtdx1", "edtdx1", "mudtdx2", "edtdx2",
    "muRx", "eRx","muRx1", "eRx1", "muRx2", "eRx2",
]


def procbeam(chain, name, maxTime = 300, maxTCalo = 100):
    dPrint = 100
    nPrint = 100
    
    hist = dict()

    cutflow = [
        'all',
        'upstream',
        'tracker0',
        'tracker1',
        'beamTime',
        'fiducial',
        'atar0',
        'atarT',
        'caloE',
        ]
    
    hist['cutflow'] = cutflow_hist(cutflow, name)

    hist['ene_fivepix'] = ROOT.TH1F(uuid.uuid4().hex, '', 100, 0, 5)
    hist['ene_fivepix_pimudif_0_1'] = ROOT.TH1F(uuid.uuid4().hex, '', 100, 0, 5)
    hist['ene_fivepix_pimudif_2_5'] = ROOT.TH1F(uuid.uuid4().hex, '', 100, 0, 5)
    hist['ene_fivepix_pimudif_6_inf'] = ROOT.TH1F(uuid.uuid4().hex, '', 100, 0, 5)

    bins = [100, 0 , 2 * timeWin, 150, 0, 75 ]
    hist["etmudif"] = ROOT.TH2F(name + "etmudif", "Edep vs Time;t (ns);Edep (MeV)", *bins)
    hist["etmudar"] = ROOT.TH2F(name + "etmudar", "Edep vs Time;t (ns);Edep (MeV)", *bins)
    hist["etpienu"] = ROOT.TH2F(name + "etpienu", "Edep vs Time;t (ns);Edep (MeV)", *bins)
    hist["etother"] = ROOT.TH2F(name + "etother", "Edep vs Time;t (ns);Edep (MeV)", *bins)
    
    bins2 = [100, 0, 2, 150, 0, 75]
    hist["edEmudif"] = ROOT.TH2F(name + "edEmudif", "Edep vs dE;dE (MeV);Edep (MeV)", *bins2)
    hist["edEmudar"] = ROOT.TH2F(name + "edEmudar", "Edep vs dE;dE (MeV);Edep (MeV)", *bins2)
    hist["edEpienu"] = ROOT.TH2F(name + "edEpienu", "Edep vs dE;dE (MeV);Edep (MeV)", *bins2)
    hist["edEother"] = ROOT.TH2F(name + "edEother", "Edep vs dE;dE (MeV);Edep (MeV)", *bins2)
    hist["edE1mudif"] = ROOT.TH2F(name + "edE1mudif", "Edep vs dE;dE (MeV);Edep (MeV)", *bins2)
    hist["edE1mudar"] = ROOT.TH2F(name + "edE1mudar", "Edep vs dE;dE (MeV);Edep (MeV)", *bins2)
    hist["edE1pienu"] = ROOT.TH2F(name + "edE1pienu", "Edep vs dE;dE (MeV);Edep (MeV)", *bins2)
    hist["edE1other"] = ROOT.TH2F(name + "edE1other", "Edep vs dE;dE (MeV);Edep (MeV)", *bins2)
    hist["edE2mudif"] = ROOT.TH2F(name + "edE2mudif", "Edep vs dE;dE (MeV);Edep (MeV)", *bins2)
    hist["edE2mudar"] = ROOT.TH2F(name + "edE2mudar", "Edep vs dE;dE (MeV);Edep (MeV)", *bins2)
    hist["edE2pienu"] = ROOT.TH2F(name + "edE2pienu", "Edep vs dE;dE (MeV);Edep (MeV)", *bins2)
    hist["edE2other"] = ROOT.TH2F(name + "edE2other", "Edep vs dE;dE (MeV);Edep (MeV)", *bins2)
    hist["edE12mudif"] = ROOT.TH2F(name + "edE12mudif", "Edep vs dE;dE (MeV);Edep (MeV)", *bins2)
    hist["edE12mudar"] = ROOT.TH2F(name + "edE12mudar", "Edep vs dE;dE (MeV);Edep (MeV)", *bins2)
    hist["edE12pienu"] = ROOT.TH2F(name + "edE12pienu", "Edep vs dE;dE (MeV);Edep (MeV)", *bins2)
    hist["edE12other"] = ROOT.TH2F(name + "edE12other", "Edep vs dE;dE (MeV);Edep (MeV)", *bins2)
    
    bins3 = [200, 0, 6, 150, 0, 75]
    hist["eRmudif"] = ROOT.TH2F(name + "eRmudif", "Edep vs dE/dx;dE/dx (MeV/mm);Edep (MeV)", *bins3)
    hist["eRmudar"] = ROOT.TH2F(name + "eRmudar", "Edep vs dE/dx;dE/dx (MeV/mm);Edep (MeV)", *bins3)
    hist["eRpienu"] = ROOT.TH2F(name + "eRpienu", "Edep vs dE/dx;dE/dx (MeV/mm);Edep (MeV)", *bins3)
    hist["eRother"] = ROOT.TH2F(name + "eRother", "Edep vs dE/dx;dE/dx (MeV/mm);Edep (MeV)", *bins3)
    hist["eR1mudif"] = ROOT.TH2F(name + "eR1mudif", "Edep vs dE/dx;dE/dx (MeV/mm);Edep (MeV)", *bins3)
    hist["eR1mudar"] = ROOT.TH2F(name + "eR1mudar", "Edep vs dE/dx;dE/dx (MeV/mm);Edep (MeV)", *bins3)
    hist["eR1pienu"] = ROOT.TH2F(name + "eR1pienu", "Edep vs dE/dx;dE/dx (MeV/mm);Edep (MeV)", *bins3)
    hist["eR1other"] = ROOT.TH2F(name + "eR1other", "Edep vs dE/dx;dE/dx (MeV/mm);Edep (MeV)", *bins3)
    hist["eR2mudif"] = ROOT.TH2F(name + "eR2mudif", "Edep vs dE/dx;dE/dx (MeV/mm);Edep (MeV)", *bins3)
    hist["eR2mudar"] = ROOT.TH2F(name + "eR2mudar", "Edep vs dE/dx;dE/dx (MeV/mm);Edep (MeV)", *bins3)
    hist["eR2pienu"] = ROOT.TH2F(name + "eR2pienu", "Edep vs dE/dx;dE/dx (MeV/mm);Edep (MeV)", *bins3)
    hist["eR2other"] = ROOT.TH2F(name + "eR2other", "Edep vs dE/dx;dE/dx (MeV/mm);Edep (MeV)", *bins3)
    
    piBins = [200, 0, 0.5, 200, 0, 10]
    hist["pidEdx1"] = ROOT.TH2F(name + "pidEdx1", "pidEdx;dx (mm);dE (MeV)", *piBins)
    
    bins4 = [200, 0, 2, 100, 0, 1.7]
    hist["mudEdx"] = ROOT.TH2F(name + "mudEdx", "mudEdx;dx (mm);dE (MeV)", *bins4)
    hist["edEdx"] = ROOT.TH2F(name + "edEdx", "edEdx;dx (mm);dE (MeV)", *bins4)  
    hist["mudEdx1"] = ROOT.TH2F(name + "mudEdx1", "mudEdx;dx (mm);dE (MeV)", *bins4)
    hist["edEdx1"] = ROOT.TH2F(name + "edEdx1", "edEdx;dx (mm);dE (MeV)", *bins4)  
    hist["mudEdx2"] = ROOT.TH2F(name + "mudEdx2", "mudEdx;dx (mm);dE (MeV)", *bins4)
    hist["edEdx2"] = ROOT.TH2F(name + "edEdx2", "edEdx;dx (mm);dE (MeV)", *bins4)
    
    bins5 = [100, 0, 10, 170, 0, 1.7]
    hist["mudtdx"] = ROOT.TH2F(name + "mudtdx", "mudtdx;dt (ps);dx (mm)", *bins5)
    hist["edtdx"] = ROOT.TH2F(name + "edtdx", "edtdx;dt (ps);dx (mm)", *bins5)  
    hist["mudtdx1"] = ROOT.TH2F(name + "mudtdx1", "mudtdx;dt (ps);dx (mm)", *bins5)
    hist["edtdx1"] = ROOT.TH2F(name + "edtdx1", "edtdx;dt (ps);dx (mm)", *bins5)  
    hist["mudtdx2"] = ROOT.TH2F(name + "mudtdx2", "mudtdx;dt (ps);dx (mm)", *bins5)
    hist["edtdx2"] = ROOT.TH2F(name + "edtdx2", "edtdx;dt (ps);dx (mm)", *bins5)
    
    
    bins6 = [170, 0, 1.7, 170, 0, 6]
    hist["muRx"] = ROOT.TH2F(name + "muRx", "muRx;dx (mm);dE/dx (MeV/mm)", *bins6)
    hist["eRx"] = ROOT.TH2F(name + "eRx", "eRx;dx (mm);dE/dx (MeV/mm)", *bins6)  
    hist["muRx1"] = ROOT.TH2F(name + "muRx1", "muRx;dx (mm);dE/dx (MeV/mm)", *bins6)
    hist["eRx1"] = ROOT.TH2F(name + "eRx1", "eRx;dx (mm);dE/dx (MeV/mm)", *bins6)  
    hist["muRx2"] = ROOT.TH2F(name + "muRx2", "muRx;dx (mm);dE/dx (MeV/mm)", *bins6)
    hist["eRx2"] = ROOT.TH2F(name + "eRx2", "eRx;dx (mm);dE/dx (MeV/mm)", *bins6)
    
    bins7 = [200, 0, 2, 200, 0, 2]
    hist["mudE12"] = ROOT.TH2F(name + "mudE12", "mudE12;dE1 (MeV);dE2 (MeV)", *bins7)
    hist["edE12"] = ROOT.TH2F(name + "edE12", "edE12;dE1 (MeV);dE2 (MeV)", *bins7)
    
    for evID, event in enumerate(chain):
        if evID > maxEvents:
            break
        if (evID == nPrint):
            print("processing ", evID)
            nPrint += dPrint
            if nPrint == 10 * dPrint:
                dPrint *= 10
            
        evType = event.info.GetType()
        eventKey = type2key(evType)
        if name == '50ps':
            if eventKey != 'mudif':
                continue
        elif name == 'pienu':
            if eventKey != 'pienu':
                continue
        else:
            raise ValueError
        
        update_cutflow(hist['cutflow'], 1)

        #first check upstream beam entrance
        us = [u for u in event.upstream if u.GetUpstreamID() == 600000 ]
        if len(us) == 0:
            #no upstream
            continue

        update_cutflow(hist['cutflow'], 2)

        T0 = us[0].GetTime()[0]
        
        #find a tracker hit in the 400 ns after
        tracker = FindTracker(event.tracker, time = T0 + timeWin, timeLimit = timeWin)
        if len(tracker) == 0:
            # no tracker hits
            continue
        update_cutflow(hist['cutflow'], 3)
            
        TT_q = -trMergeTime
        merged_tracker = []
        for tr in tracker:
            if tr.GetTime() - T0 < beamRejection:
                # beam rejection window
                continue
            if tr.GetTime() < TT_q + trMergeTime:
                # simultaneous with previous tracker hit
                continue
            merged_tracker += [tr]
            TT_q = tr.GetTime()
        if len(merged_tracker) == 0:
            continue
        update_cutflow(hist['cutflow'], 4)

        if len(merged_tracker) > 1:
            print ('evt \t x0 \t Y0 \t Z0 \t Time \t Edep')
            for track in merged_tracker:
                print('{0:d}, {1:1.3f}, {2:1.3f}, {3:1.3f}, {4:1.3f}, {5:1.3f}'.format(
                    evID, track.GetX0(), track.GetY0(), track.GetZ0(), track.GetTime(), track.GetEdep()))
            continue

        update_cutflow(hist['cutflow'], 5)

        tracker_hit = merged_tracker[0]

        
        # Fiducial Volume cut
        if not InRange(tracker_hit.GetPostPos().Theta(), theta_limits):
            continue
        update_cutflow(hist['cutflow'], 6)
        
        # Find T0 Atar hits
        atar0 = FindAtar(event.atar, time = T0 + atarWin, timeLimit= atarWin)
        if len(atar0) == 0:
            # no arrival in ATAR is seen
            continue
        update_cutflow(hist['cutflow'], 7)

        TT = tracker_hit.GetTime()
        atarT = FindAtar(event.atar, time = TT - atarWin, timeLimit = atarWin)
        if len(atarT) == 0:
            continue
        update_cutflow(hist['cutflow'], 8)
        
        caloT = FindCalo(event.calo, time = TT, timeLimit = caloWin)
        if len(caloT['edep']) == 0:
            continue
        update_cutflow(hist['cutflow'], 9)
        
        usT = FindUpstream(us[0], time = TT - atarWin, timeLimit = atarWin)
        
        etUs = sum(usT['edep'])
        etAtar = sum([a.GetEdep() for a in atarT])
        etCalo = sum(caloT['edep'])
        
        esum = etUs + etAtar + etCalo
        merged_atarT = build_merged_atar(atarT)
        ene_five = get_energy_n_nonpion_pixels(merged_atarT)
        hist['ene_fivepix'].Fill(ene_five)

        true_muon_pixel_ids = get_true_muonpix_ids(event.atar)
        if len(true_muon_pixel_ids) < 2:
            hist['ene_fivepix_pimudif_0_1'].Fill(ene_five)
        elif len(true_muon_pixel_ids) < 6:
            hist['ene_fivepix_pimudif_2_5'].Fill(ene_five)
        else:
            hist['ene_fivepix_pimudif_6_inf'].Fill(ene_five)

        hist["et" + eventKey].Fill(TT, esum)
    return hist

def get(chain, name, forceRead = False, writeHist = True):
    hist = None
    if not forceRead:
        hist = load(name)
        
    if not hist:
        hist = procbeam(chain, name)
        if writeHist:
            write(name, hist)
    return hist





if __name__ == '__main__':

    beam50ps = ROOT.TChain("sim", "sim")
    beam50ps.Add("/Users/quentin/decay_new/sim/beam_50ps*.root")
    print(beam50ps.GetEntries())

    beamPienu = ROOT.TChain("sim", "sim")
    beamPienu.Add("/Users/quentin/decay_new/sim/beam_pienu*.root")
    print(beamPienu.GetEntries())

    hist50ps = get(beam50ps, "50ps")
    histPienu = get(beamPienu, "pienu")

    for i in range(hist50ps['cutflow'].GetNbinsX()):
        print (hist50ps['cutflow'].GetXaxis().GetBinLabel(i+1), hist50ps['cutflow'].GetBinContent(i+1))
    
    
    for k, v in hist50ps.items():
        v.Scale(f_50ps)

    for j, v in histPienu.items():
        v.Scale(r_e_mu)
 
    from plotting.basics import plot_stack
    c, lg = plot_stack(
        hist50ps['ene_fivepix_pimudif_0_1'],
        hist50ps['ene_fivepix_pimudif_2_5'],
        hist50ps['ene_fivepix_pimudif_6_inf'],
        hist50ps['ene_fivepix'],
        histPienu['ene_fivepix'])

    c.Update()
    c
    
