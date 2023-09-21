import uuid
import numpy as np
import ROOT
import os
import array
from core.reco import *
from core.variables import *
from core.io import load, write
from analysis.variables import VARIABLES, VARIABLES_2D
import ROOT
import tabulate

ROOT.gROOT.SetBatch(True)
from atar.utils import geometry
atar_geo = geometry('/Users/quentin/decay_new/default.gdml')

# These should become class members
# Fiducial Volume and Time cuts:

# theta_limits = [-1, np.pi /3.]
theta_limits = [np.pi/3., 2 * np.pi / 3.]
# theta_limits = [0., np.pi/3.]


timeWin = 200 #half width

atarWin = 1 # ATAR search time window
caloWin = 5 # Calo search time window

trMergeTime = 5
beamRejection = 5 # reject 5 ns after beam arrival

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
    'dep',
    'reject_mudar',
    'ene_fivepix',
    'de_dx',
    'cos_theta',
]
def get_scaling_factors():
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


    mu_frac = sum(N_mu[beamRejection:2 * timeWin]) / sum(N_mu)
    pi_frac = np.exp(-beamRejection / tau_pi) - np.exp(-2 * timeWin / tau_pi)
    print("Integration Window (ns): ", t[beamRejection], t[2 * timeWin])
    print("Muon Lifetime Correction (numeric): ", mu_frac)
    print("Pion Lifetime Correction: ", pi_frac)

    return r_e_mu, f_50ps, f_500ns
    

def print_cutflow(hists, key='cutflow'):
    cutflow_yields = []
    for i in range(hists['cutflow'].GetNbinsX()):
        cutflow_yields += [[
            i + 1,
            hists['cutflow'].GetXaxis().GetBinLabel(i+1),
            hists['cutflow'].GetBinContent(i+1)
        ]]
    table = tabulate.tabulate(cutflow_yields, headers=['Step', 'Name', 'Entries'])
    print (table)
        




def cutflow_hist(cutflow_list, name):
    cutflow_bins = list(range(len(cutflow_list) + 1))
    hist = ROOT.TH1F(
        name + "cutflow", "",
        len(cutflow_bins) - 1,
        array.array('f', cutflow_bins))
    for i_label, label in enumerate(cutflow_list):
        hist.GetXaxis().SetBinLabel(i_label + 1, label)
    return hist
    
                        
                        

def update_cutflow(hists, bin_number, step_name='', ene=0, weight=1):
    old = hists['cutflow'].GetBinContent(bin_number)
    hists['cutflow'].SetBinContent(bin_number, old + weight)
    hists['total_energy_' + step_name].Fill(ene, weight)

    
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


def procbeam(chain, name, maxTime = 300, maxTCalo = 100, maxEvents=-1):
    dPrint = 100
    nPrint = 100
    
    hist = dict()
    hist_2d = dict()
    
    hist['cutflow'] = cutflow_hist(cutflow, name)
    for step in cutflow:
        hist['total_energy_' + step] = ROOT.TH1F(
            uuid.uuid4().hex, '', *VARIABLES['total_energy']['binning'])
        hist['total_energy_' + step].GetXaxis().SetTitle('Energy [MeV]')
        
    # 1D vars
    for _, var in VARIABLES.items():
        for _ext in [
                '',
                '_pimudif_0_1',
                '_pimudif_2_5',
                '_pimudif_6_inf']:
            hist[var['name'] + _ext] = ROOT.TH1F(name + var['name'] + _ext, '', *var['binning'])
            hist[var['name'] + _ext].GetXaxis().SetTitle(var['title'])
            if 'unit' in var.keys():
                hist[var['name'] + _ext].GetXaxis().SetTitle(
                    hist[var['name'] + _ext].GetXaxis().GetTitle() + ' [{}]'.format(
                        var['unit']))
        
    # 2D vars:
    for _, var in VARIABLES_2D.items():
        for _ext in [
                '',
                '_pimudif_0_1',
                '_pimudif_2_5',
                '_pimudif_6_inf']:
            hist_2d[var['name'] + _ext] = ROOT.TH2F(
                name + var['name'] + _ext,
                '',
                *var['binning'])
            hist_2d[var['name'] + _ext].GetXaxis().SetTitle(var['xlabel'])
            hist_2d[var['name'] + _ext].GetYaxis().SetTitle(var['ylabel'])
            if 'xunit' in var.keys():
                hist_2d[var['name'] + _ext].GetXaxis().SetTitle(
                    hist_2d[var['name'] + _ext].GetXaxis().GetTitle() + ' [{}]'.format(
                        var['xunit']))
                
            if 'yunit' in var.keys():
                hist_2d[var['name'] + _ext].GetYaxis().SetTitle(
                    hist_2d[var['name'] + _ext].GetYaxis().GetTitle() + ' [{}]'.format(
                        var['yunit']))

    n_tot = 0
    n_last_012 = 0
    n_previous_to_last_012 = 0
    for evID, event in enumerate(chain):
        if maxEvents > 0 and evID > maxEvents:
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
        elif name == '500ns':
            if eventKey != 'mudar':
                continue
        else:
            raise ValueError

        de_last, de_previous_to_last, dr_last, dr_previous_to_last, theta_last, theta_previous_to_last = get_truth_pion_travel(event.atar, verbose=False)

        dx_last_estimated = get_dx_from_de(de_last)
        dx_last_plus_previous_to_last_estimated = get_dx_from_de(
            de_last + de_previous_to_last, last_plus_prev_to_last=True)
        dx_previous_to_last = max(dx_last_plus_previous_to_last_estimated - dx_last_estimated, 0)
        n_tot += 1

        hist['delta_dx_last'].Fill(dx_last_estimated - dr_last)
        hist['delta_dx_previous_to_last'].Fill(dx_previous_to_last_estimated - dr_previous_to_last)
        hist_2d['pion_ene_last_travel'].Fill(de_last, dr_last)
        hist_2d['pion_ene_last_and_prev_travel'].Fill(de_last + de_previous_to_last, dr_last + dr_previous_to_last)
        # hist_2d['pion_ene_previous_to_last_travel'].Fill(de_last + de_previous_to_last, dr_last + dr_previous_to_last)
        hist_2d['pion_ene_previous_to_last_travel'].Fill(de_previous_to_last, dr_previous_to_last)
        hist_2d['pion_ene_previous_to_last_angle'].Fill(de_previous_to_last, theta_previous_to_last)
        # if 0.120 < dr_previous_to_last < 0.122:
        hist_2d['pion_ene_last_previous_to_last'].Fill(de_last, de_previous_to_last)
        hist_2d['pion_dr_last_previous_to_last'].Fill(dr_last, dr_previous_to_last)
        hist_2d['pion_theta_last_previous_to_last'].Fill(theta_last, theta_previous_to_last)
        if dr_last > 0.12:
            n_last_012 += 1
        if dr_previous_to_last > 0.12:
            n_previous_to_last_012 += 1
        
        update_cutflow(hist, 1, step_name=cutflow[0])

        #first check upstream beam entrance
        us = [u for u in event.upstream if u.GetUpstreamID() == 600000 ]
        if len(us) == 0:
            #no upstream
            continue

        update_cutflow(hist, 2, step_name=cutflow[1])

        T0 = us[0].GetTime()[0]
        
        #find a tracker hit in the 400 ns after
        tracker = FindTracker(event.tracker, time = T0 + timeWin, timeLimit = timeWin)
        if len(tracker) == 0:
            # no tracker hits
            continue
        update_cutflow(hist, 3, step_name=cutflow[2])

        
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
        update_cutflow(hist, 4, step_name=cutflow[3])

        if len(merged_tracker) > 1:
            print ('evt \t x0 \t Y0 \t Z0 \t Time \t Edep')
            for track in merged_tracker:
                print('{0:d}, {1:1.3f}, {2:1.3f}, {3:1.3f}, {4:1.3f}, {5:1.3f}'.format(
                    evID, track.GetX0(), track.GetY0(), track.GetZ0(), track.GetTime(), track.GetEdep()))
            continue
        update_cutflow(hist, 5, step_name=cutflow[4])

        tracker_hit = merged_tracker[0]

        
        # Fiducial Volume cut
        if not InRange(tracker_hit.GetPostPos().Theta(), theta_limits):
            continue
        update_cutflow(hist, 6, step_name=cutflow[5])
        
        # Find T0 Atar hits
        atar0 = FindAtar(event.atar, time = T0 + atarWin, timeLimit= atarWin)
        if len(atar0) == 0:
            # no arrival in ATAR is seen
            continue
        setattr(event, 'atar0', atar0)
        update_cutflow(hist, 7, step_name=cutflow[6])

        # print()
        # print(atar0[-1].GetX1(), atar0[-1].GetY1())
        
        TT = tracker_hit.GetTime()
        atarT = FindAtar(event.atar, time = TT - atarWin, timeLimit = atarWin)
        if len(atarT) == 0:
            continue
        setattr(event, 'atarT', atarT)
        update_cutflow(hist, 8, step_name=cutflow[7])
        
        caloT = FindCalo(event.calo, time = TT, timeLimit = caloWin)
        if len(caloT['edep']) == 0:
            continue
        update_cutflow(hist, 9, step_name=cutflow[8])
        
        usT = FindUpstream(us[0], time = TT - atarWin, timeLimit = atarWin)
        
        etUs = sum(usT['edep'])
        etAtar = sum([a.GetEdep() for a in atarT])
        etCalo = sum(caloT['edep'])
        
        esum = etUs + etAtar + etCalo
        merged_atarT = build_merged_atar(atarT)
        setattr(event, 'merged_atarT', merged_atarT)


        _, travel_pion_last_pixel = get_last_hit_pion_depth(atar0, atar_geo)
        ene_five = get_energy_n_nonpion_pixels(merged_atarT)

        
        
        # if travel_pion_last_pixel > 0.120 / 2.:
        #     continue
        update_cutflow(hist, 10, step_name=cutflow[9], ene=esum)

        _dE_dx, cos_theta, dz_0 = dE_dx(usT, atar0, atarT, tracker_hit, atar_geo)
        dz_0_true, dx_0_true, dE_0 = get_true_muon_first_pixel_dz(event.atar)

        # print(get_true_distance_pion_muon(event.atar))
        
        if dx_0_true != 0:
            hist['true_dE_dx'].Fill(dE_0 / dx_0_true)

        if dz_0_true != 0:
            hist['true_dE_dx_zonly'].Fill(dE_0 / dz_0_true)

        atar_mudar = FindAtar_mudar(event.atar, atar0, atarT)
        if len(atar_mudar) != 0:
            continue
        update_cutflow(hist, 11, step_name=cutflow[10], ene=esum)


        _, _, z_pion_reco, pixid_pion_reco = get_pion_stopping_point(atar0, atar_geo)
        _, _, z_pion_truth, pixid_pion_truth = get_true_pion_stopping_position(event.atar, atar_geo)
        delta_z = z_pion_reco - z_pion_truth
        delta_id = pixid_pion_reco - pixid_pion_truth
        delta_dz = dz_0 - dz_0_true
        
        hist['ene_dE_dx'].Fill(_dE_dx)
        hist['ene_fivepix'].Fill(ene_five)
        hist['delta_z_pion'].Fill(delta_z)
        hist['delta_id_pion'].Fill(delta_id)
        hist['total_energy'].Fill(esum)
        hist['delta_dz_first_muon'].Fill(delta_dz)
        hist['truth_dz_first_muon'].Fill(dz_0_true)
        hist['reco_dz_first_muon'].Fill(dz_0)
        
        true_muon_pixel_ids = get_true_muonpix_ids(event.atar)
        if len(true_muon_pixel_ids) < 2:
            hist['ene_fivepix_pimudif_0_1'].Fill(ene_five)
            hist['ene_dE_dx_pimudif_0_1'].Fill(_dE_dx)
            hist['delta_z_pion_pimudif_0_1'].Fill(delta_z)
            hist['delta_id_pion_pimudif_0_1'].Fill(delta_id)
            hist['delta_dz_first_muon_pimudif_0_1'].Fill(delta_dz)
            hist['truth_dz_first_muon_pimudif_0_1'].Fill(dz_0_true)
            hist['reco_dz_first_muon_pimudif_0_1'].Fill(dz_0)
        elif len(true_muon_pixel_ids) < 6:
            hist['ene_fivepix_pimudif_2_5'].Fill(ene_five)
            hist['ene_dE_dx_pimudif_2_5'].Fill(_dE_dx)
            hist['delta_z_pion_pimudif_2_5'].Fill(delta_z)
            hist['delta_id_pion_pimudif_2_5'].Fill(delta_id)
            hist['delta_dz_first_muon_pimudif_2_5'].Fill(delta_dz)
            hist['truth_dz_first_muon_pimudif_2_5'].Fill(dz_0_true)
            hist['reco_dz_first_muon_pimudif_2_5'].Fill(dz_0)
        else:
            hist['ene_fivepix_pimudif_6_inf'].Fill(ene_five)
            hist['ene_dE_dx_pimudif_6_inf'].Fill(_dE_dx)
            hist['delta_z_pion_pimudif_6_inf'].Fill(delta_z)
            hist['delta_id_pion_pimudif_6_inf'].Fill(delta_id)
            hist['delta_dz_first_muon_pimudif_6_inf'].Fill(delta_dz)
            hist['truth_dz_first_muon_pimudif_6_inf'].Fill(dz_0_true)
            hist['reco_dz_first_muon_pimudif_6_inf'].Fill(dz_0)

        if ene_five > 0.5: # MeV
            continue
        update_cutflow(hist, 12, step_name=cutflow[11], ene=esum)

        if _dE_dx > 2.: #Mev/mm
            continue
        update_cutflow(hist, 13, step_name=cutflow[12], ene=esum)
            
        if cos_theta < 0:
            continue
        update_cutflow(hist, 14, step_name=cutflow[13], ene=esum)

    print (n_tot, n_last_012, n_previous_to_last_012)
    return hist, hist_2d

def get(
        chain,
        name, forceRead = True, writeHist = True, maxEvents=-1):
    hist = None
    # if not forceRead:
    #     hist = load(name)
        
    if not hist:
        hist, hist_2d = procbeam(chain, name, maxEvents=maxEvents)
        # if writeHist:
        #     write(name, hist)
    return hist, hist_2d





if __name__ == '__main__':

    beam50ps = ROOT.TChain("sim", "sim")
    beam50ps.Add("/Users/quentin/decay_new/sim/beam_50ps*.root")
    print(beam50ps.GetEntries())

    #beam500ns = ROOT.TChain("sim", "sim")
    #beam500ns.Add("/Users/quentin/decay_new/sim/beam_500ns*.root")
    #print(beam500ns.GetEntries())

    
    beamPienu = ROOT.TChain("sim", "sim")
    beamPienu.Add("/Users/quentin/decay_new/sim/beam_pienu*.root")
    print(beamPienu.GetEntries())

    maxEvents = 1e6

    r_e_mu, f_50ps, f_500ns = get_scaling_factors()
    
    hist50ps, hist_2d_50ps = get(beam50ps, "50ps", maxEvents=maxEvents)
    print(20 * '-')
    print('pi-mu-e 50ps')
    print_cutflow(hist50ps)
    for k, v in hist50ps.items():
        v.Scale(f_50ps)

    # hist500ns = get(beam500ns, "500ns", maxEvents=maxEvents)
    # print(20 * '-')
    # print('pi-mu-e 500ns')
    # print_cutflow(hist500ns)
    # for k, v in hist500ns.items():
    #     v.Scale(f_500ns)

    # histPienu, hist2d_Pienu = get(beamPienu, "pienu", maxEvents=maxEvents)
    # print(20 * '-')
    # print('pi-e ')
    # print_cutflow(histPienu)
    # for j, v in histPienu.items():
    #     v.Scale(r_e_mu)



    # ====== plotting
 
    from plotting.basics import *
    
    # for _, var in VARIABLES.items():
    #     c, lg = plot_stack(
    #         hist50ps[var['name'] + '_pimudif_0_1'],
    #         hist50ps[var['name'] + '_pimudif_2_5'],
    #         hist50ps[var['name'] + '_pimudif_6_inf'],
    #         hist50ps[var['name']],
    #         histPienu[var['name']],
    #         max_events=maxEvents)
    #     c.Update()
    #     c.SaveAs('plots/{}.pdf'.format(
    #         var['name']))
    #     c.Close()
        
    #     c, lg, plotables = plot_eff(
    #         hist50ps[var['name'] + '_pimudif_0_1'],
    #         hist50ps[var['name'] + '_pimudif_2_5'],
    #         hist50ps[var['name'] + '_pimudif_6_inf'],
    #         hist50ps[var['name']],
    #         histPienu[var['name']])
    #     c.Update()
    #     c.SaveAs('plots/eff_{}.pdf'.format(
    #         var['name']))
    #     c.Close()

    # for i_cut, cut in enumerate(cutflow):
    #     if i_cut < 9:
    #         continue

    #     c_ene, lg, p = plot_ene(
    #         histPienu['total_energy_' + cut],
    #         None,
    #         # hist500ns['total_energy_' + cut],
    #         hist50ps['total_energy_' + cut],
    #         max_events=maxEvents,
    #         norm_strategy='branching_ratio')
    #     c_ene.Update()
    #     c_ene.SaveAs('plots/pr_plot_branching_ratio_{}.pdf'.format(cut))
    #     c_ene.Close()


    # c_ene, lg, p = plot_ene(
    #     histPienu['total_energy'],
    #     None,
    #     # hist500ns['total_energy'],
    #     hist50ps['total_energy'],
    #     max_events=maxEvents,
    #     norm_strategy='branching_ratio')
    # c_ene.Update()
    # c_ene.SaveAs('plots/pr_plot_branching_ratio.pdf')
    # c_ene.Close()

    for k, v in VARIABLES_2D.items():
        c_2d = ROOT.TCanvas()
        hist_2d_50ps[k].Draw('colz')
        c_2d.Update()
        if k in ('pion_ene_last_and_prev_travel', 'pion_ene_last_travel'):
            fit = ROOT.TF1('fit', 'TMath::Power(x/[0], [1])', 0, 3)
            fit.SetParameter(0, 5.5)
            fit.SetParameter(1, 1.7)
            hist_2d_50ps[k].Fit(fit, 'RSV')
            fit.Draw('same')
        c_2d.Update()
        c_2d.SaveAs('plots/{}.pdf'.format(v['name']))
