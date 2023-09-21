import os
import ROOT
import multiprocessing
from plotting.basics import plot_ene
from analysis.utils import cutflow_hist, print_cutflow
from core.parallel import FuncWorker, run_pool


def hit_dr(hit_0, hit_1):
    _dx = hit_0.GetPos().X() - hit_1.GetPos().X()
    _dy = hit_0.GetPos().Y() - hit_1.GetPos().Y()
    _dz = hit_0.GetPos().Z() - hit_1.GetPos().Z()
    _dr = ROOT.TMath.Power(_dx, 2)
    _dr += ROOT.TMath.Power(_dy, 2)
    _dr += ROOT.TMath.Power(_dz, 2)
    _dr = ROOT.TMath.Sqrt(_dr)
    return _dr

def min_hit_dr(track_0, track_1):
    dr_ref = 1e6
    for hit_0 in track_0.GetAllHits():
        for hit_1 in track_1.GetAllHits():
            _dr = hit_dr(hit_0, hit_1)
            if _dr < dr_ref:
                dr_ref = _dr
    return dr_ref

def min_hit_dt(track_0, track_1):
    dt_ref = 1e6
    for hit_0 in track_0.GetAllHits():
        for hit_1 in track_1.GetAllHits():
            _dt = abs(hit_1.GetTime() - hit_0.GetTime())
            if _dt < dt_ref:
                dt_ref = _dt
    return dt_ref

def get_summary(evt, verbose=False):
    """
    """
    # Return the first summary object found with upstream time = 0
    for i_s, summary in enumerate(evt.summary):
        upstream_times = summary.GetUSTime()
        if len(upstream_times) != 0:
            # if len(upstream_times) > 1:
                # if verbose:
                #     print ('more than 1 upstream time for this summary!')
                #     for t in upstream_times:
                #         print(t)
                # raise ValueError('More than 1 Upstream time for a given summary!')
            _time_of_interest = upstream_times[0]
            if _time_of_interest == 0:
                return summary
    return None

def histbook(sample_name, cutflow_list):
    """
    """
    h_cutflow = cutflow_hist(cutflow_list, sample_name)
    _histbook = {
        'cutflow': h_cutflow
        }

    for cut in cutflow_list:
        _histbook['energy_' + cut] = ROOT.TH1F(
            'energy' + '_' + cut + '_'  + sample_name ,
            '',
            100, 0, 100)
    _histbook['delta_r_bhabha'] = ROOT.TH1F('delta_r_bhabha', '', 10, 0, 2)
    _histbook['delta_t_bhabha'] = ROOT.TH1F('delta_t_bhabha', '', 10, 0, 5)
    
    return _histbook

def flush_histbook_to_rootfile(histbook, name):
    rfile = ROOT.TFile(name, 'recreate')
    for k, v in histbook.items():
        v.Write(k)
    rfile.Close()
    
def fill_hist(rfile_name, histbook, evt_type, verbose=False, debug=False):
    rfile = ROOT.TFile(rfile_name)
    tree = rfile.Get('rec')


    counter_nevts_more_than1pattern = 0
    counter_nevts_more_than2pattern = 0
    counter_nevts_additionalpattern_is_electron = 0
    counter_nevts_additionalpattern_is_electron_with_onesinglehit = 0
    counter_nevts_additionalpattern_is_electron_veryclosetothepositrontrack = 0
    counter_nevts_additionalpattern_is_electron_with_onesinglehit_thatistracker = 0
    counter_nevts_additionalpattern_is_positron = 0
    
    for ievt, evt in enumerate(tree):

        histbook['cutflow'].Fill('all', 1)
        if ievt % 10000 == 0:
            print('{}: {} / {} ({:1.2f})%'.format(
                rfile_name,
                ievt, tree.GetEntries(), float(ievt)/tree.GetEntries() * 100))
        
        _info = getattr(evt, 'in')[0]
        if not _info.Has(evt_type):
            continue
        histbook['cutflow'].Fill('PIMCInfoFlag', 1)
        
        if debug:
            if ievt > 5000:
                break

            
        summary = get_summary(evt, verbose=verbose)
        if summary == None:
            continue
        histbook['cutflow'].Fill('findsummary', 1)
        histbook['energy_findsummary'].Fill(summary.GetCaloEnergy() + summary.GetAtarEnergy())

        upstream_time = summary.GetUSTime()[0]
        calo_time = summary.GetCaloTime()
        if abs(upstream_time - calo_time) < 5:
            continue
        histbook['cutflow'].Fill('upstream_calo_match', 1)
        histbook['energy_upstream_calo_match'].Fill(summary.GetCaloEnergy() + summary.GetAtarEnergy())
        
        pos = summary.GetTrackerPosition()
        if len(pos) == 0:
            continue
        histbook['cutflow'].Fill('found_tracker_position', 1)
        histbook['energy_found_tracker_position'].Fill(summary.GetCaloEnergy() + summary.GetAtarEnergy())

        if pos[0].Theta() > ROOT.TMath.Pi() / 3.:
            continue

        histbook['cutflow'].Fill('pass_theta_cut', 1)
        histbook['energy_pass_theta_cut'].Fill(summary.GetCaloEnergy() + summary.GetAtarEnergy())

        if verbose:
            atar_pattern = summary.GetPattern()
            if (len(atar_pattern) > 1):
                counter_nevts_more_than1pattern += 1
                if (len(atar_pattern) > 2):
                    counter_nevts_more_than2pattern += 1
                print (20 * '-')
                print ('event', ievt, 'has len(summary) = ', len(evt.summary), len(evt.cl)) 
                print('\t summary has N(patterns) =', len(atar_pattern))

                ref_positron_X, ref_positron_Y, ref_positron_Z, ref_time = -1000, -1000, -1000, -1000
                ref_track = None
                for pat in atar_pattern:
                    if len(pat.GetTrackIndices()) == 2:
                        _track_0 = evt.track[pat.GetTrackIndices()[0]]
                        _track_1 = evt.track[pat.GetTrackIndices()[1]]
                        if _track_0.GetPID() == 211 and _track_1.GetPID() == -11:
                            ref_track = _track_1
                            ref_positron_X = _track_1.GetStopPoint().X()
                            ref_positron_Y = _track_1.GetStopPoint().Y()
                            ref_positron_Z = _track_1.GetStopPoint().Z()
                            ref_positron_time = _track_1.GetStopTime()
                
                one_pattern_just_an_electron = False
                one_pattern_just_an_electron_with_one_hit = False
                one_pattern_just_an_electron_with_one_hit_that_is_tracker = False
                one_pattern_just_an_electron_very_close_to_one_of_the_positron_hit = False
                one_pattern_just_a_positron = False
                for pat in atar_pattern:
                    if len(pat.GetTrackIndices()) == 2:
                        _pid_0 = evt.track[pat.GetTrackIndices()[0]].GetPID()
                        _pid_1 = evt.track[pat.GetTrackIndices()[1]].GetPID()
                        if _pid_0 == 211 and _pid_1 == -11:
                            continue

                    print ('\t Track indices: ', pat.GetTrackIndices())
                    for idx in pat.GetTrackIndices():
                        trk = evt.track[idx]
                        if verbose:
                            print ('\t\t Track', idx, 'has PID = ', trk.GetPID())
                        if len(pat.GetTrackIndices()) == 1 and trk.GetPID() == 11:
                            this_track_X = trk.GetStartPoint().X()
                            this_track_Y = trk.GetStartPoint().Y()
                            this_track_Z = trk.GetStartPoint().Z()
                            this_track_time = trk.GetStartTime()
                            print ('\t\t Number of hits of this track =  ', len(trk.GetAllHits()))
                            print('\t\t        Ref Stop  X, Y, Z, time = ', ref_positron_X, ref_positron_Y, ref_positron_Z, ref_positron_time)
                            print('\t\t This track Start X, Y, Z, time = ', this_track_X, this_track_Y, this_track_Z, this_track_time)
                        
                    if len(pat.GetTrackIndices()) == 1:
                        track_0 = evt.track[pat.GetTrackIndices()[0]]
                        if track_0.GetPID() == 11:
                            one_pattern_just_an_electron = True
                            if len(track_0.GetAllHits()) == 1:
                                one_pattern_just_an_electron_with_one_hit = True
                                if track_0.GetAllHits()[0].GetType() == 't':
                                    one_pattern_just_an_electron_with_one_hit_that_is_tracker = True
                            else:
                                _dr = 1e6
                                _dt = 1e6
                                if ref_track != None:
                                    _dr = min_hit_dr(track_0, ref_track)
                                    _dt = min_hit_dt(track_0, ref_track)
                                    histbook['delta_r_bhabha'].Fill(_dr)
                                    histbook['delta_t_bhabha'].Fill(_dt)
                                if _dr < 0.120:
                                    one_pattern_just_an_electron_very_close_to_one_of_the_positron_hit = True
                        if track_0.GetPID() == -11:
                            one_pattern_just_a_positron = True

                if one_pattern_just_an_electron:
                    counter_nevts_additionalpattern_is_electron += 1
                    if one_pattern_just_an_electron_with_one_hit:
                        counter_nevts_additionalpattern_is_electron_with_onesinglehit += 1
                        if one_pattern_just_an_electron_with_one_hit_that_is_tracker:
                            counter_nevts_additionalpattern_is_electron_with_onesinglehit_thatistracker += 1
                    if one_pattern_just_an_electron_very_close_to_one_of_the_positron_hit:
                        counter_nevts_additionalpattern_is_electron_veryclosetothepositrontrack += 1
                if one_pattern_just_a_positron:
                    counter_nevts_additionalpattern_is_positron += 1
                            

    output_name = os.path.join(
        'cache',
        'hist_' + rfile_name.split('/')[-1])

    print('Events with >1 pattern:', counter_nevts_more_than1pattern)
    print('Events with >2 pattern:', counter_nevts_more_than2pattern)
    print('Events at least one single electron as pattern:', counter_nevts_additionalpattern_is_electron)
    print('Events with at least one single electron with a hit less than 0.120mm from a positron hit: ', counter_nevts_additionalpattern_is_electron_veryclosetothepositrontrack)
    print('Events at least one single electron as pattern with a single hit:', counter_nevts_additionalpattern_is_electron_with_onesinglehit)
    print('Events at least one single electron as pattern with a single hit that comes from the tracker:', counter_nevts_additionalpattern_is_electron_with_onesinglehit_thatistracker)
    print('Events at least one single positron as pattern:', counter_nevts_additionalpattern_is_positron)

    flush_histbook_to_rootfile(histbook, output_name)
    
            
def plot(pienu, mudif):
    pienu.SetLineColor(ROOT.kRed)
    mudif.SetLineColor(ROOT.kGreen)
    c = ROOT.TCanvas()
    pienu.Draw('HIST')
    mudif.Draw('sameHIST')
    c.Update()
    return c


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--debug', default=False, action='store_true')
    parser.add_argument('--verbose', default=False, action='store_true')
    args = parser.parse_args()

    cutflow_list = [
        'all',
        'PIMCInfoFlag',
        'findsummary',
        'upstream_calo_match',
        'found_tracker_position',
        'pass_theta_cut',
        ]

    filename_mudif = '/Users/quentin/pioneer/patricks_data/mudif_rec.root'
    
    
    filename_pienu = '/Users/quentin/pioneer/patricks_data/pienu_rec.root'
    hists_mudif = histbook('mudif', cutflow_list)
    hists_pienu = histbook('pienu', cutflow_list)

    if args.verbose:
        fill_hist(filename_pienu, hists_pienu, ROOT.PIMCInfo.kPienu, debug=args.debug, verbose=args.verbose)
    #     fill_hist(filename_mudif, hists_mudif, ROOT.PIMCInfo.kMudif, debug=args.debug, verbose=args.verbose)

    # else:
        
    #     jobs = [
    #         FuncWorker(fill_hist, filename_mudif, hists_mudif, ROOT.PIMCInfo.kMudif, debug=args.debug, verbose=args.verbose),
    #         FuncWorker(fill_hist, filename_pienu, hists_pienu, ROOT.PIMCInfo.kPienu, debug=args.debug, verbose=args.verbose),
    #     ]
    #     run_pool(jobs, n_jobs=multiprocessing.cpu_count() - 1)
        
    # from analysis.utils import get_scaling_factor
    # mudif_scaling = get_scaling_factor()
    # r_e_mu = 1.23e-4


    # histfile_mudif = ROOT.TFile.Open(os.path.join(
    #     'cache', 'hist_' + filename_mudif.split('/')[-1]))
    # histfile_pienu = ROOT.TFile.Open(os.path.join(
    #     'cache', 'hist_' + filename_pienu.split('/')[-1]))

    # cutflow_mudif = histfile_mudif.Get('cutflow')
    # cutflow_pienu = histfile_pienu.Get('cutflow')
    # mudif_total_entries = cutflow_mudif.GetBinContent(cutflow_mudif.GetXaxis().FindBin('all'))
    # pienu_total_entries = cutflow_pienu.GetBinContent(cutflow_pienu.GetXaxis().FindBin('all'))

    # keynames = [k.GetName() for k in histfile_pienu.GetListOfKeys()]
    # store_stuff = []
    # for k in keynames:
    #     if k in ('cutflow', 'energy_all', 'energy_upstream'):
    #         continue
    #     print (k)
    #     key_mudif = histfile_mudif.FindKey(k)
    #     key_pienu = histfile_pienu.FindKey(k)
    #     hist_mudif = key_mudif.ReadObj()
    #     hist_pienu = key_pienu.ReadObj()
    #     c, lg, plotables = plot_ene(
    #         hist_pienu,
    #         None,
    #         hist_mudif,
    #         max_events=1000)
    #     c.Update()
    #     store_stuff += [c, lg, plotables]
        
    # # for k in hists_pienu.keys():
    # #     if k in ('cutflow', 'energy_all', 'energy_upstream'):
    # #         continue
    # #     print(k)
    # #     hists_mudif[k].Scale(mudif_scaling / mudif_total_entries)
    # #     hists_pienu[k].Scale(r_e_mu / pienu_total_entries)

    # #     c, lg, plotables = plot_ene(
    # #         hists_pienu[k],
    # #         None,
    # #         hists_mudif[k],
    # #         max_events=1000)
    # #     c.Update()
    # #     store_stuff += [c, lg, plotables]
    #     # c = plot(hists_pienu[k], hists_mudif[k])
    

    
    # print_cutflow(cutflow_pienu)
    # print_cutflow(cutflow_mudif)
