import ROOT
import tabulate
import array

def plot_stack(
        pimudif_0_2,
        pimudif_3_5,
        pimudif_6_inf,
        pimudif,
        pienu,
        max_events=-1,
        style='HIST'):


    pimudif_0_2.SetLineColor(ROOT.kGreen + 2)
    pimudif_3_5.SetLineColor(ROOT.kBlue + 2)
    pimudif_6_inf.SetLineColor(ROOT.kAzure + 2)
    pienu.SetLineColor(ROOT.kGray + 2)

    pimudif_0_2.SetFillColorAlpha(ROOT.kGreen + 2, 0.05)
    pimudif_3_5.SetFillColorAlpha(ROOT.kBlue + 2, 0.05)
    pimudif_6_inf.SetFillColorAlpha(ROOT.kAzure + 2, 0.05)

    pienu.SetLineColor(ROOT.kRed)
    pienu.SetLineWidth(2)


    plotables = [
        pienu,
        pimudif,
        pimudif_0_2,
        pimudif_3_5,
        pimudif_6_inf
    ]
    for p in plotables:
        # if max_events != -1:
        #     p.Scale(1. / max_events)
        #     p.SetMinimum(1e-11)
        #     p.GetYaxis().SetTitle('Branching Ratio')
        # else:
        p.SetMinimum(1e-6)
        p.GetYaxis().SetTitle('Branching Ratio #times ({})'.format(max_events))
            
    plotables = sorted(
        plotables,
        key=lambda h: h.GetBinContent(h.GetMaximumBin()),
        reverse=True)

    c = ROOT.TCanvas()
    c.SetLogy()
    lg = ROOT.TLegend(
        0.6,
        0.8,
        1 - c.GetRightMargin(),
        1 - c.GetTopMargin())
    lg.AddEntry(pimudif_0_2, "0-1 Atar Strips")
    lg.AddEntry(pimudif_3_5, "2-5 Atar Strips")
    lg.AddEntry(pimudif_6_inf, "6+ Atar Strips")
    lg.AddEntry(pimudif,'all')
    lg.AddEntry(pienu, 'pi->e')
    lg.SetFillStyle(0)
    lg.SetNColumns(2)
    plotables[0].Draw(style)
    for p in plotables[1:]:
        p.Draw('same' + style)
    lg.Draw('same')
    return c, lg

def plot_eff(
        pimudif_0_2,
        pimudif_3_5,
        pimudif_6_inf,
        pimudif,
        pienu,
        max_events=-1,
        style='HIST'):

    def fill_eff_rej(h, is_signal=True):
        eff = h.Clone("eff_" + h.GetName())
        eff.Reset()
        eff.GetXaxis().SetTitle("Upper cut on {}".format(
            h.GetXaxis().GetTitle()))
        eff.SetLineColor(h.GetLineColor())
        eff.SetFillColorAlpha(h.GetFillColor(), 0.05)
        for ibin in range(h.GetNbinsX() + 1):
            integral = h.Integral(0, h.GetNbinsX()+1)
            if is_signal:
                if integral != 0:
                    eff.SetBinContent(
                        ibin,
                        h.Integral(0, ibin) / integral)
                else:
                    eff.SetBinContent(ibin, 0)

            else:
                if h.Integral(0, ibin) != 0: 
                    eff.SetBinContent(ibin, integral / h.Integral(0, ibin))
                else:
                    eff.SetBinContent(ibin, 0.)
        return eff
    
    # return only effs
    eff_0_2 = fill_eff_rej(pimudif_0_2, is_signal=True)
    eff_3_5 = fill_eff_rej(pimudif_3_5, is_signal=True)
    eff_6_inf = fill_eff_rej(pimudif_6_inf, is_signal=True)
    eff_mudif = fill_eff_rej(pimudif, is_signal=True)
    eff_pienu = fill_eff_rej(pienu, is_signal=True)


    eff_pienu.SetLineColor(ROOT.kRed)
    eff_pienu.SetLineWidth(2)
    eff_pienu.SetFillStyle(0)

    plotables_rej = [
        eff_mudif,
        eff_0_2,
        eff_3_5,
        eff_6_inf,
        eff_pienu,
    ]
    
    plotables_rej = sorted(
        plotables_rej,
        key=lambda h: h.GetBinContent(h.GetMaximumBin()),
        reverse=True)


    c = ROOT.TCanvas()
    c.SetBottomMargin(0.2)
    c.Divide(1, 2, 0, 0)
    c.cd(1)
    
    for _p in plotables_rej:
        _p.GetYaxis().SetTitle('Cut Efficiency')
        _p.GetYaxis().SetRangeUser(0, 1.1)
        
    plotables_rej[0].Draw(style)
    for p in plotables_rej[1:]:
        p.Draw('same' + style)

    def _ratio(eff, den):
        r = eff.Clone('ratio_' + eff.GetName())
        r.GetYaxis().SetTitle('Eff(S) / Eff(B)')
        r.GetXaxis().SetTitle(eff.GetXaxis().GetTitle())
        r.Divide(den, eff)
        return r
    r_mudif = _ratio(eff_mudif, eff_pienu)
    r_0_2   = _ratio(eff_0_2, eff_pienu)
    r_3_5   = _ratio(eff_3_5, eff_pienu)
    r_6_inf = _ratio(eff_6_inf, eff_pienu)
        
    plotables_ratios = [
        r_mudif,
        r_0_2,
        r_3_5,
        r_6_inf,
        ]
    
    plotables_ratios = sorted(
        plotables_ratios,
        key=lambda h: h.GetBinContent(h.GetMaximumBin()),
        reverse=True)
    p = c.GetPad(2)
    p.cd()
    p.SetLogy()
    plotables_ratios[0].Draw('HIST')
    # plotables_ratios[0].Draw('HIST')
    for _r in plotables_ratios[1:]:
        _r.Draw('sameHIST')

    c.cd(0)
    lg = ROOT.TLegend(
        0.7,
        0.43,
        1 - c.GetRightMargin(),
        0.53)
    lg.SetNColumns(2)
    lg.AddEntry(pimudif_0_2, "0-1 Atar Strips")
    lg.AddEntry(pimudif_3_5, "2-5 Atar Strips")
    lg.AddEntry(pimudif_6_inf, "6+ Atar Strips")
    lg.AddEntry(pimudif,'all')
    lg.AddEntry(pienu, 'pi->e')
    lg.SetFillStyle(0)
    lg.Draw('same')
    c.Update()
    plotables = plotables_rej + plotables_ratios
        
    return c, lg, plotables

def plot_ene(
        pienu,
        pidar_mudar,
        pidar_mudif,
        max_events=-1,
        norm_strategy='branching_ratio'):

    pienu = pienu.Clone()
    if pidar_mudar != None:
        pidar_mudar = pidar_mudar.Clone()
    pidar_mudif = pidar_mudif.Clone()

    pienu.SetLineColor(ROOT.kRed)
    if pidar_mudar != None:
        pidar_mudar.SetLineColor(ROOT.kBlue)
    pidar_mudif.SetLineColor(ROOT.kViolet)

    pienu.SetLineWidth(2)
    if pidar_mudar != None:
        pidar_mudar.SetLineWidth(2)
    pidar_mudif.SetLineWidth(2)

    plotables = [
        pienu,
        # pidar_mudar,
        pidar_mudif
    ]

    sig_low = pienu.Integral(0, pienu.FindBin(58))
    sig_high = pienu.Integral(pienu.FindBin(58), pienu.GetNbinsX()+1)

    bkg_low = pidar_mudif.Integral(0, pidar_mudif.FindBin(58))
    bkg_high = pidar_mudif.Integral(pidar_mudif.FindBin(58), pidar_mudif.GetNbinsX()+1)

    tail_fraction = sig_low / sig_high if sig_high !=0 else 0.
    contamination_high = bkg_high / (bkg_high + sig_high)
    contamination_low = bkg_low / (bkg_low + sig_low)
    
    if norm_strategy == 'branching_ratio':
        if max_events != -1:
            for h in plotables:
                h.Scale(1. / max_events)
                h.SetMinimum(1e-11)
                h.GetYaxis().SetTitle('Branching Ratio')
    
    
    plotables = sorted(
        plotables,
        key=lambda h: h.GetBinContent(h.GetMaximumBin()),
        reverse=True)

    
    c = ROOT.TCanvas()
    c.SetLogy()
    if norm_strategy == 'branching_ratio':
        plotables[0].SetMinimum(1e-11)
        
    plotables[0].Draw('HIST')
    for p in plotables[1:]:
        p.Draw('sameHIST')
    lg = ROOT.TLegend(
        c.GetLeftMargin(),
        1 - c.GetTopMargin(),
        1,
        1)
    lg.AddEntry(pienu, 'pi#rightarrow enu')
    lg.AddEntry(pidar_mudif, 'pi [dar] - mu [dif] - e')
    # lg.AddEntry(pidar_mudar, 'pi [dar] - mu [dar] - e')
    lg.SetNColumns(len(plotables))
    lg.SetFillStyle(0)
    lg.SetBorderSize(0)
    # lg.SetTextSize(0.3 * lg.GetTextSize())
    lg.Draw()

    info = ROOT.TLatex()
    info.SetNDC(True)
    info.SetTextAlign(12)
    info.SetTextSize(0.8*info.GetTextSize())
    info.DrawLatex(
        c.GetLeftMargin() + 0.03,
        1 - c.GetTopMargin() - 0.05,
        'Tail Fraction = {0:1.4f}'.format(tail_fraction))
    info.DrawLatex(
        c.GetLeftMargin() + 0.03,
        1 - c.GetTopMargin() - 0.1,
        'c_{{high}} = {0:1.6f}'.format(contamination_high))
    info.DrawLatex(
        c.GetLeftMargin() + 0.03,
        1 - c.GetTopMargin() - 0.15,
        'c_{{low}} = {0:1.4f}'.format(contamination_low))
    return c, lg, plotables
    
# def plot_2bin(
#         pienu,
#         pidar_mudar,
#         pidar_mudif,
#         max_events=-1,
#         norm_strategy='branching_ratio'):

#     new_bins = array.array('d', [0, 58, 80])
#     pienu_2bin = pienu.Rebin(
#         2, pienu.GetName()+'_2bins', new_bins)
#     if pidar_mudar != None:
#         pidar_mudar_2bin = pidar_mudar.Rebin(
#             2, pidar_mudar.GetName()+'_2bins', new_bins)
#     pidar_mudif_2bin = pidar_mudif.Rebin(
#         2, pidar_mudif.GetName()+'_2bins', new_bins)
    
#     pienu_2bin.SetLineColor(ROOT.kRed)
#     if pidar_mudar != None:
#         pidar_mudar_2bin.SetLineColor(ROOT.kBlue)
#     pidar_mudif_2bin.SetLineColor(ROOT.kViolet)

#     pienu_2bin.SetFillColor(ROOT.kRed)
#     if pidar_mudar != None:
#         pidar_mudar_2bin.SetFillColor(ROOT.kBlue)
#     pidar_mudif_2bin.SetFillColor(ROOT.kViolet)

#     plotables = [
#         # pidar_mudar_2bin,
#         pidar_mudif_2bin,
#         pienu_2bin,
#     ]

#     stack = ROOT.THStack()
#     if norm_strategy == 'branching_ratio':
#         if max_events != -1:
#             for h in plotables:
#                 h.Scale(1. / max_events)
#                 h.SetMinimum(1e-11)
#                 h.GetYaxis().SetTitle('Branching Ratio')
    
    
#     for p in plotables:
#         stack.Add(p)
    
#     c = ROOT.TCanvas()
#     c.SetLogy()
#     stack.Draw("HIST")
#     stack.GetYaxis().SetTitle('Branching Ratio')
#     stack.SetMinimum(1e-11)
#     stack.GetXaxis().SetTitle(pienu.GetXaxis().GetTitle())
#     c.Update()
#     lg = ROOT.TLegend(
#         c.GetLeftMargin(),
#         1 - c.GetTopMargin(),
#         1,
#         1)
#     lg.AddEntry(pienu_2bin, 'pi#rightarrow enu')
#     lg.AddEntry(pidar_mudif_2bin, 'pi [dar] - mu [dif] - e')
#     # lg.AddEntry(pidar_mudar_2bin, 'pi [dar] - mu [dar] - e')
#     lg.SetNColumns(len(plotables))
#     lg.SetFillStyle(0)
#     lg.SetBorderSize(0)
#     # lg.SetTextSize(0.3 * lg.GetTextSize())
#     lg.Draw()

#     # contamination_table
    
#     return c, lg, plotables, stack


def plot_2bin(
        pienu,
        pidar_mudar,
        pidar_mudif,
        max_events=-1,
        norm_strategy='branching_ratio'):

    new_bins = array.array('d', [0, 58, 80])
    pienu_2bin = pienu.Rebin(
        2, pienu.GetName()+'_2bins', new_bins)
    if pidar_mudar != None:
        pidar_mudar_2bin = pidar_mudar.Rebin(
            2, pidar_mudar.GetName()+'_2bins', new_bins)
    pidar_mudif_2bin = pidar_mudif.Rebin(
        2, pidar_mudif.GetName()+'_2bins', new_bins)
    
    pienu_2bin.SetLineColor(ROOT.kRed)
    if pidar_mudar != None:
        pidar_mudar_2bin.SetLineColor(ROOT.kBlue)
    pidar_mudif_2bin.SetLineColor(ROOT.kViolet)

    pienu_2bin.SetFillColor(ROOT.kRed)
    if pidar_mudar != None:
        pidar_mudar_2bin.SetFillColor(ROOT.kBlue)
    pidar_mudif_2bin.SetFillColor(ROOT.kViolet)

    plotables = [
        # pidar_mudar_2bin,
        pidar_mudif_2bin,
        pienu_2bin,
    ]

    contamination_fraction = pidar_mudif_2bin.Clone('contamination' + pidar_mudif_2bin.GetName())
    total_hist = pienu_2bin.Clone('Total'+pienu_2bin.GetName())
    for ibin in [1, 2]:
        total_hist.SetBinContent(ibin, 1)
        contamination_fraction.SetBinContent(
            ibin,
            pienu_2bin.GetBinContent(ibin) / (pienu_2bin.GetBinContent(ibin) + pidar_mudif_2bin.GetBinContent(ibin)))
                                         
    total_hist.Print()

    total_hist.GetXaxis().Print()
    c = ROOT.TCanvas()
    total_hist.Draw('HIST')
    # contamination_fraction.Draw('sameHIST')
    lg = ROOT.TLegend(
        c.GetLeftMargin(),
        1 - c.GetTopMargin(),
        1,
        1)
    lg.AddEntry(pienu_2bin, 'pi#rightarrow enu')
    lg.AddEntry(pidar_mudif_2bin, 'pi [dar] - mu [dif] - e')
    # lg.AddEntry(pidar_mudar_2bin, 'pi [dar] - mu [dar] - e')
    lg.SetNColumns(len(plotables))
    lg.SetFillStyle(0)
    lg.SetBorderSize(0)
    # lg.SetTextSize(0.3 * lg.GetTextSize())
    lg.Draw()

    # contamination_table
    
    return c, lg
