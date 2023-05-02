import ROOT
def plot_stack(
        pimudif_0_2,
        pimudif_3_5,
        pimudif_6_inf,
        pimudif,
        pienu,
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


    c = ROOT.TCanvas()
    c.SetLogy()
    lg = ROOT.TLegend(0.5, 0.6, 0.8, 0.85)
    lg.AddEntry(pimudif_0_2, "0-1 Atar Strips")
    lg.AddEntry(pimudif_3_5, "2-5 Atar Strips")
    lg.AddEntry(pimudif_6_inf, "6+ Atar Strips")
    lg.AddEntry(pimudif,'all')
    lg.AddEntry(pienu, 'pi->e')
    #lg.SetFillStyle(0)

    pienu.Draw(style)
    pimudif_0_2.Draw('SAME' + style)
    pimudif_3_5.Draw('SAME' + style)
    pimudif_6_inf.Draw('SAME' + style)
    pimudif.Draw('SAME' + style)
    lg.Draw('same')
    return c, lg

