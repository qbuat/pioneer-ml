import ROOT
from core.reconstruction import reco_0, reconstruct_event



if __name__ == '__main__':

    beam50ps = ROOT.TChain("sim", "sim")
    beam50ps.Add("/Users/quentin/decay_new/sim/beam_50ps*.root")
    print(beam50ps.GetEntries())
    maxEvents = 10

    chain = beam50ps
    reconstructor = reco_0()
    for evID, event in enumerate(chain):
        if maxEvents > 0 and evID > maxEvents:
            break
        reconstruct_event(event, reconstructor)
        




