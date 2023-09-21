# -- helper functions
def FindUpstream(usArray, time = -1, timeLimit = 1, trID = -1):
    retus = {"time" : [], "edep" : []}
    
    for itrid, itime, iedep in zip(
            usArray.GetTrackID()[1:], usArray.GetTime()[1:], usArray.GetEdep()[1:]):
        if trID > 0 and trID != itrid:
            continue
        if time > 0 and abs(itime - time) > timeLimit:
            continue
        retus['time'].append(itime)
        retus['edep'].append(iedep)
    return retus


def FindTracker(trackerArray, time = -1, timeLimit = 1, trID = -1):
    TA = [t for t in trackerArray if t.GetTrackerID() == 200002]
    if trID < 0:
        tr = [t for t in TA if abs(t.GetTime() - time) < timeLimit ]
    elif time < 0:
        tr = [t for t in TA if t.GetTrackID() == trID ]
    else:
        tr = [t for t in TA if abs(t.GetTime() - time) < timeLimit and t.GetTrackID() == trID ]
    return tr


def FindCalo(caloArray, time = -1, timeLimit = 1, trID = -1):
    retcal = {"time" : [], "edep" : [], "theta" : []}
    
    for calo in caloArray:
        for itrid, itime, iedep, itheta in zip(
                calo.GetTrackID(), calo.GetTime(), calo.GetEdep(), calo.GetTheta()):
            if trID > 0 and trID != itrid:
                continue
            if time > 0 and abs(itime - time) > timeLimit:
                continue
            retcal['time'].append(itime)
            retcal['edep'].append(iedep)
            retcal['theta'].append(itheta)
    
    return retcal
            
def FindAtar(atarArray, time = -1, timeLimit = 1, trID = -1):
    if trID < 0:
        atar = [a for a in atarArray if (abs(a.GetTime() - time) < timeLimit)]
    elif time < 0:
        atar = [a for a in atarArray if a.GetTrackID() == trID ]
    else:
        atar = [a for a in atarArray if a.GetTrackID() == trID and abs(a.GetTime() - time < timeLimit)]
    return atar

# reconstruction class
class reco_0(object):

    def __init__(self):
        self.upstream = None
        self.tracker = None
        self.tracker_hit = None
        self.atar0 = None
        self.atarT = None
        self.caloT = None
        self.usT = None
        self.T0 = None
        
        # reco params
        self.timeWin = 200 # ns
        self.beamRejection = 5 # ns
        self.trMergeTime = 5 #ns
        self.atarWin = 1 # ATAR search time window
        self.caloWin = 5 # Calo search time window

    def initialize_event(self):
        self.upstream = None
        self.tracker = None
        self.tracker_hit = None
        self.atar0 = None
        self.atarT = None
        self.caloT = None
        self.usT = None
        self.T0 = None
        return True
    
    def find_upstream(self, event):
        us = [u for u in event.upstream if u.GetUpstreamID() == 600000 ]
        if len(us) == 0:
            return False
        self.upstream = us[0]
        self.T0 = us[0].GetTime()[0]
        return True
            
    def reconstruct_tracker(self, event):
        if self.T0 == None:
            return False

        self.tracker = FindTracker(
            event.tracker,
            time = self.T0 + self.timeWin,
            timeLimit = self.timeWin)
        if len(self.tracker) == 0:
            return False
        
        TT_q = -self.trMergeTime
        merged_tracker = []
        for tr in self.tracker:
            if tr.GetTime() - self.T0 < self.beamRejection:
                # beam rejection window
                continue
            if tr.GetTime() < TT_q + self.trMergeTime:
                # simultaneous with previous tracker hit
                continue
            merged_tracker += [tr]
            TT_q = tr.GetTime()
        if len(merged_tracker) == 0:
            return False
            
        if len(merged_tracker) > 1:
            print ('evt \t x0 \t Y0 \t Z0 \t Time \t Edep')
            for track in merged_tracker:
                print('{0:d}, {1:1.3f}, {2:1.3f}, {3:1.3f}, {4:1.3f}, {5:1.3f}'.format(
                    evID, track.GetX0(), track.GetY0(), track.GetZ0(), track.GetTime(), track.GetEdep()))

        self.tracker_hit = merged_tracker[0]
        return True
    
    def reconstruct_atar(self, event):
    
        # Find T0 Atar hits
        self.atar0 = FindAtar(
            event.atar,
            time = self.T0 + self.atarWin,
            timeLimit= self.atarWin)

        if self.tracker_hit == None:
            return False
        
        TT = self.tracker_hit.GetTime()
        self.atarT = FindAtar(
            event.atar,
            time = TT - self.atarWin,
            timeLimit = self.atarWin)
        if len(self.atarT) == 0:
            return False
        return True

    def reconstruct_calo(self, event):
        if self.tracker_hit == None:
            return False

        TT = self.tracker_hit.GetTime()
        self.caloT = FindCalo(
            event.calo,
            time = TT,
            timeLimit = self.caloWin)
        if len(self.caloT) == 0:
            return False
        return True

    def reconstruct_upstream(self, event):
        if self.upstream == None:
            return False
        
        if self.tracker_hit == None:
            return False

        TT = self.tracker_hit.GetTime()
        self.usT = FindUpstream(
            self.upstream,
            time = TT - self.atarWin,
            timeLimit = self.atarWin)
        if len(self.usT) == 0:
            return False
        return True

# run reco on an event
def reconstruct_event(event, reco_obj):
    """
    """
    status = reco_obj.initialize_event()
    if status == False:
        return False

    status = reco_obj.find_upstream(event)
    if status == False:
        return False

    status = reco_obj.reconstruct_tracker(event)
    if status == False:
        return False

    status = reco_obj.reconstruct_atar(event)
    if status == False:
        return False

    status = reco_obj.reconstruct_calo(event)
    if status == False:
        return False

    status = reco_obj.reconstruct_upstream(event)
    if status == False:
        return False

    return True
