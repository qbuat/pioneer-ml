def get_weighted_time(pix):
    if sum(pix['edep']) != 0:
        return sum(pix['ene_times_time']) / sum(pix['edep'])
    else:
        return 0.
    
def get_ene(pix):
    return sum(pix['edep'])
    

def get_true_muonpix_ids(atar):
    true_muon_pixel_ids = []
    for a in atar:
        if a.GetPDGID() == -13:
            if not a.GetPixelID() in true_muon_pixel_ids:
                true_muon_pixel_ids += [a.GetPixelID()]
    return true_muon_pixel_ids

def build_merged_atar(atar, merging_time=5):
    pixels_5ns_merged = []
    current_pixel = -999
    for a in atar:
        if a.GetPixelID() != current_pixel:
            current_pixel = a.GetPixelID()
            pixels_5ns_merged += [{
                    'pixelid': [a.GetPixelID()],
                    'pdgid': [a.GetPDGID()],
                    'edep': [a.GetEdep()],
                    'time': [a.GetTime()],                    
                    'ene_times_time': [a.GetEdep() * a.GetTime()],
            }]
        else:
            last_time_pix = pixels_5ns_merged[-1]['time'][-1]
            if abs(a.GetTime() - last_time_pix) < merging_time:
                pixels_5ns_merged[-1]['pixelid'].append(a.GetPixelID()) 
                pixels_5ns_merged[-1]['pdgid'].append(a.GetPDGID()) 
                pixels_5ns_merged[-1]['edep'].append(a.GetEdep())  
                pixels_5ns_merged[-1]['time'].append(a.GetTime())
                pixels_5ns_merged[-1]['ene_times_time'].append(a.GetEdep() * a.GetTime())  
            else:
                pixels_5ns_merged += [{
                    'pixelid': [a.GetPixelID()],
                    'pdgid': [a.GetPDGID()],
                    'edep': [a.GetEdep()],
                    'time': [a.GetTime()],
                    'ene_times_time': [a.GetEdep() * a.GetTime()],
                }]
    return pixels_5ns_merged



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
