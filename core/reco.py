import ROOT

def get_weighted_time(pix):
    if sum(pix['edep']) != 0:
        return sum(pix['ene_times_time']) / sum(pix['edep'])
    else:
        return 0.
    
def get_ene(pix):
    return sum(pix['edep'])

def get_atar_distance(a1, a2):
    dx = abs(a1.GetX1() - a2.GetX1())
    dy = abs(a1.GetY1() - a2.GetY1())
    dz = abs(a1.GetZ1() - a2.GetZ1())
    dr = ROOT.TMath.Sqrt(dx*dx + dy*dy + dz*dz)
    return dr, dx, dy, dz

def get_atar_single_hit_distance(a1):
    if not isinstance(a1, (list, tuple)):
        a1 = [a1]
    dx = abs(a1[-1].GetX1() - a1[0].GetX0())
    dy = abs(a1[-1].GetY1() - a1[0].GetY0())
    dz = abs(a1[-1].GetZ1() - a1[0].GetZ0())
    dr = ROOT.TMath.Sqrt(dx*dx + dy*dy + dz*dz)
    return dr, dx, dy, dz

    
def get_truth_pion_travel(atar, verbose=False):
    atar_pion = list(filter(lambda a: a.GetPDGID() == 211, atar))
    merged_atar_hits = [[]]

    previous_hit = atar_pion[0]
    for hit in atar_pion:
        if hit.GetPixelID() != previous_hit.GetPixelID():
            merged_atar_hits.append([])
        if abs(hit.GetTime() - previous_hit.GetTime()) > 5:
            merged_atar_hits.append([])
        merged_atar_hits[-1] += [hit]
        previous_hit = hit

    if len(atar_pion) == 0:
        return -1., -1, -1

    dr_last, dx_last, dy_last, dz_last= get_atar_single_hit_distance(merged_atar_hits[-1])
    theta_last = ROOT.TVector3(dx_last, dy_last, dz_last).Theta()
    de_last = sum(a.GetEdep() for a in merged_atar_hits[-1])
    if len(merged_atar_hits) > 1:
        de_previous_to_last = sum(a.GetEdep() for a in merged_atar_hits[-2])
        dr_previous_to_last, dx_previous_to_last, dy_previous_to_last, dz_previous_to_last = get_atar_single_hit_distance(merged_atar_hits[-2])
        theta_previous_to_last = ROOT.TVector3(
            dx_previous_to_last, dy_previous_to_last, dz_previous_to_last).Theta()
    else:
        de_previous_to_last = -1.
        dr_previous_to_last = -1
        dz_previous_to_last = -1
        theta_previous_to_last = -1.
        
    if verbose:
        if dr_previous_to_last > 0.12:
            print('last: {0:1.3f}, {1:1.3f}, {2:1.3f}, {3:1.3f}, {4:1.3f}'.format(
                dr_last, dx_last, dy_last, dz_last, theta_last))
            print('prev: {0:1.3f}, {1:1.3f}, {2:1.3f}, {3:1.3f}, {4:1.3f}'.format(
                dr_previous_to_last, dx_previous_to_last, dy_previous_to_last, dz_previous_to_last, theta_previous_to_last))
            print()
    return de_last, de_previous_to_last, dr_last, dr_previous_to_last, theta_last, theta_previous_to_last

def get_dx_from_de(de, last_plus_prev_to_last=False):
    """
    """
    if last_plus_prev_to_last:
        return pow(pion_ene_last_hit / 5.42, 1.73)
    return pow(pion_ene_last_hit / 5.55, 1.71)

def get_true_muon_first_pixel_dz(atar):
    for a in atar:
        if a.GetPDGID() == -13:
            dz = abs(a.GetZ1() - a.GetZ0())
            dx = ROOT.TMath.Sqrt(
                ROOT.TMath.Power(a.GetX1() - a.GetX0(), 2) + 
                ROOT.TMath.Power(a.GetY1() - a.GetY0(), 2) + 
                ROOT.TMath.Power(a.GetZ1() - a.GetZ0(), 2))
            dE = a.GetEdep() 
            return dz, dx, dE
    return 0, 0, 0

def get_true_distance_pion_muon(atar):
    atar_pion = list(filter(lambda a: a.GetPDGID() == 211, atar))
    if len(atar_pion) == 0:
        return -1.
    atar_muon = list(filter(lambda a: a.GetPDGID() == -13, atar))
    if len(atar_muon) == 0:
        return -2.
    delta_z = abs(atar_muon[0].GetZ0() - atar_pion[-1].GetZ1())
    return delta_z

def get_true_pion_stopping_position(atar, atar_geo):
    atar_pion = []
    for a in atar:
        if a.GetPDGID() == 211:
            atar_pion += [a]
    # atar_pion = list(filter(lambda a: a.GetPDGID() == 211, atar))
    x, y, z = atar_pion[-1].GetX1(), atar_pion[-1].GetY1(), atar_pion[-1].GetZ1()
    pixel_pos = atar_geo[atar_pion[-1].GetPixelID()]
    return (
        10 * pixel_pos['x'] + x,
        10 * pixel_pos['y'] + y,
        10 * pixel_pos['z'] + z,
        atar_pion[-1].GetPixelID())

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


def get_pion_stopping_plane(atar, atar_geo, debug=False):
    atar_odd = []
    atar_even = []
    x_z_truth = ROOT.TGraphErrors()
    y_z_truth = ROOT.TGraphErrors()
    x_z = ROOT.TGraphErrors()
    y_z = ROOT.TGraphErrors()
    for pix in atar:
        planeid = get_planeid(pix)
        if planeid % 2 == 0:
            # print(pix.GetPDGID(), pix.GetPixelID(), atar_geo[pix.GetPixelID()]['z'], pix.GetZ1())
            atar_even += [pix]
            x_z.AddPoint(
                atar_geo[pix.GetPixelID()]['z'],
                atar_geo[pix.GetPixelID()]['x'])
            x_z.SetPointError(
                x_z.GetN() - 1,
                0.120 / 2.,
                0.2 / 2.)
            x_z_truth.AddPoint(
                atar_geo[pix.GetPixelID()]['z'] + pix.GetZ1(),
                atar_geo[pix.GetPixelID()]['x'] + pix.GetX1())
        else:
            atar_odd += [pix]
            y_z.AddPoint(
                atar_geo[pix.GetPixelID()]['z'],
                atar_geo[pix.GetPixelID()]['y'])
            y_z.SetPointError(
                y_z.GetN() - 1,
                0.120 / 2.,
                0.2 / 2.)
            y_z_truth.AddPoint(
                atar_geo[pix.GetPixelID()]['z'] + pix.GetZ1(),
                atar_geo[pix.GetPixelID()]['y'] + pix.GetX1())

    z_min = atar_geo[atar[0].GetPixelID()]['z']
    z_max = atar_geo[atar[-1].GetPixelID()]['z']
    # print(z_min, z_max)
    fit_func_x_z = ROOT.TF1('fit_x_z', "[0]", z_min, z_max)
    # fit_func_x_z.SetParameter(0, 0)
    fit_func_x_z.SetParameter(0,  atar_geo[atar[0].GetPixelID()]['x'])
    x_z.Fit(fit_func_x_z, 'RQ')

    fit_func_y_z = ROOT.TF1('fit_y_z', "[0]", z_min, z_max)
    y_z.Fit(fit_func_y_z, 'RQ')

    if debug:
        c = ROOT.TCanvas()
        x_z.SetMarkerStyle(20)
        x_z.SetMarkerSize(2)
        x_z.Draw('AP')
        x_z.GetYaxis().SetRangeUser(-1, 1)
        x_z_truth.SetMarkerColor(2)
        x_z_truth.SetMarkerStyle(21)
        x_z_truth.Draw('sameP')
        fit_func_x_z.Draw('same')
        c.SaveAs('plots/fit_pixel_x_z_{}.pdf'.format(pix.GetPixelID()))

        c1 = ROOT.TCanvas()
        y_z.SetMarkerStyle(20)
        y_z.SetMarkerSize(2)
        y_z.Draw('AP')
        y_z.GetYaxis().SetRangeUser(-1, 1)
        y_z_truth.SetMarkerColor(2)
        y_z_truth.SetMarkerStyle(21)
        y_z_truth.Draw('sameP')
        fit_func_y_z.Draw('same')
        c1.SaveAs('plots/fit_pixel_y_z_{}.pdf'.format(pix.GetPixelID()))

    return (fit_func_x_z.GetParameter(0), fit_func_y_z.GetParameter(0))


def get_last_hit_pion_depth(atar0, atar_geo):
    pion_last_hit = atar0[-1]
    pion_ene_last_hit = 0
    for pixel in atar0:
        if pixel.GetPixelID() == pion_last_hit.GetPixelID():
            pion_ene_last_hit += pixel.GetEdep()
    # from Patrick's study
    travel_pion_last_pixel = pow(pion_ene_last_hit / 5.55, 1.71) #mm
    return pion_last_hit, travel_pion_last_pixel

def get_pion_stopping_point(atar0, atar_geo, debug=False):
    thickness = 0.120 #mm
    plane_size = 20.0 #mm
    # each ATAR plane is 2.0X2.0x0.0120 cm^3

    pion_last_hit, travel_pion_last_pixel = get_last_hit_pion_depth(atar0, atar_geo)
    z_last_pixel = atar_geo[pion_last_hit.GetPixelID()]['z'] #cm
    z_pion_stop = z_last_pixel * 10. - thickness / 2. + travel_pion_last_pixel
    
    # need to get them from the geom and also previous pixel (since each pixel gives a 2D pos)
    # x_pion_stop = pion_ene_last_hist.GetX1()
    # y_pion_stop = pion_ene_last_hist.GetY1()
    x_pion_stop, y_pion_stop = get_pion_stopping_plane(atar0, atar_geo, debug=debug)

    return x_pion_stop, y_pion_stop, z_pion_stop, pion_last_hit.GetPixelID()
    
            
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

def FindAtar_mudar(atarArray, atar0, atarT):
    atar_mudar  = []
    for a in atarArray:
        if a in atar0:
            continue
        if a in atarT:
            continue

        if a.GetTime() < atar0[0].GetTime():
            continue

        if a.GetTime() > atarT[-1].GetTime():
            continue
        atar_mudar += [a]
    return atar_mudar


def get_planeid(pix):
    ATAR_BASE_ID = 100001
    numPixels = 100
    planeid = int((pix.GetPixelID() - ATAR_BASE_ID) / numPixels)
    return planeid

def atar_ene_per_plane(atar):
    ene_per_plane = {}
    for pix in atar:
        planeid = get_planeid(pix)
        if planeid in ene_per_plane.keys():
            ene_per_plane[planeid] += pix.GetEdep()
        else:
            ene_per_plane[planeid] = pix.GetEdep()
    return ene_per_plane


def dE_dx(upstream, atar0, atarT, tracker_hit, atar_geo):
    thickness = 0.120 #mm
    plane_size = 20.0 #mm
    # each ATAR plane is 2.0X2.0x0.0120 cm^3

    pion_last_hit = atar0[-1]
    z_last_pixel = atar_geo[pion_last_hit.GetPixelID()]['z'] #cm
    # z_pion_stop = z_last_pixel * 10. - thickness / 2. + travel_pion_last_pixel
    
    x_pion_stop, y_pion_stop, z_pion_stop, _ =  get_pion_stopping_point(atar0, atar_geo, debug=False)

    # z_pion_stop = z_last_pixel * 10. - thickness / 2. + travel_pion_last_pixel
    
    travel_pion_last_pixel = z_pion_stop - (z_last_pixel * 10. - thickness / 2)
    # if travel_pion_last_pixel > thickness:
    #     print (travel_pion_last_pixel, z_pion_stop, z_last_pixel, thickness)
    # print()
    # print (x_pion_stop, y_pion_stop)
    # print (pion_ene_last_hist.GetX1(), pion_ene_last_hist.GetY1())
    ene_per_plane = atar_ene_per_plane(atarT)

    ups_vec = ROOT.TVector3(
        atar_geo[atar0[0].GetPixelID()]['x'],
        atar_geo[atar0[0].GetPixelID()]['y'],
        atar_geo[atar0[0].GetPixelID()]['z'])


    # pion_stop_vec = ROOT.TVector3(
    #     10 * x_pion_stop,
    #     10 * y_pion_stop,
    #     z_pion_stop)

    tracker_vec = ROOT.TVector3(
        tracker_hit.GetX1(),
        tracker_hit.GetY1(),
        tracker_hit.GetZ1())
    
    pion_stop_vec = ROOT.TVector3(
        10 * atar_geo[atar0[-1].GetPixelID()]['x'],
        10 * atar_geo[atar0[-1].GetPixelID()]['y'],
        z_pion_stop)
    # tracker_vec = ROOT.TVector3(
    #     10 * atar_geo[atarT[-1].GetPixelID()]['x'] + atarT[-1].GetX1(),
    #     10 * atar_geo[atarT[-1].GetPixelID()]['y'] + atarT[-1].GetY1(),
    #     10 * atar_geo[atarT[-1].GetPixelID()]['z'] + atarT[-1].GetZ1())
    # if get_planeid(atarT[-1])  % 2 != 0:
    #     tracker_vec.SetX(10 * atar_geo[atarT[-1].GetPixelID()]['x'] + atarT[-1].GetY1())
    #     tracker_vec.SetY(10 * atar_geo[atarT[-1].GetPixelID()]['y'] + atarT[-1].GetX1())

    
    
    cos_theta = (tracker_vec - pion_stop_vec).CosTheta()
    dE = atar_ene_per_plane(atarT)
    dE_dx_dict = {}
    dE_dx = 0
    dz_0 = 0

    
    for i_h, hit in enumerate(atarT):
        planeid = get_planeid(hit)
        if planeid in dE_dx_dict.keys():
            continue
        if i_h == 0:
            d_z = (thickness - travel_pion_last_pixel) #/ (cos_theta)
            if cos_theta > 0:
                if travel_pion_last_pixel > thickness:
                    d_z = 2* thickness - travel_pion_last_pixel 
                else:
                    d_z = (thickness - travel_pion_last_pixel)
            else:
                if travel_pion_last_pixel > thickness:
                    d_z = travel_pion_last_pixel - thickness
                else:
                    d_z = travel_pion_last_pixel
            dz_0 = d_z
        else:
            # d_z = thickness / abs(cos_theta)
            d_z = thickness #/ abs(cos_theta)

        if d_z < 0:
            print('negative dz!', d_z, thickness, travel_pion_last_pixel, cos_theta)
            
        dE_dx_dict[planeid] = dE[planeid] / d_z
        if dE[planeid] / d_z > dE_dx:
            dE_dx = dE[planeid] / d_z
        
    return dE_dx, cos_theta, dz_0
