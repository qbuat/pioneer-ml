import re
import xml.etree.ElementTree as ET
try:
    import uproot
    import awkward as ak
except ModuleNotFoundError:
    print("No Uproot, it's awkward")
import numpy as np
import itertools

class atar_hit(object):
    def __init__(self, edep, time, pdgid, pxlID, planeid, pixelid):
        self._edep = [edep]
        self._time = [time]
        self._pdgid = [pdgid]
        self._pxlID = pxlID
        self._planeid = planeid
        self._pixelid = pixelid
        self._merged = False
        
    def should_merge(self, other):
        if self._pxlID == other._pxlID:
            if abs(self.merged_time - other._time[0]) < 1:
                self._merged = True
                return True
        return False

    def fill(self, other):
        # if not isinstance(self._edep, list):
        #     self._edep = [self._edep]
        #     self._time = [self._time]
        #     self._pdgid = [self._pdgid]

        self._edep.append(other._edep[0])
        self._time.append(other._time[0])
        self._pdgid.append(other._pdgid[0])

        if self._pxlID != None:
            if self._pxlID != other._pxlID:
                raise ValueError
        else:
            self._pxlID = other._pxlID
        
        if self._planeid != None:
            if self._planeid != other._planeid:
                raise ValueError
        else:
            self._planeid = other._planeid


        if self._pixelid != None:
            if self._pixelid != other._pixelid:
                raise ValueError
        else:
            self._planeid = other._pixelid
            
    @property
    def merged_edep(self):
        return sum(self._edep)

    @property
    def merged_time(self):
        return sum(self._time) / len(self._time) if len(self._time)!=0 else 0.

                
            
def geometry(gdml_file):
    """
    """
    # file = "/Users/qbuat/pioneer/main/MonteCarlo/geometry/generator/test_output.gdml"
    # file = "/Users/quentin/pioneer/PIONEER/MonteCarlo/geometry/generator/test_output.gdml"
    gdmlTree = ET.parse(gdml_file)
    gdmlRoot = gdmlTree.getroot()

    atar_geo = {}
    pattern = r'target_pixel_(?P<pixid>10[0-9][0-9][0-9][0-9])intarget_containerpos'
    prog = re.compile(pattern)
    for child in gdmlRoot:
        if child.tag != 'define':
            continue
        for gc in child:
            if gc.tag != 'position':
                continue
            if not 'name' in gc.attrib.keys():
                raise ValueError
            if not 'intarget' in gc.attrib['name']:
                continue
            result = prog.match(gc.attrib['name'])
            atar_geo[int(result.group('pixid'))] = {
                'x': float(gc.attrib['x']),
                'y': float(gc.attrib['y']),
                'z': float(gc.attrib['z']),
            }
    # print (atar_geo)
    return atar_geo

def get_atar_pos(pix, variable='x'):
    if pix == -999:
        return -999
    return _atar_geo[pix][variable]


def compute_atar_ids(arr):
    ATAR_BASE_ID = 100001
    numPixels = 100
    arr['atar.planeid'] = ak.values_astype((arr['atar.pxlID'] - ATAR_BASE_ID) / numPixels, 'int64')
    arr['atar.pixelid'] = (arr['atar.pxlID'] - ATAR_BASE_ID) % numPixels


def energy_centrality(atar_arrays):
    _first_atar_hit_time = ak.min(atar_arrays['atar.time'], axis=1)
    _atar_intime = atar_arrays[np.fabs(atar_arrays['atar.time'] - _first_atar_hit_time) < 1]
    _atar_intime_central =  _atar_intime[_atar_intime['atar.pixelid'] > 4]
    _atar_intime_central =  _atar_intime_central[_atar_intime_central['atar.pixelid'] < 95]
    _atar_intime_central =  _atar_intime_central[_atar_intime_central['atar.planeid'] > 7]
    _atar_intime_central =  _atar_intime_central[_atar_intime_central['atar.planeid'] < 40]
    print (ak.num(_atar_intime_central['atar.edep']))
    _ecentral =  ak.sum(_atar_intime_central['atar.edep'], axis=1)
    _etot =  ak.sum(_atar_intime['atar.edep'], axis=1)
    return _ecentral, _etot
    # atar_arrays['etot'] = _etot
    # atar_arrays['ecentral'] = _ecentral


def unique_atar(arr):
    """
    """

    arr = ak.zip({
        _f: arr[_f] for _f in filter(lambda n: 'atar.' in n, arr.fields)})

    unique_pxl_id, counts = np.unique(arr['atar.pxlID'], return_counts=True)
    non_uniques = counts[counts!=1]
    if len(non_uniques) == 0:
        return arr

    # grouped_pxl = itertools.groupby(arr['atar.pxlID'])
    # print ([(a, [i for i in b]) for a, b in grouped_pxl])


    hits = []    
    for a in arr:
        hit = atar_hit(
            a['atar.edep'],
            a['atar.time'],
            a['atar.pdgid'],
            a['atar.pxlID'],
            a['atar.planeid'],
            a['atar.pixelid'])

        mergers = list(filter(lambda h: h.should_merge(hit), hits))
        if len(mergers) != 0:
            mergers[0].fill(hit)
        else:
            hits.append(hit)

    for hit in hits:
        if hit._merged:
            print(hit._pxlID, hit._time, hit._edep)

    arr['merged_atart.edep'] = np.array([hit.merged_edep for hit in hits])        

    print (arr)
    # for pxl, cnt in zip(unique_pxl_id, counts):
    #     if cnt != 1:
    #         _ene = arr['atar.edep'][arr['atar.pxlID'] == pxl]
    #         _time = arr['atar.time'][arr['atar.pxlID'] == pxl]
    #         _pdgid = arr['atar.pdgid'][arr['atar.pxlID'] == pxl]
    #         print (pxl, _ene, _time, _pdgid)
    # print (unique_pxl_id)
    # print (counts)
