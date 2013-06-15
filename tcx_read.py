import xml.etree.ElementTree as ET
from xml.etree.ElementTree import ElementTree
import re
from datetime import datetime
import math as m
from collections import OrderedDict
import matplotlib.pyplot as plt

def startup():
    john = athlete(78)
    john.set_bike_mass(9)
    john.set_roll_res(0.004)
    eyn = tcx_import('eynsham_cycle_2013.tcx', john)
    eyn.create_tp_list()
    eyn.list_gap_check(2)
    #eyn.create_tpdv_list()
    return(eyn)




class session(object):
    def __init__(self, athlete):
        self.athlete_cfg = athlete_cfg
        self.sumtime = []
        self.elapsedtime = []
        self.laps = []
        self.datapts = []

class athlete(object):
    def __init__(self, mass):
        self.mass = mass

    def set_bike_mass(self, mass):
        self.bike_mass = mass

    def set_roll_res(self, rr = 0.004):
        self.bike_roll_res = rr

class trackpoint:
    """This class deals with sample points from a .tcx file,
        referred to as trackpoints in the .tcx file"""
    def __init__(self, athlete, _raw_tp, _id, _tpcount, _lapnumber):
        self.raw_tp = _raw_tp
        self.id = _id
        self.tpcount = _tpcount
        self.lapnumber = _lapnumber
        self.time_rd = time_parse_garmin(self.raw_tp.find('Time').text)
        self.time_st = self.time_rd.timestamp()
        _tp_dict = self.__tp_dict__()
        self.tpchan = {}
        for key in _tp_dict:
            if self.raw_tp.find('.//' + _tp_dict[key]) != None: #check if it exists
                setattr(self, key, float(self.raw_tp.find('.//' + _tp_dict[key]).text))
                self.tpchan[key] = 'file'
            else:
                self.tpchan[key] = 'absent'


    def reverse_count(self, _tp_sum):
        self.tpcount_rev = _tp_sum - self.tpcount -1
    
  
    def __tp_dict__(self):
        """This is a dictionary of the terms used in the Garmin .tcx file"""
        _tp_dict = {
                    "lat":"LatitudeDegrees",
                    "lon":"LongitudeDegrees",
                    "alt":"AltitudeMeters",
                    "dist_vsen":"DistanceMeters",
                    "hr_bpm":"Value",
                    "speed_vsen":"Speed",
                    "power":"Watts",
                    }
        return(_tp_dict)

    def calc_tp_der(self, _dtp, direc = -1):
        """Derivative calculator
        direc is either -1 for previous or +1 for next
        """
        #check if end point
        if self.tpcount + direc < 0:
            return
        if self.tpcount_rev + direc < 0:
            return
        #calculate gps distances
        gp = global_point(self)
        gp_d = global_point(_dtp)
        self.dist_geo_step = {direc:gp.dist_geo(gp_d)}
        self.vect_geo_step = {direc:gp.vect_geo(gp_d)}
        _dtp_dict = self.__dtp_dict__()
        _tp_der_dict = self.__tp_der_dict__()
        for key in _dtp_dict:
            try:
                dict_delta = getattr(self, key)
            except:
                dict_delta = {}
            else:
                if isinstance(dict_delta, dict) != True:
                    dict_delta = {None:dict_delta}
            #These are now dictionary entries
            delta_left = (getattr(self, _dtp_dict[key]) * (-1 * direc)) 
            delta_right = (getattr(_dtp, _dtp_dict[key]) * direc)
            delta = delta_left + delta_right
            dict_delta_add = {direc:delta}
            dict_delta.update(dict_delta_add)
            setattr(self, key, dict_delta)
        for der in _tp_der_dict:
            dtp = _tp_der_dict[der]
            for key in dtp:
                try:
                    dict_delta = getattr(self, key)
                except:
                    dict_delta = {}
                else:
                    if isinstance(dict_delta, dict) != True:
                        dict_delta = {None:dict_delta}
                top = (getattr(self, dtp[key])[direc] * (-1 * direc))
                delta_bottom = getattr(self, der)[direc] * (-1 * direc)
                delta = top / delta_bottom
                dict_delta_add = {direc:delta}
                dict_delta.update(dict_delta_add)
                setattr(self, key, dict_delta)
        #check if end point neighbour
        if self.tpcount + (direc * 2) < 0:
            return
        if self.tpcount_rev - (direc * 2) < 0:
            return
        _dtp2_dict = self.__dtp2_dict__()
        _tp_der2_dict = self.__tp_der2_dict__()
        #print(dir(self), self.tpcount)
        for key in _dtp2_dict:
            try:
                dict_delta = getattr(self, key)
            except:
                dict_delta = {}
            else:
                if isinstance(dict_delta, dict) != True:
                    dict_delta = {None:dict_delta}
            delta = (getattr(self, _dtp2_dict[key])[direc] * (-1 * direc)) + (getattr(_dtp, _dtp2_dict[key])[direc] * direc)
            dict_delta_add = {direc:delta}
            dict_delta.update(dict_delta_add)
            setattr(self, key, dict_delta)
        for der in _tp_der2_dict:
            dtp = _tp_der2_dict[der]
            for key in dtp:
                try:
                    dict_delta = getattr(self, key)
                except:
                    dict_delta = {}
                else:
                    if isinstance(dict_delta, dict) != True:
                        dict_delta = {None:dict_delta}
                top = (getattr(self, dtp[key])[direc] * (-1 * direc))
                delta_bottom = getattr(self, der)[direc] * (-1 * direc)
                delta = top / delta_bottom
                dict_delta_add = {direc:delta}
                dict_delta.update(dict_delta_add)
                setattr(self, key, dict_delta)

    def __aero_drag__(self, key = -1):
        pass


    def __dtp_dict__(self):
        """Creates an ordered dictionary of tp steps"""
        _dtp_dict = OrderedDict()
        _dtp_list = [
                    'dist_vsen',
                    'time_st',
                    'speed_vsen',
                    'alt',
                    ]
        step_add = '_step'
        for sa in _dtp_list:
            _dtp_dict[(sa + step_add)] = sa
        return(_dtp_dict)

    def __tp_der_dict__(self):
        """Creates an ordered dictionary of derivative terms"""
        _tp_ddist = OrderedDict()
        _tp_ddist['gradient'] = 'alt_step'
        _tp_dtime = OrderedDict()           
        _tp_dtime['speed_ddist_geo'] = 'dist_geo_step'
        _tp_dtime['speed_ddist_vsen'] = 'dist_vsen_step'
        _tp_dtime['acc_dspeed_vsen'] = 'speed_vsen_step'
        _tp_der_dict_ = OrderedDict()
        _tp_der_dict_['dist_vsen_step'] = _tp_ddist
        _tp_der_dict_['time_st_step'] = _tp_dtime
        return(_tp_der_dict_)

    def __dtp2_dict__(self):
        """Creates an ordered dictionary of 2ndry tp steps"""
        _dtp_dict = OrderedDict()
        _dtp_list = [
                    'speed_ddist_vsen',
                    'speed_ddist_geo',
                    ]
        step_add = '_step'
        for sa in _dtp_list:
            _dtp_dict[(sa + step_add)] = sa
        return(_dtp_dict)

    def __tp_der2_dict__(self):
        """Creates an ordered dictionary of 2nd derivative terms"""
        _tp_dtime = OrderedDict()           
        _tp_dtime['acc_ddist2_vsen'] = 'speed_ddist_vsen_step'
        _tp_dtime['acc_ddist2_geo'] = 'speed_ddist_geo_step'
        _tp_der_dict_ = OrderedDict()
        _tp_der_dict_['time_st_step'] = _tp_dtime
        return(_tp_der_dict_)



class global_point:
    """Processes geo information of a trackpoint"""
    def __init__(self, _trackpoint):
        self.lat = _trackpoint.lat
        self.lon = _trackpoint.lon
        self.alt = _trackpoint.alt

    def global_rad(self):
        """Computes the radius of the earth for a given latitude and elevation
            See Wikipedia for eqns"""
        _eqt_rad = 6378137.0 #meters
        _pole_rad = 6356752.3 #meters
        _lat_rad = self.lat * m.pi/180
        _eqn_top = ((_eqt_rad**2) * m.cos(_lat_rad))**2 + ((_pole_rad**2) * m.sin(_pole_rad))**2
        _eqn_bot = ((_eqt_rad) * m.cos(_lat_rad))**2 + ((_pole_rad) * m.sin(_pole_rad))**2
        _radius = (_eqn_top/_eqn_bot)**0.5 + self.alt
        return(_radius)

    def dist_geo(self, _gpoint2):
        """Calculates the geo distance between two points, 
            default radius is close for GB,
            see http://williams.best.vwh.net/avform.htm for eqns"""
        _lat1_rads = self.lat * m.pi/180
        _lon1_rads = self.lon * m.pi/180
        _lat2_rads = _gpoint2.lat * m.pi/180
        _lon2_rads = _gpoint2.lon * m.pi/180
        _left_eqn = (m.sin((_lat1_rads-_lat2_rads)/2))**2
        _right_eqn = m.cos(_lat1_rads)*m.cos(_lat2_rads)*(m.sin((_lon1_rads-_lon2_rads)/2))**2
        _dist = self.global_rad() * 2 * m.asin(m.sqrt(_left_eqn + _right_eqn))
        return(_dist)

    def vect_geo(self, _gpoint2, heading = False):
        """Calculates the vector between two geo points
        and returns it as a tuple with """
        _lat1_rads = self.lat * m.pi/180
        _lon1_rads = self.lon * m.pi/180
        _lat2_rads = _gpoint2.lat * m.pi/180
        _lon2_rads = _gpoint2.lon * m.pi/180
        _radius = self.global_rad()
        _dist_ns = _radius * 2 * m.asin(m.sqrt((m.sin((_lat1_rads-_lat2_rads)/2))**2))
        _dist_ew = _radius * 2 * m.asin(m.sqrt(m.cos(_lat1_rads)*m.cos(_lat2_rads)*(m.sin((_lon1_rads-_lon2_rads)/2))**2))
        _dist_from_vect = m.sqrt(_dist_ns**2 + _dist_ew**2)
        _vect_ns = _dist_ns / _dist_from_vect
        _vect_ew = _dist_ew / _dist_from_vect
        if heading == False:
            return((_vect_ew, _vect_ns))
        else:
            _heading = (180 / m.pi) * m.atan((-_vect_ew)/_vect_ns)
            if _heading < 0:
                _heading = 360 + _heading
            return(_heading)

class tcx_import(ElementTree):
    '''Imports a tcx file, fname and then gives it several properties'''
    def __init__(self, fname, athlete):
        ElementTree.__init__(self)
        self.root = self.__initialize__(fname)
        self.activity = self.root.findall('.//Activity')
        self.activitytype = []
        self.id = []
        self.laps = []
        self.athlete = athlete
        for act in self.activity:
            self.activitytype.append(act.attrib)
            self.id.append(act.find('Id').text)
            for lap in range(1,len(act)-1):
                self.laps.append(act[lap])

    def __initialize__(self, _fname):
        _xmlstring = open(_fname).read()
        #Ditch the name spaces
        _xmlstring = re.sub(' xmlns="[^"]+"', '', _xmlstring)
        _root = ET.fromstring(_xmlstring)
        return(_root)

    def create_tp_list(self):
        """Create a list of the trackpoints and creates the method tplist and tplist_len"""
        _lapnumber = 1
        _tp_sum = 0
        _trackpoint_list = []
        for _lap in self.laps:
            _trackpoints_raw = _lap.findall('.//Trackpoint')
            for _raw_tp in _trackpoints_raw:
                _trackpoint_list.append(trackpoint(self.athlete, _raw_tp, self.id, _tp_sum, _lapnumber))
                _tp_sum += 1
            _lapnumber += 1
        for _trackpoint in _trackpoint_list:
            _trackpoint.reverse_count(_tp_sum)
        self.tplist = _trackpoint_list
        self.tplist_len = len(self.tplist)

    def smooth_value(self, param, length, direc = None):
        """Creates a moving average of the value over the length parameter"""
        if direc == None:
            doff = 0
        else:
            doff = direc
        if doff <= 0:
            begin = doff * -1
            end = self.tplist_len
        else:
            begin = 0
            end = self.tplist_len - doff
        if length > self.tplist_len:
            length = self.tplist_len - 1
        #blen must be odd, elen must be even 
        blen = (length // 2) +1 + doff
        elen = (length // 2)
        tlen = blen + elen
        keyout = 'sm' + str(length)
        for i in range(begin, end):
            init_vals = getattr(self.tplist[i], param)
            if isinstance(init_vals, dict) != True:
                init_vals = {None:init_vals}
            p_sum = []
            _tpcount = self.tplist[i].tpcount
            _tpcount_rev = self.tplist[i].tpcount_rev
            psum = 0
            count = 0
            err_count = 0
            if _tpcount < blen:
                lb = 0
                ub = i + elen
            elif _tpcount_rev < elen:
                lb = i - blen
                ub = self.tplist_len
            else:
                lb = i - blen
                ub = i + elen
            for j in range(lb,ub):
                try:    
                    val= getattr(self.tplist[j], param)
                except:
                    print('value does not exist')
                    err_count += 1
                    if err_count > (elen-1):
                        return
                else:
                    if isinstance(val, dict):
                        if direc != None:
                            valout = val[direc]
                        else:
                            keymin = 100000
                            count2 = 0
                            for key in val.keys():
                                if key == None:
                                   # print(i, 'key = ', key, val, count)
                                   valout = val[key]
                                   break
                                keyv = int(re.sub("\D", "", key))
                                if abs(keyv) < keymin:
                                    valout = val[key]
                                    keymin = abs(keyv)
                        val = valout
                    psum += val
                    count += 1
            dict_val = init_vals
            up_val = {keyout:(psum/count)}
            dict_val.update(up_val)
            setattr(self.tplist[i], param, dict_val)


    def create_tpdv_list(self, direc = 1):
        """Updates tplist with the derivative terms"""
        if direc <= 0:
            begin = direc * -1
            end = self.tplist_len
            for i in range(begin, end):
                self.tplist[i].calc_tp_der(self.tplist[i + direc], direc)
        else:
            begin = 0
            end = self.tplist_len - direc
            for i in range((end -1), 0, -1):
                self.tplist[i].calc_tp_der(self.tplist[i + direc], direc)

    def list_gap_check(self, gap = 1):
        #check initial values
        for key in self.tplist[0].tpchan:
            if self.tplist[0].tpchan[key] == 'absent':
                self.tplist[0].tpchan[key] == 'nonex'
        gaps = []
        for i in range(gap,(self.tplist_len - gap)):
            for key in self.tplist[i].tpchan:
                outp = self.tplist[i].tpchan[key]
                if outp == 'absent':
                    #if previous is nonexistent, continue non-existent
                    if self.tplist[i-1].tpchan[key] == 'nonex':
                        self.tplist[i].tpchan[key] = 'nonex'
                    elif self.tplist[i-1].tpchan[key] == 'absent':
                        #then the gap was too large so...
                        self.tplist[i].tpchan[key] = 'absent'
                    else:
                        for j in range(gap):
                            #if value exists on later tps
                            if self.tplist[i+j+1].tpchan[key] == 'file':
                                strval = getattr(self.tplist[i-1], key)
                                endval = getattr(self.tplist[i+j+1], key)
                                #note that interpolation is by time, not by distance
                                strdiv = self.tplist[i-1].time_st
                                enddiv = self.tplist[i+j].time_st
                                slope = ((endval - strval) / (enddiv - strdiv)) / (j + 2)
                                for k in range(j+1):
                                    val = slope * (self.tplist[i+k].time_st - strdiv) + strval
                                    setattr(self.tplist[i+k], key, val)
                                    self.tplist[i+k].tpchan[key] = 'interp'
                                    gaps.append((i+k, key, ' bridged'))
                                break
                        if self.tplist[i].tpchan[key] != 'interp':
                            for j in range(gap):
                                gaps.append((i+j,key, ' open'))
        return(gaps)

    def plot_xy(self, xchan, ychan, dx = None, dy = None):
        xval = []
        yval = []
        for i in range(2,self.tplist_len -2):
            xv = getattr(self.tplist[i], xchan)
            yv = getattr(self.tplist[i], ychan)
            #find which derivative value to use
            valouty = yv
            if isinstance(yv, dict):
                if dy != None:
                    valouty = yv[dy]
                else:
                    keymin = 100000
                    for key in yv.keys():
                        if key == None:
                            valouty = yv[key]
                            break
                        keyv = int(re.sub("\D", "", key))
                        if abs(keyv) < keymin:
                            valouty = yv[key]
                            keymin = abs(keyv)
            yv = valouty
            valoutx = xv
            if isinstance(xv, dict):
                if dx != None:
                    valoutx = xv[dx]
                else:
                    keymin = 100000
                    for key in xv.keys():
                        if key == None:
                            valoutx = xv[key]
                            break
                        keyv = int(re.sub("\D", "", key))
                        if abs(keyv) < keymin:
                            valoutx = xv[key]
                            keymin = abs(keyv)
            xv = valoutx
            xval.append(xv)
            yval.append(yv)
        plt.plot(xval, yval)
        plt.show()



def time_parse_garmin(_time_raw):
    """Parses time from garmin time stamp"""
    _time = datetime.strptime(_time_raw, '%Y-%m-%dT%H:%M:%S.%fZ')
    return(_time)


        
				
