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
        _tp_dict_file = self._tp_dict()
        _tp_dict_origin = self.tp_dict_origin()
        self.tpchan = {}
        for key in _tp_dict_file:
            #check if it exists
            if self.raw_tp.find('.//' + _tp_dict_file[key]) != None: 
                #sets value as p for primary, ie from file
                dict_val = {'p':float(self.raw_tp.find(
                            './/' + _tp_dict_file[key]).text)}
                #sets origin of value as gps, vsen, etc.
                dict_type = {_tp_dict_origin[key]:dict_val}
                setattr(self, key, dict_type)
                self.tpchan[key] = 'file'
            else:
                self.tpchan[key] = 'absent'

    def reverse_count(self, _tp_sum):
        """Adds the number until the last point,
        it is the complement of tpcount
        """
        self.tpcount_rev = _tp_sum - self.tpcount -1
  
    def add_vsen(self):
        """Updates the distance and velocity to..."""

    def _tp_dict(self):
        """This is a dictionary of the terms used in the Garmin .tcx file"""
        _tp_dict_file = {
                    "lat":"LatitudeDegrees",
                    "lon":"LongitudeDegrees",
                    "alt":"AltitudeMeters",
                    "dist":"DistanceMeters",
                    "hr_bpm":"Value",
                    "speed":"Speed",
                    "power":"Watts",
                    }
        return(_tp_dict_file)
        
    def tp_dict_origin(self):
        _tp_data_origin = {
                    "lat":"gpsm",
                    "lon":"gpsm",
                    "alt":"gpsm",
                    "dist":"vsen",
                    "hr_bpm":"strap",
                    "speed":"vsen",
                    "power":"power_tap",
                    }
        return(_tp_data_origin)

    #Notation is self.dist['vsen']['p'] would be a value
    #Delta would be self.dist_step['vsen']['p_d-1']
    #Derivative would be self.speed['vsen']['p_d-1']
    #Delta would be self.speed_step['vsen']['p_d-1_2d1']
    #Derivative would be self.acc['vsen']['p_d-1_2d1']

    def calc_tp_der(self, _dtp, direc = -1, iskey = None, debug = False):
        """Derivative calculator
        direc is either -1 for previous or +1 for next
        """
        if iskey == None:
            iskey = 'p'
        #check if end point
        if self.tpcount + direc < 0:
            return
        if self.tpcount_rev + direc < 0:
            return
        pt1_notation = iskey + '_d' + str(direc)
        if ('d' in iskey) == False:
            pt_notation = pt1_notation
        else:
            pt_notation = iskey + '_2d' + str(direc)
        #calculate gps distances
        try:
            self.dist_step['gpsm'][('p_d'+direc)]
        except:
            self.dist_step = {}
            try:
                self.dist_step['gpsm']
            except:
                gp = global_point(self, 'gpsm', 'p')
                gp_d = global_point(_dtp, 'gpsm', 'p')
                self.dist_step['gpsm'] = {pt1_notation:gp.dist_geo(gp_d)}
                self.vect_step = {}
                self.vect_step['gpsm'] = {pt1_notation:gp.vect_geo(gp_d)}
        #create time steps
        self.time_st_step = self.time_st*(-1*direc) + _dtp.time_st*direc
        #create steps of values: dist_step, speed_step, alt_step
        _step_dict = self._step_dictf()
        for key in _step_dict:
            # If key was 'dist_step' then dict would return dist
            # though dist contains two values for vsen & gpsm
            # example below
            # dist = self._step_dict['dist_step']
            dict_of_chan = getattr(self, _step_dict[key])
            dict_of_chan_d = getattr(_dtp, _step_dict[key])
            for orig in dict_of_chan:
                # orig will be something like 'vsen' which = ('p': value)
                try:
                    # check existence of value on comparison tp
                    dict_of_chan_d[orig][iskey]
                except:
                    if debug == True:
                        print('comparison value for ', self.tpcount, key, dict_of_chan_d, 
                                orig, iskey, 'does not exist')
                else:
                    try:
                        test = getattr(self, key)[orig][pt_notation]
                    except:
                        # example
                        # delta_left = dist['vsen']['p']*(-1*-1)
                        delta_left = dict_of_chan[orig][iskey] * (-1 * direc)
                        delta_right = dict_of_chan_d[orig][iskey] * direc
                        delta = delta_left + delta_right
                        dict_delta_add = {pt_notation:delta}
                        try:
                            #check if method exists from another iskey
                            #example: self.acc['vsen']
                            key_dict = getattr(self, key)[orig]
                        except:
                            #create empty dict for step values
                            value_dict = {}
                            # example
                            # full_add = {'vsen':{'p_d-1': value}}
                            full_add = {orig:dict_delta_add}
                        else:
                            # value_dict contains a previous iskey
                            key_dict.update(dict_delta_add)
                            full_add = {orig:key_dict}
                            # This value_dict will overwrite previous values
                            value_dict = getattr(self, key)
                        # This value_dict will overwrite previous values
                        value_dict.update(full_add)
                        try:
                            # example
                            # self.dist_step.update(value_dict)
                            getattr(self, key).update(value_dict) #
                        except:
                            setattr(self, key, value_dict)
                    else:
                        #value already exists
                        pass
        #Now for the actual values after the steps have been created
        _tp_deriv_dict = self._tp_deriv_dictf()
        for deriv_type in _tp_deriv_dict:
            # the two keys for der are _tp_ddist and _tp_dtime 
            # from _tp_der_dictf retreieve _tp_ddist or _tp_dtime dict
            deriv_type_dict = _tp_deriv_dict[deriv_type]
            # deriv will be an attribute: gradient, speed, or acc ~ for now
            for deriv in deriv_type_dict:
                # for example
                # speed_step = self.deriv_type_dict['acc']
                try:
                    step_type = getattr(self, deriv_type_dict[deriv])
                except:
                    if debug == True:
                        print('step value for ', self.tpcount, deriv_type_dict[deriv], 
                                    orig, iskey, 'does not exist')
                    continue
                for orig in step_type:
                    try:
                        #check if attribute already exists for this exact combo
                        getattr(self, deriv)[orig][pt_notation] 
                    except:
                        # orig and iskey combo are new
                        # numerator will be the step
                        # direc cancels out, as the step was calculated w/ direc
                        # try is needed to check if step values exist
                        # example:
                        # numerator = self.speed_step['vsen']['p_d-1']
                        try:
                            numerator = getattr(self, deriv_type_dict[deriv])[orig][pt_notation]
                        except:
                            if debug == True:
                                print('step value for ', self.tpcount, deriv_type_dict[deriv], 
                                    orig, iskey, 'does not exist')
                            continue
                        if deriv_type == 'time_st_step':
                            #the time step has no dictionary parameter
                            denom = getattr(self, deriv_type)
                        else:
                            #Currently distance is the only other parameter
                            denom = getattr(self, deriv_type)[orig][pt_notation]
                        # print('deriv =', deriv, 'orig = ', orig, 'numerator = ',
                        #       numerator,'denom = ', denom)
                        delta = numerator / denom
                        dict_delta_add = {pt_notation:delta}
                        try:
                            # check if method exists from another iskey
                            # example: self.acc['vsen']
                            existing_deriv_dict = getattr(self, deriv)[orig]
                        except:
                            # create empty dict for step values
                            value_dict = {}
                            # example
                            # full_add = {'vsen':{'p_d-1': value}}
                            full_add = {orig:dict_delta_add}
                        else:
                            # value_dict contains a previous iskey
                            existing_deriv_dict.update(dict_delta_add)
                            full_add = {orig:existing_deriv_dict}
                            # This value_dict will overwrite previous values
                            value_dict = getattr(self, deriv)
                        value_dict.update(full_add)
                        # Now check if method exists for previous combo    
                        try:
                            #example
                            # self.acc.update(value_dict)
                            getattr(self, deriv).update(value_dict) #
                        except:
                            setattr(self, deriv, value_dict)
                    

    def calc_tp_der_old(self, _dtp, direc = -1, iskey = None):
        """Derivative calculator
        direc is either -1 for previous or +1 for next
        """
        if iskey == None:
            iskey = 'p'
        #check if end point
        if self.tpcount + direc < 0:
            return
        if self.tpcount_rev + direc < 0:
            return
        #calculate gps distances
        gp = global_point(self, 'gpsm', 'p')
        gp_d = global_point(_dtp, 'gpsm', 'p')
        if ('d' in iskey) == False:
            pt_notation = iskey + '_d' + str(direc)
        else:
            pt_notation = iskey + '_2d' + str(direc)
        self.dist_step = {}
        self.dist_step['gpsm'] = {pt_notation:gp.dist_geo(gp_d)}
        self.vect_step = {}
        self.vect_step['gpsm'] = {pt_notation:gp.vect_geo(gp_d)}
        #create time steps
        self.time_st_step = self.time_st*(-1*direc) + _dtp.time_st*direc
        #create steps of values: dist_step, speed_step, alt_step
        _step_dict = self._step_dictf()
        for key in _step_dict:
            # If key was 'dist_step' then dict would return dist
            # though dist contains two values for vsen & gpsm
            # example below
            # dist = self._step_dict['dist_step']
            dict_of_chan = getattr(self, _step_dict[key])
            dict_of_chan_d = getattr(_dtp, _step_dict[key])
            # print(key, dict_of_chan)
            for orig in dict_of_chan:
                # orig will be something like 'vsen' which = ('p': value)
                # print(dict_of_chan[orig])
                try:
                    # check existence of value on comparison tp
                    dict_of_chan_d[orig][iskey]
                except:
                    print('comparison value for ', self.tpcount, key, 
                            orig, iskey, 'does not exist')
                else:
                    try:
                        test = getattr(self, key)[orig][pt_notation]
                    except:
                        # example
                        # delta_left = dist['vsen']['p']*(-1*-1)
                        delta_left = dict_of_chan[orig][iskey] * (-1 * direc)
                        delta_right = dict_of_chan_d[orig][iskey] * direc
                        delta = delta_left + delta_right
                        dict_delta_add = {pt_notation:delta}
                        try:
                            #check if method exists from another iskey
                            #example: self.acc['vsen']
                            key_dict = getattr(self, key)[orig]
                        except:
                            #create empty dict for step values
                            value_dict = {}
                            # example
                            # full_add = {'vsen':{'p_d-1': value}}
                            full_add = {orig:dict_delta_add}
                        else:
                            # value_dict contains a previous iskey
                            key_dict.update(dict_delta_add)
                            full_add = {orig:key_dict}
                        # This value_dict will overwrite previous values
                        value_dict.update(full_add)
                        try:
                            # example
                            # self.dist_step.update(value_dict)
                            getattr(self, key).update(value_dict) #
                        except:
                            setattr(self, key, value_dict)
                    else:
                        #value already exists
                        pass
        #Now for the actual values after the steps have been created
        _tp_deriv_dict = self._tp_deriv_dictf()
        for deriv in _tp_deriv_dict:
            # the two keys for der are _tp_ddist and _tp_dtime 
            # from _tp_der_dictf
            dtp = _tp_deriv_dict[deriv]
            # key will be an attribute: gradient, speed, or acc ~ for now
            for key in dtp:
                # for example
                # speed_step = self.dtp['acc']
                steps_dict = getattr(self, dtp[key])
                for orig in steps_dict:
                    try:
                        getattr(self, key)[orig][pt_notation] 
                    except:
                        # orig and iskey combo are new
                        # numerator will be the step
                        # direc cancels out, as the step was calculated w/ direc
                        # try is needed to check if step values exist
                        # example:
                        # numerator = self.speed_step['vsen']['p_d-1']
                        try:
                            numerator = (getattr(self, dtp[key]))[orig][pt_notation]
                        except:
                            print('step value for ', self.tpcount, key, 
                                    orig, iskey, 'does not exist')
                        else:
                            if deriv != 'time_st_step':
                                denom = getattr(self, deriv)[orig][pt_notation]
                            else:
                                #the time step has no dictionary parameter
                                denom = getattr(self, deriv)
                            # print('key =', key, 'orig = ', orig, 'numerator = ',
                            #       numerator,'denom = ', denom)
                            delta = numerator / denom
                            dict_delta_add = {pt_notation:delta}
                            try:
                                # check if method exists from another iskey
                                # example: self.acc['vsen']
                                key_dict = getattr(self, key)[orig]
                            except:
                                # create empty dict for step values
                                value_dict = {}
                                # example
                                # full_add = {'vsen':{'p_d-1': value}}
                                full_add = {orig:dict_delta_add}
                            else:
                                # value_dict contains a previous iskey
                                key_dict.update(dict_delta_add)
                                full_add = {orig:key_dict}
                                # This value_dict will overwrite previous values
                                value_dict = getattr(self, key)
                            value_dict.update(full_add)
                            # Now check if method exists for previous combo    
                            try:
                                #example
                                # self.acc.update(value_dict)
                                getattr(self, key).update(value_dict) #
                            except:
                                setattr(self, key, value_dict)
                    else:
                        #This combo of orig and pt_notation exist already
                        pass

    def __aero_drag__(self, key = -1):
        pass


    def _step_dictf(self):
        """Creates an ordered dictionary of tp steps"""
        _step_dict = OrderedDict()
        _step_list = [
                    'dist',
                    'speed',
                    'alt',
                    ]
        step_add = '_step'
        for sa in _step_list:
            _step_dict[(sa + step_add)] = sa
        return(_step_dict)

    def _tp_deriv_dictf(self):
        """Creates an ordered dictionary of derivative terms"""
        _tp_ddist = OrderedDict()
        _tp_ddist['gradient'] = 'alt_step'
        _tp_dtime = OrderedDict()           
        _tp_dtime['speed'] = 'dist_step'
        _tp_dtime['acc'] = 'speed_step'
        _tp_deriv_dict_ = OrderedDict()
        _tp_deriv_dict_['dist_step'] = _tp_ddist
        _tp_deriv_dict_['time_st_step'] = _tp_dtime
        return(_tp_deriv_dict_)

   


class global_point:
    """Processes geo information of a trackpoint"""
    def __init__(self, _trackpoint, origin, key):
        self.lat = _trackpoint.lat[origin][key]
        self.lon = _trackpoint.lon[origin][key]
        #self.alt = _trackpoint.alt[origin][key]

    def global_rad(self):
        """Computes the radius of the earth for a given latitude and elevation
            See Wikipedia for eqns"""
        _eqt_rad = 6378137.0 #meters
        _pole_rad = 6356752.3 #meters
        _lat_rad = self.lat * m.pi/180
        _eqn_top = ((_eqt_rad**2) * m.cos(_lat_rad))**2 + ((_pole_rad**2) * m.sin(_pole_rad))**2
        _eqn_bot = ((_eqt_rad) * m.cos(_lat_rad))**2 + ((_pole_rad) * m.sin(_pole_rad))**2
        _radius = (_eqn_top/_eqn_bot)**0.5 # + self.alt
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

    def smooth_value(self, param, length, orig = None, iskey='p'):
        """Creates a moving average of the value for the length parameter"""
        #first find begin and end of values
        #list_gap_check should close gaps beforehand
        existed = False
        begin = None
        end = None
        for i in range(0, self.tplist_len-1):
            try:
                last_succesful_i = i
                if orig == None:
                    # will fail if no orig
                    pos_orig = getattr(self.tplist[i], param).keys()
                    for poso in pos_orig:
                        orig = poso #will use random data origin
                getattr(self.tplist[i], param)[orig][iskey]
                if begin == None:
                    begin = i
            except:
                pass
        if begin == None:
            print('somethings not right')
            return
        end = last_succesful_i
        if length > (end - begin):
            length = (end - begin) - 1
        #blen must be odd, elen must be even 
        blen = (length // 2) +1 
        elen = (length // 2)
        tlen = blen + elen
        keymin = ''
        keyout = iskey + '_' + 'sm' + str(length) 
        for i in range(begin, end):
            _tpcount = self.tplist[i].tpcount
            _tpcount_rev = self.tplist[i].tpcount_rev
            psum = 0
            count = 0
            err_count = 0
            if _tpcount < blen:
                lb = begin
                ub = i + elen
            elif _tpcount_rev < elen:
                lb = i - blen
                ub = (end - begin)
            else:
                lb = i - blen
                ub = i + elen
            for j in range(lb,ub):
                try:    
                    raw_val = getattr(self.tplist[j], param)[orig][iskey]
                except:
                    print('value does not exist, thought it was already checked?')
                    err_count += 1
                    if err_count > (elen-1):
                        return
                else:       
                    psum += raw_val
                    count += 1
            up_val = {keyout:(psum/count)}
            dict_val = getattr(self.tplist[j], param)[orig]
            dict_val.update(up_val)
            getattr(self.tplist[i], param)[orig].update(dict_val)
            # print(getattr(self.tplist[i], param)[orig])

    def smooth_value_old(self, param, length, direc = None):
        """Creates a moving average of the value for the length parameter"""
        #direc is used for obtaining key as well as direction, so the none 
        #value is needed
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
        keymin = ''
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
                                   keymin = ''
                                   break
                                keyv = int(re.sub("\D", "", key))
                                if abs(keyv) < keymin:
                                    valout = val[key]
                                    keymin = abs(keyv)
                        val = valout
                    psum += val
                    count += 1
            keyout = 'sm' + str(length) + str(keymin)
            dict_val = init_vals
            up_val = {keyout:(psum/count)}
            dict_val.update(up_val)
            setattr(self.tplist[i], param, dict_val)


    def create_tpdv_list(self, direc = 1, iskey = None, debug = False):
        """Updates tplist with the derivative terms"""
        if direc <= 0:
            begin = direc * -1
            end = self.tplist_len
            for_range = range(begin, end)
        else:
            begin = 0
            end = self.tplist_len - direc
            for_range = range((end -1), 0, -1)
        for i in for_range:
            self.tplist[i].calc_tp_der(self.tplist[i + direc], direc, iskey, debug)

    def list_gap_check(self, gap = 1):
        """Interpolates across gaps in data, will only interpolate across
        a gap that is equal to or smaller than 'gap'.  If value is nonexistent
        it tags it with 'nonex'.  It is intended to be used with data in it's
        raw form, not derived data.
        """
        #check initial values
        for key in self.tplist[0].tpchan:
            if self.tplist[0].tpchan[key] == 'absent':
                self.tplist[0].tpchan[key] == 'nonex'
        gaps = []
        data_origin = self.tplist[0].tp_dict_origin()
        #data_origin_val = data
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
                                strval_dict = getattr(self.tplist[i-1], key)
                                strval = strval_dict[data_origin[key]]['p']
                                endval_dict = getattr(self.tplist[i+j+1], key)
                                endval = endval_dict[data_origin[key]]['p']
                                #note that interpolation is by time, not by distance
                                strdiv = self.tplist[i-1].time_st
                                enddiv = self.tplist[i+j].time_st
                                slope = ((endval - strval) / (enddiv - strdiv)) / (j + 2)
                                #get the data origin, should be 'gpsm'
                                d2k = next (iter (strval_dict.keys()))
                                #assign values to gap
                                for k in range(j+1):
                                    val = slope * (self.tplist[i+k].time_st - strdiv) + strval
                                    d1 = {'p':val}
                                    d2 = {d2k:d1}
                                    setattr(self.tplist[i+k], key, d2)
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


        
				
