import xml.etree.ElementTree as ET
from xml.etree.ElementTree import ElementTree
import re
from datetime import datetime
import math as m
from collections import OrderedDict
import matplotlib.pyplot as plt
import scipy.stats as scist
import numpy as np

def startup():
    john = athlete(78)
    john.set_bike_mass(9)
    john.set_roll_res(0.004)
    john.set_total_mass()
    # event = data_event('steelman_dorney.tcx', john)
    event = data_event('eynsham_cycle_2013.tcx', john)
    event.create_tp_list()
    event.gap_check_list(20)
    event.smooth_value('alt', 200)
    event.create_tpdv_list(-1)
    event.create_tpdv_list(1)
    event.create_tpi_list('acc')
    event.create_tpi_list('vect_step')
    event.update_smooth_chan()
    # event.env0 = event.environment()
    # event.average_drag(event.env0, debug=True)
    return(event)

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

    def set_total_mass(self):
        self.total_mass = self.bike_mass + self.mass

class trackpoint:
    """This class deals with sample points from a .tcx file,
        referred to as trackpoints in the .tcx file"""
    def __init__(self, athlete, _raw_tp, _id, _tpcount, _lapnumber):
        self.athlete = athlete
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
                self.tpchan[key] = {_tp_dict_origin[key]:['p']}
            else:
                self.tpchan[key] = {_tp_dict_origin[key]:['absent']}

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

    def assign_init_def_keys(self, data_key, debug=False):
        """Use data_event.data_key dict to assign default value
        key attributes"""
        for chan in data_key:
            schan = chan.replace('_def','')
            (orig, key) = data_key[chan]
            try:
                setattr(self, chan, getattr(self, schan)[orig][key])
            except:
                if debug != False:
                    print('missing channel', schan, 'at point ', self.tpcount)
                pass

    def set_def_key(self, chan, data_key):
        """Use data_event.data_key dict to assign default value
        key attributes, data_key for event must be set """
        chan_def = chan + '_def'
        dkey = data_key[chan_def]
        setattr(self, chan_def, [dkey])

    
    #Notation is self.dist['vsen']['p'] would be a value
    #Delta would be self.dist_step['vsen']['p_d-1']
    #Derivative would be self.speed['vsen']['p_d-1']
    #Delta would be self.speed_step['vsen']['p_d-1_2d1']
    #Derivative would be self.acc['vsen']['p_d-1_2d1']

    def deriv_dkey(self, initkey, direc):
        if initkey == 'p':
            dkey = 'p_d' + str(direc)
        else:
            if ('d' in initkey) == False:
                dkey = initkey + '_d' + str(direc)
            elif ('2d' in initkey) == False:
                dkey = initkey + '_2d' + str(direc)
            else:
                print('problem with keys')
                return
        return(dkey)

    def calc_tp_der(self, _data_event, direc = -1, iskey=None, debug=False):
        """Derivative calculator
        direc is either -1 for previous or +1 for next
        """
        if debug == True:
            print('\n\n', self.tpcount)
        _dtp = _data_event.tplist[self.tpcount + direc]
        data_key = _data_event.data_key
        #check if end point, if so end function
        if self.tpcount + direc < 0:
            return
        if self.tpcount_rev + direc < 0:
            return
        if iskey == None:
            iskey = 'p'
        pt1_notation = 'p_d' + str(direc) #used for lat/lon coordinates
        # <-- may not be needed with default types?
        if ('d' in iskey) == False:
            pt_notation = iskey + '_d' + str(direc)
        else:
            pt_notation = iskey + '_2d' + str(direc)
        #calculate gps distances, but first see what exists
        (geo_orig, geo_key) = data_key['lon_def']
        gp = global_point(self, geo_orig, geo_key)
        gp_d = global_point(_dtp, geo_orig, geo_key)
        _dist_geo = gp.dist_geo(gp_d)
        _vect_geo_wd = gp.vect_geo(gp_d)
        #Set direction of vector based upon order of direc
        _vect_geo = tuple(i*(-1*direc) for i in _vect_geo_wd)
        try:
            is_gpsm = self.dist_step.keys()
        except:
            #no dist_step exists
            self.dist_step = {geo_orig:{pt1_notation:_dist_geo}}
            self.vect_step = {geo_orig:{pt1_notation:_vect_geo}}
            # assign defaults
            self.vect_step_def = _vect_geo
            data_key['vect_step_def'] = (geo_orig, pt1_notation)
        else:
            if geo_orig in is_gpsm:
                # will either append ['gpsm'] or overwrite existing 'pt1_note'
                self.dist_step[geo_orig].update({pt1_notation:_dist_geo})
                self.vect_step[geo_orig].update({pt1_notation:_vect_geo})
                # assign defaults
                self.vect_step_def = _vect_geo
                data_key['vect_step_def'] = (geo_orig, pt1_notation)
            else:
                # This is an unlikely case, but for completeness
                self.dist_step.update({geo_orig:{pt1_notation:_dist_geo}})
                self.vect_step.update({geo_orig:{pt1_notation:_vect_geo}})
                # Does not assing defaults
        #create time steps
        self.time_st_step = self.time_st*(-1*direc) + _dtp.time_st*direc
        #create steps of values: dist_step, speed_step, alt_step
        _step_dict = self._step_dictf()
        for d_step in _step_dict:
            # If d_step was 'dist_step' then dict would return dist
            # though dist contains two values for vsen & gpsm
            # example below
            # dist = self._step_dict['dist_step']
            dict_of_chan = getattr(self, _step_dict[d_step])
            dict_of_chan_d = getattr(_dtp, _step_dict[d_step])
            for orig in dict_of_chan:
                # orig will be something like 'vsen' which = ('p': value)
                for ckey in dict_of_chan[orig].keys():
                    try:
                        # check existence of value on comparison tp
                        dict_of_chan_d[orig][ckey]
                    except:
                        if (debug == True) :
                            print('comparison value for ', self.tpcount, 
                                d_step, dict_of_chan_d, orig, iskey, 
                                'does not exist')
                        continue
                    if (debug == True) :
                        print('comparison value for ', self.tpcount, 
                                d_step, dict_of_chan_d, orig, ckey, ' exists') 
                    try:
                        getattr(self, d_step)[orig][ckey]
                        #will pass if already exists
                    except:
                        # generate deriv_key
                        dkey = self.deriv_dkey(ckey, direc)
                        # example
                        # delta_left = dist['vsen']['p']*(-1*-1)
                        delta_left = dict_of_chan[orig][ckey] * (-1 * direc)
                        delta_right = dict_of_chan_d[orig][ckey] * direc
                        delta = delta_left + delta_right
                        dict_delta_add = {dkey:delta}
                        try:
                            #check if method exists from another iskey
                            #example: self.acc['vsen']
                            d_step_dict = getattr(self, d_step)[orig]
                        except:
                            # orig nonexistent, create empty dict for step values
                            value_dict = {}
                            # example
                            # full_add = {'vsen':{'p_d-1': value}}
                            full_add = {orig:dict_delta_add}
                        else:
                            # orig exists, value_dict contains a previous iskey
                            d_step_dict.update(dict_delta_add)
                            full_add = {orig:d_step_dict}
                            # This value_dict will overwrite previous values
                            value_dict = getattr(self, d_step)
                        # This value_dict will overwrite previous values
                        value_dict.update(full_add)
                        try:
                            # update trackpoint.attribute
                            # self.dist_step.update(value_dict)
                            getattr(self, d_step).update(value_dict) #
                        except:
                            # create trackpoint.attribute
                            setattr(self, d_step, value_dict)
                        try:
                            # check if def_d_step exists
                            def_d_step = d_step + '_def'
                            getattr(self, def_d_step)
                            # to assign def_d_step, use assign
                        except:
                            #no def_d_step, assign it
                            try:
                                # does data key value exist?
                                dkey = data_key[def_d_step]
                            except:
                                # no, then set this one
                                data_key[def_d_step] = [orig, pt_notation]
                                dkey = [orig, pt_notation]
                            setattr(self, def_d_step, delta)
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
                    for ckey in step_type[orig].keys():
                        try:
                            #check if attribute already exists for this exact combo
                            getattr(self, deriv)[orig][ckey]
                            #will pass if already exists
                        except:
                            # orig and iskey combo are new
                            # numerator will be the step
                            # direc cancels out, as the step was calculated w/ direc
                            # try is needed to check if step values exist
                            # example:
                            # numerator = self.speed_step['vsen']['p_d-1']
                            try:
                                numerator = (getattr(self, deriv_type_dict[deriv])
                                            [orig][ckey])
                            except:
                                if debug == True:
                                    print('step value for ', self.tpcount, 
                                            deriv_type_dict[deriv], orig, 
                                            ckey, 'does not exist')
                                    print('THIS SHOULD NEVER APPEAR?')
                                continue
                            if deriv_type == 'time_st_step':
                                #the time step has no dictionary parameter
                                denom = getattr(self, deriv_type)
                            else:
                                #Currently distance is the only other parameter
                                try:
                                    denom = getattr(self, deriv_type)[orig][ckey]
                                except:
                                    #try non smoothed denominator
                                    # temp_notation = re.sub("_sm\d*_", "_", pt_notation)
                                    # denom = getattr(self, deriv_type)[orig][temp_notation]
                                    # try default denominator
                                    deriv_type_def = deriv_type + '_def'
                                    denom = getattr(self, deriv_type_def)
                            # print('deriv =', deriv, 'orig = ', orig, 'numerator = ',
                            #       numerator,'denom = ', denom)
                            if denom !=0:
                                delta = numerator / denom
                            else:
                                if debug == True:
                                    print('denominator is 0 for ', self.tpcount, 
                                            deriv_type_dict[deriv], orig, ckey)
                                continue
                            dict_delta_add = {ckey:delta}
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
                            try:
                                # check if def_d_step exists
                                def_deriv = deriv + '_def'
                                getattr(self, def_deriv)
                                # to assign def_d_step, use assign
                            except:
                                #no def_d_step, assign it
                                try:
                                # does data key value exist?
                                    dkey = data_key[deriv]
                                except:
                                # no, then set this one
                                    data_key[def_deriv] = [orig, pt_notation]
                                    dkey = [orig, pt_notation]
                                setattr(self, def_deriv, delta)
        return(data_key)


    def calc_tp_interp( self, chan, def_key, orig=None, new_def=True):
        """Calculates a new interpolated value from two iskeys"""
        # at some point a fancy iskey namer will be invented
        # the function could also run through all channels?
        iskey = 'p_di'
        iskey1='p_d-1'
        iskey2='p_d1'
        avail_keys = []
        (orig_def, key_def) = def_key
        chan_def = chan + '_def'
        if orig == None:
            orig_keys = getattr(self, chan).keys()
            for orig in orig_keys:
                p_keys = getattr(self, chan)[orig].keys()
                if (iskey1 in p_keys) and (iskey2 in p_keys):
                    avail_keys.append(orig)
        else:
            p_keys = getattr(self, chan)[orig].keys()
            if (iskey1 in p_keys) and (iskey2 in p_keys):
                avail_keys = [orig]
        if avail_keys == []:
            print('No suitable keys for', chan, orig)
            return(False)
        for orig in avail_keys:
            before = getattr(self, chan)[orig][iskey1]
            after = getattr(self, chan)[orig][iskey2]
            if isinstance(before, tuple):
                (before1, before2) = before
                (after1, after2) = after
                val1 = (before1 + after1)/2
                val2 = (before2 + after2)/2
                getattr(self, chan)[orig].update({iskey:(val1,val2)})
                if (orig == orig_def) and (new_def == True):
                    setattr(self, chan_def, (val1, val2))
            else:
                val = (before + after)/2
                getattr(self, chan)[orig].update({iskey:val})
                if (orig == orig_def) and (new_def == True):
                    setattr(self, chan_def, val)
        return(True)

    def calc_tp_interp_old( self, chan, orig=None, new_def=True, iskey1='p_d-1', iskey2='p_d1'):
        """DO NOT DELETE...yet - More advanced than newer version"""
        """Calculates a new interpolated value from two iskeys"""
        # at some point a fancy iskey namer will be invented
        # the function could also run through all channels?
        iskey = 'p_di'
        avail_keys = []
        chan_def = chan + '_def'
        if orig == None:
            orig_keys = getattr(self, chan).keys()
            for orig in orig_keys:
                p_keys = getattr(self, chan)[orig].keys()
                if (iskey1 in p_keys) and (iskey2 in p_keys):
                    avail_keys.append(orig)
        else:
            p_keys = getattr(self, chan)[orig].keys()
            if (iskey1 in p_keys) and (iskey2 in p_keys):
                avail_keys = [orig]
        if avail_keys == []:
            print('No suitable keys for', chan, orig)
            return
        for orig in avail_keys:
            before = getattr(self, chan)[orig][iskey1]
            after = getattr(self, chan)[orig][iskey2]
            if isinstance(before, tuple):
                (before1, before2) = before
                (after1, after2) = after
                val1 = (before1 + after1)/2
                val2 = (before2 + after2)/2
                getattr(self, chan)[orig].update({iskey:(val1,val2)})
            else:
                val = (before + after)/2
                getattr(self, chan)[orig].update({iskey:val})

    def aero_dragF_calc(self, cd, env):
        """Calculates the drag force for an assigned drag coefficient"""
        airspeed = self.apparent_wind(env) # headwind
        dragF = 0.5* env.rho * cd * airspeed **2
        return(dragF)

    def aero_dragF_est(self, env, disp = False):
        """Calculates the drag by summing the drive force, inertial force,
        rolling resitance, and gravitational force.  Returns a force"""
        gf = self.__grav_force__(env)
        rrf = self.__rr_force__(env)
        accf = self.__acc_force__()
        df = self.__drive_force__()
        _aero_dragF_est = df + accf + rrf + gf
        if disp != False:
            print('grav f=',round(gf,2),' rr f=',round(rrf,2),
                    ' acc f=',round(accf,2),' df=',round(df,2))
        return(_aero_dragF_est)

    def __grav_force__(self, env):
        """Calculates the force acting on the rider due to gravity tangential
        to the road"""
        gradient = self.gradient_def
        slope_comp = m.sin(m.atan(self.gradient_def))
        _grav_force = -env.g_alt * self.athlete.total_mass *slope_comp
        return(_grav_force)

    def __rr_force__(self, env):
        """Calculates the force from rolling resistance using the force
        normal to the road and the assigned rolling resistance coefficient"""
        gradient = self.gradient_def
        norm_comp = m.cos(m.atan(self.gradient_def))
        norm_force = env.g_rr*self.athlete.total_mass*norm_comp
        _rr_force = norm_force*self.athlete.bike_roll_res*-1
        return(_rr_force)

    def __acc_force__(self):
        """Calculates the inertial force from acceleration"""
        _acc_force = -1*self.athlete.total_mass*self.acc_def
        return(_acc_force)

    def __drive_force__(self):
        """Calculates the drive force from the wheel, by dividing the power by
        the road speed.  Allows for different origins for power and velocity as
        well as different data sets"""
        power = self.power_def
        speed = self.speed_def
        if speed <= 0:
            return('speed <= 0')
        _drive_force = power / speed
        return(_drive_force)

    def aero_dragC_est(self, env):
        """Calculates a drag coefficient using drag force and airspeed"""
        dragF = self.aero_dragF_est(env)
        airspeed = self.apparent_wind(env)
        #dragF = .5*rho*cd*speed**2
        if airspeed > 0:
            coeff_d = (2*dragF)/(env.rho*(airspeed**2))
        else:
            coeff_d = 1
            print('airspeed < 0')
        return(coeff_d)

    def apparent_wind(self, env, xwind=False):
        """Determines the head wind and cross wind components for the rider
        returning a tuple with head wind first.  Wind direction is specified
        as a heading with 0° being North"""
        (bk_head_ew, bk_head_ns) = self.vect_step_def
        wind_ew = env.wind_speed * m.sin(env.wind_dir * m.pi/180)
        wind_ns = env.wind_speed * m.cos(env.wind_dir * m.pi/180)
        road_speed_ew = self.speed_def * bk_head_ew
        road_speed_ns = self.speed_def * bk_head_ns
        app_wind_ew = road_speed_ew + wind_ew
        app_wind_ns = road_speed_ns + wind_ns
        app_wind = m.sqrt(app_wind_ns**2 + app_wind_ew**2)
        app_heading_rad = comp_to_heading(app_wind_ew, app_wind_ns, rad=True)
        bk_head_rad = comp_to_heading(bk_head_ew, bk_head_ns, rad=True)
        app_wind_rad = bk_head_rad - app_heading_rad
        headwind = app_wind * m.cos(app_wind_rad)
        crosswind = app_wind * m.sin(app_wind_rad)
        if xwind != False:
            return(headwind, crosswind)
        return(headwind)

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
        _eqn_top = (((_eqt_rad**2) * m.cos(_lat_rad))**2 + 
                    ((_pole_rad**2) * m.sin(_pole_rad))**2)
        _eqn_bot = (((_eqt_rad) * m.cos(_lat_rad))**2 +
                    ((_pole_rad) * m.sin(_pole_rad))**2)
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
        _right_eqn = (m.cos(_lat1_rads)*m.cos(_lat2_rads)*
                    (m.sin((_lon1_rads-_lon2_rads)/2))**2)
        _dist = self.global_rad() * 2 * m.asin(m.sqrt(_left_eqn + _right_eqn))
        return(_dist)

    def vect_geo(self, _gpoint2, heading=False):
        """Calculates the vector between two geo points
        and returns it as a tuple with components (East-West, North-South)
        alternatively if heading is set to True, output will be a heading
        in degrees"""
        _lat1_rads = self.lat * m.pi/180
        _lon1_rads = self.lon * m.pi/180
        _lat2_rads = _gpoint2.lat * m.pi/180
        _lon2_rads = _gpoint2.lon * m.pi/180
        _radius = self.global_rad()
        #Note that direction is lost with the ()**2
        _dist_ns = (_radius * 2 * 
            m.asin(m.sqrt((m.sin((_lat1_rads-_lat2_rads)/2))**2)))
        _dist_ew = (_radius * 2 * m.asin(m.sqrt(m.cos(_lat1_rads)*
            m.cos(_lat2_rads)*(m.sin((_lon1_rads-_lon2_rads)/2))**2)))
        _dist_from_vect = m.sqrt(_dist_ns**2 + _dist_ew**2)
        if _dist_from_vect == 0:
            return((0,0))
        #Find the direction
        ns_diff = self.lat - _gpoint2.lat
        if ns_diff != 0:
            north_is_pos = ns_diff / m.fabs(ns_diff)
        else:
            north_is_pos = 0
        ew_diff = self.lon - _gpoint2.lon
        if ew_diff != 0:
            east_is_pos = ew_diff / m.fabs(ew_diff)
        else:
            east_is_pos = 0
        # print(_dist_from_vect) # for error check
        _vect_ns = north_is_pos * _dist_ns / _dist_from_vect
        _vect_ew = east_is_pos * _dist_ew / _dist_from_vect
        if heading != True:
            return((_vect_ew, _vect_ns))
        else:            
            heading_deg = comp_to_heading(_vect_ew, _vect_ns, rad=False)
            return(heading_deg)


class data_event(ElementTree):
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

    def __init_data_key_dict__(self, tpnum=0):
        """defines default data keys for data_event, requires tplist
        uses the tp_dict_origin function from trackpoint class"""
        _tp_dict_origin = self.tplist[tpnum].tp_dict_origin()
        data_key_dict = {}
        for chan_name in _tp_dict_origin:
            str_name = chan_name + '_def'
            val_name = [_tp_dict_origin[chan_name], 'p']
            data_key_dict[str_name] = val_name
        self.data_key = data_key_dict

    def update_data_key(self, chan=None, orig=None, key=None, 
                                debug = True):
        if chan and orig and key == None:
            pass
        elif (chan==None) or (orig==None) or (key== None):
            print("'chan' 'orig' and 'key' must all have values to make a change,"
                    " leave blank to update tplist")
            return
        else:
            str_name = chan + '_def'
            self.data_key[str_name] = [orig, key]
        for i in range(self.tplist_len):
            self.tplist[i].assign_init_def_keys(self.data_key, debug)

    def update_smooth_chan(self, chan=None, smooth=None, 
                        def_big=True, def_neg=True):
        """Updates data_key and tplist to largest channel unless def_big is
        specified, individual channels, smoothnesses can be specified"""
        midpoint = self.tplist_len // 2
        #keys are two levels deep
        attribs = self.tplist[midpoint].__dict__.items()
        key_list = []
        val_list = []
        for attrib in attribs:
            key_list.append(attrib[0])
            val_list.append(attrib[1])
        # remove tpchan
        loc = key_list.index('tpchan')
        key_list.pop(loc)
        val_list.pop(loc)
        chan_list = []
        for i in range(len(key_list)):
            if isinstance(val_list[i], dict):
                for k in val_list[i].keys():
                    klist = val_list[i][k]
                    for m in klist.keys():
                        if re.search('sm', m):
                            # get value of smoothness
                            sm_val = int(re.findall('_sm\d+', m)[0][3:])
                            if re.search('d', m):                              
                                d_val = int(re.findall('_d-*\d+', m)[0][2:])
                            else:
                                d_val = 0
                            # chan = chan name, orig, key
                            chan = [key_list[i], k, m, sm_val, d_val] # chan name
                            chan_list.append(chan)
        # prioritise smooth channels
        new_chan_list = [chan_list[0]]
        for chan in chan_list:
            n_chan_name = []
            for n_chan in new_chan_list:
                n_chan_name.append(n_chan[0])
            if (chan[0] in n_chan_name) == False:
                new_chan_list.append(chan)
                continue
            for i in range(len(new_chan_list)):
                if chan[0] in new_chan_list[i]:
                    # which one stays?
                    if chan[3] == new_chan_list[i][3]:
                        if chan[4] == -1:
                            if def_neg == True:
                                #replace new_chan
                                new_chan_list[i] = chan
                                break
                            else:
                                break
                    elif chan[3] > new_chan_list[i][3]:
                        if def_big == True:
                            # replace new_chan
                            new_chan_list[i] = chan
                            break
                        else:
                            break
                    else:
                        if def_big == True:
                            break
                        else:
                            new_chan_list[i] = chan
                            break
                else:
                    continue
        for chan in new_chan_list:
            self.update_data_key(chan[0], chan[1], chan[2], debug=False)
        

    def create_tp_list(self):
        """Create a list of the trackpoints and creates the 
        method tplist, tplist_len, and sets default data types"""
        _lapnumber = 1
        _tp_sum = 0
        _trackpoint_list = []
        for _lap in self.laps:
            _trackpoints_raw = _lap.findall('.//Trackpoint')
            for _raw_tp in _trackpoints_raw:
                _trackpoint_list.append(trackpoint(self.athlete, 
                        _raw_tp, self.id, _tp_sum, _lapnumber))
                _tp_sum += 1
            _lapnumber += 1
        for _trackpoint in _trackpoint_list:
            _trackpoint.reverse_count(_tp_sum)
        self.tplist = _trackpoint_list
        self.tplist_len = len(self.tplist)
        self.__init_data_key_dict__() #set default data types
        for _trackpoint in _trackpoint_list:
            # assign default data type for trackpoint
            _trackpoint.assign_init_def_keys(self.data_key)

    def smooth_value(self, chan, length, origkey=(None, None), new_def=True):
        """Creates a moving average of the value for the length specified"""
        #first find begin and end of values
        #list_gap_check should close gaps beforehand
        existed = False
        begin = None
        end = None
        (orig, iskey) = origkey
        if orig == None:
            chan_def = chan + '_def'
            def_keys = self.data_key[chan_def]
            orig = def_keys[0]
            iskey = def_keys[1]
        # Look for continuous data section
        for i in range(0, self.tplist_len-1):
            try:
                last_succesful_i = i
                getattr(self.tplist[i], chan)[orig][iskey]
                if begin == None:
                    begin = i
            except:
                pass
        if begin == None:
            print('No origin value was found')
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
        if new_def == True:
            # assign new value to chan_def
            # this is essential for __set_def_key__
            def_keys_n = [orig, keyout]
            self.data_key[chan_def] = def_keys_n
        for i in range(begin, end):
            _tpcount = self.tplist[i].tpcount
            _tpcount_rev = self.tplist[i].tpcount_rev
            psum = 0
            count = 0
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
                raw_val = getattr(self.tplist[j], chan)[orig][iskey]
                psum += raw_val
                count += 1
            up_val = {keyout:(psum/count)}
            dict_val = getattr(self.tplist[j], chan)[orig]
            dict_val.update(up_val)
            # create smoothed channel
            getattr(self.tplist[i], chan)[orig].update(dict_val)
            if new_def == True:
                self.tplist[i].set_def_key(chan, self.data_key)
            # print(getattr(self.tplist[i], chan)[orig])

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
            p_data_key = self.tplist[i].calc_tp_der(self, direc, 
                                        iskey, debug)
            if p_data_key == None:
                continue
            if len(p_data_key) > len(self.data_key):
                self.data_key = p_data_key

    def create_tpi_list(self, chan, orig=None, new_def=True):
        """Interpolates between direc -1 and +1 values, 
        writes a new channel p_di"""
        for_range = range(1, self.tplist_len-1)
        chan_def = chan + '_def'
        def_key = self.data_key[chan_def]
        (korig, junk) = def_key
        true_success = False
        for i in for_range:
            success = self.tplist[i].calc_tp_interp(chan, def_key, orig,
                                    new_def=True)
            if success == True:     # determines whether it has interpolated
                true_success = True # at least once, if so sets as default
        if true_success == True:
            self.data_key[chan_def] = [korig, 'p_di']

    def average_drag_old(self, env, lspeed_bound = 8,
                    bound_pc=20, debug=False):
        wind_speed = env.wind_speed
        wind_dir = env.wind_dir
        rho = env.rho
        lower_pc = bound_pc
        upper_pc = 100 - bound_pc
        wdrag_array = []
        tp_used = []
        for i in range(400,self.tplist_len-1):
        # for i in range(325, 350):
            try:
                speed = self.tplist[i].speed_def
                print(i)
                coeff_d = self.tplist[i].aero_dragC_est(env)
                print(speed, coeff_d)
                if speed > lspeed_bound:
                    tp_used.append(self.tplist[i].tpcount)
                    #for j in range(int((speed**2)//1)):
                    wdrag_array.append(coeff_d)
            except:
                if debug == True:
                    print('point ', i, ' is incomplete')
        wlow_bound = scist.scoreatpercentile(wdrag_array, lower_pc)
        whigh_bound = scist.scoreatpercentile(wdrag_array, upper_pc)
        wmean_bd = scist.tmean(wdrag_array, limits=(wlow_bound, whigh_bound))
        print('bounded = ', wlow_bound, wmean_bd, whigh_bound)
        coeff_d1 = wmean_bd
        wdrag_array = []
        tp_used2 = []
        x = []
        y = []
        #calculate wind
        wind_segments = 8
        segment_div = 360/wind_segments
        segment_bound = np.linspace(0,360 -segment_div, wind_segments)
        segment_mid = []
        for i in segment_bound:
            segment_mid.append(i + segment_div / 2) # segment midpoint
        wind_section = OrderedDict()
        for i in range(wind_segments):
            wind_section[segment_mid[i]] = []
        for i in tp_used:
            try:
                airspeed = self.tplist[i].apparent_wind(env)
                dragFE = self.tplist[i].aero_dragF_est(env)#disp=True)
            except:
                if debug == True:
                    print('2nd try point ', i, ' is incomplete')
            else:
                if dragFE > 0:
                    headwind = (m.sqrt(dragFE / (.5*rho*coeff_d1))) - airspeed
                else:
                    headwind = -1 * (m.sqrt(-dragFE / (.5*rho*coeff_d1))) - airspeed
                (hx, hy) = self.tplist[i].vect_step_def
                bk_head_deg = comp_to_heading(hx, hy, rad=False)
                for j in range(wind_segments):
                    # if 23deg < (0 + 45deg)
                    if bk_head_deg < (segment_bound[j] + segment_div):
                        wind_section[segment_mid[j]].append(headwind)
                        # print(i, segment_mid[j], headwind, bk_head_deg)
                        break
        wind_section_avg = OrderedDict()
        for i in wind_section:
            low_bound = scist.scoreatpercentile(wind_section[i], lower_pc)
            high_bound = scist.scoreatpercentile(wind_section[i], upper_pc)
            wind_section_avg[i] = scist.tmean(wind_section[i], 
                        limits=(low_bound, high_bound))
        print(wind_section_avg)
        x1 = [0]
        y1 = [0]
        wind_speed_dir = {}
        for i in wind_section_avg:
            wind_comp_x = wind_section_avg[i]*m.sin(i*m.pi/180)
            wind_comp_y = wind_section_avg[i]*m.cos(i*m.pi/180)
            j = heading_complement(i)
            wind_speed_dir[i] = (wind_section_avg[i] - wind_section_avg[j])/2
            x1.append(wind_comp_x)
            y1.append(wind_comp_y)
            x1.append(0)
            y1.append(0)
        wind_speed_dir = OrderedDict(sorted(wind_speed_dir.items(), key=lambda k: k[1]))
        print(wind_speed_dir)
        splot_xy(x1, y1)
        # return((wind_head_high, wind_head_low, xmean_high, ymean_high, xmean_low, ymean_low))

    def average_drag(self, env, lspeed_bound = 8,
                    bound_pc=20, debug=False):
        wind_speed = env.wind_speed
        wind_dir = env.wind_dir
        rho = env.rho
        lower_pc = bound_pc
        upper_pc = 100 - bound_pc
        wdrag_array = []
        tp_used = []
        for i in range(400,self.tplist_len-1):
        # for i in range(325, 350):
            try:
                speed = self.tplist[i].speed_def
                coeff_d = self.tplist[i].aero_dragC_est(env)
                if speed > lspeed_bound:
                    tp_used.append(self.tplist[i].tpcount)
                    #for j in range(int((speed**2)//1)):
                    wdrag_array.append(coeff_d)
            except:
                if debug == True:
                    print('point ', i, ' is incomplete')
        wlow_bound = scist.scoreatpercentile(wdrag_array, lower_pc)
        whigh_bound = scist.scoreatpercentile(wdrag_array, upper_pc)
        wmean_bd = scist.tmean(wdrag_array, limits=(wlow_bound, whigh_bound))
        print('bounded = ', wlow_bound, wmean_bd, whigh_bound)
        return(wdrag_array)


    def gap_check_list(self, gap = 1):
        """Interpolates across gaps in data, will only interpolate across
        a gap that is equal to or smaller than 'gap'.  If value is nonexistent
        it tags it with 'nonex'.  It is intended to be used with data in it's
        raw form, not derived data.
        """
        #check initial values
        data_origin = self.tplist[0].tp_dict_origin()
        for key in self.tplist[0].tpchan:
            for orig in self.tplist[0].tpchan[key]:
                if self.tplist[0].tpchan[key][orig][0] == 'absent':
                    self.tplist[0].tpchan[key][orig] = ['nonex']
        gaps = []
        
        #data_origin_val = data
        for i in range(gap,(self.tplist_len - gap)):
            for key in self.tplist[i].tpchan:
                for orig in self.tplist[i].tpchan[key]:
                    outp = self.tplist[i].tpchan[key][orig][0]
                    if outp == 'absent':
                        #if previous is nonexistent, continue non-existent
                        if self.tplist[i-1].tpchan[key][orig][0] == 'nonex':
                            self.tplist[i].tpchan[key][orig][0] = 'nonex'
                        elif self.tplist[i-1].tpchan[key][orig][0] == 'absent':
                            #then the gap was too large so...
                            self.tplist[i].tpchan[key][orig][0] = 'absent'
                        else:
                            for j in range(gap):
                                #if value exists on later tps
                                if self.tplist[i+j+1].tpchan[key][orig][0] == 'p':
                                    strval_dict = getattr(self.tplist[i-1], key)
                                    strval = strval_dict[data_origin[key]]['p']
                                    endval_dict = getattr(self.tplist[i+j+1], key)
                                    endval = endval_dict[data_origin[key]]['p']
                                    # note that interpolation is by time, 
                                    # not by distance
                                    strdiv = self.tplist[i-1].time_st
                                    enddiv = self.tplist[i+j].time_st
                                    slope = (((endval - strval)/(enddiv - strdiv))
                                            / (j + 2))
                                    #get the data origin, should be 'gpsm'
                                    d2k = next (iter (strval_dict.keys()))
                                    #assign values to gap
                                    for k in range(j+1):
                                        val = slope * (self.tplist[i+k].time_st
                                                     - strdiv) + strval
                                        d1 = {'p':val}
                                        d2 = {d2k:d1}
                                        setattr(self.tplist[i+k], key, d2)
                                        self.tplist[i+k].tpchan[key][orig] = 'interp'
                                        gaps.append((i+k, key, ' bridged'))
                                    break
                            if self.tplist[i].tpchan[key] != 'interp':
                                for j in range(gap):
                                    gaps.append((i+j,key, ' open'))
        return(gaps)

    class environment(object):
        def __init__(self):
            self.temp = 20
            self.lat = 51.5
            self.lon = 0
            self.alt = 0
            self.bar_p = 102400
            self.g = self.__grav_from_lat__()
            self.g_rr = self.g
            self.g_alt = self.g
            self.rho = self.__density_pvnrt__()
            self.wind_dir = 0
            self.wind_speed = 0

        def __grav_from_lat__(self):
          """Derived from Helmert's Eqn"""
          g = (9.8061999 - 0.0259296 * m.cos(2 * self.lat) 
            + 0.0000567 * (m.cos (2 * self.lat))**2)
          #alt_g = (9.780327*(1+ 0.0053024 * (m.sin (self.lat))**2 -
                    # 0.0000058 * (m.sin (2*self.lat))**2))
          #print(alt_g)
          return(g)

        def __density_pvnrt__(self):
            temp_k = self.temp + 273
            r_gas = 287.058
            rho = self.bar_p / (temp_k * r_gas)
            return(rho)

        def set_wind(self, wind_speed, wind_dir):
            self.wind_speed = wind_speed
            self.wind_dir = wind_dir

    def plot_xy(self, xchan, ychan, xorig='vsen', yorig='gpsm', 
                xkey='p', ykey='p'):
        xval = []
        yval = []
        for i in range(self.tplist_len):
            try:
                xv = getattr(self.tplist[i], xchan)[xorig][xkey]
                yv = getattr(self.tplist[i], ychan)[yorig][ykey]
            except:
                continue
            else:            
                xval.append(xv)
                yval.append(yv)
        plt.plot(xval, yval)
        plt.show()

def splot_xy( xchan, ychan):
    
        plt.plot(xchan, ychan)
        plt.show()


def time_parse_garmin(_time_raw):
    """Parses time from garmin time stamp"""
    _time = datetime.strptime(_time_raw, '%Y-%m-%dT%H:%M:%S.%fZ')
    return(_time)

def comp_to_heading(ew, ns, rad=False):
    """Returns a heading in degrees from components, if rad is
    not False, returns heading in radians"""
    if ns == 0:
        heading = m.pi / 2 # 90 deg
    else:
        heading =  m.atan(ew / ns)
    if ns < 0:
        heading += m.pi # add 180 deg due to atan error 
    if rad == False:
        heading = (180 / m.pi) * heading
        if heading < 0:
                heading = 360 + heading
    return(heading)

def heading_to_comp(heading, rad=False):
    """Converts a heading into components East-West, North-South
    if rad=True heading is in radians """
    if rad == False:
        heading = heading * m.pi/180
    ew = m.sin(heading)
    ns = m.cos(heading)
    return(ew, ns)
				
def heading_complement(heading, rad=False):
    """Gives 180 degrees from the current heading"""
    if rad == False:
        heading += 180
        if heading >= 360:
            heading -= 360
    else:
        heading += m.pi
        if heading >= 2*m.pi:
            heading -= 2*m.pi
    return(heading)