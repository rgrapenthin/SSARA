#! /usr/bin/env python3
###############################################################################
# isce2hdf5.py
#
#  Project:  Seamless SAR Archive
#  Purpose:  Create HDF5 interferogram product 
#  Author:   Scott Baker
#  Created:  April 2015
#
###############################################################################
#  Copyright (c) 2015, Scott Baker
#
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.
###############################################################################

import os
import sys
import glob
import argparse
import datetime
import xml.etree.ElementTree as ET

import numpy as np
import h5py

from iscesys.Parsers.FileParserFactory import createFileParser
import isce
import isceobj
import pickle
from mroipac.geolocate.Geolocate import Geolocate

def read_float32(infile,length,width):
  '''Reads roi_pac unw, cor, or hgt data.

  Requires the file path and returns amplitude and phase
  Usage:
    amplitude, phase, rscDictionary = readUnw('/Users/sbaker/Desktop/geo_070603-070721_0048_00018.unw')
  '''
  oddindices = np.where(np.arange(length*2)&1)[0]
  data = np.fromfile(infile,np.float32,length*2*width).reshape(length*2,width)
  a = np.array([data.take(oddindices-1,axis=0)]).reshape(length,width)
  p = np.array([data.take(oddindices,axis=0)]).reshape(length,width)
  return a, p

def read_complex64(infile,length,width):
  '''Reads roi_pac int or slc data.

  Requires the file path and returns amplitude and phase
  Usage:
    amp, phase, rscDictionary = readInt('/Users/sbaker/Desktop/geo_070603-070721_0048_00018.int')
  '''
  data = np.fromfile(infile,np.complex64,length*2*width).reshape(length,width)
  a = np.array([np.hypot(data.real,data.imag)]).reshape(length,width)
  p = np.array([np.arctan2(data.imag,data.real)]).reshape(length,width)
  return a, p

def read_dem(infile,length,width):
  '''Read a roipac dem file.

  Input:
    roi_pac format dem file
  '''
  d=np.fromfile(infile,dtype=np.int16).reshape(length,width)
  return d

def footprintFromPickle():
    insar = pickle.load(open('PICKLE/preprocess','rb'))
    planet = insar.masterFrame.getInstrument().getPlatform().getPlanet()
    earlySquint = insar.masterFrame._squintAngle
    orbit = insar.masterFrame.getOrbit()
    lookSide = int(insar.masterFrame.getInstrument().getPlatform().pointingDirection)
    geolocate = Geolocate()
    geolocate.wireInputPort(name='planet',object=planet)
    nearRange = insar.masterFrame.getStartingRange()
    farRange = insar.masterFrame.getFarRange()
    earlyStateVector = orbit.interpolateOrbit(insar.masterFrame.getSensingStart())
    lateStateVector = orbit.interpolateOrbit(insar.masterFrame.getSensingStop())
    nearEarlyCorner,nearEarlyLookAngle,nearEarlyIncAngle =geolocate.geolocate(earlyStateVector.getPosition(),earlyStateVector.getVelocity(),nearRange,earlySquint,lookSide)
    farEarlyCorner,farEarlyLookAngle,farEarlyIncAngle = geolocate.geolocate(earlyStateVector.getPosition(),earlyStateVector.getVelocity(),farRange,earlySquint,lookSide)
    nearLateCorner,nearLateLookAngle,nearLateIncAngle = geolocate.geolocate(lateStateVector.getPosition(),lateStateVector.getVelocity(),nearRange,earlySquint,lookSide)
    farLateCorner,farLateLookAngle,farLateIncAngle = geolocate.geolocate(lateStateVector.getPosition(),lateStateVector.getVelocity(),farRange,earlySquint,lookSide)
    wkt = "POLYGON((%f %f, %f %f, %f %f, %f %f, %f %f))" % ( nearEarlyCorner.getLongitude(),nearEarlyCorner.getLatitude(),farEarlyCorner.getLongitude(),farEarlyCorner.getLatitude(),farLateCorner.getLongitude(),farLateCorner.getLatitude(), nearLateCorner.getLongitude(),nearLateCorner.getLatitude(),nearEarlyCorner.getLongitude(),nearEarlyCorner.getLatitude() )
    return wkt

def footprintFromLogFile():
    """WARNING:  MASSIVE KLUDGE AHEAD"""
    f = open('isce.log')
    log = f.readlines()
    f.close()
    lats = []
    lons = []
    for line_number, line in enumerate(log, 1):
        if 'contrib.frameUtils.FrameInfoExtractor' in line and 'Corner' in line:
            lats.append(log[line_number-1].strip().split(":")[-1])
            lons.append(log[line_number].strip().split(":")[-1])
    poly_lats = [lats[0],lats[1],lats[3],lats[2],lats[0]]
    poly_lons =  [lons[0],lons[1],lons[3],lons[2],lons[0]]
    wkt = "POLYGON((" + ",".join([lon+' '+lat for lat,lon in zip(poly_lats,poly_lons)]) + "))"
    return wkt

def parse():
    '''Command line parser.

    You can change/add defaults for any of these if you work a lot with the same mission for example.   You should
    also update the defaults for -software, -software_version, and -institution.
    '''
    parser = argparse.ArgumentParser(description='Create HDF5 interferogram product from ISCE output')
    ## REQUIRED METADATA FOR ARCHIVING PRODUCTS ##
#    parser.add_argument('-mission', dest='mission', action='store', help='Name of the mission', type=str, required=True)
#    parser.add_argument('-relative_orbit', dest='relative_orbit', action='store', help='Relative orbit/Track/Path number', type=int, required=True)
    parser.add_argument('-processing_type', dest='processing_type', action='store', help='Type of processing: INTERFEROGRAM, LOS_VELOCITY,...', type=str, default='INTERFEROGRAM')
#    parser.add_argument('-footprint', dest='scene_footprint', action='store', help='WKT Polygon for the area covered by the swath', type=str)
    parser.add_argument('-swath', dest='beam_swath', action='store', help='Swath name without underscores', type=str )
    ## RECOMMENDED METADATA ##
    parser.add_argument('-beam_mode', dest='beam_mode', action='store', help='', type=str)
    parser.add_argument('-frame', dest='frame', action='store', help='', type=int)
#    parser.add_argument('-flight', dest='flight_direction', action='store', help='', type=str, default='A')
#    parser.add_argument('-look', dest='look_direction', action='store', help='', type=str, default='R')
#    parser.add_argument('-polarization', dest='polarization', action='store', help='', type=str)
    parser.add_argument('-software', dest='processing_software', action='store', help='', type=str, default='ISCE')
    parser.add_argument('-software_version', dest='processing_software_version', action='store', help='', type=str,default='2.0.0')
    parser.add_argument('-atmos_correct_method', dest='processing_atmos_correct_method', action='store', help='', type=str)
    parser.add_argument('-institution', dest='processing_facility', action='store', help='', type=str, default='')
#    parser.add_argument('-master_platform', dest='master_platform', action='store', help='', type=str)
#    parser.add_argument('-master_orbit', dest='master_absolute_orbit', action='store', help='', type=int)
#    parser.add_argument('-slave_platform', dest='slave_platform', action='store', help='', type=str)
#    parser.add_argument('-slave_orbit', dest='slave_absolute_orbit', action='store', help='', type=int)
#    parser.add_argument('-', dest='', action='store', help='', type=str)
    clos = parser.parse_args()
    return clos

def main(argv):
    # GET THE COMMAND LINE OPTIONS
    clos = parse()

    ### READ GEOCODE DATASETS ###
    # these are hardwired in here, change if you have different naming conventions or want to include different products
    int_file = "filt_topophase.flat.geo"
    xmlfile = int_file+".xml"
    cor_file = "phsig.cor.geo"
    unw_file = "filt_topophase.flat.unw.geo"
    rdr_file = 'los.rdr.geo'

    PA = createFileParser('xml')
    dictOut,dictFact,dictMisc = PA.parse(xmlfile)
    width = dictOut['width']
    length = dictOut['length']
    xstep = dictOut['Coordinate1']['delta']
    ystep = dictOut['Coordinate2']['delta']
    north = dictOut['Coordinate2']['startingvalue']
    south = north + length*ystep
    west = dictOut['Coordinate1']['startingvalue']
    east = west + width*xstep

    inta,intp = read_complex64(int_file,length,width)
    unwa,unwp = read_float32(unw_file,length,width)
    rdra,rdrp = read_float32(rdr_file,length,width)
    corp = np.fromfile(cor_file,dtype=np.float32).reshape(length,width)

    #################################
    ###  METADATA
    #################################
    # we will grab most of the metadata from the insarProc.xml file
    tree = ET.parse('insarProc.xml')
    root = tree.getroot()
    first_date = datetime.datetime.strptime(root.find('master/frame/SENSING_START').text,'%Y-%m-%d %H:%M:%S.%f') 
    last_date = datetime.datetime.strptime(root.find('slave/frame/SENSING_START').text,'%Y-%m-%d %H:%M:%S.%f') 
    meta_dict = {}
    ## MANDATORY METADATA ##
    meta_dict['mission'] = root.find('master/platform/MISSION').text.replace("'","").replace("b","")
    meta_dict['beam_swath'] = clos.beam_swath 
    if root.find('master/frame/TRACK_NUMBER').text == 'None':
        if 'CSK' in meta_dict['mission']:
            if meta_dict['mission']=='CSKS4':
                meta_dict['relative_orbit']  = (int(root.find('master/frame/ORBIT_NUMBER').text)-193) % 237
            else:
                meta_dict['relative_orbit']  = int(root.find('master/frame/ORBIT_NUMBER').text) % 237
            meta_dict['mission'] = meta_dict['mission'][:3] # only take the CSK part
    else:
        meta_dict['relative_orbit'] = int(root.find('master/frame/TRACK_NUMBER').text)
    meta_dict['first_date'] = first_date.strftime("%Y%m%d") 
    meta_dict['last_date'] = last_date.strftime("%Y%m%d")
    meta_dict['processing_type'] = clos.processing_type # SET AS A DEFAULT IN parse()
#    meta_dict['scene_footprint'] = footprintFromPickle()
    meta_dict['scene_footprint'] = footprintFromLogFile()

    ## RECOMMENDED METADATA ##
    if clos.beam_mode:
        meta_dict['beam_mode'] = clos.beam_mode
    meta_dict['frame'] = 0000
    meta_dict['flight_direction'] = root.find('master/frame/PASS_DIRECTION').text.replace("'","").replace("b","")
    if int(root.find('master/lookSide').text) == -1:
        meta_dict['look_direction'] = 'R' 
    else:
        meta_dict['look_direction'] = 'L'
    meta_dict['wavelength'] = float(root.find('master/wavelength').text)
    meta_dict['polarization'] = root.find('master/frame/POLARIZATION').text.replace("'","").replace("b","")
    meta_dict['prf'] = float(root.find('master/prf').text)
    meta_dict['master_platform'] = root.find('master/platform/MISSION').text.replace("'","").replace("b","")
    meta_dict['master_absolute_orbit'] = int(root.find('master/frame/ORBIT_NUMBER').text)
    meta_dict['master_sensing_start'] = root.find('master/frame/SENSING_START').text 
    meta_dict['master_sensing_stop'] = root.find('master/frame/SENSING_STOP').text
    meta_dict['slave_platform'] = root.find('slave/platform/MISSION').text.replace("'","").replace("b","")
    meta_dict['slave_absolute_orbit'] = int(root.find('slave/frame/ORBIT_NUMBER').text)
    meta_dict['slave_sensing_start'] = root.find('slave/frame/SENSING_START').text
    meta_dict['slave_sensing_stop'] = root.find('slave/frame/SENSING_STOP').text

    meta_dict['width'] = width
    meta_dict['length'] = length
    meta_dict['xstep'] = xstep
    meta_dict['ystep'] = ystep
    meta_dict['north'] = north
    meta_dict['south'] = south
    meta_dict['west'] = west
    meta_dict['east'] = east
    meta_dict['ellipsoid'] = 'WGS84'


    meta_dict['incidence_angle'] = ''
    meta_dict['producer_names'] = ''

    meta_dict['processing_facility'] = clos.processing_facility
    meta_dict['processing_software'] = clos.processing_software
    meta_dict['processing_software_version'] = clos.processing_software_version 
    if clos.processing_atmos_correct_method:
        meta_dict['processing_atmos_correct_method'] = clos.processing_atmos_correct_method
    meta_dict['processing_dem'] = 'SRTM1'
    meta_dict['history'] = 'H5 file created: %s' % datetime.datetime.utcnow()

    meta_dict['average_coherence'] = np.mean(corp)
    meta_dict['max_coherence'] = np.nanmax(corp)
#    meta_dict['percent_unwrapped'] = ''
#    meta_dict['percent_atmos'] = ''
    meta_dict['baseline_perp'] = float(root.find('baseline/perp_baseline_top').text) 
    meta_dict['temporal_baseline'] = abs((first_date-last_date).days)
 
    ## add any other metadata to dictionary, for example if you parsed some other XML file into a dictionary
#    for key,value in wraprsc.iteritems():
#        meta_dict[key] = value

    print( 'Creating HDF5 file containing the geocoded *int, *unw, *cor, los, and dem ' )
    filename_root = '%s_%s_%03d_%04d_%s-%s_%04d_%05d' % (meta_dict['mission'],meta_dict['beam_swath'],meta_dict['relative_orbit'],meta_dict['frame'],first_date.strftime("%Y%m%d"),last_date.strftime("%Y%m%d"),meta_dict['temporal_baseline'],meta_dict['baseline_perp']) 
    h5file = os.getcwd() + '/'+filename_root+'.h5' 
    ## OPEN HDF5 FILE ##
    f = h5py.File(h5file)
    ## CREATE GEOCODE GROUP ##
    group = f.create_group('GEOCODE')
    ## CREATE GEOCODE DATASETS ##
    if not os.path.basename('unwrapped_interferogram') in group:
        dset = group.create_dataset('unwrapped_interferogram', data=unwp, compression='gzip')
    if not os.path.basename('wrapped_interferogram') in group:
        dset = group.create_dataset('wrapped_interferogram', data=intp, compression='gzip')
    if not os.path.basename('correlation') in group:
        dset = group.create_dataset('correlation', data=corp, compression='gzip')
    if not os.path.basename('incidence_angle') in group:
        dest = group.create_dataset('incidence_angle', data=rdrp, compression='gzip')
#    if not os.path.basename('digital_elevatino_model') in group:
#        dest = group.create_dataset('digital_elevation_model',data=dem,compression='gzip')

    ## WRITE ATTRIBUTES TO THE HDF ##
    for key,value in meta_dict.items():
        f.attrs[key] = value

    f.close()

if __name__ == '__main__':
    main(sys.argv[:])

