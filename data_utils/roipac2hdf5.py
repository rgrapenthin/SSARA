#! /usr/bin/env python
###############################################################################
# roipac2hdf5.py
#
#  Project:  Seamless SAR Archive
#  Purpose:  Create HDF5 interferogram product 
#  Author:   Scott Baker
#  Created:  July 2014
#
###############################################################################
#  Copyright (c) 2014, Scott Baker
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

import numpy as np
import h5py

def read_rsc_file(rscfile):
  '''Read the .rsc file into a python dictionary structure.

  '''
  rsc_dict = dict(np.loadtxt(rscfile,dtype=str))
  return rsc_dict

def read_float32(floatfile):
  '''Reads roi_pac unw, cor, or hgt data.

  Requires the file path and returns amplitude and phase
  Usage:
    amplitude, phase, rscDictionary = readUnw('geo_070603-070721_0048_00018.unw')
  '''
  rscContents = read_rsc_file(floatfile + '.rsc')
  width = int(rscContents['WIDTH'])
  length = int(rscContents['FILE_LENGTH'])
  oddindices = np.where(np.arange(length*2)&1)[0]
  data = np.fromfile(floatfile,np.float32,length*2*width).reshape(length*2,width)
  a = np.array([data.take(oddindices-1,axis=0)]).reshape(length,width)
  p = np.array([data.take(oddindices,axis=0)]).reshape(length,width)
  return a, p, rscContents

def read_complex64(complexfile):
  '''Reads roi_pac int or slc data.

  Requires the file path and returns amplitude and phase
  Usage:
    amp, phase, rscDictionary = readInt('geo_070603-070721_0048_00018.int')
  '''
  rscContents = read_rsc_file(complexfile + '.rsc')
  width = int(rscContents['WIDTH'])
  length = int(rscContents['FILE_LENGTH'])
  data = np.fromfile(complexfile,np.complex64,length*2*width).reshape(length,width)
  a = np.array([np.hypot(data.real,data.imag)]).reshape(length,width)
  p = np.array([np.arctan2(data.imag,data.real)]).reshape(length,width)
  return a, p, rscContents

def read_dem(demfile):
  '''Read a roipac dem file.

  Input:
    roi_pac format dem file
  '''
  rscContents = read_rsc_file(demfile + '.rsc')
  width = int(rscContents['WIDTH'])
  length = int(rscContents['FILE_LENGTH'])
  d=np.fromfile(demfile,dtype=np.int16).reshape(length,width)
  return d, rscContents

def parse():
    '''Command line parser.

    You can change/add defaults for any of these if you work a lot with the same mission for example.   You should
    also update the defaults for -software, -software_version, and -institution.
    '''
    parser = argparse.ArgumentParser(description='Create HDF5 interferogram product from GMTSAR output')
    parser.add_argument('-rsc1', dest='rsc1', action='store', help='SLC rsc file for the first date or master', type=str,required=True)
    parser.add_argument('-rsc2', dest='rsc2', action='store', help='SLC rsc file for the second date or slave', type=str,required=True)
    ## REQUIRED METADATA FOR ARCHIVING PRODUCTS ##
#    parser.add_argument('-mission', dest='mission', action='store', help='Name of the mission', type=str, required=True)
#    parser.add_argument('-relative_orbit', dest='relative_orbit', action='store', help='Relative orbit/Track/Path number', type=int, required=True)
    parser.add_argument('-processing_type', dest='processing_type', action='store', help='Type of processing: INTERFEROGRAM, LOS_VELOCITY,...', type=str, default='INTERFEROGRAM')
#    parser.add_argument('-footprint', dest='scene_footprint', action='store', help='WKT Polygon for the area covered by the swath', required=True, type=str)
    parser.add_argument('-swath', dest='beam_swath', action='store', help='Swath name without underscores', type=str, required=True)
    ## RECOMMENDED METADATA ##
    parser.add_argument('-beam_mode', dest='beam_mode', action='store', help='', type=str)
    parser.add_argument('-frame', dest='frame', action='store', help='', type=int)
#    parser.add_argument('-flight', dest='flight_direction', action='store', help='', type=str, default='A')
#    parser.add_argument('-look', dest='look_direction', action='store', help='', type=str, default='R')
#    parser.add_argument('-polarization', dest='polarization', action='store', help='', type=str)
    parser.add_argument('-software', dest='processing_software', action='store', help='', type=str, default='ROI_PAC')
    parser.add_argument('-software_version', dest='processing_software_version', action='store', help='', type=str,default='3.1')
    parser.add_argument('-atmos_correct_method', dest='processing_atmos_correct_method', action='store', help='', type=str)
    parser.add_argument('-institution', dest='processing_facility', action='store', help='', type=str, default='UNAVCO')
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

    rsc_master = read_rsc_file(clos.rsc1)
    rsc_slave = read_rsc_file(clos.rsc2)
    first_date = datetime.datetime.strptime(rsc_master['DATE'],'%y%m%d')
    last_date = datetime.datetime.strptime(rsc_slave['DATE'],'%y%m%d')
    rsc_baseline = read_rsc_file('%s_%s_baseline.rsc' % (first_date.strftime("%y%m%d"),last_date.strftime("%y%m%d")))

    ### READ GEOCODE DATASETS ###
    geo_root = 'geo_%s-%s' % (first_date.strftime("%y%m%d"),last_date.strftime("%y%m%d"))
    # the file root is hardcoded based on the standard processing produces same/standard filenames each time
    # if you have modified output names, you might need to update this to match your configuration
    ampunw,unw,unwrsc = read_float32(geo_root+'.unw')
    ampwrap,wrapped,wraprsc = read_complex64(geo_root+'.int')
    ampcor,cor,corrsc = read_float32(geo_root+'.cor')
    ampinc, phaseinc, incrsc = read_float32('geo_incidence.unw')
    dem,demrsc = read_dem('../DEM/roipac.dem')

    # these define the footprint of the scene and are used to create the WKT POLYGON below
    lats = [wraprsc['LAT_REF1'],wraprsc['LAT_REF3'],wraprsc['LAT_REF4'],wraprsc['LAT_REF2'],wraprsc['LAT_REF1']]
    lons = [wraprsc['LON_REF1'],wraprsc['LON_REF3'],wraprsc['LON_REF4'],wraprsc['LON_REF2'],wraprsc['LON_REF1']]

    #################################
    ###  METADATA
    #################################
    meta_dict = {}
    ## MANDATORY METADATA ##
    meta_dict['mission'] = rsc_master['PLATFORM']
    meta_dict['beam_swath'] = clos.beam_swath 
    meta_dict['relative_orbit'] = int(rsc_master['TRACK'])
    meta_dict['scene_footprint'] = "POLYGON((" + ",".join([lon+' '+lat for lat,lon in zip(lats,lons)]) + "))" 
    meta_dict['first_date'] = first_date.strftime("%Y%m%d")
    meta_dict['last_date'] = last_date.strftime("%Y%m%d")
    meta_dict['processing_type'] = clos.processing_type # SET AS A DEFAULT IN parse()

    ## RECOMMENDED METADATA ##
    if clos.beam_mode:
        meta_dict['beam_mode'] = clos.beam_mode
    meta_dict['frame'] = int(rsc_master['FIRST_FRAME'])
    meta_dict['flight_direction'] = rsc_master['ORBIT_DIRECTION'].upper()
    if int(rsc_master['ANTENNA_SIDE']) == -1:
        meta_dict['look_direction'] = 'R' 
    else:
        meta_dict['look_direction'] = 'L'
    meta_dict['polarization'] = rsc_master['POLARIZATION']
    meta_dict['prf'] = float(rsc_master['PRF'])
    meta_dict['master_platform'] = rsc_master['PLATFORM']
    meta_dict['master_absolute_orbit'] = int(rsc_master['ORBIT_NUMBER'])
    meta_dict['master_doppler'] = rsc_master['DOPPLER_RANGE0']+', '+rsc_master['DOPPLER_RANGE1']+', '+rsc_master['DOPPLER_RANGE2']+', '+rsc_master['DOPPLER_RANGE3']
    meta_dict['slave_platform'] = rsc_slave['PLATFORM']
    meta_dict['slave_absolute_orbit'] = int(rsc_slave['ORBIT_NUMBER'])
    meta_dict['slave_doppler'] = rsc_slave['DOPPLER_RANGE0']+', '+rsc_slave['DOPPLER_RANGE1']+', '+rsc_slave['DOPPLER_RANGE2']+', '+rsc_slave['DOPPLER_RANGE3']

    meta_dict['processing_facility'] = clos.processing_facility
    meta_dict['processing_software'] = clos.processing_software
    meta_dict['processing_software_version'] = clos.processing_software_version 
    if clos.processing_atmos_correct_method:
        meta_dict['processing_atmos_correct_method'] = clos.processing_atmos_correct_method
    meta_dict['processing_dem'] = 'SRTM'
    meta_dict['history'] = 'H5 file created: %s' % datetime.datetime.utcnow()

    meta_dict['average_coherence'] = np.mean(cor)
    meta_dict['max_coherence'] = np.nanmax(cor)
#    meta_dict['percent_unwrapped'] = ''
#    meta_dict['percent_atmos'] = ''
    meta_dict['baseline_perp'] = np.mean([float(rsc_baseline['P_BASELINE_TOP_HDR']),float(rsc_baseline['P_BASELINE_BOTTOM_HDR'])])
    meta_dict['temporal_baseline'] = abs((first_date-last_date).days)
 
    ## APPEND baseline.rsc contents to the metadata dictionary
    for key,value in wraprsc.iteritems():
        meta_dict[key] = value
    for key,value in rsc_baseline.iteritems():
        meta_dict[key] = value

    print 'Creating HDF5 file containing geo*int, geo*unw, geo*cor, and geo_incidence.unw ' 
    filename_root = '%s_%s_%03d_%04d_%s-%s_%04d_%05d' % (meta_dict['mission'],meta_dict['beam_swath'],meta_dict['relative_orbit'],meta_dict['frame'],meta_dict['first_date'],meta_dict['last_date'],meta_dict['temporal_baseline'],meta_dict['baseline_perp']) 
    h5file = os.getcwd() + '/'+filename_root+'.h5' 
    ## OPEN HDF5 FILE ##
    f = h5py.File(h5file)
    ## CREATE GEOCODE GROUP ##
    group = f.create_group('GEOCODE')
    ## CREATE GEOCODE DATASETS ##
    if not os.path.basename('unwrapped_interferogram') in group:
        dset = group.create_dataset('unwrapped_interferogram', data=unw, compression='gzip')
    if not os.path.basename('wrapped_interferogram') in group:
        dset = group.create_dataset('wrapped_interferogram', data=wrapped, compression='gzip')
    if not os.path.basename('correlation') in group:
        dset = group.create_dataset('correlation', data=cor, compression='gzip')
    if not os.path.basename('incidence_angle') in group:
        dest = group.create_dataset('incidence_angle', data=ampinc, compression='gzip')
    if not os.path.basename('digital_elevatino_model') in group:
        dest = group.create_dataset('digital_elevation_model',data=dem,compression='gzip')

    ## WRITE ATTRIBUTES TO THE HDF ##
    for key,value in meta_dict.iteritems():
        f.attrs[key] = value

    f.close()

if __name__ == '__main__':
    main(sys.argv[:])
