#! /usr/bin/env python
###############################################################################
# gmtsar2hdf5.py
#
#  Project:  Seamless SAR Archive (SSARA)
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
import datetime
import argparse

import h5py
from osgeo import gdal

def read_prm(prm_file):
    prm_dict = {}
    for line in open(prm_file):
        c = line.split("=")
        prm_dict[c[0].strip()] = str.replace(c[1], '\n', '').strip()
    return prm_dict

def parse():
    '''Command line parser.

    You can change/add defaults for any of these if you work a lot with the same mission for example.   You should 
    also update the defaults for -software, -software_version, and -institution.
    '''
    parser = argparse.ArgumentParser(description='Create HDF5 interferogram product from GMTSAR output')
    parser.add_argument('-prm1', dest='prm1', action='store', help='PRM file for the first date or master', type=str, default='master.PRM')
    parser.add_argument('-prm2', dest='prm2', action='store', help='PRM file for the second date or slave', type=str, default='slave.PRM')
    ## REQUIRED METADATA FOR ARCHIVING PRODUCTS ##
    parser.add_argument('-mission', dest='mission', action='store', help='Name of the mission. Will need to use this for ALOS2', type=str ) 
    parser.add_argument('-relative_orbit', dest='relative_orbit', action='store', help='Relative orbit/Track/Path number', type=int, required=True) 
    parser.add_argument('-processing_type', dest='processing_type', action='store', help='Type of processing: INTERFEROGRAM, LOS_VELOCITY,...', type=str, default='INTERFEROGRAM') 
    parser.add_argument('-footprint', dest='scene_footprint', action='store', help='WKT Polygon for the area covered by the swath', required=True, type=str) 
    parser.add_argument('-swath', dest='beam_swath', action='store', help='Swath name without underscores', type=str)
    ## RECOMMENDED METADATA ##
    parser.add_argument('-beam_mode', dest='beam_mode', action='store', help='', type=str)
    parser.add_argument('-frame', dest='frame', action='store', help='', type=int)
    parser.add_argument('-flight', dest='flight_direction', action='store', help='', type=str, default='A')
    parser.add_argument('-look', dest='look_direction', action='store', help='', type=str, default='R')
    parser.add_argument('-polarization', dest='polarization', action='store', help='', type=str)
    parser.add_argument('-software', dest='processing_software', action='store', help='', type=str, default='GMTSAR')
    parser.add_argument('-software_version', dest='processing_software_version', action='store', help='', type=str,default='9.4')
    parser.add_argument('-atmos_correct_method', dest='processing_atmos_correct_method', action='store', help='', type=str)
    parser.add_argument('-institution', dest='processing_facility', action='store', help='', type=str, default='UNAVCO')
    parser.add_argument('-master_platform', dest='master_platform', action='store', help='', type=str)
    parser.add_argument('-master_orbit', dest='master_absolute_orbit', action='store', help='', type=int)
    parser.add_argument('-slave_platform', dest='slave_platform', action='store', help='', type=str)
    parser.add_argument('-slave_orbit', dest='slave_absolute_orbit', action='store', help='', type=int)
#    parser.add_argument('-', dest='', action='store', help='', type=str)
    clos = parser.parse_args()
    return clos

def main(argv):
    # GET THE COMMAND LINE OPTIONS 
    clos = parse()

    print 'Creating HDF5 file containing correlation, wrapped, and unwrapped datasets'
    prm_master = read_prm(clos.prm1) # SET AS A DEFAULT IN parse() 
    prm_slave = read_prm(clos.prm2) # SET AS A DEFAULT IN parse()

    first_date = datetime.datetime.strptime(prm_master['SC_clock_start'].split(".")[0],'%Y%j') 
    last_date = datetime.datetime.strptime(prm_slave['SC_clock_start'].split(".")[0],'%Y%j')

    #################################
    ###  METADATA
    #################################
    meta_dict = {}
    ## MANDATORY METADATA ##
    if clos.mission:
        meta_dict['mission'] = clos.mission
    else:
        mission2sc_id = {'1': 'ERS','2':'ERS','3':'RS1','4':'ENV1','5':'ALOS','6':'','7':'TSX','8':'CSK','9':'RS2'}
        meta_dict['mission'] = mission2sc_id[prm_master['SC_identity']]
    meta_dict['beam_swath'] = clos.beam_swath
    meta_dict['relative_orbit'] = clos.relative_orbit
    meta_dict['first_date'] = datetime.datetime.strptime(prm_master['SC_clock_start'].split(".")[0],'%Y%j').strftime("%Y%m%d")
    meta_dict['last_date'] = datetime.datetime.strptime(prm_slave['SC_clock_start'].split(".")[0],'%Y%j').strftime("%Y%m%d") 
    meta_dict['scene_footprint'] = clos.scene_footprint 
    meta_dict['processing_type'] = clos.processing_type # SET AS A DEFAULT IN parse()

    ## RECOMMENDED METADATA ##
    if clos.beam_mode:
        meta_dict['beam_mode'] = clos.beam_mode
    if clos.frame:
        meta_dict['frame'] = clos.frame
    if clos.polarization:
        meta_dict['polarization'] = clos.polarization
    if 'orbdir' in prm_master:
        meta_dict['flight_direction'] = prm_master['orbdir'] 
    else:
        meta_dict['flight_direction'] = clos.flight_direction # SET AS A DEFAULT IN parse()
    meta_dict['look_direction'] = clos.look_direction # SET AS A DEFAULT IN parse()
    meta_dict['prf'] = prm_master['PRF']
    meta_dict['wavelength'] = prm_master['radar_wavelength']

    meta_dict['processing_facility'] = clos.processing_facility  # SET AS A DEFAULT IN parse()
    meta_dict['processing_software'] = clos.processing_software # SET AS A DEFAULT IN parse() 
    meta_dict['processing_software_version'] = clos.processing_software_version # SET AS A DEFAULT IN parse()
    if clos.processing_atmos_correct_method:
        meta_dict['processing_atmos_correct_method'] = clos.processing_atmos_correct_method
    meta_dict['processing_dem'] = 'SRTM1'
    meta_dict['history'] = 'H5 file created: %s' % datetime.datetime.utcnow()
    meta_dict['description'] = """Interferogram generated with %s verion %s by %s""" % (clos.processing_software,clos.processing_software_version,clos.processing_facility) 


    # CHECK IF MASTER AND SLAVE PLATFORMS ARE GIVEN, OTHERWISE USE MISSION
    if not clos.master_platform:
        clos.master_platform = meta_dict['mission']
    if not clos.slave_platform:
        clos.slave_platform = meta_dict['mission']

    meta_dict['master_platform'] = clos.master_platform
    if clos.master_absolute_orbit:
        meta_dict['master_absolute_orbit'] = clos.master_absolute_orbit
    meta_dict['master_doppler'] = prm_master['fd1']
#    meta_dict['master_scene'] = ''
    meta_dict['slave_platform'] = clos.slave_platform
    if clos.slave_absolute_orbit:
        meta_dict['slave_absolute_orbit'] = clos.slave_absolute_orbit
    meta_dict['slave_doppler'] = prm_slave['fd1']
#    meta_dict['slave_scene'] = ''

#    meta_dict['average_coherence'] = 0.0
#    meta_dict['max_coherence'] = 0.0
#    meta_dict['percent_unwrapped'] = 0.0
#    meta_dict['percent_atmos'] = 0.0
    meta_dict['baseline_perp'] = float(prm_slave['baseline_center'])
    meta_dict['temporal_baseline'] = abs((first_date-last_date).days)


    #################################
    ###  DATASETS: 
    ###  amplitude, correlation, wrapped phase, unwrapped phase, incidence angle, troposphere, dem, model
    #################################
    h5file = os.getcwd() + '/%s_%s_%03d_%04d_%s-%s_%04d_%05d.h5' % (meta_dict['mission'],meta_dict['beam_swath'],meta_dict['relative_orbit'],meta_dict['frame'],meta_dict['first_date'],meta_dict['last_date'],meta_dict['temporal_baseline'],meta_dict['baseline_perp'])
    f = h5py.File(h5file)
    group = f.create_group("GEOCODE")

    dset = gdal.Open('phase_ll.grd')
    meta_dict['X_FIRST'] = dset.GetGeoTransform()[0]
    meta_dict['X_STEP'] = dset.GetGeoTransform()[1]
    meta_dict['X_UNIT'] = 'degrees'
    meta_dict['Y_FIRST'] = dset.GetGeoTransform()[3]
    meta_dict['Y_STEP'] = dset.GetGeoTransform()[5]
    meta_dict['Y_UNIT'] = 'degrees'
    meta_dict['FILE_LENGTH'] = dset.RasterYSize
    meta_dict['WIDTH'] = dset.RasterXSize
    meta_dict['north'] = dset.GetGeoTransform()[3]
    meta_dict['west'] = dset.GetGeoTransform()[0]
    meta_dict['south'] = meta_dict['north'] + meta_dict['FILE_LENGTH']*meta_dict['Y_STEP']
    meta_dict['east'] = meta_dict['west'] + meta_dict['WIDTH']*meta_dict['X_STEP']
    if not os.path.basename('wrapped_interferogram') in group:
        group.create_dataset('wrapped_interferogram', data=dset.ReadAsArray(), compression='gzip')
    if not os.path.basename('unwrapped_interferogram') in group:
        group.create_dataset('unwrapped_interferogram', data=gdal.Open('unwrap_ll.grd').ReadAsArray(), compression='gzip')
    if not os.path.basename('wrapped_filtered_interferogram') in group:
        group.create_dataset('wrapped_filtered_interferogram', data= gdal.Open('phasefilt_ll.grd').ReadAsArray(), compression='gzip')
    if not os.path.basename('correlation') in group:
        group.create_dataset('correlation', data=gdal.Open('corr_ll.grd').ReadAsArray(), compression='gzip')  
#    if not os.path.basename('incidence') in group:
#        group.create_dataset('incidence', data=gdal.Open('look_ll.grd').ReadAsArray(), compression='gzip')
    for key,value in sorted(meta_dict.iteritems()):
        f.attrs[key] = value
    f.close()

if __name__ == '__main__':
    main(sys.argv[:])
    pass
