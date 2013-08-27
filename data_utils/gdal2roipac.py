#! /usr/bin/env python
###############################################################################
#  gdal2roipac.py
#
#  Project:  Seamless SAR Archive
#  Purpose:  Convert any format supported by GDAL to ROI_PAC format
#  Author:   Scott Baker
#  Created:  August 2013
#
###############################################################################
#  Copyright (c) 2013, Scott Baker 
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

import sys
from osgeo import gdal

if len(sys.argv)<3:
    print "USAGE: %s IN_FILE OUT_FILE" % sys.argv[0]
    print "\nExample: %s dem.tif roipac.dem" % sys.argv[0]
    exit()
#######################################
##### CREATE A ROI_PAC FORMAT DEM #####
#######################################
indataset = gdal.Open(sys.argv[1])
data = indataset.ReadAsArray()
data.flatten()
data.tofile(sys.argv[2])
adfGeoTransform = indataset.GetGeoTransform(can_return_null = True)
x_first = '%.12f' % adfGeoTransform[0]
x_step = '%.12f' % adfGeoTransform[1]
y_first = '%.12f' % adfGeoTransform[3]
y_step = '%.12f' % adfGeoTransform[5]
with open('roipac.dem.rsc','w') as RSC:
    RSC.write('WIDTH          '+str(indataset.RasterXSize)+'\n')
    RSC.write('FILE_LENGTH    '+str(indataset.RasterYSize)+'\n')
    RSC.write('X_FIRST        '+x_first+'\n')
    RSC.write('Y_FIRST        '+y_first+'\n')
    RSC.write('X_STEP         '+x_step+'\n')
    RSC.write('Y_STEP         '+y_step+'\n')
    RSC.write('Z_SCALE        1\n')
    RSC.write('Z_OFFSET       0\n')
    RSC.write('X_UNIT         degrees\n')
    RSC.write('Y_UNIT         degrees\n')
    RSC.write('PROJECTION     LATLON')
