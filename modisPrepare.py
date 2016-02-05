import sys
import os
import gdal
import glob
import re
import numpy as np
import osr
import subprocess
from numpy import inf
from numpy import nan

"""process MODIS surface reflectance images of California
"""

def checkDataLen(year, hdfpath):
	# get all dowladed fies, take the 006 version 
	alldownload = glob.glob(hdfpath + '\*006*.hdf')
	
	# check if file # reasonable 
	h08v04 = []
	h08v05 = []
	h09v04 = []
	for x in alldownload:
		if 'h08v04' in x:
			h08v04.append(x)
		if 'h08v05' in x:
			h08v05.append(x)
		if 'h09v04' in x:
			h09v04.append(x)

	h08v04 = sorted(h08v04, key=lambda x: re.search('A20\d+', x).group(0))
	h08v05 = sorted(h08v05, key=lambda x: re.search('A20\d+', x).group(0))
	h09v04 = sorted(h09v04, key=lambda x: re.search('A20\d+', x).group(0))

	print len(h08v04)  #h08v04
	print len(h08v05)  #h08v05
	print len(h09v04)  #h09v04

	# # further examine the file names 
	# checklist=[]
	# print 'data for h09v04'
	# for i in range(len(h09v04)):
	# 	if 'A'+str(year) in h09v04[i]:
	# 		# print h08v04[i]
	# 		checklist.append(re.search('M\w+\d+A1.A20\d+', h09v04[i]).group(0))

	# conlist = []
	# for i in range(1, 366, 8):
	# 	conlist.append("%03d"%i)
	# print year, '***************', sorted(checklist)
	# print year, '***************', conlist
	# print year, '***************', len(checklist)
	# print year, '***************', len(conlist)

def convert2tif(inhdf):
	# convert hdf4 to tif format
	try: 
		ds = gdal.Open(inhdf) 
		bands = ds.GetSubDatasets()
	except:
		print "cannot open hdf file"

	# read in teh subdatesets of hdf file--b1:red; b2:nir
	b1 = gdal.Open(ds.GetSubDatasets()[0][0])
	b2 = gdal.Open(ds.GetSubDatasets()[1][0])

	# get project, geotransform, rows, cols of the data
	projection = b1.GetProjection()
	geotransform = b1.GetGeoTransform()
	rows = b1.RasterYSize
	cols = b1.RasterXSize
	# driver = b1.GetDriver()

	# set output driver to be geotiff
	driver = gdal.GetDriverByName( 'GTiff' )  
	outDs = driver.Create(inhdf[:-4] + '.tif', cols, rows, 2, gdal.GDT_Int32)

	b1 = np.array(b1.ReadAsArray())
	b2 = np.array(b2.ReadAsArray())

	outBand1 = outDs.GetRasterBand(1)
	outBand2 = outDs.GetRasterBand(2)
		
	# write gdd data
	outBand1.WriteArray(b1)
	outBand2.WriteArray(b2)

	# flush data to disk, set the NoData value and calculate stats
	outBand1.FlushCache()
	outBand1.SetNoDataValue(-28672)
	outBand2.FlushCache()
	outBand2.SetNoDataValue(-28672)

	# georeference the image and set the projection
	outDs.SetGeoTransform(geotransform)
	outDs.SetProjection(projection)

	# reproejct 
	# transform = osr.CoordinateTransformation(projection, osr.SpatialReference().ImportFromEPSG(4269))

def mosaicTifs(alldownload):
	
	h08v04 = []
	h08v05 = []
	h09v04 = []
	for x in alldownload:
		print x
		if 'h08v04' in x:
			h08v04.append(x)
		if 'h08v05' in x:
			h08v05.append(x)
		if 'h09v04' in x:
			h09v04.append(x)

	h08v04 = sorted(h08v04, key=lambda x: re.search('A20\d+', x).group(0))
	h08v05 = sorted(h08v05, key=lambda x: re.search('A20\d+', x).group(0))
	h09v04 = sorted(h09v04, key=lambda x: re.search('A20\d+', x).group(0))

	# # https://jgomezdans.github.io/stitching-together-modis-data.html
	crsString = '-t_srs EPSG:4269'
	inshp =  r'G:\database\pur\PURGIS\shapefiles\cnty24k09_1_state_poly.shp'

	for i in range(len(h08v05)):
		# assert input data are from the same date
		if re.search('A2\d+', h08v05[i]).group(0) == re.search('A2\d+', h08v04[i]).group(0) == re.search('A2\d+', h09v04[i]).group(0):
			
			file1_b1 = 'HDF4_EOS:EOS_GRID:"'+ h08v05[i] +'":MOD_Grid_500m_Surface_Reflectance:sur_refl_b01'
			file2_b1 = 'HDF4_EOS:EOS_GRID:"'+ h08v04[i] +'":MOD_Grid_500m_Surface_Reflectance:sur_refl_b01'
			file3_b1 = 'HDF4_EOS:EOS_GRID:"'+ h09v04[i] +'":MOD_Grid_500m_Surface_Reflectance:sur_refl_b01'
			outVRT1 = os.path.dirname( h08v05[i] ) + '\mergedB1_' + re.search('A2\d+', h08v05[i]).group(0) +'.vrt'
			outTIF1 = outVRT1[:-4] + '.tif'
			print outVRT1
			
			mergeCom1 = ' '.join(['C:\\OSGeo4W\\bin\\gdalbuildvrt.exe', outVRT1, file1_b1, file2_b1, file3_b1])
			subprocess.Popen(mergeCom1, shell=True).wait()
			
			reprojCom1 = ' '.join(['C:\\OSGeo4W\\bin\\gdalwarp.exe',  crsString, outVRT1, outTIF1])
			subprocess.Popen(reprojCom1, shell=True).wait()

			clipCom1 = ' '.join(['C:\\OSGeo4W\\bin\\gdalwarp.exe','-cutline', inshp, '-crop_to_cutline', outTIF1, outTIF1[:-4]+'_ca.tif'])
			subprocess.Popen(clipCom1, shell=True).wait()

			# processing band2 

			file1_b2 = 'HDF4_EOS:EOS_GRID:"'+ h08v05[i] +'":MOD_Grid_500m_Surface_Reflectance:sur_refl_b02'
			file2_b2 = 'HDF4_EOS:EOS_GRID:"'+ h08v04[i] +'":MOD_Grid_500m_Surface_Reflectance:sur_refl_b02'
			file3_b2 = 'HDF4_EOS:EOS_GRID:"'+ h09v04[i] +'":MOD_Grid_500m_Surface_Reflectance:sur_refl_b02'
			outVRT2 = os.path.dirname( h08v05[i] ) + '\mergedB2_' + re.search('A2\d+', h08v05[i]).group(0) +'.vrt'
			outTIF2 = outVRT2[:-4] + '.tif'
			print outVRT2
			
			mergeCom2 = ' '.join(['C:\\OSGeo4W\\bin\\gdalbuildvrt.exe', outVRT2, file1_b2, file2_b2, file3_b2])
			subprocess.Popen(mergeCom2, shell=True).wait()

			reprojCom2 = ' '.join(['C:\\OSGeo4W\\bin\\gdalwarp.exe',  crsString, outVRT2, outTIF2])
			subprocess.Popen(reprojCom2, shell=True).wait()

			clipCom2 = ' '.join(['C:\\OSGeo4W\\bin\\gdalwarp.exe','-cutline', inshp, '-crop_to_cutline', outTIF2, outTIF2[:-4]+'_ca.tif' ]) #-dstnodata "-28672"
			subprocess.Popen(clipCom2, shell=True).wait()

			# processing band3: all in 500 meters

			file1_b3 = 'HDF4_EOS:EOS_GRID:"'+ h08v05[i] +'":MOD_Grid_500m_Surface_Reflectance:sur_refl_b03'
			file2_b3 = 'HDF4_EOS:EOS_GRID:"'+ h08v04[i] +'":MOD_Grid_500m_Surface_Reflectance:sur_refl_b03'
			file3_b3 = 'HDF4_EOS:EOS_GRID:"'+ h09v04[i] +'":MOD_Grid_500m_Surface_Reflectance:sur_refl_b03'
			outVRT3 = os.path.dirname( h08v05[i] ) + '\mergedB3_' + re.search('A2\d+', h08v05[i]).group(0) +'.vrt'
			outTIF3 = outVRT3[:-4] + '.tif'
			print outVRT3
			
			mergeCom3 = ' '.join(['C:\\OSGeo4W\\bin\\gdalbuildvrt.exe', outVRT3, file1_b3, file2_b3, file3_b3])
			subprocess.Popen(mergeCom3, shell=True).wait()

			reprojCom3 = ' '.join(['C:\\OSGeo4W\\bin\\gdalwarp.exe',  crsString, outVRT3, outTIF3])
			subprocess.Popen(reprojCom3, shell=True).wait()

			clipCom3 = ' '.join(['C:\\OSGeo4W\\bin\\gdalwarp.exe','-cutline', inshp, '-crop_to_cutline', outTIF3, outTIF3[:-4]+'_ca.tif' ]) #-dstnodata "-28672"
			subprocess.Popen(clipCom3, shell=True).wait()

			os.remove(outVRT1)
			os.remove(outVRT2)
			os.remove(outVRT3)
			os.remove(outTIF1)
			os.remove(outTIF2)
			os.remove(outTIF3)

def viCal(band1, band2, band3 = None):
	"""	calculate NDVI and EVI, when band3 == None, export NDVI only, 
		when band3 is given, export EVI as the 2nd layer
		MOD09Q1 products contains 250-m band1 and band2 (R and NIR)
		MOD09A1 products contains 500-m band1, band2, and band3(BLUE)
	"""
	# open tif files bands using gdal 
	b1ds = gdal.Open(band1)
	b2ds = gdal.Open(band2)

	# get info from one of the data source
	rows = b1ds.RasterYSize
	cols = b1ds.RasterXSize
	driver = b1ds.GetDriver()
	geotransform = b1ds.GetGeoTransform()
	resolutuion  = geotransform[1]
	projection = b1ds.GetProjection()
	nodata =  b1ds.GetRasterBand(1).GetNoDataValue()

	# get pixel value as numpy array 
	b1 = np.array(b1ds.GetRasterBand(1).ReadAsArray())
	b2 = np.array(b2ds.GetRasterBand(1).ReadAsArray())
	b1 = b1.astype(float)
	b2 = b2.astype(float)

	# ndvi calculation NDVI =  ( NIR - RED ) / ( NIR + RED )
	with np.errstate(divide='ignore', invalid='ignore'):
		ndvi = (b2*0.0001-b1*0.0001)/(b2*0.0001+b1*0.0001)
	print'sampel some nvdi value', ndvi[500][300:320]
	ndvi[ndvi== -inf] = nan 
	ndvi[ndvi== inf] = nan 

	# setup the output file info 
	outFileName = os.path.join(os.path.dirname(band1), 'vi_' + re.search('A\d+_ca.tif', band1).group(0))
	print 'output file name %s' %outFileName

	if band3 is None: 
		outDs = driver.Create(outFileName, cols, rows, 1, gdal.GDT_Float32)
		outBand1 = outDs.GetRasterBand(1)

		# write gdd data
		outBand1.WriteArray(ndvi)

		# flush data to disk, set the NoData value and calculate stats
		outBand1.FlushCache()
		outBand1.SetNoDataValue(nodata)

	else:
			
		b3ds = gdal.Open(band3)
		b3 = np.array(b3ds.GetRasterBand(1).ReadAsArray()) # resampel needed 
		b3 = b3.astype(float)

		# EVI = 2.5 * [( NIR - RED ) / ( NIR + 6.0 * RED - 7.5 * BLUE+ 1.0 )]
		with np.errstate(divide='ignore', invalid='ignore'):
			evi = 2.5 * ((b2*0.0001-b1*0.0001) / (b2*0.0001 + 6*b1*0.0001 -7.5*b3*0.0001 + 1))
			evi[evi== -inf] = nan 
			evi[evi== inf] = nan 
		print evi[500][300:320]
	
		# set output datasource with 2 bands 
		outDs = driver.Create(outFileName, cols, rows, 2, gdal.GDT_Float32)
		outBand1 = outDs.GetRasterBand(1)
		outBand2 = outDs.GetRasterBand(2)

		# write gdd data
		outBand1.WriteArray(ndvi)
		outBand2.WriteArray(evi)

		# flush data to disk, set the NoData value and calculate stats
		outBand1.FlushCache()
		outBand1.SetNoDataValue(nodata)
		outBand2.FlushCache()
		outBand2.SetNoDataValue(nodata)
	
	# georeference the image and set the projection
	outDs.SetGeoTransform(geotransform)
	outDs.SetProjection(projection)

	del ndvi, rows, cols, driver, projection, nodata


if __name__== "__mian__":
	print "tested on a windows desktop"


