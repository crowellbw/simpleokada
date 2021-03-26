#!/usr/bin/python

import math
import coord_tools
import numpy

def fault_CMT(lat,lon,depth,strikeF,dipF,nstr,ndip,LEN,WID):
	[x0,y0]=coord_tools.ll2utm(lon, lat, lon, lat)
	x0 = x0/1000
	y0 = y0/1000
	
	fault_alt = numpy.zeros([nstr*ndip, 1])
	fault_lon = numpy.zeros([nstr*ndip, 1])
	fault_lat = numpy.zeros([nstr*ndip, 1])
	fault_lon1 = numpy.zeros([nstr*ndip, 1])
	fault_lat1 = numpy.zeros([nstr*ndip, 1])
	fault_lon2 = numpy.zeros([nstr*ndip, 1])
	fault_lat2 = numpy.zeros([nstr*ndip, 1])
	fault_lon3 = numpy.zeros([nstr*ndip, 1])
	fault_lat3 = numpy.zeros([nstr*ndip, 1])
	fault_lon4 = numpy.zeros([nstr*ndip, 1])
	fault_lat4 = numpy.zeros([nstr*ndip, 1])
	strike = strikeF*numpy.ones([nstr*ndip, 1]) 
	dip = dipF*numpy.ones([nstr*ndip, 1]) 


	if (WID/2*math.sin(dipF*math.pi/180) > depth):#sets the initial top depth  - either depth-width/2*sin(dip) or 0 depending on how wide and close to surface fault is
		z0 = 0
		x0 = x0 - LEN/2*math.sin(strikeF*math.pi/180) - depth*math.cos(dipF*math.pi/180)*math.sin((strikeF+90)*math.pi/180)
		y0 = y0 - LEN/2*math.cos(strikeF*math.pi/180) - depth*math.cos(dipF*math.pi/180)*math.cos((strikeF+90)*math.pi/180)
	else:
		z0 = depth - WID/2*math.sin(dipF*math.pi/180)
		x0 = x0 - LEN/2*math.sin(strikeF*math.pi/180) - WID/2*math.cos(dipF*math.pi/180)*math.sin((strikeF+90)*math.pi/180)
		y0 = y0 - LEN/2*math.cos(strikeF*math.pi/180) - WID/2*math.cos(dipF*math.pi/180)*math.cos((strikeF+90)*math.pi/180)

	DLEN = LEN/nstr
	DWID = WID/ndip

	xdoff = DWID*math.cos(dipF*math.pi/180)*math.sin((strikeF+90)*math.pi/180)
	ydoff = DWID*math.cos(dipF*math.pi/180)*math.cos((strikeF+90)*math.pi/180)

	xsoff = DLEN*math.sin(strikeF*math.pi/180)
	ysoff = DLEN*math.cos(strikeF*math.pi/180)


	k=0
	for j in range (0, ndip):
		for i in range (0, nstr):
			fault_X = x0 + (0.5+i)*xsoff + (0.5+j)*xdoff
			fault_Y = y0 + (0.5+i)*ysoff + (0.5+j)*ydoff

			fault_X1 = x0 + (i)*xsoff + (j)*xdoff
			fault_Y1 = y0 + (i)*ysoff + (j)*ydoff
	
			fault_X2 = x0 + (i+1)*xsoff + (j)*xdoff
			fault_Y2 = y0 + (i+1)*ysoff + (j)*ydoff

			fault_X3 = x0 + (i+1)*xsoff + (j+1)*xdoff
			fault_Y3 = y0 + (i+1)*ysoff + (j+1)*ydoff

			fault_X4 = x0 + (i)*xsoff + (j+1)*xdoff
			fault_Y4 = y0 + (i)*ysoff + (j+1)*ydoff

			fault_alt[k] = z0 + (0.5+j)*DWID*math.sin(dipF*math.pi/180)
			(fault_lon[k], fault_lat[k]) = coord_tools.utm2ll(fault_X*1000,fault_Y*1000,lon,lat)

			(fault_lon1[k], fault_lat1[k]) = coord_tools.utm2ll(fault_X1*1000,fault_Y1*1000,lon,lat)
			(fault_lon2[k], fault_lat2[k]) = coord_tools.utm2ll(fault_X2*1000,fault_Y2*1000,lon,lat)
			(fault_lon3[k], fault_lat3[k]) = coord_tools.utm2ll(fault_X3*1000,fault_Y3*1000,lon,lat)
			(fault_lon4[k], fault_lat4[k]) = coord_tools.utm2ll(fault_X4*1000,fault_Y4*1000,lon,lat)
			k = k+1

	return (fault_lon, fault_lat, fault_alt, strike, dip, DLEN*numpy.ones([nstr*ndip,1]), DWID*numpy.ones([nstr*ndip,1]))
			
