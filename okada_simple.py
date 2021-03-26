#!/usr/bin/python
import math
import numpy
import okadagreen
from coord_tools import ll2utm
from location_init import location
from fault_plane import fault_CMT



stafile = 'sites.txt' #list of stations used 
[sta_lat,sta_lon,sta_alt] = location(stafile) #Location initialization


eqlat = 43.01 #earthquake latitude
eqlon = -125  #earthquake longitude
eqdep = 10 #depth in km, positive down
faultlen = 50 #fault length in km
faultwid = 20 #fault width in km
nlen = 20 #number of along length segments
nwid = 10 #number of along width segments
strikeang = 270 #strike angle, clockwise from north
dipang = 90 #dip angle, following Okada's definition (strike of zero has dip towards the east)
SS = 1 #strike slip, assigns this to all fault segments 
DS = 0 #dip slip 



[fault_lon,fault_lat,fault_alt,strike,dip,LEN,WID]=fault_CMT(eqlat,eqlon,eqdep,strikeang,dipang,nlen,nwid,faultlen,faultwid)

l1=len(sta_lat)
l2=len(fault_lat)

LEN = LEN*1000 
WID = WID*1000

xrs = numpy.zeros([l2,l1])
yrs = numpy.zeros([l2,l1])
zrs = numpy.zeros([l2,l1])

for i in range (0, l2):
        for j in range (0, l1):
                (x1,y1) = ll2utm(sta_lon[j],sta_lat[j], eqlon,eqlat)
                (x2,y2) = ll2utm(fault_lon[i],fault_lat[i], eqlon,eqlat)
                xrs[i,j] = (x1-x2)
                yrs[i,j] = (y1-y2)
                zrs[i,j] = sta_alt[j]+fault_alt[i]*1000



G = okadagreen.greenF(xrs, yrs, zrs, strike, dip, WID, LEN) #Compute Green's functions

S = numpy.zeros([l2*2,1])
for i in range (0,l2):
        S[2*i,0] = SS
        S[2*i+1,0] = DS

U = numpy.dot(G,S)

EN = numpy.zeros([l1,1])
NN = numpy.zeros([l1,1])
UN = numpy.zeros([l1,1])
for i in range (0, l1):
        EN[i,0]=U[3*i,0]
        NN[i,0]=U[3*i+1,0]
        UN[i,0]=U[3*i+2,0]



ffo = open('output.txt','w') #output file

for j in range(0, l1):
        vert = "{0:.4f}".format(float(UN[j]*100)) #I have it output cm since I assume slip is meters
        north = "{0:.4f}".format(float(NN[j]*100))
        east = "{0:.4f}".format(float(EN[j]*100))
        lonsta = "{0:.4f}".format(float(sta_lon[j]))
        latsta = "{0:.4f}".format(float(sta_lat[j]))
        ffo.write(lonsta+' '+latsta+' '+north+' '+east+' '+vert+'\n')
ffo.close()


