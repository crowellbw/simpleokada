#!/usr/bin/python
import time
import math
import numpy
import csv

def location(stafile):
        with open(stafile) as f:
                for i, l in enumerate(f):
                        pass
                sta_length = i+1

        sta_lat = numpy.zeros([sta_length,1])
        sta_lon = numpy.zeros([sta_length,1])
        sta_alt = numpy.zeros([sta_length,1])

        k=0
        with open(stafile, 'rt') as f:
                reader = csv.reader(f, delimiter=',')
                for row in reader:
                        sta_lat[k,0] = row[2]
                        sta_lon[k,0] = row[3]
                        sta_alt[k,0] = row[4]
                        k=k+1

        return (sta_lat, sta_lon, sta_alt)







