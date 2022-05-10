#!/usr/bin/env python

#############################################
#                   Usage                   #
#!!! argv[1]: xml filename of DPC or POSP   #
#!!! argv[2]: output aeronet text filename  #
#############################################

#############################################
#           Data Format Confugure           #
#!!! All points   : AVG=10                  #
#!!! Daily average: AVG=20                  #
#############################################
avg="AVG=10"

import requests as req
import datetime as dt
import numpy as np
import sys

hourup=dt.timedelta(hours=1)
urlpre="https://aeronet.gsfc.nasa.gov/cgi-bin/print_web_data_v3?"

def usage():
	print("Usage:")
	print("\t%s xml output"%(sys.argv[0]))
	sys.exit(5)

def xmldt(filename):
	r=""
	with open(filename) as fo:
		for line in fo:
			if -1!=line.find("StartTime"):
				idx=line.find(">")+1
				r+=dt.datetime.strptime(line[idx:idx+13],"%Y-%m-%d %H").strftime("year=%Y&month=%m&day=%d&hour=%H&")
			elif -1!=line.find("EndTime"):
				idx=line.find(">")+1
				r+=(dt.datetime.strptime(line[idx:idx+13],"%Y-%m-%d %H")+hourup).strftime("year2=%Y&month2=%m&day2=%d&hour2=%H&")
			elif -1!=line.find("NadirLong"):
				lonarr=np.array(
				line[line.find(">")+1:line.find("</")].split(","),dtype=float)
				r+="lon1=%07.2f&"%(lonarr.min())
				r+="lon2=%07.2f&"%(lonarr.max())
			elif -1!=line.find("NadirLat"):
				latarr=np.array(
				line[line.find(">")+1:line.find("</")].split(","),dtype=float)
				r+="lat1=%06.2f&"%(latarr.min())
				r+="lat2=%06.2f"%(latarr.max())
	return r

def savetolocal(url, local):
	with open(local, "wb") as fo:
		fo.write(req.get(url).content)

#main
if 3!=len(sys.argv):usage()
sat=xmldt(sys.argv[1])
if not sat:
	print("%s open failed!\nAborted.")
	sys.exit(1)
prod=""
with open("aeronet_datatype.conf") as fo:prod=fo.readline()[0:7]
if not prod:
	print("aeronet_datatype.conf NOT found. Using AOD15=1 instead")
	prod="AOD15=1"
savetolocal(urlpre+"&".join([avg,prod,sat]),sys.argv[2])

