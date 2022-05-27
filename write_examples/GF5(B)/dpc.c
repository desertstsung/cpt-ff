/*
 *file: write_examples/GF5(B)/dpc.c
 *descreption:
 *  write cpt format file from DPC sensor
 *syntax:
 *  a.out hdf5 [ptxt, [cpt]]
 *init date: May/27/2022
 *last modify: May/27/2022
 *
 */

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <fcntl.h>
#include <math.h>

#include <curl/curl.h>

#include "hdf5.h"
#include "../../read/readcpt.h"


/*  Site info settings  */
#define WR_CPT_INVSITEBYTES 180
#define WR_CPT_DATELEN      11
#define WR_CPT_TIMELEN      9
#define WR_CPT_DTLEN        (WR_CPT_DATELEN+WR_CPT_TIMELEN)
#define WR_CPT_NPARAM       8


/*  Satellite info settings  */
#define WR_CPT_DPCNBANDS  8
#define WR_CPT_DPCDSNAME  "Data_Fields"
#define WR_CPT_DPCGLNAME  "Geolocation_Fields"
#define WR_CPT_DPCDTLEN   20
#define WR_CPT_DPCINAME   "I"
#define WR_CPT_DPCQNAME   "Q"
#define WR_CPT_DPCUNAME   "U"
#define WR_CPT_DPCLONNAME "Longitude"
#define WR_CPT_DPCLATNAME "Latitude"
#define WR_CPT_DPCSLNAME  "Sea_Land_Flags"
#define WR_CPT_DPCALTNAME "Surface_Altitude"
#define WR_CPT_DPCSZNAME  "Sol_Zen_Ang"
#define WR_CPT_DPCVZNAME  "View_Zen_Ang"
#define WR_CPT_DPCSANAME  "Sol_Azim_Ang"
#define WR_CPT_DPCVANAME  "View_Azim_Ang"
#define WR_CPT_DPCCNTRWV  (int16_t[WR_CPT_DPCNBANDS]) \
                           {443, -490, 565, -670, 763, 765, -865, 910}

#define WR_CPT_LONLIM_MIN  -180
#define WR_CPT_LONLIM_MAX  180
#define WR_CPT_LATLIM_MIN  -90
#define WR_CPT_LATLIM_MAX  90

#define WR_CPT_INVALIDLON(lon) ((lon < WR_CPT_LONLIM_MIN) || (lon > WR_CPT_LONLIM_MAX))
#define WR_CPT_INVALIDLAT(lat) ((lat < WR_CPT_LATLIM_MIN) || (lat > WR_CPT_LATLIM_MAX))
