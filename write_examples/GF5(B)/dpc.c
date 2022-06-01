/*
 *file: write_examples/GF5(B)/dpc.c
 *descreption:
 *  write cpt format file from DPC sensor
 *syntax:
 *  a.out DPC_prefix [ptxt, [cpt]]
 *init date: May/27/2022
 *last modify: May/30/2022
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
#define WR_CPT_NPARAM       4


/*  Satellite info settings  */
#define WR_CPT_DPCNBANDS  8
#define WR_CPT_DPCDSNAME  "Data_Fields"
#define WR_CPT_DPCGLNAME  "Geolocation_Fields"
#define WR_CPT_DPCLONNAME "Longitude"
#define WR_CPT_DPCLATNAME "Latitude"
#define WR_CPT_DPCSLNAME  "Sea_Land_Flags"
#define WR_CPT_DPCALTNAME "Surface_Altitude"
#define WR_CPT_DPCSZNAME  "Sol_Zen_Ang"
#define WR_CPT_DPCVZNAME  "View_Zen_Ang"
#define WR_CPT_DPCSANAME  "Sol_Azim_Ang"
#define WR_CPT_DPCVANAME  "View_Azim_Ang"
#define WR_CPT_DPCSCLNAME "Scale_Factor"
#define WR_CPT_DPC443SUF  "_B443.h5"
#define WR_CPT_DPC490SUF  "_B490.h5"
#define WR_CPT_DPC565SUF  "_B565.h5"
#define WR_CPT_DPC670SUF  "_B670.h5"
#define WR_CPT_DPC763SUF  "_B763.h5"
#define WR_CPT_DPC765SUF  "_B765.h5"
#define WR_CPT_DPC865SUF  "_B865.h5"
#define WR_CPT_DPC910SUF  "_B910.h5"
#define WR_CPT_DPCCNTRWV  (int16_t[WR_CPT_DPCNBANDS]) \
                           {443, -490, 565, -670, 763, 765, -865, 910}

#define WR_CPT_XMLSUFFIX ".xml"
#define WR_CPT_XMLENDTAG "</ProductMetaData>"
#define WR_CPT_XMLSTTAG  "StartTime"
#define WR_CPT_XMLETTAG  "EndTime"
#define WR_CPT_XMLNCTAG  "SamplesCount"
#define WR_CPT_XMLNRTAG  "LineCount"
#define WR_CPT_XMLNLTAG  "MaxObservingAngles"

#define WR_CPT_LONLIM_MIN  -180
#define WR_CPT_LONLIM_MAX  180
#define WR_CPT_LATLIM_MIN  -90
#define WR_CPT_LATLIM_MAX  90

#define WR_CPT_INVALIDLON(lon)      ((lon < WR_CPT_LONLIM_MIN) || (lon > WR_CPT_LONLIM_MAX))
#define WR_CPT_INVALIDLAT(lat)      ((lat < WR_CPT_LATLIM_MIN) || (lat > WR_CPT_LATLIM_MAX))
#define WR_CPT_INVALIDGEO(lon, lat) (WR_CPT_INVALIDLON(lon) || WR_CPT_INVALIDLAT(lat))

#define WR_CPT_VALIDLON(lon)      ((lon > WR_CPT_LONLIM_MIN) && (lon < WR_CPT_LONLIM_MAX))
#define WR_CPT_VALIDLAT(lat)      ((lat > WR_CPT_LATLIM_MIN) && (lat < WR_CPT_LATLIM_MAX))
#define WR_CPT_VALIDGEO(lon, lat) (WR_CPT_VALIDLON(lon) && WR_CPT_VALIDLAT(lat))

#define WR_CPT_VALIDMASK(mask) (mask <= 100)
#define WR_CPT_VALIDALT(alt) ((alt < 10000) && (alt > -15000))

#define WR_CPT_SECDIFFMAX ((uint64_t) 60*30)  /*  30 min  */

#define WR_CPT_URLPREFIX "https://aeronet.gsfc.nasa.gov/cgi-bin/print_web_data_v3?"
#define WR_CPT_DTFNAME   "dpc_datatype.conf"
#define WR_CPT_URLDTLEN  7
#define WR_CPT_URLAVG    "AVG=10"
#define WR_CPT_PTSUFFIX  "ptxt"
#define WR_CPT_PTSUFLEN  strlen(WR_CPT_PTSUFFIX)
#define WR_CPT_URLDT1LEN 33
#define WR_CPT_URLDT2LEN 37
#define WR_CPT_URLGEOLEN 49

#define WR_CPT_CPTSUFFIX  "cpt"

struct wr_cpt_dpcband {
	hid_t iid;   /*  entrance of intensity    */
	hid_t szid;  /*  entrance of sol zen ang  */
	hid_t vzid;  /*  entrance of sat zen ang  */
	hid_t said;  /*  entrance of sol azi ang  */
	hid_t vaid;  /*  entrance of sat azi ang  */
	hid_t aid;   /*  entrance of altitude     */
	hid_t sid;   /*  entrance of sea-land     */
	
	hid_t gdid;  /*  dateset group id  */
	hid_t ggid;  /*  geo-loc group id  */
	hid_t fid;   /*  file id           */
	float *lon, *lat;
};

struct wr_cpt_dpcbandp {
	hid_t iid;   /*  entrance of intensity    */
	hid_t qid;   /*  entrance of polar q      */
	hid_t uid;   /*  entrance of polar u      */
	hid_t szid;  /*  entrance of sol zen ang  */
	hid_t vzid;  /*  entrance of sat zen ang  */
	hid_t said;  /*  entrance of sol azi ang  */
	hid_t vaid;  /*  entrance of sat azi ang  */
	hid_t aid;   /*  entrance of altitude     */
	hid_t sid;   /*  entrance of sea-land     */
	
	hid_t gdid;  /*  dateset group id  */
	hid_t ggid;  /*  geo-loc group id  */
	hid_t fid;   /*  file id           */
	float *lon, *lat;
};

struct wr_cpt_dpc {
	struct wr_cpt_dpcband  b443,
	                       b565,
	                       b763,
	                       b765,
	                       b910;
	struct wr_cpt_dpcbandp b490,
	                       b670,
	                       b865;
	
	hid_t h2id;  /*  2-d hyperslab  */
	hid_t h3id;  /*  3-d hyperslab  */
	hid_t m2id;  /*  2-d mem id     */
	hid_t m3id;  /*  3-d mem id     */
	
	hsize_t l2id[2];  /*  2-d hyper len  */
	hsize_t l3id[3];  /*  3-d hyper len  */
	
	uint8_t  nlayer;        /*  count of layer      */
	uint16_t ncol, nrow;    /*  count of row/col    */
	uint32_t nelements;     /*  full elements       */
	float    secsperline;   /*  timelapse per line  */
	double   scaleobs;      /*  scale for obs       */
	double   scaleang;      /*  scale for ang       */
	uint64_t secswhenscan;  /*  timestamp of begin  */
	uint64_t secswhenend;   /*  timestamp of end    */
	float   *lon, *lat;     /*  full lon/lat        */
	
	float lonmin, latmin;   /*  left bottom pixel   */
	float lonmax, latmax;   /*  right top pixel     */
};


/*  Error numbers  */
enum WR_CPT_ERR {
WR_CPT_EINVARG = 1,
WR_CPT_EOPEN,
WR_CPT_EINVSITE,
WR_CPT_EMEM,
WR_CPT_ENORCD,
WR_CPT_E
};


/*  XXX  */
const static size_t _cpt_1byte = sizeof(int8_t );
const static size_t _cpt_2byte = sizeof(int16_t);
const static size_t _cpt_4byte = sizeof(int32_t);
const static size_t _cpt_8byte = sizeof(int64_t);
const static size_t _cpt_parsz = sizeof(double[WR_CPT_NPARAM]);


/*  Declaration of main function  */
int wrcpt(const char *prefix, const char *ptxtfname, const char *cptfname);

/*  Program entrance  */
int main(int argc, char *argv[]) {
	if (3 == argc) {
		CPT_ECHOWITHTIME("Mode: download remote ptxt without creating cpt");
		return wrcpt(argv[1], argv[2], NULL);
	} else if (4 == argc) {
		return wrcpt(argv[1], argv[2], argv[3]);
	} else {
		CPT_ERRECHOWITHTIME("Usage: %s DPC_prefix [ptxt, [cpt]]", argv[0]);
		return WR_CPT_EINVARG;
	}
}

/*TODO move this fn to read/ dir
 *  Free several allocated space
 *  e.g.
 *      void *ptr1 = malloc(size1),
             *ptr2 = malloc(size2);
 *      cpt_freethemall(2, &ptr1, &ptr2);
 */
int cpt_freethemall(uint8_t n, ...)
{
	va_list ap;
	void **p;
	
	va_start(ap, n);
	while (n-- > 0) {
		p = va_arg(ap, void **);
		CPT_FREE(*p);
	}
	va_end(ap);
	
	return 0;
}

/*TODO move this fn to read/ dir
 *  Free one or more cpt_point st pointer
 */
int cpt_freepointall(struct cpt_point **p, uint16_t n)
{
	if (*p) {
		while (n-- > 0)
			CPT_FREE((*p+n)->params);
		CPT_FREE(*p);
	}
	
	return 0;
}

/*TODO move this fn to read/ dir
 *  Free one or more cpt_pt st pointer
 */
int cpt_freeptall(struct cpt_pt **p, uint32_t n)
{
	if (*p) {
		while (n-- > 0)
			cpt_freepointall(&(*p+n)->points, (*p+n)->nt);
		CPT_FREE(*p);
	}
	
	return 0;
}

/*TODO move this fn to read/ dir
 *  Free one or more cpt_channel st pointer
 */
int cpt_freechannelall(struct cpt_channel **p, uint16_t n)
{
	if (*p) {
		while (n-- > 0)
			cpt_freethemall(2, &(*p+n)->obs, &(*p+n)->ang);
		CPT_FREE(*p);
	}
	
	return 0;
}

/*TODO move this fn to read/ dir
 *  Free one or more cpt_pixel st pointer
 */
int cpt_freepixelall(struct cpt_pixel **p, uint16_t n)
{
	if (*p) {
		while (n-- > 0)
			cpt_freechannelall(&(*p+n)->channels, (*p+n)->nchannel);
		CPT_FREE(*p);
	}
	
	return 0;
}

/*TODO move this fn to read/ dir
 *  Free one or more cpt_px st pointer
 */
int cpt_freepxall(struct cpt_px **p, uint32_t n)
{
	if (*p) {
		while (n-- > 0) {
			cpt_freepixelall(&(*p+n)->centrepixel, 1);
			cpt_freepixelall(&(*p+n)->vicinity, (*p+n)->nvicinity);
		}
		CPT_FREE(*p);
	}
	
	return 0;
}

/*  Init essential info from xml  */
static int initfromxml(const char *fname, struct wr_cpt_dpc *st)
{
	int      fd;
	char    *buffer, *pbuf, *pprev, *line, *pline;
	char     dt[WR_CPT_DTLEN];
	uint32_t flen, llen, maxlen;
	struct tm stm;
	
	/*  Try openning text file  */
	if ((fd = open(fname, O_RDONLY)) < 0) {
		CPT_ERROPEN(fname);
		return WR_CPT_EOPEN;
	}
	
	/*  Buffer all content  */
	buffer = malloc(sizeof(char[flen = lseek(fd, 0, SEEK_END)]));
	if (!buffer) {
		CPT_ERRMEM(buffer);
		return WR_CPT_EMEM;
	}
	pread(fd, buffer, flen, 0);
	pprev = pbuf = buffer;
	close(fd);
	
	line = NULL;
	maxlen = 1;
	
	/*  Handle each line  */
	do {
		while (*++pbuf != '\n') ;
		
		llen = pbuf-pprev;
		if (llen > maxlen)
			line = realloc(line, sizeof(char[maxlen = llen]));
		/*
		 *  Set the rest of memory to 0 to avoid `strstr` return wrong value
		 *  Wrong case without this fix:
		 *    1. current line len(llen) is not greater than max line line(maxlen)
		 *    2. buffer(line) remains the same as previous
		 *    3. the front part of buffer(line) is replaced by newer line contents
		 *    4. the end part of buffer(line) remains as previous
		 *    5. needle substring happens to appear in the end part of previous one
		 *    6. strstr return true, but needle don't really appear in this line of text
		 */
		else
			memset(line+llen, 0, maxlen-llen);
		memcpy(line, pprev, llen);
		
		/*  End of XML  */
		if ((pline = strstr(line, WR_CPT_XMLENDTAG))) break;
		
		/*  Start time  */
		if ((pline = strstr(line, WR_CPT_XMLSTTAG))) {
			sscanf(pline, "%*[^'>']>%[^'<']", dt);
			strptime(dt, "%Y-%m-%d %H:%M:%S", &stm);
			st->secswhenscan = timegm(&stm);
			goto next_line;
		}
		
		/*  Ending time  */
		if ((pline = strstr(line, WR_CPT_XMLETTAG))) {
			sscanf(pline, "%*[^'>']>%[^'<']", dt);
			strptime(dt, "%Y-%m-%d %H:%M:%S", &stm);
			st->secswhenend = timegm(&stm);
			goto next_line;
		}
		
		/*  N of row  */
		if ((pline = strstr(line, WR_CPT_XMLNRTAG))) {
			sscanf(pline, "%*[^'>']>%hu", &st->nrow);
			goto next_line;
		}
		
		/*  N of col  */
		if ((pline = strstr(line, WR_CPT_XMLNCTAG))) {
			sscanf(pline, "%*[^'>']>%hu", &st->ncol);
			goto next_line;
		}
		
		/*  N of layer  */
		if ((pline = strstr(line, WR_CPT_XMLNLTAG))) {
			sscanf(pline, "%*[^'>']>%hhu", &st->nlayer);
			goto next_line;
		}
		
		next_line:
		pprev = ++pbuf;
	} while (*pbuf != '\0');
	
	/*  N of elements per layer  */
	st->nelements = (uint32_t) st->ncol * (uint32_t) st->nrow;
	
	/*  Timelapse  */
	st->secsperline = st->secswhenend - st->secswhenscan;
	if (0 == st->secsperline) {
		CPT_ERRECHOWITHTIME("BAD time formatting in XML!");
		return WR_CPT_E;
	}
	st->secsperline /= st->nrow;
	
	/*  Set pointers to NULL  */
	pbuf = pprev = pline = NULL;
	
	/*  Free allocated buffer  */
	line = realloc(line, sizeof(char[maxlen]));
	cpt_freethemall(2, &buffer, &line);
	
	return 0;
}

/*  Init band st from individual DPC h5  */
static int _initfromh5(const char *fname, const char *name, struct wr_cpt_dpcband *st, const size_t sbuf)
{
	hid_t id;
	
	/*  Open h5 entrance  */
	st->fid  = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
	st->gdid = H5Gopen(st->fid, WR_CPT_DPCDSNAME, H5P_DEFAULT);
	st->ggid = H5Gopen(st->fid, WR_CPT_DPCGLNAME, H5P_DEFAULT);
	
	st->iid  = H5Dopen(st->gdid, name, H5P_DEFAULT);
	st->szid = H5Dopen(st->gdid, WR_CPT_DPCSZNAME , H5P_DEFAULT);
	st->vzid = H5Dopen(st->gdid, WR_CPT_DPCVZNAME , H5P_DEFAULT);
	st->said = H5Dopen(st->gdid, WR_CPT_DPCSANAME , H5P_DEFAULT);
	st->vaid = H5Dopen(st->gdid, WR_CPT_DPCVANAME , H5P_DEFAULT);
	
	st->aid  = H5Dopen(st->ggid, WR_CPT_DPCALTNAME, H5P_DEFAULT);
	st->sid  = H5Dopen(st->ggid, WR_CPT_DPCSLNAME , H5P_DEFAULT);
	
	/*  Buffer all lat/lon  */
	id = H5Dopen(st->ggid, WR_CPT_DPCLATNAME, H5P_DEFAULT);
	st->lat = malloc(sbuf);
	H5Dread(id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, st->lat);
	H5Dclose(id);
	
	id = H5Dopen(st->ggid, WR_CPT_DPCLONNAME, H5P_DEFAULT);
	st->lon = malloc(sbuf);
	H5Dread(id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, st->lon);
	H5Dclose(id);
	
	return 0;
}

/*  Init polar band st from individual DPC h5  */
static int _initfromh5p(const char *fname, const char *namei, const char *nameq, const char *nameu,
                        struct wr_cpt_dpcbandp *st, const size_t sbuf)
{
	hid_t id;
	
	/*  Open h5 entrance  */
	st->fid  = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
	st->gdid = H5Gopen(st->fid, WR_CPT_DPCDSNAME, H5P_DEFAULT);
	st->ggid = H5Gopen(st->fid, WR_CPT_DPCGLNAME, H5P_DEFAULT);
	
	st->iid  = H5Dopen(st->gdid, namei, H5P_DEFAULT);
	st->qid  = H5Dopen(st->gdid, nameq, H5P_DEFAULT);
	st->uid  = H5Dopen(st->gdid, nameu, H5P_DEFAULT);
	st->szid = H5Dopen(st->gdid, WR_CPT_DPCSZNAME , H5P_DEFAULT);
	st->vzid = H5Dopen(st->gdid, WR_CPT_DPCVZNAME , H5P_DEFAULT);
	st->said = H5Dopen(st->gdid, WR_CPT_DPCSANAME , H5P_DEFAULT);
	st->vaid = H5Dopen(st->gdid, WR_CPT_DPCVANAME , H5P_DEFAULT);
	
	st->aid  = H5Dopen(st->ggid, WR_CPT_DPCALTNAME, H5P_DEFAULT);
	st->sid  = H5Dopen(st->ggid, WR_CPT_DPCSLNAME , H5P_DEFAULT);
	
	/*  Buffer all lat/lon  */
	id = H5Dopen(st->ggid, WR_CPT_DPCLATNAME, H5P_DEFAULT);
	st->lat = malloc(sbuf);
	H5Dread(id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, st->lat);
	H5Dclose(id);
	
	id = H5Dopen(st->ggid, WR_CPT_DPCLONNAME, H5P_DEFAULT);
	st->lon = malloc(sbuf);
	H5Dread(id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, st->lon);
	H5Dclose(id);
	
	return 0;
}

/*  Init st from DPC h5  */
static int initfromh5(const char *prefix, struct wr_cpt_dpc *st)
{
	hid_t attr;
	char *fname;
	float scale;
	const size_t bufsize = sizeof(float[st->nelements]);
	
	/*  443  */
	asprintf(&fname, "%s%s", prefix, WR_CPT_DPC443SUF);
	_initfromh5(fname, "I443", &st->b443, bufsize);
	
	/*  565  */
	asprintf(&fname, "%s%s", prefix, WR_CPT_DPC565SUF);
	_initfromh5(fname, "I565", &st->b565, bufsize);
	
	/*  763  */
	asprintf(&fname, "%s%s", prefix, WR_CPT_DPC763SUF);
	_initfromh5(fname, "I763", &st->b763, bufsize);
	
	/*  765  */
	asprintf(&fname, "%s%s", prefix, WR_CPT_DPC765SUF);
	_initfromh5(fname, "I765", &st->b765, bufsize);
	
	/*  910  */
	asprintf(&fname, "%s%s", prefix, WR_CPT_DPC910SUF);
	_initfromh5(fname, "I910", &st->b910, bufsize);
	
	/*  490  */
	asprintf(&fname, "%s%s", prefix, WR_CPT_DPC490SUF);
	_initfromh5p(fname, "I490P", "Q490P", "U490P", &st->b490, bufsize);
	
	/*  670  */
	asprintf(&fname, "%s%s", prefix, WR_CPT_DPC670SUF);
	_initfromh5p(fname, "I670P", "Q670P", "U670P", &st->b670, bufsize);
	
	/*  865  */
	asprintf(&fname, "%s%s", prefix, WR_CPT_DPC865SUF);
	_initfromh5p(fname, "I865P", "Q865P", "U865P", &st->b865, bufsize);
	
	/*  Scale Factor  */
	attr = H5Aopen(st->b443.iid, WR_CPT_DPCSCLNAME, H5P_DEFAULT);
	H5Aread(attr, H5T_NATIVE_FLOAT, &scale);
	st->scaleobs = scale;
	H5Aclose(attr);
	
	attr = H5Aopen(st->b443.szid, WR_CPT_DPCSCLNAME, H5P_DEFAULT);
	H5Aread(attr, H5T_NATIVE_FLOAT, &scale);
	st->scaleang = scale;
	H5Aclose(attr);
	
	/*  Space for hyper-reading  */
	st->h2id = H5Screate_simple(2, (hsize_t[2]) {st->nrow, st->ncol}, NULL);
	st->h3id = H5Screate_simple(3, (hsize_t[3]) {st->nlayer, st->nrow, st->ncol}, NULL);
	st->m2id = H5Screate_simple(1, (hsize_t[1]) {1}, NULL);
	st->m3id = H5Screate_simple(1, (hsize_t[1]) {st->nlayer}, NULL);
	
	st->l3id[0] = st->nlayer;
	st->l3id[1] = st->l3id[2] = 1;
	st->l2id[0] = st->l2id[1] = 1;
	
	/*  Find boundery corner  */
	st->lon = malloc(bufsize);
	st->lat = malloc(bufsize);
	st->lonmin = WR_CPT_LONLIM_MAX;
	st->lonmax = WR_CPT_LONLIM_MIN;
	st->latmin = WR_CPT_LATLIM_MAX;
	st->latmax = WR_CPT_LATLIM_MIN;
	for (uint32_t index = 0; index < st->nelements; ++index) {
		if (WR_CPT_VALIDGEO(st->b443.lon[index], st->b443.lat[index])) {
			st->lon[index] = st->b443.lon[index];
			st->lat[index] = st->b443.lat[index];
		} else if (WR_CPT_VALIDGEO(st->b490.lon[index], st->b490.lat[index])) {
			st->lon[index] = st->b490.lon[index];
			st->lat[index] = st->b490.lat[index];
		} else if (WR_CPT_VALIDGEO(st->b565.lon[index], st->b565.lat[index])) {
			st->lon[index] = st->b565.lon[index];
			st->lat[index] = st->b565.lat[index];
		} else if (WR_CPT_VALIDGEO(st->b670.lon[index], st->b670.lat[index])) {
			st->lon[index] = st->b670.lon[index];
			st->lat[index] = st->b670.lat[index];
		} else if (WR_CPT_VALIDGEO(st->b763.lon[index], st->b763.lat[index])) {
			st->lon[index] = st->b763.lon[index];
			st->lat[index] = st->b763.lat[index];
		} else if (WR_CPT_VALIDGEO(st->b765.lon[index], st->b765.lat[index])) {
			st->lon[index] = st->b765.lon[index];
			st->lat[index] = st->b765.lat[index];
		} else if (WR_CPT_VALIDGEO(st->b865.lon[index], st->b865.lat[index])) {
			st->lon[index] = st->b865.lon[index];
			st->lat[index] = st->b865.lat[index];
		} else if (WR_CPT_VALIDGEO(st->b910.lon[index], st->b910.lat[index])) {
			st->lon[index] = st->b910.lon[index];
			st->lat[index] = st->b910.lat[index];
		} else {
			st->lon[index] = st->b910.lon[index];
			st->lat[index] = st->b910.lat[index];
			continue;
		}
		
		if (st->lon[index] < st->lonmin)
			st->lonmin = st->lon[index];
		if (st->lat[index] < st->latmin)
			st->latmin = st->lat[index];
		if (st->lon[index] > st->lonmax)
			st->lonmax = st->lon[index];
		if (st->lat[index] > st->latmax)
			st->latmax = st->lat[index];
	}
	
	CPT_FREE(fname);
	cpt_freethemall(WR_CPT_DPCNBANDS, &st->b443.lon, &st->b490.lon, &st->b565.lon, &st->b670.lon,
	                                  &st->b763.lon, &st->b765.lon, &st->b865.lon, &st->b910.lon);
	cpt_freethemall(WR_CPT_DPCNBANDS, &st->b443.lat, &st->b490.lat, &st->b565.lat, &st->b670.lat,
	                                  &st->b763.lat, &st->b765.lat, &st->b865.lat, &st->b910.lat);
	
	return 0;
}

/*  Init DPC st  */
static int initdpcst(const char *prefix, struct wr_cpt_dpc **st)
{
	char *xmlfname;
	
	*st = malloc(sizeof(struct wr_cpt_dpc));
	asprintf(&xmlfname, "%s%s", prefix, WR_CPT_XMLSUFFIX);
	
	initfromxml(xmlfname, *st);
	initfromh5(prefix, *st);
	
	CPT_FREE(xmlfname);

	return 0;
}

/*  Destroy DPC st  */
static int cleandpcst(struct wr_cpt_dpc **st)
{
	/*  443  */
	H5Dclose((*st)->b443.iid);
	H5Dclose((*st)->b443.aid);
	H5Dclose((*st)->b443.sid);
	H5Dclose((*st)->b443.szid);
	H5Dclose((*st)->b443.vzid);
	H5Dclose((*st)->b443.said);
	H5Dclose((*st)->b443.vaid);
	H5Gclose((*st)->b443.gdid);
	H5Gclose((*st)->b443.ggid);
	H5Fclose((*st)->b443.fid);
	
	/*  565  */
	H5Dclose((*st)->b565.iid);
	H5Dclose((*st)->b565.aid);
	H5Dclose((*st)->b565.sid);
	H5Dclose((*st)->b565.szid);
	H5Dclose((*st)->b565.vzid);
	H5Dclose((*st)->b565.said);
	H5Dclose((*st)->b565.vaid);
	H5Gclose((*st)->b565.gdid);
	H5Gclose((*st)->b565.ggid);
	H5Fclose((*st)->b565.fid);
	
	/*  763  */
	H5Dclose((*st)->b763.iid);
	H5Dclose((*st)->b763.aid);
	H5Dclose((*st)->b763.sid);
	H5Dclose((*st)->b763.szid);
	H5Dclose((*st)->b763.vzid);
	H5Dclose((*st)->b763.said);
	H5Dclose((*st)->b763.vaid);
	H5Gclose((*st)->b763.gdid);
	H5Gclose((*st)->b763.ggid);
	H5Fclose((*st)->b763.fid);
	
	/*  765  */
	H5Dclose((*st)->b765.iid);
	H5Dclose((*st)->b765.aid);
	H5Dclose((*st)->b765.sid);
	H5Dclose((*st)->b765.szid);
	H5Dclose((*st)->b765.vzid);
	H5Dclose((*st)->b765.said);
	H5Dclose((*st)->b765.vaid);
	H5Gclose((*st)->b765.gdid);
	H5Gclose((*st)->b765.ggid);
	H5Fclose((*st)->b765.fid);
	
	/*  910  */
	H5Dclose((*st)->b910.iid);
	H5Dclose((*st)->b910.aid);
	H5Dclose((*st)->b910.sid);
	H5Dclose((*st)->b910.szid);
	H5Dclose((*st)->b910.vzid);
	H5Dclose((*st)->b910.said);
	H5Dclose((*st)->b910.vaid);
	H5Gclose((*st)->b910.gdid);
	H5Gclose((*st)->b910.ggid);
	H5Fclose((*st)->b910.fid);
	
	/*  490  */
	H5Dclose((*st)->b490.iid);
	H5Dclose((*st)->b490.qid);
	H5Dclose((*st)->b490.uid);
	H5Dclose((*st)->b490.aid);
	H5Dclose((*st)->b490.sid);
	H5Dclose((*st)->b490.szid);
	H5Dclose((*st)->b490.vzid);
	H5Dclose((*st)->b490.said);
	H5Dclose((*st)->b490.vaid);
	H5Gclose((*st)->b490.gdid);
	H5Gclose((*st)->b490.ggid);
	H5Fclose((*st)->b490.fid);
	
	/*  670  */
	H5Dclose((*st)->b670.iid);
	H5Dclose((*st)->b670.qid);
	H5Dclose((*st)->b670.uid);
	H5Dclose((*st)->b670.aid);
	H5Dclose((*st)->b670.sid);
	H5Dclose((*st)->b670.szid);
	H5Dclose((*st)->b670.vzid);
	H5Dclose((*st)->b670.said);
	H5Dclose((*st)->b670.vaid);
	H5Gclose((*st)->b670.gdid);
	H5Gclose((*st)->b670.ggid);
	H5Fclose((*st)->b670.fid);
	
	/*  865  */
	H5Dclose((*st)->b865.iid);
	H5Dclose((*st)->b865.qid);
	H5Dclose((*st)->b865.uid);
	H5Dclose((*st)->b865.aid);
	H5Dclose((*st)->b865.sid);
	H5Dclose((*st)->b865.szid);
	H5Dclose((*st)->b865.vzid);
	H5Dclose((*st)->b865.said);
	H5Dclose((*st)->b865.vaid);
	H5Gclose((*st)->b865.gdid);
	H5Gclose((*st)->b865.ggid);
	H5Fclose((*st)->b865.fid);
	
	H5Sclose((*st)->h2id);
	H5Sclose((*st)->h3id);
	H5Sclose((*st)->m2id);
	H5Sclose((*st)->m3id);
	
	cpt_freethemall(3, &(*st)->lat, &(*st)->lon, st);
	
	return 0;
}

/*
 *  CURL part aims at downloading aeronet site file automatically.
 *  fn write_data: callback of write fn when there is data received.
 */
static size_t write_data(char *ptr, size_t size, size_t nmemb, void *userdata)
{
	return fwrite(ptr, size, nmemb, (FILE *)userdata);
}

/*
 *  CURL part aims at downloading aeronet site file automatically.
 *  fn downpt: download text from aeronet website.
 */
static int downpt(struct wr_cpt_dpc *dpcst, const char *ptfname)
{
	int   fd;
	char *url, *dt, *dt1, *dt2, *bound;
	FILE *fp;
	CURL *curl_handle;
	time_t sec;
	struct tm *stm;
	
	/*  Get URL  */
	if ((fd = open(WR_CPT_DTFNAME, O_RDONLY)) < 0) {
		CPT_ERROPEN(WR_CPT_DTFNAME);
		return WR_CPT_EOPEN;
	}
	dt = malloc(sizeof(int8_t[WR_CPT_URLDTLEN]));
	read(fd, dt, WR_CPT_URLDTLEN);
	close(fd);
	
	sec = dpcst->secswhenscan - (time_t)60*30;
	stm = gmtime(&sec);
	dt1 = malloc(sizeof(int8_t[WR_CPT_URLDT1LEN]));
	strftime(dt1, UINT8_MAX, "year=%Y&month=%m&day=%d&hour=%H", stm);
	
	sec = dpcst->secswhenend  + (time_t)60*90;
	stm = gmtime(&sec);
	dt2 = malloc(sizeof(int8_t[WR_CPT_URLDT2LEN]));
	strftime(dt2, UINT8_MAX, "year2=%Y&month2=%m&day2=%d&hour2=%H", stm);
	
	asprintf(&bound, "lon1=%07.2f&lon2=%07.2f&lat1=%06.2f&lat2=%06.2f",
	         dpcst->lonmin, dpcst->lonmax, dpcst->latmin, dpcst->latmax);
	
	asprintf(&url, WR_CPT_URLPREFIX"%s&%s&%s&%s&%s",
	         WR_CPT_URLAVG, dt, dt1, dt2, bound);
	
	/*  CURL part  */
	curl_global_init(CURL_GLOBAL_ALL);
	
	/*  Init the curl session  */
	curl_handle = curl_easy_init();
	
	/*  Set URL to get here  */
	curl_easy_setopt(curl_handle, CURLOPT_URL, url);

#ifdef CPT_DEBUG
	CPT_ECHOWITHTIME("download %s from %s", ptfname, url);
#endif
	
	/*  Switch on/off full protocol/debug output while testing  */
	curl_easy_setopt(curl_handle, CURLOPT_VERBOSE, 0L);
	
	/*  Disable/enable progress meter, set to 0L to enable it  */
#ifdef CPT_DEBUG
	curl_easy_setopt(curl_handle, CURLOPT_NOPROGRESS, 0L);
#else
	curl_easy_setopt(curl_handle, CURLOPT_NOPROGRESS, 1L);
#endif
	
	/*  Send all data to this function  */
	curl_easy_setopt(curl_handle, CURLOPT_WRITEFUNCTION, write_data);
	
	/*  File already exist ?  */
	fp = fopen(ptfname, "wb");
	
	/*  Write the page body to this file handle  */
	curl_easy_setopt(curl_handle, CURLOPT_WRITEDATA, fp);
	
	/*  Get it!  */
	curl_easy_perform(curl_handle);
	
	/*  Close the header file  */
	fclose(fp);
	
	/*  Cleanup of heap  */
	cpt_freethemall(5, &dt, &url, &dt1, &dt2, &bound);
	
	/*  Cleanup curl stuff  */
	curl_easy_cleanup(curl_handle);
	curl_global_cleanup();
	
	return 0;
}

/*
 *  Load site info into cpt_pt st
 */
static int querysda(const char *fname, struct cpt_pt **allpt, uint32_t *ptcount)
{
	int      fd;
	char    *buffer, *pbuf, *pprev, *line;
	char     rcddate[WR_CPT_DATELEN], rcdtime[WR_CPT_TIMELEN], rcddt[WR_CPT_DTLEN];
	float    lon, lat, prevlon, prevlat;
	double   tparams[WR_CPT_NPARAM];
	int16_t  alt;
	uint32_t flen, llen, maxlen;
	struct cpt_point *ppoint;
	struct cpt_pt *ppt;
	struct tm stm;
	
	/*  Try openning text file  */
	if ((fd = open(fname, O_RDONLY)) < 0) {
		CPT_ERROPEN(fname);
		return WR_CPT_EOPEN;
	}
	
	/*  Bad content  */
	if ((flen = lseek(fd, 0, SEEK_END)) < WR_CPT_INVSITEBYTES) {
		CPT_ERRECHOWITHTIME("%s contains NO valid site info", fname);
		close(fd);
		return WR_CPT_EINVSITE;
	}
	
	/*  Buffer all content  */
	buffer = malloc(sizeof(char[flen]));
	if (!buffer) {
		CPT_ERRMEM(buffer);
		return WR_CPT_EMEM;
	}
	pread(fd, buffer, flen, 0);
	pbuf = buffer;
	close(fd);
	
	/*  Points to line below title  */
	do {
		while (*++pbuf != '\n');
	} while (*++pbuf != 'A');
	while (*++pbuf != '\n');
	pprev = ++pbuf;
	
	line   = NULL;
	maxlen = 1;  /*  max length of text line  */
	
	prevlon = WR_CPT_LONLIM_MAX+1;
	prevlat = WR_CPT_LATLIM_MAX+1;
	
	*allpt   = NULL;
	*ptcount = 0;
	
	/*  Handle each line  */
	do {
		while (*++pbuf != '\n') ;
		
		llen = pbuf-pprev;
		if (llen > maxlen)
			line = realloc(line, sizeof(char[maxlen = llen]));
		memcpy(line, pprev, llen);
		
		/*  End of records  */
		if ('<' == *line) break;
		
		/*  Read info into stack  */
		sscanf(line, "%*[^','],"
		             "%[^','],%[^','],"  /*  A.read date time  */
		             "%*[^','],%*[^','],"
		             "%lf,%lf,"          /*  B.read total and fine AOD500  */
		             "%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],"
		             "%lf,"              /*  C.read AE  */
		             "%*[^','],"
		             "%lf,"              /*  D.read fine AE  */
		             
		             /*  Many skip  */
		             "%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],"
		             "%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],"
		             "%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],"
		             "%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],"
		             "%*[^','],"
		             
		             "%f,%f,%hd",        /*  E.read geolocation  */
		             
		             rcddate, rcdtime,               /*  A.receive  */
		             tparams, tparams+1,             /*  B.receive  */
		             tparams+2,                      /*  C.receive  */
		             tparams+3,                      /*  D.receive  */
		             &lat, &lon, &alt);              /*  E.receive  */
		
		/*  Site of current line first appear  */
		if ((lon != prevlon) && (lat != prevlat)) {
			prevlon = lon;
			prevlat = lat;
			
			*allpt = realloc(*allpt, sizeof(struct cpt_pt[++*ptcount]));
			ppt = *allpt+*ptcount-1;
			ppt->nt  = 1;
			ppt->lon = lon;
			ppt->lat = lat;
			ppt->alt = alt;
			
			ppt->points = malloc(CPT_POINTSIZE);
			ppoint = ppt->points;
		} else {  /*  site appeared before  */
			ppt->points = realloc(ppt->points, CPT_POINTSIZE*(++ppt->nt));
			ppoint = ppt->points+ppt->nt-1;
		}
		
		ppoint->params = malloc(_cpt_parsz);
		memcpy(ppoint->params, tparams, _cpt_parsz);
		
		snprintf(rcddt, WR_CPT_DTLEN, "%s,%s", rcddate, rcdtime);
		strptime(rcddt, "%d:%m:%Y,%H:%M:%S", &stm);
		ppoint->seconds = timegm(&stm);
		
		pprev = ++pbuf;
	} while (*pbuf != '\0');
	
	/*  Set pointers to NULL  */
	pbuf = pprev = NULL;
	ppoint = NULL;
	ppt = NULL;
	
	/*  Free allocated buffer  */
	line = realloc(line, sizeof(char[maxlen]));
	cpt_freethemall(2, &buffer, &line);
	
	return 0;
}

/*  Set hyper space before read partial DPC data  */
static int setdpchyper(struct wr_cpt_dpc *st, uint16_t ir, uint16_t ic)
{
	return
	H5Sselect_hyperslab(st->h3id, H5S_SELECT_SET, (hsize_t[3]) {0,ir,ic}, NULL, st->l3id, NULL) < 0
	||
	H5Sselect_hyperslab(st->h2id, H5S_SELECT_SET, (hsize_t[2]) {ir,ic}, NULL, st->l2id, NULL) < 0;
}

/*  Load certain channel  */
static int loadchannel(struct wr_cpt_dpc *st, struct wr_cpt_dpcband *pb, struct cpt_channel *pc,
                       uint8_t *mask, int16_t *alt, uint8_t *masknset, uint8_t *altnset)
{
	uint8_t   ilayer;
	int16_t  *obs;
	uint16_t *ang;
	
	if (*masknset) {
		H5Dread(pb->sid, H5T_NATIVE_UINT8, st->m2id,
		        st->h2id, H5P_DEFAULT, mask);
		if (WR_CPT_VALIDMASK(*mask))
			*masknset = 0;
	}
	if (*altnset) {
		H5Dread(pb->aid, H5T_NATIVE_INT16, st->m2id,
		        st->h2id, H5P_DEFAULT, alt);
		if (WR_CPT_VALIDALT(*alt))
			*altnset = 0;
	}
	
	ang = malloc(sizeof(int16_t[st->nlayer]));
	
	H5Dread(pb->szid, H5T_NATIVE_UINT16, st->m3id, st->h3id, H5P_DEFAULT, ang);
	for (ilayer = 0; ilayer < st->nlayer; ++ilayer)
		pc->ang[ilayer] = st->scaleang * ang[ilayer];
	
	H5Dread(pb->vzid, H5T_NATIVE_UINT16, st->m3id, st->h3id, H5P_DEFAULT, ang);
	for (ilayer = 0; ilayer < st->nlayer; ++ilayer)
		pc->ang[ilayer + st->nlayer] = st->scaleang * ang[ilayer];
	
	H5Dread(pb->said, H5T_NATIVE_UINT16, st->m3id, st->h3id, H5P_DEFAULT, ang);
	for (ilayer = 0; ilayer < st->nlayer; ++ilayer)
		pc->ang[ilayer + 2*st->nlayer] = st->scaleang * ang[ilayer];
	
	H5Dread(pb->vaid, H5T_NATIVE_UINT16, st->m3id, st->h3id, H5P_DEFAULT, ang);
	for (ilayer = 0; ilayer < st->nlayer; ++ilayer)
		pc->ang[ilayer + 3*st->nlayer] = st->scaleang * ang[ilayer];
	
	CPT_FREE(ang);
	obs = malloc(sizeof(int16_t[st->nlayer]));
	
	H5Dread(pb->iid, H5T_NATIVE_INT16, st->m3id, st->h3id, H5P_DEFAULT, obs);
	for (ilayer = 0; ilayer < st->nlayer; ++ilayer)
		pc->obs[ilayer] = st->scaleobs * obs[ilayer];
	
	CPT_FREE(obs);
	
	return 0;
}

/*  Load certain polar channel  */
static int loadchannelp(struct wr_cpt_dpc *st, struct wr_cpt_dpcbandp *pb, struct cpt_channel *pc,
                        uint8_t *mask, int16_t *alt, uint8_t *masknset, uint8_t *altnset)
{
	uint8_t   ilayer;
	int16_t  *obs;
	uint16_t *ang;
	
	if (*masknset) {
		H5Dread(pb->sid, H5T_NATIVE_UINT8, st->m2id,
		        st->h2id, H5P_DEFAULT, mask);
		if (WR_CPT_VALIDMASK(*mask))
			*masknset = 0;
	}
	if (*altnset) {
		H5Dread(pb->aid, H5T_NATIVE_INT16, st->m2id,
		        st->h2id, H5P_DEFAULT, alt);
		if (WR_CPT_VALIDALT(*alt))
			*altnset = 0;
	}
	
	ang = malloc(sizeof(int16_t[st->nlayer]));
	
	H5Dread(pb->szid, H5T_NATIVE_UINT16, st->m3id, st->h3id, H5P_DEFAULT, ang);
	for (ilayer = 0; ilayer < st->nlayer; ++ilayer)
		pc->ang[ilayer] = st->scaleang * ang[ilayer];
	
	H5Dread(pb->vzid, H5T_NATIVE_UINT16, st->m3id, st->h3id, H5P_DEFAULT, ang);
	for (ilayer = 0; ilayer < st->nlayer; ++ilayer)
		pc->ang[ilayer + st->nlayer] = st->scaleang * ang[ilayer];
	
	H5Dread(pb->said, H5T_NATIVE_UINT16, st->m3id, st->h3id, H5P_DEFAULT, ang);
	for (ilayer = 0; ilayer < st->nlayer; ++ilayer)
		pc->ang[ilayer + 2*st->nlayer] = st->scaleang * ang[ilayer];
	
	H5Dread(pb->vaid, H5T_NATIVE_UINT16, st->m3id, st->h3id, H5P_DEFAULT, ang);
	for (ilayer = 0; ilayer < st->nlayer; ++ilayer)
		pc->ang[ilayer + 3*st->nlayer] = st->scaleang * ang[ilayer];
	
	CPT_FREE(ang);
	obs = malloc(sizeof(int16_t[st->nlayer]));
	
	H5Dread(pb->iid, H5T_NATIVE_INT16, st->m3id, st->h3id, H5P_DEFAULT, obs);
	for (ilayer = 0; ilayer < st->nlayer; ++ilayer)
		pc->obs[ilayer] = st->scaleobs * obs[ilayer];
	
	H5Dread(pb->qid, H5T_NATIVE_INT16, st->m3id, st->h3id, H5P_DEFAULT, obs);
	for (ilayer = 0; ilayer < st->nlayer; ++ilayer)
		pc->obs[ilayer + st->nlayer] = st->scaleobs * obs[ilayer];
	
	H5Dread(pb->uid, H5T_NATIVE_INT16, st->m3id, st->h3id, H5P_DEFAULT, obs);
	for (ilayer = 0; ilayer < st->nlayer; ++ilayer)
		pc->obs[ilayer + 2*st->nlayer] = st->scaleobs * obs[ilayer];
	
	CPT_FREE(obs);
	
	return 0;
}

/*  Load certain location DPC pixel  */
static int loadpxfromst(struct cpt_pixel *pixel, struct wr_cpt_dpc *st, uint16_t ir, uint16_t ic)
{
	uint8_t  ispolar, masknset, altnset;
	struct cpt_channel *pchannel;
	
	const uint32_t idx = (uint32_t) ir*st->ncol+ic;
	const void *bands[WR_CPT_DPCNBANDS] = {&st->b443, &st->b490, &st->b565, &st->b670,
	                                       &st->b763, &st->b765, &st->b865, &st->b910};
	
	pixel->lon = st->lon[idx];
	pixel->lat = st->lat[idx];
	
	pixel->nlayer   = st->nlayer;
	pixel->nchannel = WR_CPT_DPCNBANDS;
	
	/*  Set hyperslab  */
	setdpchyper(st, ir, ic);
	masknset = altnset = 1;
	
	/*  Load channels  */
	pixel->channels = malloc(sizeof(struct cpt_channel[pixel->nchannel]));
	for (uint16_t channel = 0; channel < pixel->nchannel; ++channel) {
		pchannel = pixel->channels+channel;
		pchannel->centrewv = WR_CPT_DPCCNTRWV[channel];
		ispolar = pchannel->centrewv < 0;
		
		/*
		 *  Allocate memory
		 *  4 stands for sz, vz, sa and va
		 *  3 stands for I/Q/U
		 */
		pchannel->ang = malloc(sizeof(double[pixel->nlayer][4]));
		pchannel->obs = malloc(sizeof(double[pixel->nlayer][ispolar ? 3 : 1]));
		
		if (ispolar) {
			loadchannelp(st, (struct wr_cpt_dpcbandp *) bands[channel],
			             pchannel, &pixel->mask, &pixel->alt, &masknset, &altnset);
		} else {
			loadchannel(st, (struct wr_cpt_dpcband *) bands[channel],
			            pchannel, &pixel->mask, &pixel->alt, &masknset, &altnset);
		}
	}
	
	pchannel = NULL;

	return 0;
}

/*  Pairing Pt and Px  */
static uint32_t pairdpc(struct cpt_pt *allpt, uint32_t ptcount, struct wr_cpt_dpc *dpcst,
                        struct cpt_px **pairpx, struct cpt_pt **pairpt)
{
	uint8_t  ipoint, npoint, pointsta, ivicinity,
	         rowntop, rownbottom, colnleft, colnright;
	int16_t  rowv, colv;
	uint16_t row, col;
	uint32_t ipt, idx, ptxcount;
	uint64_t sec;
	
	struct cpt_pt    *ppt,    *ppairpt;
	struct cpt_px    *ppx;
	struct cpt_point *ppoint, *ppairpoint;
	
	const uint16_t rowlimit = dpcst->nrow-1,
	               collimit = dpcst->ncol-1;
	const float rowcoef = rowlimit / 2.f / WR_CPT_LATLIM_MAX,
	            colcoef = collimit / 2.f / WR_CPT_LONLIM_MAX;
	
	ptxcount = 0;
	*pairpx  = NULL;
	*pairpt  = NULL;
	for (ipt = 0; ipt < ptcount; ++ipt) {
		ppt = allpt+ipt;
		row = rowcoef * (WR_CPT_LATLIM_MAX-ppt->lat);
		col = colcoef * (WR_CPT_LONLIM_MAX+ppt->lon);
		idx = (uint32_t) row * dpcst->ncol + col;
		
		if (WR_CPT_INVALIDGEO(dpcst->lon[idx], dpcst->lat[idx]))
			goto next_pt;
		
		sec = dpcst->secswhenscan + (rowlimit-row)*dpcst->secsperline;
		
		npoint   = 0;
		pointsta = UINT8_MAX;
		for (ipoint = 0; ipoint < ppt->nt; ++ipoint) {
			ppoint = ppt->points+ipoint;
			if (((ppoint->seconds > sec) ?
			(ppoint->seconds-sec) : (sec-ppoint->seconds)) < WR_CPT_SECDIFFMAX) {
				++npoint;
				if (UINT8_MAX == pointsta)
					pointsta = ipoint;
			}
		}
		
		if (!npoint)
			goto next_pt;
		
		/*  Pt  */
		*pairpt = realloc(*pairpt, sizeof(struct cpt_pt[++ptxcount]));
		ppairpt = *pairpt + ptxcount-1;
		ppairpt->nt  = npoint;
		ppairpt->alt = ppt->alt;
		ppairpt->lon = ppt->lon;
		ppairpt->lat = ppt->lat;
		ppairpt->points = malloc(sizeof(struct cpt_point[ppairpt->nt]));
		for (ipoint = 0; ipoint < npoint; ++ipoint) {
			ppoint = ppt->points + pointsta + ipoint;
			ppairpoint = ppairpt->points + ipoint;
			ppairpoint->seconds = ppoint->seconds;
			ppairpoint->params  = malloc(_cpt_parsz);
			memcpy(ppairpoint->params, ppoint->params, _cpt_parsz);
		}
		
		/*  Px  */
		*pairpx = realloc(*pairpx, sizeof(struct cpt_px[ptxcount]));
		ppx = *pairpx + ptxcount-1;
		ppx->centrepixel = malloc(CPT_PIXELSIZE);
		ppx->seconds     = sec;
		
		/*  Centre pixel  */
		loadpxfromst(ppx->centrepixel, dpcst, row, col);
		
		/*  Vicinity count  */
		rowntop    = (row != 0);
		rownbottom = (row != rowlimit);
		colnleft   = (col != 0);
		colnright  = (col != collimit);
		switch (rowntop+rownbottom+colnleft+colnright) {
		case 4: ppx->nvicinity = 8; break;
		case 3: ppx->nvicinity = 5; break;
		case 2: ppx->nvicinity = 3; break;
		default: ppx->nvicinity = 0;
		}
		
		/*  Vicinity memory manage  */
		if (ppx->nvicinity)
			ppx->vicinity = malloc(sizeof(struct cpt_pixel[ppx->nvicinity]));
		else
			ppx->vicinity = NULL;
		
		/*  Vicinity assignment  */
		ivicinity = 0;
		for (rowv = -1; rowv < 2; ++rowv) {
			if (((-1 == rowv) && !rowntop) || ((1 == rowv) && !rownbottom))
				continue;
			for (colv = -1; colv < 2; ++colv) {
				if (((-1 == colv) && !colnleft) || ((1 == colv) && !colnright))
					continue;
				/*  Centre pixel already load before  */
				if ((0 == rowv) && (0 == colv))
					continue;
				loadpxfromst(ppx->vicinity+ivicinity++, dpcst, row+rowv, col+colv);
			}
		}
		
		next_pt:
		continue;
	}
	
	ppx = NULL;
	ppt = ppairpt = NULL;
	ppoint = ppairpoint = NULL;
	
	return ptxcount;
}

/*  Safer implementation of write fn  */
static ssize_t safewrite(int filedes, const void *buffer, size_t size)
{
	ssize_t ret;
	
	while ((ret = write(filedes, buffer, size)) > 0) {
		if (ret >= size) {
			return ret-size;
		} else {
			buffer += ret;
			size   -= ret;
		}
	}
	
	return ret;
}

/*  Write pixel individual  */
static int writepixel(int filedes, struct cpt_pixel *pixel)
{
	struct cpt_channel *pchannel;

	/*  Geolocation  */
	safewrite(filedes, &pixel->lon, _cpt_4byte);
	safewrite(filedes, &pixel->lat, _cpt_4byte);
	safewrite(filedes, &pixel->alt, _cpt_2byte);
	safewrite(filedes, &pixel->mask, _cpt_1byte);
	
	/*  Dimensions  */
	safewrite(filedes, &pixel->nchannel, _cpt_1byte);
	safewrite(filedes, &pixel->nlayer, _cpt_1byte);
	
	/*  Channel  */
	for (uint8_t ichannel = 0; ichannel < pixel->nchannel; ++ichannel) {
		pchannel = pixel->channels+ichannel;
		safewrite(filedes, &pchannel->centrewv, _cpt_2byte);
		safewrite(filedes, pchannel->obs,
		          sizeof(double[pixel->nlayer][(pchannel->centrewv < 0) ? 3 : 1]));
		safewrite(filedes, pchannel->ang, sizeof(double[pixel->nlayer][4]));
	}
	pchannel = NULL;

	return 0;
}

/*  Export struct to file  */
static int writecpttofile(const char *fname, struct cpt_ff *st)
{
	int fd;
	uint8_t  ipoint, ivicinity;
	uint32_t iptx;
	struct cpt_pt *ppt;
	struct cpt_px *ppx;
	struct cpt_point *ppoint;
	
	/*  File already exist ?  */
	fd = open(fname, O_PATH);
	if (fd > 0) {
		close(fd);
#ifdef CPT_DEBUG
		fd = open(fname, O_WRONLY|O_CREAT, S_IRUSR|S_IWUSR);
#else
		return EEXIST;
#endif
	} else {
		fd = open(fname, O_WRONLY|O_CREAT, S_IRUSR|S_IWUSR);
	}
	
	/*  Header  */
	safewrite(fd, st->hdr->magic_number, CPT_MAGICLEN);
	safewrite(fd, &st->hdr->ver, _cpt_1byte);
	safewrite(fd, &st->hdr->nptx, _cpt_4byte);
	safewrite(fd, &st->hdr->nparam, _cpt_1byte);
	
	/*  Data/Ptx  */
	for (iptx = 0; iptx < st->hdr->nptx; ++iptx) {
		
		/*  Pt  */
		ppt = st->data->pt+iptx;
		safewrite(fd, &ppt->lon, _cpt_4byte);
		safewrite(fd, &ppt->lat, _cpt_4byte);
		safewrite(fd, &ppt->alt, _cpt_2byte);
		safewrite(fd, &ppt->nt , _cpt_1byte);
		for (ipoint = 0; ipoint < ppt->nt; ++ipoint) {
			ppoint = ppt->points+ipoint;
			safewrite(fd, &ppoint->seconds, _cpt_8byte);
			safewrite(fd, ppoint->params, _cpt_parsz);
		}
		
		/*  Px  */
		ppx = st->data->px+iptx;
		safewrite(fd, &ppx->seconds, _cpt_8byte);
		writepixel(fd, ppx->centrepixel);
		safewrite(fd, &ppx->nvicinity, _cpt_1byte);
		for (ivicinity = 0; ivicinity < ppx->nvicinity; ++ivicinity) {
			writepixel(fd, ppx->vicinity+ivicinity);
		}
		
	}
	
	/*  Ending  */
	safewrite(fd, st->ending, CPT_ENDINGLEN);
	close(fd);
	
	ppt = NULL;
	ppx = NULL;
	ppoint = NULL;

	return 0;
}

/*  Definition of main function  */
int wrcpt(const char *prefix, const char *ptxtfname, const char *cptfname)
{
	int ret;
	uint32_t ptcount, ptxcount;
	struct cpt_pt *allpt  = NULL,
	              *pairpt = NULL;
	struct cpt_px *pairpx = NULL;
	struct cpt_ff  cptout;
	struct cpt_ptx ptx;
	struct cpt_header  hdr;
	struct wr_cpt_dpc *dpcst;
	
	/*  Px prepare  */
	initdpcst(prefix, &dpcst);
	
	/*  Whether download site info locally  */
	if (cptfname) {
		int fd = open(ptxtfname, O_PATH);
		if (fd > 0) {
			close(fd);
			CPT_ECHOWITHTIME("Mode: using local ptxt %s to create cpt %s",
			                 ptxtfname, cptfname);
		} else {
			CPT_ECHOWITHTIME("Mode: download remote ptxt %s to create cpt %s",
			                 ptxtfname, cptfname);
			downpt(dpcst, ptxtfname);
		}
	} else {
		downpt(dpcst, ptxtfname);
		goto cleanup;
	}
	
	/*  All Pt load  */
	if ((ret = querysda(ptxtfname, &allpt, &ptcount))) {
		cpt_freeptall(&allpt, ptcount);
		return ret;
	}
	
	/*  Exit when there's no record in Pt  */
	if (!ptcount) {
		cpt_freeptall(&allpt, ptcount);
		cleandpcst(&dpcst);
		return WR_CPT_ENORCD;
	}
	
	/*  Get paired Px and Pt  */
	ptxcount = pairdpc(allpt, ptcount, dpcst, &pairpx, &pairpt);
#ifdef CPT_DEBUG
	for (uint32_t i = 0; i < ptxcount; ++i) {
		CPT_ECHOWITHTIME("No.%03d: lon %9.4f lat %8.4f with %2d points",
		                 i+1, (pairpx+i)->centrepixel->lon,
		                 (pairpx+i)->centrepixel->lat, (pairpt+i)->nt);
	}
#endif
	
	/*  Construct cpt  */
	hdr.ver    = CPT_VERSION;
	hdr.nptx   = ptxcount;
	hdr.nparam = WR_CPT_NPARAM;
	hdr.magic_number = CPT_MAGIC;
	ptx.pt = pairpt;
	ptx.px = pairpx;
	cptout.hdr    = &hdr;
	cptout.data   = &ptx;
	cptout.ending = CPT_ENDING;
	
	/*  Write to file  */
	writecpttofile(cptfname, &cptout);
	
	/*  Cleanup  */
	cleanup:
	ret = cpt_freeptall(&allpt, ptcount);
	ret = cpt_freepxall(&pairpx, ptxcount);
	ret = cpt_freeptall(&pairpt, ptxcount);
	ret = cleandpcst(&dpcst);
	
	ptx.pt = NULL;
	ptx.px = NULL;
	cptout.hdr  = NULL;
	cptout.data = NULL;
	
	return 0;
}

