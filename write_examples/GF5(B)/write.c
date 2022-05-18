/*
 *file: write_examples/GF5(B)/write.c
 *descreption:
 *  write cpt format file
 *arguements:
 *  [1]: downloaded in-situ aeronet file
 *  [2]: DPC/POSP data file
 *  [3]: output cpt
 *init date: May/10/2022
 *last modify: May/18/2022
 *
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <fcntl.h>
#include <math.h>

#include "hdf5.h"
#include "../../read/readcpt.h"


/*  Site info settings  */
#define WR_CPT_INVSITEBYTES 180
#define WR_CPT_DATELEN      11
#define WR_CPT_TIMELEN      9
#define WR_CPT_DTLEN        (WR_CPT_DATELEN+WR_CPT_TIMELEN)
#define WR_CPT_NPARAM       8


/*  Satellite info settings  */
#define WR_CPT_POSPNBANDS  9
#define WR_CPT_POSPDSNAME  "Data_Fields"
#define WR_CPT_POSPGLNAME  "Geolocation_Fields"
#define WR_CPT_POSPINAME   "I"
#define WR_CPT_POSPQNAME   "Q"
#define WR_CPT_POSPUNAME   "U"
#define WR_CPT_POSPLONNAME "Longitude"
#define WR_CPT_POSPLATNAME "Latitude"
#define WR_CPT_POSPSLNAME  "Sea_Land_Flags"
#define WR_CPT_POSPALTNAME "Surface_Altitude"
#define WR_CPT_POSPSZNAME  "Sol_Zen_Ang"
#define WR_CPT_POSPVZNAME  "View_Zen_Ang"
#define WR_CPT_POSPSANAME  "Sol_Azim_Ang"
#define WR_CPT_POSPVANAME  "View_Azim_Ang"
#define WR_CPT_POSPCNTRWV  (int16_t[WR_CPT_POSPNBANDS]) \
                           {-380, -410, -443, -490, -670, -865, -1380, -1610, -2250}

#define WR_CPT_XMLSUFFIX   "xml"
#define WR_CPT_XMLSUFLEN   strlen(WR_CPT_XMLSUFFIX)
#define WR_CPT_XMLENDTAG   "</ProductMetaData>"
#define WR_CPT_XMLSTTAG    "StartTime"
#define WR_CPT_XMLETTAG    "EndTime"
#define WR_CPT_XMLLONTAG   "NadirLong"
#define WR_CPT_XMLLATTAG   "NadirLat"

#define WR_CPT_LONLIM_MIN  -180
#define WR_CPT_LONLIM_MAX  180
#define WR_CPT_LATLIM_MIN  -90
#define WR_CPT_LATLIM_MAX  90

#define WR_CPT_INVALIDLON(lon) ((lon < WR_CPT_LONLIM_MIN) || (lon > WR_CPT_LONLIM_MAX))
#define WR_CPT_INVALIDLAT(lat) ((lat < WR_CPT_LATLIM_MIN) || (lat > WR_CPT_LATLIM_MAX))

#define WR_CPT_SECDIFFMAX ((uint64_t) 60*30)  /*  30 min  */

struct wr_cpt_posp {
	hid_t iid;   /*  entrance of intensity    */
	hid_t qid;   /*  entrance of polar q      */
	hid_t uid;   /*  entrance of polar u      */
	hid_t szid;  /*  entrance of sol zen ang  */
	hid_t vzid;  /*  entrance of sat zen ang  */
	hid_t said;  /*  entrance of sol azi ang  */
	hid_t vaid;  /*  entrance of sat azi ang  */
	hid_t aid;   /*  entrance of altitude     */
	hid_t sid;   /*  entrance of sea-land     */
	
	hid_t h2id;  /*  2-d hyperslab  */
	hid_t h3id;  /*  3-d hyperslab  */
	hid_t msid;  /*  mem space id   */
	
	hid_t gdid;  /*  dateset group id  */
	hid_t ggid;  /*  geo-loc group id  */
	hid_t fid;   /*  file id           */
	
	hsize_t l2id[2];  /*  2-d hyper len  */
	hsize_t l3id[3];  /*  3-d hyper len  */
	
	uint16_t ncol, nrow;    /*  count of row/col    */
	float    secsperline;   /*  timelapse per line  */
	uint64_t secswhenscan;  /*  timestamp of begin  */
	uint64_t secswhenend;   /*  timestamp of end    */
	double  *lon, *lat;     /*  full lon/lat        */
	
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
int wrcpt(const char *ptname, const char *pxname, const char *cptname);

/*  Program entrance  */
int main(int argc, char *argv[]) {
	if (4 == argc) {
		return wrcpt(argv[1], argv[2], argv[3]);
	} else {
		CPT_ERRECHOWITHTIME("Usage: %s sitefile satefile output", argv[0]);
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
static int cpt_freethemall(uint8_t n, ...)
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
static int cpt_freepointall(struct cpt_point **p, uint16_t n)
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
static int cpt_freeptall(struct cpt_pt **p, uint16_t n)
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
static int cpt_freechannelall(struct cpt_channel **p, uint16_t n)
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
static int cpt_freepixelall(struct cpt_pixel **p, uint16_t n)
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
static int cpt_freepxall(struct cpt_px **p, uint16_t n)
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

/*
 *  Load site info into cpt_pt st
 */
static int querypt(const char *fname, struct cpt_pt **allpt, uint16_t *ptcount)
{
	int      fd;
	char    *buffer, *pbuf, *pprev, *line;
	char     rcddate[WR_CPT_DATELEN], rcdtime[WR_CPT_TIMELEN], rcddt[WR_CPT_DTLEN];
	float    lon, lat, alt, prevlon, prevlat;
	double   tparams[WR_CPT_NPARAM];
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
		CPT_ERRECHOWITHTIME("%s contains NO valid site info\n", fname);
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
		             "%[^','],%[^','],"    /*  A.read date time  */
		             "%*[^','],%*[^','],"
		             "%lf,%lf,%lf,"        /*  B.read 1640, 1020 and 870  */
		             "%*[^','],%*[^','],"
		             "%lf,"                /*  C.read 675  */
		             "%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],"
		             "%lf,"                /*  D.read 500  */
		             "%*[^','],%*[^','],"
		             "%lf,"                /*  E.read 440  */
		             "%*[^','],%*[^','],"
		             "%lf,%lf,"            /*  F.read 380 and 340  */
		             
		             /*  Many skip  */
		             "%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],"
		             "%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],"
		             "%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],"
		             "%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],"
		             "%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],"
		             "%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],%*[^','],"
		             
		             "%f,%f,%hd",          /*  G.read geolocation  */
		             
		             rcddate, rcdtime,               /*  A.receive  */
		             tparams, tparams+1, tparams+2,  /*  B.receive  */
		             tparams+3,                      /*  C.receive  */
		             tparams+4,                      /*  D.receive  */
		             tparams+5,                      /*  E.receive  */
		             tparams+6, tparams+7,           /*  F.receive  */
		             &lat, &lon, &alt);              /*  G.receive  */
		
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
		ppoint->seconds = mktime(&stm);
		
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

/*  Init st from POSP h5  */
static int pospopenall(const char *fname, struct wr_cpt_posp *st)
{
	hid_t    latid, lonid, space;
	size_t   size;
	hsize_t  dim[2];
	uint32_t index;
	
	/*  Open h5 entrance  */
	st->fid  = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
	st->gdid = H5Gopen(st->fid, WR_CPT_POSPDSNAME, H5P_DEFAULT);
	st->ggid = H5Gopen(st->fid, WR_CPT_POSPGLNAME, H5P_DEFAULT);
	
	st->iid = H5Dopen(st->gdid, WR_CPT_POSPINAME, H5P_DEFAULT);
	st->qid = H5Dopen(st->gdid, WR_CPT_POSPQNAME, H5P_DEFAULT);
	st->uid = H5Dopen(st->gdid, WR_CPT_POSPUNAME, H5P_DEFAULT);
	
	st->aid  = H5Dopen(st->ggid, WR_CPT_POSPALTNAME, H5P_DEFAULT);
	st->sid  = H5Dopen(st->ggid, WR_CPT_POSPSLNAME , H5P_DEFAULT);
	st->szid = H5Dopen(st->ggid, WR_CPT_POSPSZNAME , H5P_DEFAULT);
	st->vzid = H5Dopen(st->ggid, WR_CPT_POSPVZNAME , H5P_DEFAULT);
	st->said = H5Dopen(st->ggid, WR_CPT_POSPSANAME , H5P_DEFAULT);
	st->vaid = H5Dopen(st->ggid, WR_CPT_POSPVANAME , H5P_DEFAULT);
	
	/*  Buffer all lat/lon  */
	latid = H5Dopen(st->ggid, WR_CPT_POSPLATNAME, H5P_DEFAULT);
	lonid = H5Dopen(st->ggid, WR_CPT_POSPLONNAME, H5P_DEFAULT);
	space = H5Dget_space(lonid);
	H5Sget_simple_extent_dims(space, dim, NULL);
	H5Sclose(space);
	st->nrow = dim[0];
	st->ncol = dim[1];
	st->lat = malloc(size = sizeof(double[st->nrow][st->ncol]));
	st->lon = malloc(size);
	H5Dread(latid, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, st->lat);
	H5Dread(lonid, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, st->lon);
	H5Dclose(latid);
	H5Dclose(lonid);
	
	/*  Additional step  */
	st->secsperline /= st->nrow;
	
	/*  Space for hyper-reading  */
	st->h2id = H5Screate_simple(2, (hsize_t[2]) {st->nrow, st->ncol}, NULL);
	st->h3id = H5Screate_simple(3, (hsize_t[3]) {WR_CPT_POSPNBANDS, st->nrow, st->ncol}, NULL);
	st->msid = H5Screate_simple(1, (hsize_t[1]) {1}, NULL);
	
	st->l3id[0] = st->l3id[1] = st->l3id[2] = 1;
	st->l2id[0] = st->l2id[1] = 1;
	
	/*  Find boundery corner  */
	st->lonmin = WR_CPT_LONLIM_MAX;
	st->lonmax = WR_CPT_LONLIM_MIN;
	st->latmin = WR_CPT_LATLIM_MAX;
	st->latmax = WR_CPT_LATLIM_MIN;
	for (index = 0; index < (uint32_t) st->nrow*st->ncol; ++index) {
		if (WR_CPT_INVALIDLON(st->lon[index]) || WR_CPT_INVALIDLAT(st->lat[index]))
			continue;
		
		if ((st->lon[index] < st->lonmin) && (st->lat[index] < st->latmin)) {
			st->lonmin = st->lon[index];
			st->latmin = st->lat[index];
		}
		if ((st->lon[index] > st->lonmax) && (st->lat[index] > st->latmax)) {
			st->lonmax = st->lon[index];
			st->latmax = st->lat[index];
		}
	}
	
	return 0;
}

/*  Set hyper space before read partial data  */
static int pospsethyper(struct wr_cpt_posp *st, uint16_t ib, uint16_t ir, uint16_t ic)
{
	return
	H5Sselect_hyperslab(st->h3id, H5S_SELECT_SET, (hsize_t[3]) {ib,ir,ic}, NULL, st->l3id, NULL) < 0
	||
	H5Sselect_hyperslab(st->h2id, H5S_SELECT_SET, (hsize_t[2]) {ir,ic}, NULL, st->l2id, NULL) < 0;
}

/*  Load certain location POSP pixel  */
static int loadpxfromst(struct cpt_pixel *pixel, struct wr_cpt_posp *st,
                        uint16_t ir, uint16_t ic, uint8_t allocate)
{
	double   tmp;
	uint8_t  ispolar;
	uint32_t idx = (uint32_t) ir*st->ncol+ic;
	struct cpt_channel *pchannel;
	
	pixel->lon = st->lon[idx];
	pixel->lat = st->lat[idx];
	
	pixel->nlayer   = 1;
	pixel->nchannel = WR_CPT_POSPNBANDS;
	
	/*  Load channel-unrelated data  */
	pospsethyper(st, 0, ir, ic);
	H5Dread(st->sid, H5T_NATIVE_DOUBLE, st->msid, st->h2id, H5P_DEFAULT, &tmp);
	pixel->mask = tmp;
	H5Dread(st->aid, H5T_NATIVE_DOUBLE, st->msid, st->h2id, H5P_DEFAULT, &tmp);
	pixel->alt  = tmp;
	
	/*  Load channels  */
	if (allocate)
		pixel->channels = malloc(sizeof(struct cpt_channel[pixel->nchannel]));
	for (uint16_t channel = 0; channel < pixel->nchannel; ++channel) {
		pchannel = pixel->channels+channel;
		pchannel->centrewv = WR_CPT_POSPCNTRWV[channel];
		ispolar = pchannel->centrewv < 0;
		
		/*
		 *  Allocate memory
		 *  4 stands for sz, vz, sa and va
		 *  3 stands for I/Q/U
		 */
		if (allocate) {
			pchannel->ang = malloc(sizeof(double[pixel->nlayer][4]));
			pchannel->obs = malloc(sizeof(double[pixel->nlayer][ispolar ? 3 : 1]));
		}
		
		/*  set hyperslab  */
		pospsethyper(st, channel, ir, ic);
		
		/*  Load scanning angles  */
		H5Dread(st->szid, H5T_NATIVE_DOUBLE, st->msid, st->h2id, H5P_DEFAULT, pchannel->ang);
		H5Dread(st->vzid, H5T_NATIVE_DOUBLE, st->msid, st->h2id, H5P_DEFAULT, pchannel->ang+1);
		H5Dread(st->said, H5T_NATIVE_DOUBLE, st->msid, st->h2id, H5P_DEFAULT, pchannel->ang+2);
		H5Dread(st->vaid, H5T_NATIVE_DOUBLE, st->msid, st->h2id, H5P_DEFAULT, pchannel->ang+3);
		
		/*  Load I/Q/U  */
		H5Dread(st->iid, H5T_NATIVE_DOUBLE, st->msid, st->h3id, H5P_DEFAULT, pchannel->obs);
		if (ispolar)
		{
			H5Dread(st->qid, H5T_NATIVE_DOUBLE, st->msid,
			        st->h3id, H5P_DEFAULT, pchannel->obs+1);
			H5Dread(st->uid, H5T_NATIVE_DOUBLE, st->msid,
			        st->h3id, H5P_DEFAULT, pchannel->obs+2);
		}
	}
	
	pchannel = NULL;

	return 0;
}

/*  Init essential info from xml  */
static int pospinfoinit(const char *fname, struct wr_cpt_posp **pospst)
{
	int      fd;
	char    *buffer, *pbuf, *pprev, *line, *pline;
	char    *xmlfname;
	char     dtbeg[WR_CPT_DTLEN], dtend[WR_CPT_DTLEN];
	size_t   xmlnlen;
	uint32_t flen, llen, maxlen;
	struct tm stm;
	
	xmlfname = malloc(xmlnlen = strlen(fname)-1);
	xmlnlen -= WR_CPT_XMLSUFLEN;
	
	memcpy(xmlfname, fname, xmlnlen);
	memcpy(xmlfname+xmlnlen, WR_CPT_XMLSUFFIX, WR_CPT_XMLSUFLEN);
	
	/*  Try openning text file  */
	if ((fd = open(xmlfname, O_RDONLY)) < 0) {
		CPT_ERROPEN(xmlfname);
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
	maxlen  = 1;
	*pospst = malloc(sizeof(struct wr_cpt_posp));
	
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
		if (pline = strstr(line, WR_CPT_XMLENDTAG)) break;
		
		/*  Start time  */
		if (pline = strstr(line, WR_CPT_XMLSTTAG)) {
			sscanf(pline, "%*[^'>']>%[^'<']", dtbeg);
			strptime(dtbeg, "%Y-%m-%d %H:%M:%S", &stm);
			(*pospst)->secswhenscan = mktime(&stm);
			goto next_line;
		}
		
		/*  Ending time  */
		if (pline = strstr(line, WR_CPT_XMLETTAG)) {
			sscanf(pline, "%*[^'>']>%[^'<']", dtbeg);
			strptime(dtbeg, "%Y-%m-%d %H:%M:%S", &stm);
			(*pospst)->secswhenend = mktime(&stm);
			(*pospst)->secsperline = (*pospst)->secswhenend - (*pospst)->secswhenscan;
			goto next_line;
		}
		
		/*  Longitude bounds  */
		if (pline = strstr(line, WR_CPT_XMLLONTAG)) {
			//TODO
			goto next_line;
		}
		
		/*  Latitude bounds  */
		if (pline = strstr(line, WR_CPT_XMLLATTAG)) {
			//TODO
			goto next_line;
		}
		
		/*
		 *  Number of col/row inside POSP XML does NOT corresponding
		 *  the self describing dimension inside its hdf5, sometimes.
		 *  So this part of code is deprecated.
		 */
		#if 0
		/*  N of row  */
		if (pline = strstr(line, WR_CPT_XMLNRTAG)) {
			sscanf(pline, "%*[^'>']>%hu", &pospst->nrow);
			pospst->secsperline /= pospst->nrow;
			if (!pospst->secsperline) pospst->secsperline = 1;
			goto next_line;
		}
		
		/*  N of col  */
		if (pline = strstr(line, WR_CPT_XMLNCTAG)) {
			sscanf(pline, "%*[^'>']>%hu", &pospst->ncol);
			goto next_line;
		}
		#endif
		next_line:
		pprev = ++pbuf;
	} while (*pbuf != '\0');
	
	/*  Set pointers to NULL  */
	pbuf = pprev = pline = NULL;
	
	/*  Free allocated buffer  */
	line = realloc(line, sizeof(char[maxlen]));
	cpt_freethemall(3, &buffer, &line, &xmlfname);
	
	return 0;
}

/*  Destroy POSP st  */
static int pospcleanst(struct wr_cpt_posp **st)
{
	H5Dclose((*st)->iid);
	H5Dclose((*st)->qid);
	H5Dclose((*st)->uid);
	H5Dclose((*st)->aid);
	H5Dclose((*st)->sid);
	H5Dclose((*st)->szid);
	H5Dclose((*st)->vzid);
	H5Dclose((*st)->said);
	H5Dclose((*st)->vaid);
	
	H5Sclose((*st)->h2id);
	H5Sclose((*st)->h3id);
	H5Sclose((*st)->msid);
	
	H5Gclose((*st)->gdid);
	H5Gclose((*st)->ggid);
	
	H5Fclose((*st)->fid);
	
	cpt_freethemall(3, &(*st)->lat, &(*st)->lon, st);
	
	return 0;
}

/*  Pairing Pt and Px  */
static uint16_t posppair(struct cpt_pt *allpt, uint16_t ptcount, struct wr_cpt_posp *pospst,
                         struct cpt_px **pairpx, struct cpt_pt **pairpt)
{
	float    lonres, latres, diff, diffmin,
	         ptgeodiff[ptcount], londiff, latdiff;
	uint8_t  appendpx, ipoint, npoint, npointmax, ivicinity,
	         rownottop, rownotbottom, colnotleft, colnotright;
	uint8_t *pointloc;
	int16_t  rowv, colv;
	uint16_t row, col, rowlimit, collimit, ipt, iptnear, ptxcount;
	uint32_t idx, idxv;
	uint64_t linesec;
	struct cpt_pt *ppt,
	              *ppairpt;
	struct cpt_px *ppx;
	struct cpt_point *ppoint,
	                 *ppairpoint;
	
	ptxcount = 0;
	rowlimit = pospst->nrow-1;
	collimit = pospst->ncol-1;
	*pairpx  = NULL;
	*pairpt  = NULL;
	
	npointmax = 0;
	for (ipt = 0; ipt < ptcount; ++ipt) {
		npoint = (allpt+ipt)->nt;
		if (npoint > npointmax)
			npointmax = npoint;
	}
	pointloc = malloc(sizeof(int8_t[npointmax]));
	
	/*  Init with 0/flase  */
	for (ipt = 0; ipt < ptcount; ++ipt) {
		ptgeodiff[ipt] = 0;
	}
	
	for (row = 0; row < pospst->nrow; ++row) {
		/*  Vars can be put outside of column loop  */
		linesec = pospst->secswhenscan + (rowlimit-row)*pospst->secsperline;
		if (linesec > pospst->secswhenend)
			linesec = pospst->secswhenend;
		
		idx = (uint32_t) row * pospst->ncol;
		rownottop = (row != 0);
		rownotbottom = (row != rowlimit);
	for (col = 0; col < pospst->ncol; ++col, ++idx) {
		/*  Nearest Pt of current pixel  */
		diffmin = UINT64_MAX;
		for (ipt = 0; ipt < ptcount; ++ipt) {
			ppt  = allpt+ipt;
			diff = fabsf(ppt->lon - pospst->lon[idx]) + fabsf(ppt->lat - pospst->lat[idx]);
			if (diff < diffmin) {
				diffmin = diff;
				iptnear = ipt;
			}
		}
		
		/*
		 *  Although the chosen Pt is the nearest site among all Pt,
		 *  it may also locates at the outside of image.
		 *  lonres and latres are estimated from two nearby pixels
		 *  as lon/lat difference minimium threshold.
		 */
		lonres = 0;
		if (colnotleft  = (col != 0))
			lonres = fabsf(pospst->lon[idx-1] - pospst->lon[idx]);
		if (colnotright = (col != collimit)) {
			if (lonres) {
				lonres += fabsf(pospst->lon[idx+1] - pospst->lon[idx]);
				lonres /= 2;
			} else {
				lonres = fabsf(pospst->lon[idx+1] - pospst->lon[idx]);
			}
		}
		
		latres = 0;
		if (rownottop)
			latres = fabsf(pospst->lat[idx-pospst->ncol] - pospst->lat[idx]);
		if (rownotbottom) {
			if (latres) {
				latres += fabsf(pospst->lat[idx+pospst->ncol] - pospst->lat[idx]);
				latres /= 2;
			} else {
				latres = fabsf(pospst->lat[idx+pospst->ncol] - pospst->lat[idx]);
			}
		}
		
		ppt = allpt+iptnear;
		londiff = fabsf(ppt->lon - pospst->lon[idx]);
		latdiff = fabsf(ppt->lat - pospst->lat[idx]);
		if ((londiff > lonres) || (latdiff > latres))
			goto next_pixel;
		
		/*
		 *  At most 9 pixels may match one site,
		 *  use ptgeodiff to determine which is the closest.
		 */
		diff = londiff+latdiff;
		appendpx = 1;
		if (ptgeodiff[iptnear]) {
			if (diff > ptgeodiff[iptnear]) {
				goto next_pixel;
			} else {
				appendpx = 0;
				ptgeodiff[iptnear] = diff;
			}
		} else {
			ptgeodiff[iptnear] = diff;
		}
		
		/*
		 *  The pixel matches this closest site, spatially.
		 *  We would also check them temporally, based on the datetime of satellite pixel.
		 */
		npoint  = 0;
		for (ipoint = 0; ipoint < ppt->nt; ++ipoint) {
			ppoint = ppt->points + ipoint;
			if (((ppoint->seconds > linesec) ?
			(ppoint->seconds-linesec) : (linesec-ppoint->seconds)) < WR_CPT_SECDIFFMAX) {
				pointloc[npoint++] = ipoint;
			}
		}
		
		if (!npoint)
			goto next_pixel;
		
		if (appendpx) {
			/*
			 *  Finally, the Points and the Pixels get paired.
			 *  Points is multiple as its several obs from certian location.
			 *  Pixels is multiple as its several location from one time.
			 */
			*pairpt = realloc(*pairpt, sizeof(struct cpt_pt[++ptxcount]));
			ppairpt = *pairpt + ptxcount-1;
			ppairpt->nt  = npoint;
			ppairpt->alt = ppt->alt;
			ppairpt->lon = ppt->lon;
			ppairpt->lat = ppt->lat;
			ppairpt->points = malloc(sizeof(struct cpt_point[ppairpt->nt]));
			for (ipoint = 0; ipoint < npoint; ++ipoint) {
				ppoint = ppt->points + pointloc[ipoint];
				ppairpoint = ppairpt->points + ipoint;
				ppairpoint->seconds = ppoint->seconds;
				ppairpoint->params = malloc(_cpt_parsz);
				memcpy(ppairpoint->params, ppoint->params, _cpt_parsz);
			}
			
			/*  Px init or reload  */
			*pairpx = realloc(*pairpx, sizeof(struct cpt_px[ptxcount]));
			ppx = *pairpx + ptxcount-1;
			ppx->centrepixel = malloc(CPT_PIXELSIZE);
			ppx->nvicinity   = 0;
			ppx->vicinity    = NULL;
		}
		
		/*  Centre pixel  */
		loadpxfromst(ppx->centrepixel, pospst, row, col, appendpx);
		
		/*  Vicinity count  */
		switch (rownottop+rownotbottom+colnotleft+colnotright) {
		case 4: ivicinity = 8; break;
		case 3: ivicinity = 5; break;
		case 2: ivicinity = 3; break;
		default: ivicinity = 0;
		}
		
		/*  Vicinity memory manage  */
		if (ivicinity) {
			ppx->vicinity = realloc(ppx->vicinity,
			                        sizeof(struct cpt_pixel[ppx->nvicinity = ivicinity]));
		}

		
		/*  Vicinity assignment  */
		ivicinity = 0;
		for (rowv = -1; rowv < 2; ++rowv) {
			if (((-1 == rowv) && !rownottop) || ((1 == rowv) && !rownotbottom))
				continue;
			for (colv = -1; colv < 2; ++colv) {
				if (((-1 == colv) && !colnotleft) || ((1 == colv) && !colnotright))
					continue;
				/*  Centre pixel already load before  */
				if ((0 == rowv) && (0 == colv))
					continue;
				loadpxfromst(ppx->vicinity+ivicinity++, pospst,
				             row+rowv, col+colv, appendpx);
			}
		}
		
		next_pixel:
		continue;
	}
	}
	
	ppx = NULL;
	ppt = ppairpt = NULL;
	ppoint = ppairpoint = NULL;
	CPT_FREE(pointloc);
	
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
	uint8_t  ipoint, ivicinity;
	uint16_t iptx;
	int fd;
	
	/*  File already exist ?  */
	fd = open(fname, O_PATH);
	if (fd > 0) {
		close(fd);
#ifdef CPT_DEBUG
		return EEXIST;
#else
		fd = open(fname, O_WRONLY|O_CREAT, S_IRUSR|S_IWUSR);
#endif
	} else {
		fd = open(fname, O_WRONLY|O_CREAT, S_IRUSR|S_IWUSR);
	}
	
	/*  Header  */
	safewrite(fd, st->hdr->magic_number, CPT_MAGICLEN);
	safewrite(fd, &st->hdr->ver, _cpt_1byte);
	safewrite(fd, &st->hdr->nptx, _cpt_2byte);
	safewrite(fd, &st->hdr->nparam, _cpt_1byte);
	
	/*  Data/Ptx  */
	for (iptx = 0; iptx < st->hdr->nptx; ++iptx) {
		
		/*  Pt  */
		safewrite(fd, &(st->data->pt+iptx)->lon, _cpt_4byte);
		safewrite(fd, &(st->data->pt+iptx)->lat, _cpt_4byte);
		safewrite(fd, &(st->data->pt+iptx)->alt, _cpt_2byte);
		safewrite(fd, &(st->data->pt+iptx)->nt , _cpt_1byte);
		for (ipoint = 0; ipoint < (st->data->pt+iptx)->nt; ++ ipoint) {
			safewrite(fd, &((st->data->pt+iptx)->points+ipoint)->seconds, _cpt_8byte);
			safewrite(fd, ((st->data->pt+iptx)->points+ipoint)->params, _cpt_parsz);
		}
		
		/*  Px  */
		safewrite(fd, &(st->data->px+iptx)->seconds, _cpt_8byte);
		writepixel(fd, (st->data->px+iptx)->centrepixel);
		safewrite(fd, &(st->data->px+iptx)->nvicinity, _cpt_1byte);
		for (ivicinity = 0; ivicinity < (st->data->px+iptx)->nvicinity; ++ivicinity) {
			writepixel(fd, (st->data->px+iptx)->vicinity+ivicinity);
		}
		
	}
	
	/*  Ending  */
	safewrite(fd, st->ending, CPT_ENDINGLEN);
	close(fd);

	return 0;
}

/*  Definition of main function  */
int wrcpt(const char *ptname, const char *pxname, const char *cptname)
{
	int ret;
	uint16_t  ptcount, ptxcount;
	struct cpt_pt *allpt;
	struct cpt_pt *pairpt;
	struct cpt_px *pairpx;
	struct wr_cpt_posp *pospst;
	struct cpt_header hdr;
	struct cpt_ptx ptx;
	struct cpt_ff cptout;
	
	/*  Px prepare  */
	pospinfoinit(pxname, &pospst);
	pospopenall(pxname, pospst);
	
	/*  All Pt load  */
	if (ret = querypt(ptname, &allpt, &ptcount)) {
		cpt_freeptall(&allpt, ptcount);
		return ret;
	}
	
	/*  Exit when there's no record in Pt  */
	if (!ptcount) {
		cpt_freeptall(&allpt, ptcount);
		pospcleanst(&pospst);
		return WR_CPT_ENORCD;
	}
	
	/*  Get paired Px and Pt  */
	ptxcount = posppair(allpt, ptcount, pospst, &pairpx, &pairpt);
#ifdef CPT_DEBUG
	for (uint16_t i = 0; i < ptxcount; ++i) {
		printf("No.%03d: lon %9.4f lat %8.4f with %2d points\n",
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
	writecpttofile(cptname, &cptout);
	
	/*  Cleanup  */
	ret = cpt_freeptall(&allpt, ptcount);
	ret = cpt_freepxall(&pairpx, ptxcount);
	ret = cpt_freeptall(&pairpt, ptxcount);
	ret = pospcleanst(&pospst);
	
	ptx.pt = NULL;
	ptx.px = NULL;
	cptout.hdr  = NULL;
	cptout.data = NULL;
	
	return 0;
}

