/*
 *file: write_examples/GF5(B)/write.c
 *descreption:
 *  write cpt format file
 *arguements:
 *  [1]: downloaded in-situ aeronet file
 *  [2]: DPC/POSP data file
 *  [3]: output cpt
 *init date: May/10/2022
 *last modify: May/12/2022
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

#include "hdf5.h"
#include "../../read/readcpt.h"


/*  Site info settings  */
#define WR_CPT_INVSITEBYTES 180
#define WR_CPT_DATELEN      11
#define WR_CPT_TIMELEN      9
#define WR_CPT_DTLEN        (WR_CPT_DATELEN+WR_CPT_TIMELEN)


/*  Satellite info settings  */
#define WR_CPT_POSPNBANDS  9
#define WR_CPT_POSPDSNAME  "Data_Fields"
#define WR_CPT_POSPGLNAME  "Geolocation_Fields"
#define WR_CPT_POSPINAME   "I"
#define WR_CPT_POSPQNAME   "Q"
#define WR_CPT_POSPUNAME   "U"
#define WR_CPT_POSPLATNAME "Latitude"
#define WR_CPT_POSPLONNAME "Longitude"
#define WR_CPT_POSPSLNAME  "Sea_Land_Flags"
#define WR_CPT_POSPALTNAME "Surface_Altitude"
#define WR_CPT_POSPSZNAME  "Sol_Zen_Ang"
#define WR_CPT_POSPVZNAME  "View_Zen_Ang"
#define WR_CPT_POSPSANAME  "Sol_Azim_Ang"
#define WR_CPT_POSPVANAME  "View_Azim_Ang"

#define WR_CPT_XMLSUFFIX   "xml"
#define WR_CPT_XMLSUFLEN   strlen(WR_CPT_XMLSUFFIX)
#define WR_CPT_XMLENDTAG   "</ProductMetaData>"
#define WR_CPT_XMLSTTAG    "StartTime"
#define WR_CPT_XMLETTAG    "EndTime"

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
	hid_t m2id;  /*  2-d mem id     */
	hid_t m3id;  /*  3-d mem id     */  
	
	hid_t gdid;  /*  dateset group id  */
	hid_t ggid;  /*  geo-loc group id  */
	hid_t fid;   /*  file id           */
	
	hsize_t l2id[2];  /*  2-d hyper len  */
	hsize_t l3id[3];  /*  3-d hyper len  */
	
	uint16_t ncol;  /*  count of col  */
	uint16_t nrow;  /*  count of row  */
	uint16_t secsperline;   /*  timelapse per line  */
	uint64_t secswhenscan;  /*  timestamp of begin  */
	double *lat;  /*  full lat  */
	double *lon;  /*  full lon  */
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


/*  Declaration of main function  */
int wrcpt(const char *ptname, const char *psname, const char *cptname);

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
 *  Free one or more cpt_pt st pointer
 */
static int cpt_freeptall(struct cpt_pt **p, uint16_t n)
{
	if (*p) {
		while (n-- > 0)
			CPT_FREE((*p+n)->params);
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
 *  Free one or more cpt_ps st pointer
 */
static int cpt_freepsall(struct cpt_ps **p, uint16_t n)
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
	uint32_t flen, llen, maxlen;
	struct cpt_pt *ppt;
	struct tm     *stm;
	
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
	
	stm    = malloc(sizeof(struct tm));
	line   = malloc(1);
	maxlen = 1;
	*allpt   = malloc(1);
	*ptcount = 0;
	
	/*  Handle each line  */
	do {
		while (*++pbuf != '\n') ;
		
		llen = pbuf-pprev;
		if (llen > maxlen) {
			maxlen = llen;
			line = realloc(line, sizeof(char[llen]));
		}
		memcpy(line, pprev, llen);
		
		/*  End of records  */
		if ('<' == *line) break;
		
		*allpt = realloc(*allpt, CPT_PTSIZE*(++*ptcount));
		ppt = *allpt+*ptcount-1;
		ppt->params = malloc(sizeof(double[8]));
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
		             
		             rcddate, rcdtime,                           /*  A.receive  */
		             ppt->params, ppt->params+1, ppt->params+2,  /*  B.receive  */
		             ppt->params+3,                              /*  C.receive  */
		             ppt->params+4,                              /*  D.receive  */
		             ppt->params+5,                              /*  E.receive  */
		             ppt->params+6, ppt->params+7,               /*  F.receive  */
		             &ppt->lat, &ppt->lon, &ppt->alt);           /*  G.receive  */
		
		snprintf(rcddt, WR_CPT_DTLEN, "%s,%s", rcddate, rcdtime);
		strptime(rcddt, "%d:%m:%Y,%H:%M:%S", stm);
		ppt->seconds = mktime(stm);
		printf("Time: %s Sec: %lld Lat: %f Lon: %f Alt:%hd AOD500: %lf \n",
		       rcddt, ppt->seconds, ppt->lat,ppt->lon,ppt->alt,ppt->params[4]);
		
		pprev = ++pbuf;
	} while (*pbuf != '\0');
	
	/*  Set pointers to NULL  */
	pbuf = pprev = NULL;
	ppt  = NULL;
	
	/*  Free allocated buffer  */
	line = realloc(line, sizeof(char[maxlen]));
	cpt_freethemall(3, &buffer, &line, &stm);
	
	return 0;
}

/*  Init st from POSP h5  */
static int pospopenall(const char *fname, struct wr_cpt_posp *st)
{
	hid_t   latid, lonid;
	hid_t   space;
	size_t  size;
	
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
	
	/*  Space for hyper-reading  */
	st->h2id = H5Screate_simple(2, (hsize_t[2]) {st->nrow, st->ncol}, NULL);
	st->h3id = H5Screate_simple(3, (hsize_t[3]) {WR_CPT_POSPNBANDS, st->nrow, st->ncol}, NULL);
	st->m2id = H5Screate_simple(1, (hsize_t[1]) {1}, NULL);
	st->m3id = H5Screate_simple(1, (hsize_t[1]) {WR_CPT_POSPNBANDS}, NULL);
	
	st->l3id[0] = WR_CPT_POSPNBANDS;
	st->l3id[1] = st->l3id[2] = 1;
	st->l2id[0] = st->l2id[1] = 1;
	
	/*  Buffer all lat/lon  */
	latid = H5Dopen(st->ggid, WR_CPT_POSPLATNAME, H5P_DEFAULT);
	lonid = H5Dopen(st->ggid, WR_CPT_POSPLONNAME, H5P_DEFAULT);
	st->lat = malloc(size = sizeof(double[st->nrow][st->ncol]));
	st->lon = malloc(size);
	H5Dread(latid, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, st->lat);
	H5Dread(lonid, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, st->lon);
	H5Dclose(latid);
	H5Dclose(lonid);
	
	return 0;
}

/*  Set hyper space before read partial data  */
static int pospsethyper(struct wr_cpt_posp *st, uint16_t ir, uint16_t ic)
{
	return
	H5Sselect_hyperslab(st->h3id, H5S_SELECT_SET, (hsize_t[3]) {0, ir, ic}, NULL, st->l3id, NULL) < 0
	||
	H5Sselect_hyperslab(st->h2id, H5S_SELECT_SET, (hsize_t[2]) {ir, ic}, NULL, st->l2id, NULL) < 0;
}

static int pospinfoinit(const char *fname, struct wr_cpt_posp *pospst)
{
	char  *xmlfname;
	size_t len;
	int      fd;
	char    *buffer, *pbuf, *pprev, *line, *pline;
	char     dtbeg[WR_CPT_DTLEN], dtend[WR_CPT_DTLEN];
	uint32_t flen, llen, maxlen;
	struct tm *stm;
	
	xmlfname = malloc(len = strlen(fname)-1);
	len -= WR_CPT_XMLSUFLEN;
	
	memcpy(xmlfname, fname, len);
	memcpy(xmlfname+len, WR_CPT_XMLSUFFIX, WR_CPT_XMLSUFLEN);
	
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
	
	stm    = malloc(sizeof(struct tm));
	line   = malloc(1);
	maxlen = 1;
	
	/*  Handle each line  */
	do {
		while (*++pbuf != '\n') ;
		
		llen = pbuf-pprev;
		if (llen > maxlen) {
			maxlen = llen;
			line = realloc(line, sizeof(char[llen]));
		}
		memcpy(line, pprev, llen);
		
		/*  End of records  */
		if (pline = strstr(line, WR_CPT_XMLENDTAG)) break;
		
		snprintf(rcddt, WR_CPT_DTLEN, "%s,%s", rcddate, rcdtime);
		strptime(rcddt, "%d:%m:%Y,%H:%M:%S", stm);
		ppt->seconds = mktime(stm);
		printf("Time: %s Sec: %lld Lat: %f Lon: %f Alt:%hd AOD500: %lf \n",
		       rcddt, ppt->seconds, ppt->lat,ppt->lon,ppt->alt,ppt->params[4]);
		
		pprev = ++pbuf;
	} while (*pbuf != '\0');
	
	/*  Set pointers to NULL  */
	pbuf = pprev = pline = NULL;
	
	/*  Free allocated buffer  */
	line = realloc(line, sizeof(char[maxlen]));
	cpt_freethemall(4, &buffer, &line, &stm, &xmlfname);
	
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
	H5Sclose((*st)->m2id);
	H5Sclose((*st)->m3id);
	
	H5Gclose((*st)->gdid);
	H5Gclose((*st)->ggid);
	
	H5Fclose((*st)->fid);
	
	cpt_freethemall(3, &(*st)->lat, &(*st)->lon, st);
	
	return 0;
}

/*  Definition of main function  */
int wrcpt(const char *ptname, const char *psname, const char *cptname)
{
	int            ret;
	uint16_t       ptcount;
	struct cpt_pt *allpt;
	struct cpt_ps *allps;
	struct wr_cpt_posp *pospst;
	
	if (ret = querypt(ptname, &allpt, &ptcount)) {
		cpt_freeptall(&allpt, ptcount);
		return ret;
	}
	
	if (ptcount) {
		allps = malloc(CPT_PSSIZE*ptcount);
	} else {
		cpt_freeptall(&allpt, ptcount);
		return WR_CPT_ENORCD;
	}
	
	pospst = malloc(sizeof(struct wr_cpt_posp));
	pospinfoinit(psname, pospst);
	
	/*  Cleanup  */
	ret = cpt_freeptall(&allpt, ptcount);
	
	return 0;
}

