/*
 *file: write_examples/GF5(B)/write.c
 *descreption:
 *  write cpt format file
 *arguements:
 *  [1]: downloaded in-situ aeronet file
 *  [2]: DPC/POSP data file
 *  [3]: output cpt
 *init date: May/10/2022
 *last modify: May/11/2022
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


#define WR_CPT_INVSITEBYTES 180
#define WR_CPT_DATELEN      11
#define WR_CPT_TIMELEN      9
#define WR_CPT_DTLEN        (WR_CPT_DATELEN+WR_CPT_TIMELEN)


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

/*  Definition of main function  */
int wrcpt(const char *ptname, const char *psname, const char *cptname)
{
	int            ret;
	uint16_t       ptcount;
	struct cpt_pt *allpt;
	struct cpt_ps *allps;
	
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
	
	/*  Cleanup  */
	ret = cpt_freeptall(&allpt, ptcount);
	
	return 0;
}

