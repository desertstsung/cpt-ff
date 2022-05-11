/*
 *file: write_examples/GF5(B)/write.c
 *descreption:
 *  write cpt format file
 *arguements:
 *  [1]: downloaded in-situ aeronet file
 *  [2]: DPC/POSP data file
 *  [3]: output cpt
 *init date: May/10/2022
 *last modify: May/10/2022
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


/*  Error numbers  */
enum WR_CPT_ERR {
WR_CPT_EINVARG = 1,
WR_CPT_EOPEN,
WR_CPT_EINVSITE,
WR_CPT_EMEM,
WR_CPT_E
};

/*  Declaration of main function  */
int wrcpt(const char *ptname, const char *psname, const char *cptname);

/*  Program entrance  */
int main(int argc, char *argv[]) {
	if (4 == argc) {
		return wrcpt(argv[1], argv[2], argv[3]);
	} else {
		CPT_ERRECHOWITHTIME("Usage: %s sitefile, satefile, output", argv[0]);
		return WR_CPT_EINVARG;
	}
}

static int querypt(const char *fname, struct cpt_pt **allpt, uint16_t *ptcount)
{
	int      fd;
	char    *buffer, *pbuf, *ppos, *pprev, *line;
	char     rcddate[11], rcdtime[9], rcddt[20];
	double  *pparam;
	uint32_t flen, llen, maxlen;
	struct cpt_pt *ppt;
	struct tm     *stm;
	
	/*  Try openning text file  */
	if ((fd = open(fname, O_RDONLY)) < 0) {
		CPT_ERROPEN(fname);
		return WR_CPT_EOPEN;
	}
	
	/*  Bad content  */
	if ((flen = lseek(fd, 0, SEEK_END)) < 180) {
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
	pprev = ppos = ++pbuf;
	
	stm    = malloc(sizeof(struct tm));
	line   = malloc(1);
	maxlen = 1;
	*allpt   = malloc(1);
	*ptcount = 0;
	
	/*  Handle each line  */
	do {
		while (*++pbuf != '\n') ;
		
		llen = pbuf-pprev;
		if (llen > maxlen) maxlen = llen;
		line = realloc(line, sizeof(char[llen]));
		memcpy(line, pprev, llen);
		
		/*  End of records  */
		if ('<' == *line) break;
		
		*allpt = realloc(*allpt, CPT_PTSIZE*(++*ptcount));
		ppt = *allpt+*ptcount-1;
		ppt->params = malloc(sizeof(double[7]));
		pparam = (double *) ppt->params;
		sscanf(line, "%*[^','],"           /*  Skip site name  */
		             "%[^','],%[^','],"    /*  Read date time  */
		             "%*[^','],%*[^','],"  /*  Skip DoY  */
		             "%lf,%lf,%lf,", rcddate, rcdtime, pparam, pparam+1, pparam+2);
		snprintf(rcddt, 20, "%s,%s", rcddate, rcdtime);
		strptime(rcddt, "%d:%m:%Y,%H:%M:%S", stm);
		ppt->seconds = mktime(stm);
		printf("%s %lld\n", rcddt, ppt->seconds);
		
		pprev = ++pbuf;
	} while (*pbuf != '\0');
	
	/*  Free line buffer  */
	line = realloc(line, sizeof(char[maxlen]));
	free(line);
	
	return 0;
}

/*  Definition of main function  */
int wrcpt(const char *ptname, const char *psname, const char *cptname)
{
	struct cpt_pt *allpt;
	int            ret;
	uint16_t       ptcount;
	
	if (ret = querypt(ptname, &allpt, &ptcount))
		return ret;
	
	return 0;
}

