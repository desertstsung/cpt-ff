/*
 *file: read/readcpt.c
 *descreption:
 *  read cpt format file
 *init date: May/10/2022
 *last modify: Jul/27/2022
 *
 */

#include "readcpt.h"

const static size_t _cpt_1byte = sizeof(int8_t );
const static size_t _cpt_2byte = sizeof(int16_t);
const static size_t _cpt_4byte = sizeof(int32_t);
const static size_t _cpt_8byte = sizeof(int64_t);

#ifdef CPT_DEBUG
int main(int argc, char *argv[])
{
	if (argc != 2) {
		CPT_ERRECHOWITHTIME("Usage: %s cptfile", argv[0]);
		return 1;
	}
	
	uint8_t  nparam;
	uint32_t nptx;
	struct cpt_ptx ptx;
	
	cpt_readall(argv[1], &ptx, &nptx, &nparam);
	for (uint16_t i = 0; i < nptx; ++i) {
		printf("No.%03d: lon %9.4f lat %8.4f with %2d points (%s)\n",
		       i+1, (ptx.px+i)->centrepixel->lon,
		       (ptx.px+i)->centrepixel->lat, (ptx.pt+i)->nt,
		       (ptx.pt+i)->name);
	}
	
	cpt_freeptall(&ptx.pt, nptx);
	cpt_freepxall(&ptx.px, nptx);
	
	return 0;
}
#endif

int cpt_readall(const char *fname, struct cpt_ptx *ptx, uint32_t *nptx, uint8_t *nparam)
{
	int fd;
	size_t  sparams;
	uint8_t mgc[CPT_MAGICLEN], ending[CPT_ENDINGLEN], ver;
	uint8_t  ipoint, ivicinity;
	uint16_t iptx;
	struct cpt_pt *ppt;
	struct cpt_px *ppx;
	struct cpt_point *ppoint;
	
	if ((fd = open(fname, O_RDONLY)) < 0) {
		CPT_ERROPEN(fname);
		return 1;
	}
	
	/*  Header check  */
	read(fd, mgc, CPT_MAGICLEN);
	if (memcmp(mgc, CPT_MAGIC, CPT_MAGICLEN)) {
		close(fd);
		CPT_ERRECHOWITHTIME("%s is NOT a cpt file!", fname);
		return 2;
	}
	
	/*  Version check  */
	read(fd, &ver, _cpt_1byte);
	if (CPT_VERSION != ver) {
		close(fd);
		CPT_ERRECHOWITHTIME("%s is a cpt file in version %d.%d!\n"
		                    "while current lib is %d.%d",
		                    fname, ver>>4, ver|0b00001111, CPT_VER_MAJOR, CPT_VER_MINOR);
	}
	
	/*  Meta info  */
	read(fd, nptx, _cpt_4byte);
	read(fd, nparam, _cpt_1byte);
	sparams = sizeof(double[*nparam]);
	
	ptx->pt = malloc(sizeof(struct cpt_pt[*nptx]));
	ptx->px = malloc(sizeof(struct cpt_px[*nptx]));
	
	/*  Data  */
	char namec;
	uint8_t namelen; // Assume all name is shorter than 255 characters.
	for (iptx = 0; iptx < *nptx; ++iptx) {
		
		/*  Pt  */
		ppt = ptx->pt+iptx;
		
		ppt->name = NULL;
		namelen = 0;
		while (read(fd, &namec, _cpt_1byte) && (namec != '\0') && (++namelen)) {
			ppt->name = realloc(ppt->name, namelen+1);
			ppt->name[namelen-1] = namec;
			ppt->name[namelen] = '\0';
		}
		
		read(fd, &ppt->lon, _cpt_4byte);
		read(fd, &ppt->lat, _cpt_4byte);
		read(fd, &ppt->alt, _cpt_2byte);
		read(fd, &ppt->nt , _cpt_1byte);
		ppt->points = malloc(sizeof(struct cpt_point[ppt->nt]));
		for (ipoint = 0; ipoint < ppt->nt; ++ipoint) {
			ppoint = ppt->points+ipoint;
			read(fd, &ppoint->seconds, _cpt_8byte);
			ppoint->params = malloc(sparams);
			read(fd, ppoint->params, sparams);
		}
		
		/*  Px  */
		ppx = ptx->px+iptx;
		read(fd, &ppx->seconds, _cpt_8byte);
		ppx->centrepixel = malloc(sizeof(struct cpt_pixel));
		readpixel(fd, ppx->centrepixel);
		read(fd, &ppx->nvicinity, _cpt_1byte);
		ppx->vicinity = malloc(sizeof(struct cpt_pixel[ppx->nvicinity]));
		for (ivicinity = 0; ivicinity < ppx->nvicinity; ++ivicinity) {
			readpixel(fd, ppx->vicinity+ivicinity);
		}
		
	}
	
	/*  Ending  */
	read(fd, ending, CPT_ENDINGLEN);
	close(fd);
	if (memcmp(ending, CPT_ENDING, CPT_ENDINGLEN)) {
		CPT_ERRECHOWITHTIME("%s has NO ending, the results may be incorrect", fname);
		return 2;
	}
	
	ppt = NULL;
	ppx = NULL;
	ppoint = NULL;
	
	return 0;
}

int readpixel(int filedes, struct cpt_pixel *pixel)
{
	/*  Geolocation  */
	read(filedes, &pixel->lon, _cpt_4byte);
	read(filedes, &pixel->lat, _cpt_4byte);
	read(filedes, &pixel->alt, _cpt_2byte);
	read(filedes, &pixel->mask, _cpt_1byte);
	
	/*  Dimensions  */
	read(filedes, &pixel->nchannel, _cpt_1byte);
	read(filedes, &pixel->nlayer, _cpt_1byte);
	
	if (pixel->nchannel) {
		size_t angsize = sizeof(double[pixel->nlayer][4]);
		size_t _obssize = sizeof(double[pixel->nlayer]);
		
		/*  Channel  */
		size_t obssize;
		struct cpt_channel *pchannel;
		pixel->channels = malloc(sizeof(struct cpt_channel[pixel->nchannel]));
		for (uint8_t ichannel = 0; ichannel < pixel->nchannel; ++ichannel) {
			pchannel = pixel->channels+ichannel;
			read(filedes, &pchannel->centrewv, _cpt_2byte);
			
			pchannel->obs = malloc(obssize = _obssize*((pchannel->centrewv < 0) ? 3 : 1));
			read(filedes, pchannel->obs, obssize);
			pchannel->ang = malloc(angsize);
			read(filedes, pchannel->ang, angsize);
		}
		pchannel = NULL;
	} else {
		pixel->nlayer = 0;
		pixel->channels = NULL;
	}
	
	read(filedes, &pixel->nextra, _cpt_1byte);
	if (pixel->nextra) {
		pixel->extra = malloc(sizeof(double[pixel->nextra]));
		for (uint8_t iextra = 0; iextra < pixel->nextra; ++iextra) {
			read(filedes, pixel->extra+iextra, _cpt_8byte);
		}
	} else {
		pixel->extra = NULL;
	}

	return 0;
}

/*
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

/*
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

/*
 *  Free one or more cpt_pt st pointer
 */
int cpt_freeptall(struct cpt_pt **p, uint32_t n)
{
	if (*p) {
		while (n-- > 0) {
			cpt_freepointall(&(*p+n)->points, (*p+n)->nt);
			CPT_FREE((*p+n)->name);
		}
		CPT_FREE(*p);
	}
	
	return 0;
}

/*
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

/*
 *  Free one or more cpt_pixel st pointer
 */
int cpt_freepixelall(struct cpt_pixel **p, uint16_t n)
{
	if (*p) {
		while (n-- > 0) {
			cpt_freechannelall(&(*p+n)->channels, (*p+n)->nchannel);
			CPT_FREE((*p+n)->extra);
		}
		CPT_FREE(*p);
	}
	
	return 0;
}

/*
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

