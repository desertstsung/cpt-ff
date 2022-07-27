/*
 *file: read/readcpt.h
 *init date: May/10/2022
 *last modify: Jul/27/2022
 *
 */


#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdarg.h>
#include <errno.h>
#include <string.h>


/*  Const numbers  */
#define CPT_MAGICLEN  13
#define CPT_ENDINGLEN 16
#define CPT_MAGIC     (uint8_t[CPT_MAGICLEN]) \
                      {0x02, 'L', 'e', 'r', 'S', 'A', 'T', '@', 'c', 'p', 't', '\n', 0x03}
#define CPT_ENDING    (uint8_t[CPT_ENDINGLEN]) \
                      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}


/*  Version  */
#define CPT_VER_MAJOR (uint8_t) 0
#define CPT_VER_MINOR (uint8_t) 1
#define CPT_VERSION   ((CPT_VER_MAJOR<<4) | CPT_VER_MINOR)


/*  Structures in cpt hierarchy  */
struct cpt_header {
	uint8_t  ver;
	uint8_t  nparam;
	uint32_t nptx;
	uint8_t *magic_number;
};

struct cpt_channel {
	int16_t centrewv;
	double  *obs;
	double  *ang;
};

struct cpt_pixel {
	uint8_t mask;
	uint8_t nchannel;
	uint8_t nlayer;
	uint8_t nextra;
	int16_t alt;
	float   lat;
	float   lon;
	double *extra;
	struct  cpt_channel *channels;
};
#define CPT_PIXELSIZE (sizeof(struct cpt_pixel))

struct cpt_px {
	uint8_t  nvicinity;
	uint64_t seconds;
	struct  cpt_pixel *centrepixel;
	struct  cpt_pixel *vicinity;
};

struct cpt_point {
	uint64_t seconds;
	double  *params;
};
#define CPT_POINTSIZE  (sizeof(struct cpt_point))

struct cpt_pt {
	uint8_t nt;
	int16_t alt;
	float   lon;
	float   lat;
	struct cpt_point *points;
};

struct cpt_ptx {
	struct cpt_pt *pt;
	struct cpt_px *px;
};

struct cpt_ff {
	struct cpt_header *hdr;
	struct cpt_ptx    *data;
	void  *ending;
};


/*  Useful fn  */
#define CPT_FREE(ptr) \
	do { \
		if (ptr) { \
			free(ptr); \
			ptr = NULL; \
		} \
	} while (0)

time_t curtime;
#define __CPT_ECHOWITHTIME(stream, ...) \
	do { \
		time(&curtime); \
		fprintf(stream, "[%15.15s] ", ctime(&curtime)+4); \
		fprintf(stream, __VA_ARGS__); \
		fprintf(stream, "\n"); \
	} while (0)
#define CPT_ECHOWITHTIME(...) __CPT_ECHOWITHTIME(stdout, __VA_ARGS__)


/*  Error prompt  */
#define CPT_ERRLOC fprintf(stderr, "File: %s, Fn: %s, Ln: %d\n", __FILE__, __FUNCTION__, __LINE__)

#define CPT_ERRECHOWITHTIME(...) __CPT_ECHOWITHTIME(stderr, __VA_ARGS__)

#define CPT_ERROPEN(fname) \
	do { \
		CPT_ERRECHOWITHTIME("ERROR %d %s: %s\n", errno, strerror(errno), fname); \
		CPT_ERRLOC; \
	} while (0)

#define CPT_ERRFIO(fd) \
	do { \
		close(fd); \
		CPT_ERRECHOWITHTIME("ERROR %d %s\n", errno, strerror(errno)); \
		CPT_ERRLOC; \
	} while (0)

#define CPT_ERRMEM(ptr) \
	do { \
		CPT_FREE(ptr); \
		CPT_ERRECHOWITHTIME("MEM PANIC %d %s\n", errno, strerror(errno)); \
		CPT_ERRLOC; \
	} while (0)


/*  fn  */
int cpt_readall(const char *fname, struct cpt_ptx *ptx, uint32_t *nptx, uint8_t *nparam);
int readpixel(int filedes, struct cpt_pixel *pixel);
int cpt_freethemall(uint8_t n, ...);
int cpt_freepointall(struct cpt_point **p, uint16_t n);
int cpt_freeptall(struct cpt_pt **p, uint32_t n);
int cpt_freechannelall(struct cpt_channel **p, uint16_t n);
int cpt_freepixelall(struct cpt_pixel **p, uint16_t n);
int cpt_freepxall(struct cpt_px **p, uint32_t n);

