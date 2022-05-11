/*
 *file: read/readcpt.h
 *init date: May/10/2022
 *last modify: May/11/2022
 *
 */


#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdarg.h>
#include <errno.h>


/*  Const numbers  */
#define CPT_MAGIC  (uint8_t[13]) {0x02, 'L', 'e', 'r', 'S', 'A', 'T', '@', 'c', 'p', 't', '\n', 0x03}
#define CPT_ENDING (uint8_t[16]) {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}


/*  Structures in cpt hierarchy  */
struct cpt_header {
	uint8_t  ver;
	uint8_t  nparam;
	uint16_t nptps;
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
	int16_t alt;
	float   lat;
	float   lon;
	struct  cpt_channel *channels;
};

struct cpt_ps {
	uint8_t nvicinity;
	int8_t  minutesdiff;
	struct  cpt_pixel *centrepixel;
	struct  cpt_pixel *vicinity;
};
#define CPT_PSSIZE (sizeof(struct cpt_ps))

struct cpt_pt {
	int16_t  alt;
	float    lat;
	float    lon;
	uint64_t seconds;
	double  *params;
};
#define CPT_PTSIZE (sizeof(struct cpt_pt))

struct cpt_ptps {
	struct cpt_pt *pt;
	struct cpt_ps *ps;
};

struct cpt_ff {
	struct cpt_header *hdr;
	struct cpt_ptps   *data;
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





