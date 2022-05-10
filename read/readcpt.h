/*
 *file: read/readcpt.h
 *init date: May/10/2022
 *last modify: May/10/2022
 *
 */


#include <stdint.h>


#define CPT_MAGIC  (uint8_t[13]) {0x02, 'L', 'e', 'r', 'S', 'A', 'T', '@', 'c', 'p', 't', '\n', 0x03}
#define CPT_ENDING (uint8_t[16]) {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}


struct cpt_header {
	uint8_t  ver;
	uint8_t  nparam;
	uint16_t nptps;
	uint8_t *magic_number;
};

struct cpt_channel {
	int16_t centrewv;
	void   *obs;
	void   *ang;
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

struct cpt_pt {
	int16_t  alt;
	float    lat;
	float    lon;
	uint64_t seconds;
	void    *params;
};

struct cpt_ptps {
	struct cpt_pt *pt;
	struct cpt_ps *ps;
};

struct cpt_ff {
	struct cpt_header *hdr;
	struct cpt_ptps   *data;
	void  *ending;
};


