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
#include <time.h>
#include <unistd.h>
#include <fcntl.h>

#include "hdf5.h"
#include "../../read/readcpt.h"


/*  Error numbers  */
enum WR_CPT_ERR {
WR_CPT_EINVARG = 1,
WR_CPT_E
};

/*  Declaration of main function  */
int wrcpt(const char *ptname, const char *psname, const char *cptname);

/*  Program entrance  */
int main(int argc, char *argv[]) {
	return (3 == argc) ? wrcpt(argv[1], argv[2], argv[3]) : WR_CPT_EINVARG;
}

/*  Definition of main function  */
int wrcpt(const char *ptname, const char *psname, const char *cptname)
{
	
	return 0;
}

