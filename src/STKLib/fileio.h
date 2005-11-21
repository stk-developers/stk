/***************************************************************************
 *   copyright            : (C) 2004 by Lukas Burget,UPGM,FIT,VUT,Brno     *
 *   email                : burget@fit.vutbr.cz                            *
 ***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef FILEIO_H
#define FILEIO_H

#include <string>
#include <stdio.h>
#include "common.h"

using namespace std;

typedef struct {
    INT_32 nSamples;              /* Structure for HTK header */
    INT_32 sampPeriod;
    INT_16 sampSize;
    UINT_16 sampKind;
}
HTK_Header;


#ifdef __cplusplus
  extern "C" {
#endif


int WriteHTKHeader (FILE * fp_out, HTK_Header header, bool swap);
int WriteHTKFeature (FILE * fp_out, FLOAT *out, size_t fea_len, bool swap);
int ReadHTKHeader (FILE * fp_in, HTK_Header *header, bool swap);

int ReadHTKFeature (FILE * fp_in, FLOAT *in, size_t fea_len, bool swap,
                    bool decompress, FLOAT *A, FLOAT *B);

int Mkdir4File(const char *file_name);

typedef struct {
  char *last_file_name;
  char *last_cmn_file;
  char *last_cvn_file;
  char *last_cvg_file;
  FILE *fp;
  FLOAT *cmn;
  FLOAT *cvn;
  FLOAT *cvg;
  HTK_Header last_header;
  FLOAT *A;
  FLOAT *B;
} RHFBuffer;

FLOAT *ReadHTKFeatures(
  char *file_name,
  bool swap,
  int extLeft,
  int extRight,
  int targetKind,
  int derivOrder,
  int *derivWinLen,
  HTK_Header *header,
  const char *cmn_file,
  const char *cvn_file,
  const char *cvg_file,
  RHFBuffer *buff);

#ifdef __cplusplus
}
#endif



#endif // FILEIO_H
