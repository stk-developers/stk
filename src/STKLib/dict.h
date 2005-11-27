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

#ifndef DICT_H
#define DICT_H

#include "hmms.h"
#include "common.h"

extern const char *dict_filter;

typedef struct _Word Word;
typedef struct _Pronun Pronun;

struct _Word {
  char *mpName;
  Pronun **pronuns;
  int npronunsInDict;  // Number of word pronuns found in dictionary file
  int npronuns;        // Can be higher than npronunsInDict as new empty
                       // pronuns can be added while reading the network
};


union _Model 
{
//    Hmm  *hmm;
    char *mpName;
};

struct _Pronun {
  Word  *word;
  char  *outSymbol;
  int   variant_no;
  FLOAT prob;
  int   nmodels;
  union _Model *model;
};

#ifdef __cplusplus
  extern "C" {
#endif

void ReadDictionary(
  const char *dictFileName,
  struct my_hsearch_data *wordHash,
  struct my_hsearch_data *phoneHash);

void FreeDictionary(struct my_hsearch_data *wordHash);

#ifdef __cplusplus
}
#endif

#endif // DICT_H
