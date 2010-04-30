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

#include "dict.h"
#include "common.h"
#include "mymath.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

const char *dict_filter;

using namespace STK;

void ReadDictionary(
    const char *dictFileName,
    MyHSearchData *wordHash,
    MyHSearchData *phoneHash)
{
  char line[1024];
  int line_no = 0;
  FILE *fp;
  Word *word;
  char *word_name, *model_name, *chrptr, *tchrptr;
  Pronun *new_pronun;

  if ((fp = STK::my_fopen(dictFileName, "rt", dict_filter)) == NULL) {
    STK::Error("Cannot open file: '%s'", dictFileName);
  }
  while (fgets(line, sizeof(line), fp)) {
    ENTRY e, *ep;
    
    line_no++;
    if (strlen(line) == sizeof(line)-1) {
      STK::Error("Line is too long (%s:%d)", dictFileName, line_no);
    }
    
    line[strlen(line)-1] = '\0'; // remove final '\n' caracter
    if (STK::getHTKstr(word_name = line, &chrptr)) {
      STK::Error("%s (%s:%d)", chrptr, dictFileName, line_no);
    }
    e.key  = word_name;
    e.data = NULL;
    my_hsearch_r(e, FIND, &ep, wordHash);

    if (ep != NULL) {
      word = (Word *) ep->data;
    } else {
      e.key  = strdup(word_name);
      word = (Word *) malloc(sizeof(Word));
      if (e.key == NULL || word  == NULL) STK::Error("Insufficient memory");
      word->mpName = e.key;
      word->npronuns = 0;
      word->pronuns  = NULL;
      e.data = word;

      if (!my_hsearch_r(e, ENTER, &ep, wordHash)) {
        STK::Error("Insufficient memory");
      }
    }

    word->npronuns++;
    word->npronunsInDict = word->npronuns;

    word->pronuns = (Pronun **) realloc(word->pronuns, word->npronuns * sizeof(Pronun *));
    new_pronun = (Pronun *) malloc(sizeof(Pronun));
    if (word->pronuns == NULL || new_pronun == NULL) STK::Error("Insufficient memory");

    word->pronuns[word->npronuns-1] = new_pronun;
    new_pronun->mpWord       = word;
    new_pronun->outSymbol  = word->mpName;
    new_pronun->nmodels    = 0;
    new_pronun->model      = NULL;
    new_pronun->variant_no = word->npronuns;

    if (*chrptr == '[') {

      model_name = ++chrptr;
      while (*chrptr && *chrptr != ']') chrptr++;

      if (*chrptr != ']') STK::Error("Matching `]' is missing (%s:%d)",
                                model_name, dictFileName, line_no);
      *chrptr++ = '\0';
      new_pronun->outSymbol = *model_name != '\0' ?  strdup(model_name) : NULL;
    }

    Str2Number(chrptr, &new_pronun->prob, &tchrptr);
    
    if(tchrptr != chrptr) {
      new_pronun->prob = my_log(new_pronun->prob);
      chrptr = tchrptr;
    }

    while (*chrptr) {
      if (STK::getHTKstr(model_name = chrptr, &chrptr)) {
        STK::Error("%s (%s:%d)", chrptr, dictFileName, line_no);
      }

      e.key  = model_name;
      e.data = NULL;
      my_hsearch_r(e, FIND, &ep, phoneHash);

      if (ep == NULL) {
        e.key  = strdup(model_name);
        e.data = e.key;

        if (e.key == NULL || !my_hsearch_r(e, ENTER, &ep, phoneHash)) {
          STK::Error("Insufficient memory");
        }
        ep->data = e.data;
      }
      new_pronun->nmodels++;
      new_pronun->model=(_Model*) realloc(new_pronun->model,
                                new_pronun->nmodels*sizeof(*new_pronun->model));
      if (new_pronun->model == NULL) STK::Error("Insufficient memory");
      new_pronun->model[new_pronun->nmodels-1].mpName = static_cast <char*> (ep->data);
    }
  }
  if (ferror(fp) || STK::my_fclose(fp)) {
    STK::Error("Cannot read dictionary file %s", dictFileName);
  }
}

void FreeDictionary(MyHSearchData *wordHash) {
  size_t i;
  size_t j;
  
  for (i = 0; i < wordHash->mNEntries; i++) 
  {
    Word *word = (Word *) wordHash->mpEntry[i]->data;
    for (j = 0; j < static_cast<size_t>(word->npronuns); j++) {
      free(word->pronuns[j]);
    }
    free(word->pronuns);
    free(wordHash->mpEntry[i]->data);
    free(wordHash->mpEntry[i]->key);
  }
  my_hdestroy_r(wordHash, 0);
}
