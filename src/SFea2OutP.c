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

#define VERSION "0.1 12/19/2002"
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <malloc.h>
#include "STKLib/fileio.h"
#include "STKLib/common.h"
#include "STKLib/Models.h"
#include "STKLib/Decoder.h"

using namespace STK;


void usage(char *progname)
{
  char *tchrptr;
  if ((tchrptr = strrchr(progname, '\\')) != NULL) progname = tchrptr+1;
  if ((tchrptr = strrchr(progname, '/')) != NULL) progname = tchrptr+1;
  fprintf(stderr,
"\nUSAGE: %s [options] DataFiles...\n\n"
" Option                                       Default\n\n"
" -l s       Dir to store output prob. files                 Feature file dir\n"
" -y s       Output file extension                           lop\n"
" -A         Print command line arguments                    Off\n"
" -D [N ...] Append derivatives using window lengths N...    No change\n"
" -H mmf     Load HMM macro file mmf\n"
" -N         Natural read order for feature file             Swap bytes\n"
" -S f       Set script file to f                            None\n"
" -T N       Set trace flags to N                            0\n"
" -V         Print version information                       Off\n"
"\n"
" %s is Copyright (C) 2004 Lukas Burget et al. and\n"
" licensed under the GNU General Public License, version 2.\n"
" Bug reports, feedback, etc, to: burget@fit.vutbr.cz\n"
"\n", progname, progname);
  exit(-1);
}

int main(int argc, char *argv[]) {
  HtkHeader header;
  ModelSet hset;
  FILE *sfp, *ofp = NULL;
  FLOAT  *obsMx;
  FLOAT *opp;
  int i, j, m;
  char line[1024];
  char outFile[1024];
  char *chrptr;
  long lval, time;

  typedef struct _StrListElem StrListElem;
  struct _StrListElem {
    StrListElem *next;
    char *mpPhysical;
    char logical[1];
  };

  StrListElem *feature_files = NULL;
  int nfeature_files = 0;
  StrListElem *file_name = NULL;
  StrListElem **last = &feature_files;

  Macro *macro;

  char *out_dir   = NULL;
  char *out_ext   = "lop";

  int trace_flag      =  0;
  int print_cmdline   =  0;
  int print_version   =  0;
  int swap_features   =  1;
  int derivOrder      = -1;
  int *derivWinLengths= NULL;
  int startFrmExt     =  0;
  int endFrmExt       =  0;
  int targetKind      =  PARAMKIND_ANON;


  if (argc == 1) usage(argv[0]);

  //InitHMMSet(&hset, 0);
  hset.Init();

  for (;;) {
    int opt = getopt(argc, argv, "-l:y:ADH:NS:T:V");

    if (opt == -1) break;

    switch (opt) {
      case 'l': out_dir  = optarg; break;
////      case 'x': hmm_ext      = optarg; break;
      case 'y': out_ext  = optarg; break;

      case 'N': swap_features = 0;     break;
      case 'V': print_version = 1;     break;

      case 'H': hset.ParseMmf(optarg, NULL); break;

      case 1:   *last = (StrListElem *) malloc(sizeof(StrListElem)+strlen(optarg));
                if (!*last) Error("Insufficient memory");
                chrptr = strcpy((*last)->logical, optarg);
                for (; *chrptr; chrptr++) if (*chrptr == '\\') *chrptr = '/';
                chrptr = strchr((*last)->logical, '=');
                if (chrptr) *chrptr = '\0';
                (*last)->mpPhysical = chrptr ? chrptr+1: (*last)->logical;
                last = &(*last)->next;
                *last = NULL;
                nfeature_files++;
                break;


      case 'S': if ((sfp = fopen(optarg, "rt")) == NULL) {
                  Error("Cannot open script file %s", optarg);
                }

                while (fscanf(sfp, "%s", line) == 1) {
                  *last = (StrListElem*) malloc(sizeof(StrListElem)+strlen(line));
                  if (!*last) Error("Insufficient memory");
                  chrptr = strcpy((*last)->logical, line);
                  for (; *chrptr; chrptr++) if (*chrptr == '\\') *chrptr = '/';
                  chrptr = strchr((*last)->logical, '=');
                  if (chrptr) *chrptr = '\0';
                  (*last)->mpPhysical = chrptr ? chrptr+1: (*last)->logical;
                  last = &(*last)->next;
                  nfeature_files++;
                }
                *last = NULL;
                fclose(sfp);
                break;

      case 'T': trace_flag = strtol(optarg, &chrptr, 0);
                if (!*optarg || *chrptr) {
                  Error("Integer number is expected after option -T");
                }
                break;

      case 'D': derivOrder = 0;
                for (; optind < argc; optind++) {
                lval = strtol(argv[optind], &chrptr, 0);
                if (!*argv[optind] || *chrptr) {
                  break;
                }

                derivWinLengths = (int *)
                  realloc(derivWinLengths, ++derivOrder * sizeof(int));
                if (derivWinLengths == NULL) Error("Insufficient memory");
                derivWinLengths[derivOrder-1] = lval;
              }
              break;


      case '?': usage(argv[0]);

      default:  Error("Option -%c not suported yet.", opt);
    }
  }

  if (print_cmdline) {
    for (i=0; i < argc; i++) {
      fputs(argv[i], stdout);
      putchar(' ');
    }
    putchar('\n');
  }

  if (print_version) puts("Version: "VERSION"\n");
  if (feature_files == NULL) return 0;
  if ((opp = (FLOAT *) malloc(hset.mNStates * sizeof(float))) == NULL) {
    Error("Insufficient memory");
  }

  if (hset.outPDF_kind != KID_DiagC) {
    Error("Invalid HMM set distribution kind <%s>", gpKwds[hset.outPDF_kind]);
  }

  for (file_name = feature_files; file_name; file_name = file_name->mpNext) {
    obsMx = ReadHTKFeatures(file_name->mpPhysical, swap_features,
                            startFrmExt, endFrmExt, targetKind,
                            derivOrder, derivWinLengths, &header);

    if (hset.mInputVectorSize != header.mSampleSize / sizeof(float)) {
      Error("Vector size [%d] in '%s' is incompatible with HMM set [%d]",
            header.mSampleSize/sizeof(float), file_name->mpPhysical, hset.mInputVectorSize);
    }

    MakeFileName(outFile, file_name->logical, out_dir, out_ext);

    if ((ofp = fopen(outFile, "wb")) == NULL) {
      Error("Cannot open output probability file: '%s'", outFile);
    }

    header.mSampleKind = 9;
    header.mSampleSize = hset.mNStates * sizeof(float);

    if (WriteHTKHeader(ofp, header, 1)) {
      Error("Cannot write to output probability file: '%s'", outFile);
    }

    time = -hset.mTotalDelay;
    ResetXFormInstances(&hset);

    for (i = 0; i < header.mNSamples; i++) {
      UpdateStacks(&hset, obsMx + i * hset.mInputVectorSize, ++time, FORWARD);
      if (time <= 0) continue;

      for (m = 0; m < hset.mHmmHash.mNEntries; m++) {
        macro = (Macro *) hset.mHmmHash.mpEntry[m]->data;
        if (macro->mpData->mpMacro != macro) continue;

        for (j = 0; j < ((Hmm *) macro->mpData)->mNStates-2; j++) {
          State *state = ((Hmm *) macro->mpData)->state[j];
          if (!state->mpMacro) { // take only non-shared states
            opp[state->mID] =
              DiagCGaussianMixtureDensity(state, obsMx + i * hset.mInputVectorSize, NULL);
          }
        }
      }

      for (m = 0; m < hset.mStateHash.mNEntries; m++) {
        State * state;
        macro = (Macro *) hset.mStateHash.mpEntry[m]->data;
        if (macro->mpData->mpMacro != macro) continue;

        state = (State *) macro->mpData;
        opp[state->mID] =
          DiagCGaussianMixtureDensity(state, obsMx + i * hset.mInputVectorSize, NULL);
      }

      if (WriteHTKFeature (ofp, opp, hset.mNStates, 1)) {
        Error("Cannot write to output probability file: '%s'", outFile);
      }
    }

    fclose(ofp);
    free(obsMx);
  }
  return 0;
}
