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

#define VERSION "0.1 "__TIME__" "__DATE__
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <malloc.h>
#include <assert.h>
#include "STKLib/fileio.h"
#include "STKLib/common.h"
#include "STKLib/Models.h"
#include "STKLib/Viterbi.h"
#include "STKLib/labels.h"
#ifndef WIN32
#include <unistd.h>
#else
#include "getopt.h"
#endif

void ReadXformList(HMMSet *hmm_set, char *xformListFileName);

void usage(char *progname)
{
  char *tchrptr;
  if ((tchrptr = strrchr(progname, '\\')) != NULL) progname = tchrptr+1;
  if ((tchrptr = strrchr(progname, '/')) != NULL) progname = tchrptr+1;
  fprintf(stderr,
"\nUSAGE: %s [options] hmmList SER1.acc [weight1] SER2.acc [weight2]...\n\n"
" Option                                                     Default\n\n"
//" -d s     dir to find hmm definitions                     current
" -h s       Save combined accumulators to file s            ./SER0.acc\n"
" -l s       Preserve stats for xforms listed in file s      Off\n"
" -n         Normalize accumulators                          Off\n"
" -s s       Print statistics to file s                      Off\n"
//" -x s     Extension for hmm files                         none
" -A         Print command line arguments                    Off\n"
////" -D       Display configuration variables                 Off
" -M dir     Dir to write accumulator file                   Current\n"
" -H mmf     Load HMM macro file mmf                         \n"
//" -S f       Set script file to f                            None\n"
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
  HMMSet hset;
  double dval;
  double *weights = NULL;
  long   lval;
  char   *chrptr;
  int    i;
  FLOAT totLogLike      = 0;
  int   totFrames       = 0;
//  MyHSearchData labelHash;

  typedef struct _StrListElem StrListElem;
  struct _StrListElem {
    StrListElem *next;
    char str[1];
  };

  StrListElem *acc_files = NULL;
  int nacc_files = 0;
  StrListElem *file_name = NULL;
  StrListElem **last = &acc_files;

  char *out_acc         = NULL;
  char *out_acc_dir     = NULL;
  char *xform_list_file = NULL;
  char *hmm_list_file   = NULL;
  char *stat_file       = NULL;

  int trace_flag      =  0;
  int print_cmdline   =  0;
  int print_version   =  0;
  int normalize       =  0;

  if (argc == 1) usage(argv[0]);

  //InitHMMSet(&hset, 1);
  hset.Init(1);

  for (;;) {
    int opt = getopt(argc, argv, "-h:l:ns:AH:M:T:V");

    if (opt == -1) break;

    switch (opt) {
      case 'A': print_cmdline   = 1;       break;
      case 'V': print_version   = 1;       break;
      case 'n': normalize       = 1;       break;


      case 'l': xform_list_file = optarg;  break;
      case 'h': out_acc         = optarg;  break;
      case 's': stat_file       = optarg;  break;
      case 'M': out_acc_dir     = optarg;  break;

      case 'H': ReadHMMSet(optarg, &hset);
                break;

      case 1:   if (hmm_list_file == NULL) {
                  hmm_list_file = optarg;
                  break;
                }

                *last = (StrListElem *) malloc(sizeof(StrListElem)+strlen(optarg));
                if (!*last) Error("Insufficient memory");
                chrptr = strcpy((*last)->str, optarg);
                for (; *chrptr; chrptr++) if (*chrptr == '\\') *chrptr = '/';
                last = &(*last)->next;
                *last = NULL;
                nacc_files++;

                weights = realloc(weights, sizeof(double) * nacc_files);
                weights[nacc_files-1] = 1.0;
                if (optind < argc) {
                  dval = strtod(argv[optind], &chrptr);
                  if (!*chrptr) {
                    optind++;
                    weights[nacc_files-1] = dval;
                  }
                }
                break;

//      case 'w': dval = strtod(optarg, &chrptr);
//                if (!*optarg || *chrptr) {
//                  Error("Decimal number is expected after option -%c", opt);
//                }
//                switch (opt) {
//                  case 'w': min_mix_weight = dval; break;
//                }
//                break;

      case 'T': lval = strtol(optarg, &chrptr, 0);
                if (!*optarg || *chrptr) {
                  Error("Integer number is expected after option -%c", opt);
                }
                switch (opt) {
                  case 'T': trace_flag    = lval; break;
                }
                break;

      case '?': usage(argv[0]);

      default:  Error("Option -%c not suported yet.", opt);
    }
  }

  if (print_version) {
    puts("Version: "VERSION"\n");
    if (hmm_list_file == NULL) exit(0);
  }

//  if (hmm_list_file == NULL) usage(argv[0]);

  if (print_cmdline) {
    for (i=0; i < argc; i++) {
      fputs(argv[i], stdout);
      putchar(' ');
    }
    putchar('\n');
  }

//  labelHash = ReadHMMList(&hset, &hset, hmm_list_file);


  if (xform_list_file != NULL) {
    ReadXformList(&hset, xform_list_file);
  }

  hset.AllocateAccumulatorsForXformStats();
  //ResetAccumsForHMMSet(&hset);
  hset.ResetAccums();

  for (i=0,file_name=acc_files;file_name!=NULL;i++,file_name=file_name->mpNext) {
    if (trace_flag & 1) {
      TraceLog("Processing file '%s'", file_name->str);
    }

    FLOAT P;
    long S;

    //ReadAccums(file_name->str, weights[i], &hset, &S, &P, 0);
    hset.ReadAccums(file_name->str, weights[i], &S, &P, 0);
    totFrames  += S;
    totLogLike += P;

    if (trace_flag & 1) {
      TraceLog("Frames: %d,  Log likelihood: %e", S, P);
    }
  }
  if (trace_flag & 2) {
    TraceLog("Total number of frames: %d\nTotal log likelihood: %e",
             totFrames, totLogLike);
  }
  if (stat_file) WriteHMMStats(stat_file, &hset);
  if (normalize) NormalizeAccums(&hset);

  WriteAccums(out_acc ? out_acc : "SER0.acc", out_acc_dir, &hset,
              totFrames, totLogLike);

  //ReleaseHMMSet(&hset);
  hset.Release();

  while (acc_files) {
    file_name = acc_files;
    acc_files = acc_files->mpNext;
    free(file_name);
  }

  return 0;
}
