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
#include "STKLib/Models.h"
#include "STKLib/common.h"
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <malloc.h>
#include <assert.h>
#ifndef WIN32
#include <unistd.h>
#else
#include "getopt.h"
#endif

void usage(char *progname)
{
  char *tchrptr;
  if ((tchrptr = strrchr(progname, '\\')) != NULL) progname = tchrptr+1;
  if ((tchrptr = strrchr(progname, '/')) != NULL) progname = tchrptr+1;
  fprintf(stderr,
"\nUSAGE: %s [options] hmmList num1.acc den1.acc...          \n\n"
" Option                                                     Default\n\n"
////" -c f   Mixture pruning threshold                       10.0
//" -d s     dir to find hmm definitions                     current
" -e f [f]   Smoothing constant E [h]                        2.0 2.0\n"
" -h mmf     Save all HMMs to macro file mmf                 as src\n"
//" -l f       Update xforms listed in file f                  Off\n"
" -m N       Min examples needed per model                   3\n"
" -o s       Extension for new hmm files                     as src\n"
" -s s       Print statistics to file s                      Off\n"
" -t f       I-smoothing tau constant                        0.0\n"
" -u tmvwsxo Update t)rans m)eans v)ars w)ghts s)tats x)form tmvwsx\n"
" -v f       Set minimum variance to f                       0.0\n"
" -w f       Set mix weight floor to f*MINMIX                1.0\n"
//" -x s     Extension for hmm files                         none
" -A         Print command line arguments                    Off\n"
" -B         Save HMM macro files as binary                  Off\n"
////" -D       Display configuration variables                 Off
////" -F fmt   Set source data format to fmt                   as config
////" -G fmt   Set source label format to fmt                  as config
" -H mmf     Load HMM macro file mmf                         \n"
" -M dir     Dir to write HMM macro files                    Current\n"
" -S f       Set script file to f                            None\n"
" -T N       Set trace flags to N                            0\n"
" -V         Print version information                       Off\n"
//" -O f f   Trans and output prob scale factor  1.0 1.0     \n"
"\n"
" %s is Copyright (C) 2004 Lukas Burget et al. and\n"
" licensed under the GNU General Public License, version 2.\n"
" Bug reports, feedback, etc, to: burget@fit.vutbr.cz\n"
"\n", progname, progname);
  exit(-1);
}

int main(int argc, char *argv[]) {
  HMMSet hset;
  FILE *sfp;
  int i, fcnt = 0;
  char line[1024];
  char *chrptr;
  double dval;
  long   lval;
//  MyHSearchData labelHash;

  FLOAT totLogLike      = 0;
  FLOAT totLogPosterior = 0;
  int   totFrames       = 0;

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

  UpdateMask update_mask = 0;

  double min_variance   = 0.0;
  double min_mix_weight = 1.0;
  double E_constant     = 2.0;
  double h_constant     = 2.0;
  double I_smoothing    = 0.0;

  char *out_hmm_dir     = NULL;
  char *out_MMF         = NULL;
  char *out_hmm_ext     = NULL;
  char *xform_list_file = NULL;
  char *hmm_list_file   = NULL;
  char *stat_file       = NULL;

  int trace_flag      =  0;
  int print_cmdline   =  0;
  int print_version   =  0;
  int min_examples    =  3;
  int stats_binary    =  0;
  int hmms_binary     =  0;
  enum {UT_ML=0, UT_MMI, UT_MPE} update_type =  UT_MMI;


  if (argc == 1) usage(argv[0]);

  //InitHMMSet(&hset, 1);
  hset.Init(MODEL_SET_WITH_ACCUM);

  for (;;) {
    int opt = getopt(argc, argv, "-eh:l:m:o:s:t:u:v:w:ABH:M:S:T:U:V");
    if (opt == -1) break;

    switch (opt) {
      case 'A': print_cmdline   = 1;       break;
      case 'V': print_version   = 1;       break;
      case 'B': hmms_binary     = 1;       break;

      case 'l': xform_list_file = optarg;  break;
      case 'h': out_MMF         = optarg;  break;
      case 'M': out_hmm_dir     = optarg;  break;
      case 'o': out_hmm_ext     = optarg;  break;
      case 's': stat_file       = optarg;  break;

      case 'H': ReadHMMSet(optarg, &hset);
                break;

      case 1:   if (hmm_list_file == NULL) {
                  hmm_list_file = optarg;
                  break;
                }

                *last = (StrListElem *) malloc(sizeof(StrListElem)+strlen(optarg));
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

      case 'u': for (chrptr = optarg; *chrptr; chrptr++) {
                  switch (*chrptr) {
                    case 't': update_mask |= UM_TRANSITION; break;
                    case 'm': update_mask |= UM_MEAN;       break;
                    case 'v': update_mask |= UM_VARIANCE;   break;
                    case 'w': update_mask |= UM_WEIGHT;     break;
                    case 's': update_mask |= UM_XFSTATS;    break;
                    case 'x': update_mask |= UM_XFORM;      break;
                    case 'o': update_mask |= UM_OLDMEANVAR; break;
                    case 'c': update_mask |= UM_CWEIGHTS;   break;
                    default:
                      Error("Some of characters 'tmvwsxo' are expected after option -u");
                  }
                }
                break;

      case 't':
      case 'v':
      case 'w': dval = strtod(optarg, &chrptr);
                if (!*optarg || *chrptr) {
                  Error("Decimal number is expected after option -%c", opt);
                }
                switch (opt) {
                  case 't': I_smoothing    = dval; break;
                  case 'v': min_variance   = dval; break;
                  case 'w': min_mix_weight = dval; break;
                }
                break;


/*      case 'O': if (optind < argc) {
                  transp_scale = strtod(argv[optind], &chrptr);
                  if (!*chrptr && ++optind < argc) {
                    outprb_scale = strtod(argv[optind], &chrptr);
                    if (!*chrptr && ++optind < argc) {
                      occprb_scale = strtod(argv[optind], &chrptr);
                    }
                  }
                }

                if (++optind > argc || *chrptr) {
                  Error("Three decimal number are expected after option -%c", opt);
                }
                break;*/

      case 'e': if (optind < argc) {
                  E_constant = strtod(argv[optind], &chrptr);
                }

                if (optind++ >= argc || *chrptr) {
                  Error("One or two decimal number are expected after option -%c", opt);
                }

                if (optind >= argc) break;

                dval = strtod(argv[optind], &chrptr);

                if (!*chrptr) {
                  ++optind;
                  h_constant = dval;
                }
                break;

      case 'T':
      case 'm': lval = strtol(optarg, &chrptr, 0);
                if (!*optarg || *chrptr) {
                  Error("Integer number is expected after option -%c", opt);
                }
                switch (opt) {
                  case 'T': trace_flag    = lval; break;
                  case 'm': min_examples  = lval; break;
                }
                break;

      case 'U': for (chrptr=optarg; *chrptr; chrptr++) *chrptr = toupper(*chrptr);
                if (!strcmp(optarg, "ML"))       update_type = UT_ML;
                else if (!strcmp(optarg, "MPE")) update_type = UT_MPE;
                else if (!strcmp(optarg, "MMI")) update_type = UT_MMI;
                else Error("ML, MMI or MPE expected after option -%c but found '%s'",
                           opt, optarg);
                break;


      case '?': usage(argv[0]);

      default:  Error("Option -%c not suported yet.", opt);
    }
  }

  if (print_version) {
    puts("Version: "VERSION"\n");
    if (hmm_list_file == NULL) exit(0);
  }

  if (hmm_list_file == NULL) usage(argv[0]);

  if (print_cmdline) {
    for (i=0; i < argc; i++) {
      fputs(argv[i], stdout);
      putchar(' ');
    }
    putchar('\n');
  }

  if (nfeature_files & 1) {
    Error("MMI update requires even number of accumulator files");
    nfeature_files /= 2;
  }

//  labelHash = ReadHMMList(&hset, &hset, hmm_list_file);


  if (xform_list_file != NULL) {
    ReadXformList(&hset, xform_list_file);
  }

  hset.AllocateAccumulatorsForXformStats();

  hset.mMmiUpdate = update_type;

  //ResetAccumsForHMMSet(&hset);
  hset.ResetAccums();

  hset.MMI_E               = E_constant;
  hset.MMI_h               = h_constant;
  hset.MMI_tauI            = I_smoothing;
  hset.minOccurances       = min_examples;
  hset.mMinMixWeight        = min_mix_weight * MIN_WEGIHT;
  hset.updateMask          = update_mask ? update_mask :
                             UM_TRANSITION | UM_MEAN    | UM_VARIANCE |
                             UM_WEIGHT     | UM_XFSTATS | UM_XFORM;

  if ((hset.updateMask & (UM_MEAN | UM_VARIANCE)) &&
     !hset.mAllMixuresUpdatableFromStatAccums) {
    Warning("Statistic are estimated for Xform not being "
            "a single linear Xform on the input of a mixture. "
            "Means and variances will not be updated");
    hset.updateMask &= ~(UM_MEAN | UM_VARIANCE);
  }


  for (file_name = feature_files; file_name != NULL; file_name = file_name->mpNext->mpNext) {
    if (trace_flag & 1) {
      TraceLog("Processing file pair %d/%d '%s' <-> %s",  ++fcnt,
              nfeature_files/2, file_name->mpPhysical,file_name->mpNext->mpPhysical);
    }

    FLOAT P;
    long S;

    //ReadAccums(file_name->mpPhysical, 1.0, &hset, &S, &P, 0);
    hset.ReadAccums(file_name->mpPhysical, 1.0, &S, &P, 0);
    totFrames  += S;
    totLogLike += P;

    if (trace_flag & 1) {
      TraceLog("[%d frames] %f", S, P/S);
    }

    FLOAT mmi_P;
    // Second set of accums is for compeating models
    //ReadAccums(file_name->mpNext->mpPhysical, 1.0, &hset, &S, &mmi_P, update_type);
    hset.ReadAccums(file_name->mpNext->mpPhysical, 1.0, &S, &mmi_P, update_type);
    totLogPosterior += P - mmi_P;

    if (trace_flag & 1) {
      TraceLog("[%d frames] %f", S, mmi_P/S);
    }
  }

  if (trace_flag & 2) {
      TraceLog("Total number of frames: %d\nTotal log likelihood: %e"
               "\nTotal log posterior: %e",
               totFrames, totLogLike, totLogPosterior);
  }

  if (stat_file) WriteHMMStats(stat_file, &hset);

  Macro *macro = FindMacro(&hset.mVarianceHash, "varFloor1");

    if (macro != NULL || (float) min_variance > 0.0) {
      Variance *tmpvar = macro ? (Variance *) macro->mpData : NULL;
//      assert(!tmpvar || hset.mInputVectorSize == tmpvar->VectorSize());

      //***
      // old malloc
      //hset.varFloor = (Variance *) malloc(sizeof(Variance)+((tmpvar ? tmpvar->VectorSize() : hset.mInputVectorSize)-1)*sizeof(FLOAT));
      //if (hset.varFloor == NULL) Error("Insufficient memory");
      //
      //hset.varFloor->VectorSize() = tmpvar ? tmpvar->VectorSize() : hset.mInputVectorSize;
      
      hset.varFloor = new Variance((tmpvar ? tmpvar->VectorSize() : hset.mInputVectorSize), false);

      for (i = 0; i < hset.varFloor->VectorSize(); i++) {
        if (macro) {
          
          hset.varFloor->mVector[i] =
             tmpvar && ((float) min_variance <= 0.0 ||
                        tmpvar->mVector[i] < 1/min_variance)
             ? tmpvar->mVector[i] : 1 / min_variance;
          
          // old vector
          //hset.varFloor->mpVectorO[i] =
          //   tmpvar && ((float) min_variance <= 0.0 ||
          //              tmpvar->mpVectorO[i] < 1/min_variance)
          //   ? tmpvar->mpVectorO[i] : 1 / min_variance;
        }
      }
    }

  // Required by WriteXformStatsAndRunCommands and UpdateHMMSetFromAccums
  ScanHMMSet(&hset, MTM_MEAN|MTM_VARIANCE, NULL, NormalizeStatsForXform, 0);

  if (hset.updateMask & UM_XFSTATS) {
    WriteXformStatsAndRunCommands(out_hmm_dir, stats_binary, &hset);
  }

  if (hset.updateMask != UM_XFSTATS) 
  {
    //UpdateHMMSetFromAccums(out_hmm_dir, &hset);
    hset.UpdateFromAccums(out_hmm_dir);
    
    WriteHMMSet(out_MMF, out_hmm_dir, out_hmm_ext, hmms_binary, &hset);
  }

  //ReleaseHMMSet(&hset);
  hset.Release();

//  my_hdestroy_r(&labelHash, 0);

  delete hset.varFloor;
  while (feature_files) {
    file_name = feature_files;
    feature_files = feature_files->mpNext;
    free(file_name);
  }

  return 0;
}
