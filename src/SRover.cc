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

#define VERSION "0.5 "__TIME__" "__DATE__

// This alignment is kind of wierd - check it!!!
//six     !NULL   !NULL   !NULL   !NULL
//oh      !NULL   !NULL   !NULL   !NULL
//six     six     six     six     six
//!NULL   !NULL   !NULL   !NULL   oh

#include "STKLib/labels.h"
#include "STKLib/common.h"
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <malloc.h>
#include <assert.h>

using namespace STK;


void usage(char *progname)
{
  char *tchrptr;
  if ((tchrptr = strrchr(progname, '\\')) != NULL) progname = tchrptr+1;
  if ((tchrptr = strrchr(progname, '/')) != NULL) progname = tchrptr+1;
  fprintf(stderr,
"\nUSAGE: %s [options] labelList MLF1 [weight1] MLF2 [weight2]...\n\n"
" Option                                         Default\n\n"
"  -b       Random choice in the case of a tie   Take first\n"
"  -c i i i Substitution/insertion/deletion cost 10 7 7\n"
"  -i s     Output transcriptions to MLF s       Off\n"
"  -l s     Dir to store label files             Current\n"
"  -o  s    Output label formating NCST          None\n"
"  -f f     Word frequencyvoting weight          1.0\n"
"  -a f     Average confidence voting weight     0.0\n"
"  -m f     Maximum confidence voting weight     0.0\n"
"  -s f     Confidence sum voting weight         0.0\n"
"  -t       Output aligned trascriptions         Off\n"
"  -r i     Time-mediated alignment (1 2)        0\n"
"  -n f     Confidence of !NULL                  0.0\n"
"  -u s     Ignore labels not in labelList       Off\n"
"  -y s     Output label file extension          rec\n"
"  -A       Print command line arguments         Off\n"
//"  -I mlf   Load master label file mlf\n"
//"  -L dir   Set input label (or net) dir         Current\n"
//"  -S f     Set script file to f                 None\n"
"  -T N     Set trace flags to N                 0\n"
"  -V       Print version information            Off\n"
//"  -X ext   Set input label (or net) file ext    lab\n"
"\n"
" %s is Copyright (C) 2004 Lukas Burget et al. and\n"
" licensed under the GNU General Public License, version 2.\n"
" Bug reports, feedback, etc, to: burget@fit.vutbr.cz\n"
"\n", progname, progname);
  exit(-1);
}

//typedef struct _StrListElem StrListElem;
struct StrListElem 
{
  StrListElem *mpNext;
  char str[1];
};

StrListElem *MLF_files = NULL;
StrListElem *MLF_fn = NULL;
StrListElem **last = &MLF_files;

int main(int argc, char *argv[]) {
  int trace_flag       = 0;
  int print_cmdline    = 0;
  int print_version    = 0;
  int output_alignment = 0;
  int ignore_unknown_labels = 0;
  int rnd_from_best    = 0;

  LabelFormat in_lbl_frm = {0};
  LabelFormat out_lbl_frm = {0};
  char *out_lbl_dir    = NULL;
  char *out_lbl_ext    = "alg";
  char *out_MLF_fn     = NULL;
  int  *lab2cand       = NULL; // label index to candidate index mapping (used for voting)
//  char **labelList     = NULL;
  MyHSearchData labelHash = {0};
//  int  nlabels;
  int i;


  struct candidate {
    Label* label;
    int    votes;
    float  confsum;
    float  maxconf;
  } *candidates = NULL;

  float wrdfreq_w      = 1.0;
  float avgconf_w      = 0.0;
  float maxconf_w      = 0.0;
  float confsum_w      = 0.0;

  #define VOTING_SCORE(cand) ((cand).votes   / (float) nMLFs * wrdfreq_w + \
                              (cand).confsum / (float) nMLFs * confsum_w + \
                              (cand).confsum / (cand).votes  * avgconf_w + \
                              (cand).maxconf                 * maxconf_w)
  int  nMLFs    = 0;
  FILE **MLFfps = NULL;
  FILE *out_MLF_fp = NULL;
//  char line[1024];
  char label_file[1024];
//  char *base_name;
//  char *extension;
  char *chrptr;
//  Transcription transcription, hypothesis;
  Label *tlptr, *llptr, *labels;
  long *weights = NULL;
  double dval;

  if (argc == 1) usage(argv[0]);

  for (;;) {
    int opt = getopt(argc, argv, "-bi:l:y:o:r:tf:a:m:s:n:ucAT:V");

    if (opt == -1) break;

    switch (opt) {
      case 'A': print_cmdline    = 1;  break;
      case 'V': print_version    = 1;  break;
      case 't': output_alignment = 1;  break;
      case 'u': ignore_unknown_labels = 1; break;
      case 'b': rnd_from_best    = 1;  break;

      case 'i': out_MLF_fn   = optarg; break;
      case 'l': out_lbl_dir  = optarg; break;
      case 'y': out_lbl_ext  = optarg; break;

      case 1:  if (labelHash.mTabSize == 0) {
                  labelHash = readLabelList(optarg);
               } else {
                  *last = (StrListElem *) malloc(sizeof(StrListElem)+strlen(optarg));
                  if (!*last) Error("Insufficient memory");
                  chrptr = strcpy((*last)->str, optarg);
                  for (; *chrptr; chrptr++) if (*chrptr == '\\') *chrptr = '/';
                  last = &(*last)->mpNext;
                  *last = NULL;
                  ++nMLFs;

                  weights = (long int *) realloc(weights, sizeof(long) * nMLFs);
                  weights[nMLFs-1] = 1;

                  if (optind < argc) {
                    long weight = strtol(argv[optind], &chrptr, 0);
                    if (!*chrptr) {
                      optind++;
                      weights[nMLFs-1] = weight;
                    }
                  }
                }
                break;

      case 'c': if (optind < argc) {
                  substitution_cost = strtol(argv[optind], &chrptr, 0);
                  if (!*chrptr && ++optind < argc) {
                    insertion_cost = strtol(argv[optind], &chrptr, 0);
                    if (!*chrptr && ++optind < argc) {
                      deletion_cost = strtol(argv[optind], &chrptr, 0);
                    }
                  }
                }
                if (optind >= argc || *chrptr) {
                  Error("Three integer number are expected after option -%c", opt);
                }
                optind++;
                break;

      case 'f':
      case 'a':
      case 'm':
      case 's':
      case 'n': dval = strtod(optarg, &chrptr);
                if (!*optarg || *chrptr) {
                  Error("Decimal number is expected after option -%c", opt);
                }
                switch (opt) {
                  case 'f': wrdfreq_w = dval; break;
                  case 'a': avgconf_w = dval; break;
                  case 'm': maxconf_w = dval; break;
                  case 's': confsum_w = dval; break;
                  case 'n': null_cost = dval; break;
                }
                break;
                
      case 'o': for (chrptr = optarg; *chrptr; chrptr++) {
                  switch (*chrptr) {
                    case 'N': out_lbl_frm.SCORE_NRM = 1; break;
                    case 'S': out_lbl_frm.SCORE_OFF = 1; break;
                    case 'C': out_lbl_frm.CENTRE_TM = 1; break;
                    case 'T': out_lbl_frm.TIMES_OFF = 1; break;
                    default:
                      Error("Some of characters 'NCST' are expected after option -%c", opt);
                  }
                }
                break;


      case 'r': time_med_alg = strtol(optarg, &chrptr, 0);
                if (!*optarg || *chrptr) {
                  Error("Integer number is expected after option -%c", opt);
                }
                break;

      case 'T': trace_flag = strtol(optarg, &chrptr, 0);
                if (!*optarg || *chrptr) {
                  Error("Integer number is expected after option -%c", opt);
                }
                break;

      case '?': usage(argv[0]);

      default:  Error("Option -%c not suported yet.", opt);
    }
  }

  if (print_version) {
    puts("Version: "VERSION"\n");
    if (!nMLFs) exit(0);
  }

  if (!nMLFs) usage(argv[0]);

  if (print_cmdline) {
    for (i=0; i < argc; i++) {
      fputs(argv[i], stdout);
      putchar(' ');
    }
    putchar('\n');
  }

  MLFfps = (FILE **) malloc(sizeof(FILE *) * nMLFs);
  if (!MLFfps) Error("Insufficient memory");

  lab2cand = (int *) malloc(sizeof(int) * (labelHash.mNEntries + 1)); //one more for !NULL
  if (!lab2cand) Error("Insufficient memory");
  for (i = 0; i < labelHash.mNEntries + 1; i++) lab2cand[i] = -1;

  candidates = (struct candidate *) malloc(sizeof(struct candidate) * nMLFs);
  if (!candidates) Error("Insufficient memory");

  for (MLF_fn = MLF_files, i=0; MLF_fn; MLF_fn = MLF_fn->mpNext, i++) {
    MLFfps[i] = OpenInputMLF(MLF_fn->str);
  }

  out_MLF_fp = OpenOutputMLF(out_MLF_fn);

  ///////// For each record (label file) from MLFs////////
  for (;;) {

    //// Read one record from each MLF and align them ////

    label_file[0]='\0'; // For 1st MLF name of label file is not given ->
    // - 1st MLF is read sequentially record by record.
    // - In 2nd to nth MLF, records are searched by name obtained from 1st MLF.

    for (i=0, MLF_fn=MLF_files; i < nMLFs; i++, MLF_fn=MLF_fn->mpNext) {
      Label *tlabs;

      if (!OpenInputLabelFile(label_file, NULL, NULL, MLFfps[i], MLF_fn->str)) {
        return 0; // No more records in first MLF
      }

      tlabs = ReadLabels(MLFfps[i], &labelHash,
                           ignore_unknown_labels ? UL_IGNORE : UL_ERROR,
                           in_lbl_frm, 100000, label_file, MLF_fn->str, NULL);

      if (i == 0) labels = tlabs;
      else AlingTranscriptions(&labels, tlabs, weights);
    }

    out_MLF_fp = OpenOutputLabelFile(label_file, out_lbl_dir, out_lbl_ext,
                                     out_MLF_fp, out_MLF_fn);

    if (output_alignment) {
      WriteLabels(out_MLF_fp, labels, out_lbl_frm, 100000, label_file, out_MLF_fn);
    } else {

      ///////// For each corespondence set... ////////
      for (tlptr = labels; tlptr != NULL; tlptr=tlptr->mpNext) {

        ////////////// Do the voting //////////////

        int ncandidates = 0, candidate;
        struct candidate *best_candidate;
        float  voting_score, best_voting_score;

        // count votes
        for (llptr = tlptr; llptr !=NULL; llptr = llptr->mpNextLevel) {

          candidate = lab2cand[(int) llptr->mpData];

          if (candidate == -1) { // new candidate
            lab2cand[(int) llptr->mpData] = candidate = ncandidates++;
            candidates[candidate].label = llptr;
            candidates[candidate].votes = 0;
            candidates[candidate].maxconf = llptr->mScore;
            candidates[candidate].confsum = 0.0;
          }

          candidates[candidate].votes++;
          candidates[candidate].confsum += llptr->mScore;
          if (llptr->mScore > candidates[candidate].maxconf) {
            candidates[candidate].maxconf = llptr->mScore;
          }
        }

        // reset lab2cand[] to -1
        for (candidate=0; candidate<ncandidates; candidate++) {
          lab2cand[(int) candidates[candidate].label->mpData] = -1;
        }

        //choose the best word
        best_candidate    = &candidates[0];
        best_voting_score = VOTING_SCORE(candidates[0]);

        for (candidate=1; candidate<ncandidates; candidate++) {
          voting_score = VOTING_SCORE(candidates[candidate]);
          if (voting_score > best_voting_score) {
            best_candidate    = &candidates[candidate];
            best_voting_score = voting_score;
          }
        }

        if (rnd_from_best) {
          int bestCnt = 0;

          for (candidate=1; candidate<ncandidates; candidate++) {
            if (VOTING_SCORE(candidates[candidate]) == best_voting_score) {
              bestCnt++;
            }
          }

          if (bestCnt > 1) {
            bestCnt = 1+(int) ((float) bestCnt * rand() / (RAND_MAX+1.0));
            for (candidate=1; candidate<ncandidates; candidate++) {
              if (VOTING_SCORE(candidates[candidate]) == best_voting_score
                 && --bestCnt == 0) {
                best_candidate = &candidates[candidate];
                break;
              }
            }
          }
        }

        if (best_candidate->label->mpName) {
          Label label = *best_candidate->label;
          label.mpNext  = label.mpNextLevel = NULL;
          label.mScore = best_voting_score;
          WriteLabels(out_MLF_fp, &label, out_lbl_frm, 100000, label_file, out_MLF_fn);
        }
      }
    }
    CloseOutputLabelFile(out_MLF_fp, out_MLF_fn);
    ReleaseLabels(labels);  
  }
  return 0;
}
