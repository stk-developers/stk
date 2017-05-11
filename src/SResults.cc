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

#define SVN_DATE       "$Date$"
#define SVN_AUTHOR     "$Author$"
#define SVN_REVISION   "$Revision$"
#define SVN_ID         "$Id$"


#define MODULE_VERSION "2.0.1 "__TIME__" "__DATE__" "SVN_ID  
//#define MODULE_VERSION "0.1 "__TIME__" "__DATE__

#include "STKLib/labels.h"
#include "STKLib/common.h"
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <malloc.h>
#include <assert.h>
#include <time.h>

using namespace STK;


void usage(char *progname)
{
  char *tchrptr;
  if ((tchrptr = strrchr(progname, '\\')) != NULL) progname = tchrptr+1;
  if ((tchrptr = strrchr(progname, '/')) != NULL) progname = tchrptr+1;
  fprintf(stderr,
"\n%s version " MODULE_VERSION "\n"
"\nUSAGE: %s [options] labelList testMLF\n\n"
" Option                                                    Default\n\n"
"  -o i i i  Substitution/insertion/deletion cost            10 7 7\n"
"  -f        Enable full results                             Off\n"
"  -j f      Keyword score thresholt                         Off\n"
"  -r i      Time-mediated mAlignment (1 2)                   0\n"
"  -t        Output time aligned transcriptions              Off\n"
"  -u f      False alarm time units (hours)                  1.0\n"
"  -v i      ROC table size (time units)                     10\n"
//"  -w[l]     Enable word spotting analysis [lower bound]     Off\n"
"  -w        Enable word spotting analysis. Words that are\n"
"            hits having the same score as some fucking\n"
"            false alarms are taken into account first       Off\n"
"  -words_that_are_fucking_false_alarms_having_the_same_score_as_some_hits_are_taken_into_account_first\n"
"            Enable word spotting analysis. Words that are\n"
"            fucking false alarms having the same score\n"
"            as some hits are taken into account first       Off\n"
"  -A        Print command line arguments                    Off\n"
"  -I mlf    Load master label file mlf\n"
"  -L dir    Set input label (or net) dir                    Current\n"
//"  -S f     Set script file to f                 none\n"
"  -T N      Set trace flags to N                            0\n"
"  -V        Print version information                       Off\n"
"  -X ext    Set input label file ext                        lab\n"
"\n"
" %s is Copyright (C) 2004 Lukas Burget et al. and\n"
" licensed under the GNU General Public License, version 2.\n"
" Bug reports, feedback, etc, to: burget@fit.vutbr.cz\n"
"\n", progname, progname, progname);
  exit(-1);
}

//typedef struct _StrListElem StrListElem;
class StrListElem 
{
public:
  StrListElem * mpNext;
  char          str[1];
};

StrListElem *MLF_files = NULL;
StrListElem *file_name = NULL;
StrListElem **last = &MLF_files;

typedef struct 
{
  unsigned id:8*sizeof(unsigned)-1;
  unsigned is_hit:1;
  FLOAT score;
} t_kwd_tab;

typedef struct {
  int actual;
  int hits;
  int FA;
  float *ROC;

} t_kwd_stats;

int lower_bound = 0;

int main(int argc, char *argv[]) {
  int trace_flag    = 0;
  int print_cmdline = 0;
  int print_version = 0;
  int word_spotting = 0;
  int full_results  = 0;
  int time_alig_out = 0;
  MyHSearchData labelHash = {0};
  char rec_lab_fn[1024];
  char ref_lab_fn[1024];
  char *chrptr;
  char *in_lbl_dir = NULL;
  const char *in_lbl_ext = "lab";
  char *ref_MLF_fn   = NULL;
  FILE *ref_MLF_fp  = NULL;
  char *rec_MLF_fn = NULL;
//  LabelFormat in_lbl_frm = {0};
  STKNetworkOutputFormat in_lbl_frm = {0};

  FILE *rec_MLF_fp;
  float totalTime = 0.0;
  FLOAT score_threshold = LOG_0;
  float FATimeUnit = 1.0;
  int maxFAperHour = 10;
  int nonKeywords  = 0;
  double dval;
  long   lval;
  struct LabelStats ref_stats;
  struct LabelStats rec_stats;
  Label *reclabels, *reflabels, *tlptr;
  int i, kwd_tab_size = 1000, nkwds = 0;
  int stot = 0, ssub = 0, tot = 0, ins = 0, del = 0, sub = 0;

  int start_time_cmp(const void *a, const void *b);
  int score_cmp(const void *a, const void *b);

  t_kwd_tab   *kwd_tab   = NULL;
  t_kwd_stats *kwd_stats = NULL;

  if (argc == 1) usage(argv[0]);

  for (;;) {
    int opt = getopt(argc, argv, "-fj:or:tu:v:w::AI:X:L:T:V");

    if (opt == -1) break;

    switch (opt) {
      case 'A': print_cmdline = 1; break;
      case 'V': print_version = 1; break;
      case 'f': full_results  = 1; break;
      case 't': time_alig_out = 1; break;

      case 'I': ref_MLF_fn = optarg; break;
      case 'L': in_lbl_dir = optarg; break;
      case 'X': in_lbl_ext = optarg; break;

      case 'w': word_spotting = 1;
                if (optarg) {
                  if (!strcmp(optarg, "l") ||
                     !strcmp(optarg, "ords_that_are_fucking_false_alarms_"
                                     "having_the_same_score_as_some_hits_are_"
                                     "taken_into_account_first")) {
                    lower_bound = 1;
                  } else {
                    Error("Invalid option %c%s", opt, optarg);
                  }
                }
                break;


      case  1:  if (labelHash.mTabSize == 0) {
                  labelHash = readLabelList(optarg);
                } else if (!rec_MLF_fn) {
                  rec_MLF_fn = optarg;
                } else {
                  Error("Unexpected parameter %s", optarg);
                }
                break;

      case 'o': if (optind < argc) {
                  substitution_cost = strtol(argv[optind], &chrptr, 0);
                  if (!*chrptr && ++optind < argc) {
                    insertion_cost  = strtol(argv[optind], &chrptr, 0);
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

      case 'j':
      case 'u': dval = strtod(optarg, &chrptr);
               if (!*optarg || *chrptr) {
                  Error("Decimal number is expected after option -%c", opt);
                }
                switch (opt) {
                  case 'j': score_threshold = dval; break;
                  case 'u': FATimeUnit      = dval; break;
                }
                break;

      case 'v': lval = strtol(optarg, &chrptr, 0);
                if (!*optarg || *chrptr) {
                  Error("Integer number is expected after option -%c", opt);
                }
                switch (opt) {
                  case 'v': maxFAperHour = lval; break;
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
    puts("Version: "MODULE_VERSION"\n");
    if (labelHash.mTabSize == 0) exit(0);
  }

  if (!rec_MLF_fn) usage(argv[0]);

  if (print_cmdline) {
    int i;
    for (i=0; i < argc; i++) {
      fputs(argv[i], stdout);
      putchar(' ');
    }
    putchar('\n');
  }

  ref_MLF_fp = OpenInputMLF(ref_MLF_fn);
  rec_MLF_fp = OpenInputMLF(rec_MLF_fn);

  if (word_spotting) {
    kwd_tab   = (t_kwd_tab *)   malloc(sizeof(t_kwd_tab)   * kwd_tab_size);
    kwd_stats = (t_kwd_stats *) malloc(sizeof(t_kwd_stats) * (labelHash.mNEntries+2));
    if (kwd_tab == NULL || kwd_stats == NULL) Error("Insufficient memory");

    for (size_t n = 0; n <= (labelHash.mNEntries+1); n++) {
      kwd_stats[n].FA = kwd_stats[n].actual = kwd_stats[n].hits = 0;
    }
  }

  while (rec_lab_fn[0]='\0', // Rec. MLF is read sequentially record by record
        OpenInputLabelFile(rec_lab_fn, NULL, NULL, rec_MLF_fp, rec_MLF_fn)) {

    reclabels = ReadLabels(rec_MLF_fp, &labelHash, UL_IGNORE, in_lbl_frm,
                           /*sampleRate*/ 1, rec_lab_fn, rec_MLF_fn, &rec_stats);

    strcpy(ref_lab_fn, rec_lab_fn);
    ref_MLF_fp = OpenInputLabelFile(ref_lab_fn, in_lbl_dir, in_lbl_ext,
                                    ref_MLF_fp, ref_MLF_fn);

    reflabels = ReadLabels(ref_MLF_fp, &labelHash, UL_IGNORE, in_lbl_frm,
                           1, ref_lab_fn, ref_MLF_fn, &ref_stats);

    if (!word_spotting) {
      int err = 0;

      AlingTranscriptions(&reclabels, reflabels, NULL);

      for (tlptr=reclabels; tlptr!=NULL; tlptr=tlptr->mpNext) {
        if (tlptr->mpNextLevel->mpName != tlptr->mpName) err = 1;

        if (!tlptr->mpNextLevel->mpName) ins++;
        else {
          if (!tlptr->mpName) del++;
          else if (tlptr->mpNextLevel->mpName != tlptr->mpName) sub++;
          tot++;
        }
      }

      stot++;
      if (err) ssub++;

      if (err && time_alig_out) {
        Label *reclptr;
        int rec;
        printf("Aligned transcription: %s vs %s", ref_lab_fn, rec_lab_fn);
        for (rec=0; rec < 2; rec++) {
          fputs(rec ? "\n REC:" : "\n LAB:", stdout);
          for (reclptr=reclabels; reclptr!=NULL; reclptr=reclptr->mpNext) {
            Label *reflptr = reclptr->mpNextLevel;
            char *reclstr = (char*) (reclptr->mpName ? reclptr->mpName : "");
            char *reflstr = (char*) (reflptr->mpName ? reflptr->mpName : "");
            int lablen = HIGHER_OF(strlen(reclstr), strlen(reflstr));

            printf(" %-*s", lablen, rec ? reclstr: reflstr);
          }
        }
        puts("");
      }
    } else {
      int n = 0;
      Label **ltab;

      nonKeywords += rec_stats.nLabelsTotal - rec_stats.nLabelsRead;
      ltab = (Label **) malloc(sizeof(Label *) * rec_stats.nLabelsRead);
      if (ltab == NULL) Error("Insufficient memory");

      for (n=0, tlptr=reclabels; tlptr!=NULL; tlptr=tlptr->mpNext) {
        if (tlptr->mScore < score_threshold) continue;
        ltab[n++] = tlptr;
      }

      qsort(ltab, n, sizeof(Label *), start_time_cmp);

      if (nkwds + n > kwd_tab_size) {
        kwd_tab_size = HIGHER_OF(kwd_tab_size * 2, nkwds + n);
        kwd_tab = (t_kwd_tab *)realloc(kwd_tab, sizeof(t_kwd_tab)*kwd_tab_size);
        if (kwd_tab == NULL) Error("Insufficient memory");
      }

      tlptr = reflabels;
      for (i = 0; i < n; i++) {
        Label *tlptr2;
        while (tlptr && ((tlptr->mStart+tlptr->mStop)/2) < ltab[i]->mStart) tlptr = tlptr->mpNext;

        kwd_tab[nkwds].id     = reinterpret_cast<size_t>(ltab[i]->mpData) - 1;
        kwd_tab[nkwds].score  = ltab[i]->mScore;
        kwd_tab[nkwds].is_hit = 0;

        for (tlptr2 = tlptr;
            tlptr2 && ((tlptr2->mStart+tlptr2->mStop)/2) <= ltab[i]->mStop;
            tlptr2 = tlptr2->mpNext) {
          if (ltab[i]->mpName == tlptr2->mpName) {
            kwd_tab[nkwds].is_hit = 1;
            break;
          }
        }
        nkwds++;
      }

      for (tlptr = reflabels; tlptr != NULL; tlptr=tlptr->mpNext) {
        kwd_stats[reinterpret_cast<size_t>(tlptr->mpData) - 1].actual++;
      }

      if (ref_stats.endTime > LOG_MIN) {
        totalTime += ref_stats.endTime / (36e9 * FATimeUnit);
      }

      free(ltab);
      ReleaseLabels(reflabels);
    }
    ReleaseLabels(reclabels);
  }

  time_t tm = time(NULL);
  printf(
    "======================= STK Results Analysis ========================\n"
    "  Date: %s", ctime(&tm));
  if (ref_MLF_fn) printf("  Ref: %s\n", ref_MLF_fn);
  printf("  Rec: %s\n", rec_MLF_fn);

  if (!word_spotting) {
    printf(
      "------------------------- Overall Results ---------------------------\n"
      "SENT: %%Correct=%.2f [H=%d, S=%d, N=%d]\n"
      "WORD: %%Corr=%.2f, Acc=%.2f [H=%d, D=%d, S=%d, I=%d, N=%d]",
      100.0 * (stot-ssub)/stot, stot-ssub, ssub, stot,
      100.0 * (tot-del-sub)/tot, 100.0 * (tot-del-ins-sub)/tot,
      tot-del-sub, del, sub, ins, tot);
  } else {
    FLOAT sufficientThreshold = LOG_0;
    
    printf("  Ref. total time: %f hours\n", totalTime  * FATimeUnit);
    if (nonKeywords)printf("  %d non keywords found in test files\n",nonKeywords);

    int nROCbins = (int) (maxFAperHour * totalTime + 0.49999) + 1;
    float lastBinWght =   maxFAperHour * totalTime - nROCbins + 1;
    
    for (size_t n = 0; n <= (labelHash.mNEntries+1); n++) {
      kwd_stats[n].ROC = (float *) malloc(nROCbins * sizeof(float));
      if (kwd_stats[n].ROC == NULL) Error("Insufficient memory");
    }

    qsort(kwd_tab, nkwds, sizeof(kwd_tab[0]), score_cmp);
    
    // Fill ROC tables
    t_kwd_stats *glob_ks = &kwd_stats[labelHash.mNEntries+1];
    for (i = 0; i < nkwds; i++) {
      t_kwd_stats *ks = &kwd_stats[kwd_tab[i].id];
      if (kwd_tab[i].is_hit) {
        ks->hits++;
        glob_ks->hits++;
      } else {
        if (ks->FA < nROCbins) {
          ks->ROC[ks->FA] = ks->actual ? (float) ks->hits / ks->actual : 0.0;
          sufficientThreshold = kwd_tab[i].score;
        }
        if (glob_ks->FA < nROCbins) {
          glob_ks->ROC[glob_ks->FA] = glob_ks->hits;
        }
        ks->FA++;
        glob_ks->FA++;
      }
    }
        
    for (size_t n = 0; n < labelHash.mNEntries; n++) {
      t_kwd_stats *ks = &kwd_stats[n];
      for (i = ks->FA; i < nROCbins; i++) {
        ks->ROC[i] = ks->actual ? (float) ks->hits / ks->actual : 0.0;
        sufficientThreshold = LOG_0;
      }
    }
    
    if (score_threshold > LOG_MIN) {
      if (sufficientThreshold > LOG_MIN) {
        printf("  Sufficient score threshold was: %f\n", sufficientThreshold);
      } else {
        printf("  Score threshold was not sufficient\n");
      }
    }
    
    // Calculate global counts
    for (size_t n = 0; n < labelHash.mNEntries; n++) {
      kwd_stats[labelHash.mNEntries].actual += kwd_stats[n].actual;
      kwd_stats[labelHash.mNEntries].hits   += kwd_stats[n].hits;
      kwd_stats[labelHash.mNEntries].FA   += kwd_stats[n].FA;
    }
    
    for (i = 0; i < nROCbins; i++) {
      kwd_stats[labelHash.mNEntries].ROC[i] = 0.0;
      for (size_t n = 0; n < labelHash.mNEntries; n++) {
        kwd_stats[labelHash.mNEntries].ROC[i] += kwd_stats[n].ROC[i] * kwd_stats[n].actual;
      }
      kwd_stats[labelHash.mNEntries].ROC[i] /= kwd_stats[labelHash.mNEntries].actual;
    }

    float hits=0;
    for (i = 0; i < nROCbins; i++) {
      hits=kwd_stats[labelHash.mNEntries+1].ROC[i];
      kwd_stats[labelHash.mNEntries+1].ROC[i] = kwd_stats[labelHash.mNEntries].actual ? (float) hits / kwd_stats[labelHash.mNEntries].actual : 0.0;
    }
    for (i = kwd_stats[labelHash.mNEntries+1].FA; i < nROCbins; i++) {
      kwd_stats[labelHash.mNEntries+1].ROC[i] = kwd_stats[labelHash.mNEntries].actual ? (float) kwd_stats[labelHash.mNEntries].hits / kwd_stats[labelHash.mNEntries].actual : 0.0;
    }
    kwd_stats[labelHash.mNEntries+1].actual=kwd_stats[labelHash.mNEntries].actual;
    
                
    if (full_results) {
      fputs(
      "------------------------- ROC Information ---------------------------"
      "\n\n      KeyWord:", stdout);
      for (i = 0; i <= maxFAperHour; i++) printf(" %4d", i);
      fputs("\n", stdout);

      for (size_t n = 0; n <= (labelHash.mNEntries+1); n++) {
        t_kwd_stats *ks = &kwd_stats[n];
        printf("%13s:", n == labelHash.mNEntries ? "Overall" : n == (labelHash.mNEntries+1) ? "Global" : labelHash.mpEntry[n]->key);
        for (i = 0; i <= maxFAperHour; i++) {
          float y1, y2, hitr;
          float a = (i * totalTime + 0.5);
          int   b = (int) a;

          assert(b < nROCbins);

          a -= b;
          y1 = b > 0        ? ks->ROC[b-1] :
              nROCbins > 1 ? ks->ROC[0] * 2.0 - ks->ROC[1] : ks->ROC[0];
          y2 = ks->ROC[b];
          hitr = (y2 * a + y1 * (1.0 - a)) * 100.0;
          hitr = hitr > 0.0 ? hitr : 0.0;
          printf(hitr != 100.0 ? " %4.1f" : " 100.", hitr);
        }
        puts("");
      }
    }
    puts(
    "------------------------- Figures of Merit --------------------------"
    "\n\n      KeyWord:    #Hits     #FAs  #Actual      FOM");

    for (size_t n = 0; n <= (labelHash.mNEntries+1); n++) {
      t_kwd_stats *ks = &kwd_stats[n];
      float FOM = 0.0;

      for (i = 0; i < nROCbins-1; i++) FOM += ks->ROC[i];
      FOM += ks->ROC[nROCbins-1] * lastBinWght;
      FOM *=  100.0 / (maxFAperHour * totalTime);

      printf("%13s:", n == labelHash.mNEntries ? "Overall" : n == (labelHash.mNEntries+1) ? "Global" : labelHash.mpEntry[n]->key);
      printf(" %8d %8d %8d %8.2f\n", ks->hits, ks->FA, ks->actual, FOM);
    }
  }
  puts(
  "\n=====================================================================");

  return 0;
}

int start_time_cmp(const void *a, const void *b)
{
  return (*(Label**) a)->mStart - (*(Label**) b)->mStart;
}

int score_cmp(const void *a, const void *b)
{
  t_kwd_tab *kta = (t_kwd_tab *) a;
  t_kwd_tab *ktb = (t_kwd_tab *) b;

  return kta->score == ktb->score
         ? (lower_bound
            ? kta->is_hit - ktb->is_hit  // For lower_bound prefer false alarms
            : ktb->is_hit - kta->is_hit) // otherwise prefer hits
         : kta->score  < ktb->score
           ?  1 : -1;
}
