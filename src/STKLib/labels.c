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

#include "viterbi.h"
#include "labels.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

const Label init_label = {
  NULL, UNDEF_TIME, UNDEF_TIME, NULL, NULL, -1, 0, NULL
};

const char *transc_filter  = NULL;
const char *transc_ofilter = NULL;

int IsMLF(FILE *MLFfp)
{
  char line[] = "#!MLF!#";
  int i, c;
  while(isspace(c = getc(MLFfp)));

  for(i=0; (c == line[i] ||
           (c == 'N' && i == 3) // !!!To be removed, #!MNF!# will be no more supported
           ) && i < sizeof(line)-1; i++) {
    c = getc(MLFfp);
  }

  if(i < sizeof(line)-1) {
   ungetc(c, MLFfp);
   while(--i >= 0) ungetc(line[i], MLFfp);
   return 0;
  }

  while(c != EOF && c != '\n') c = getc(MLFfp);
  return 1;
}


FILE *OpenOutputMLF(const char *out_MLF)
{
  FILE *lfp = NULL;

  if(out_MLF) {
    if((lfp = my_fopen(out_MLF, "wt", transc_ofilter)) == NULL) {
      Error("Cannot open output MLF file %s", out_MLF);
    }
/*    if(!strcmp(out_MLF, "-")) {
      lfp = stdout;
    } else if(transc_ofilter) {
      char *filter = expandFilterCommand(transc_ofilter, out_MLF);

      if((lfp = popen(filter, "w")) == NULL) {
        Error("Cannot popen output filter %s", filter);
      }
      free(filter);
    } else if((lfp = fopen(out_MLF, "wt")) == NULL) {
      Error("Cannot open output MLF file %s", out_MLF);
    }*/
    if(fputs("#!MLF!#\n", lfp) == EOF) {
      Error("Cannot write to output MLF file %s", out_MLF);
    }
  }
  return lfp;
}

FILE *OpenOutputLabelFile(char *file_name, const char *out_lbl_dir, const char *out_lbl_ext,
                          FILE *MLFfp, const char *out_MLF)
{
  char label_file[1024];
  FILE *lfp;

  if(!out_MLF && !strcmp(file_name, "-")) {
    return stdout;
  }

  MakeFileName(label_file, file_name, out_lbl_dir, out_lbl_ext);
  strcpy(file_name, label_file);

  if(out_MLF) {
    assert(MLFfp);
    lfp = MLFfp;
    if(fprintf(lfp, "\"%s\"\n", label_file) < 0) {
      Error("Cannot write to output MLF file %s", out_MLF);
    }
  } else if((lfp = my_fopen(label_file, "wt", transc_ofilter)) == NULL) {
    Error("Cannot open output label file %s", label_file);
  }
  /*if(transc_ofilter) {
    char *filter = expandFilterCommand(transc_ofilter, label_file);

    if((lfp = popen(filter, "w")) == NULL) {
      Error("Cannot popen output filter %s", filter);
    }
    free(filter);
  } else if((lfp = fopen(label_file, "wt")) == NULL) {
    Error("Cannot open output label file %s", label_file);
  }*/

  strcpy(file_name, label_file);
  return lfp;
}

void WriteLabels(FILE* lfp, Label *label, LabelFormat labelFormat, long smpPrd,
                 const char *label_file, const char *out_MLF)
{
  Label *level, *prev = NULL;
  int ctm = labelFormat.CENTRE_TM;

  for(;label != NULL; label = label->next) {
    if(!(labelFormat.TIMES_OFF)) {
      long long time = UNDEF_TIME;
      for(level = label; level != NULL; level = level->nextLevel) {
        if(time == UNDEF_TIME || time < level->start) time = level->start;
      }
      if(time != UNDEF_TIME) {
        fprintf(lfp, "%08lld", (long long)
                smpPrd * (2 * time + ctm) / 2 - labelFormat.left_extent);

        time = UNDEF_TIME;
        for(level = label; level != NULL; level = level->nextLevel) {
          if(time == UNDEF_TIME || time > level->stop) time = level->stop;
        }
        if(time != UNDEF_TIME) {
          fprintf(lfp," %08lld", (long long)
                  smpPrd * (2 * time - ctm) / 2 + labelFormat.right_extent);
        }
      }
    }

    for(level = label; level != NULL; level = level->nextLevel) {
      if(prev == level) break;


      putc(' ', lfp);
      fprintHTKstr(lfp, level->mpName);

      if(level->id >= 0) {
        fprintf(lfp, "[%d]", level->id + 2);
      }

      if(!(labelFormat.SCORE_OFF) && level->score != 0.0) {
        // !!! only acoustic scores should be normalized !!!
        fprintf(lfp, " %g", labelFormat.SCORE_NRM
                            ? level->score / (level->stop - level->start)
                            : level->score);
      }
      if(prev) prev = prev->nextLevel;
    }
    fputc('\n', lfp);
    if(ferror(lfp)) {
      Error("Cannot write to output label file %s", out_MLF ? out_MLF : label_file);
    }
    prev = label;
  }
}

void CloseOutputLabelFile(FILE *lfp, const char *out_MLF)
{
  if(out_MLF) {
    if(fputs(".\n", lfp) == EOF) {
      Error("Cannot write to output MLF file %s", out_MLF);
    }
  } else {
    my_fclose(lfp);
//    if(lfp != stdout) fclose(lfp);
  }
}

struct my_hsearch_data readLabelList(char *labelListFileName)
{
  char line[1024];
  FILE *fp;
  int nlabels=0;
  struct my_hsearch_data hash = {0};

  if((fp = fopen(labelListFileName, "rt")) == NULL) {
      Error("Cannot open file: '%s'", labelListFileName);
  }

  if(!my_hcreate_r(100, &hash)) {
     Error("Insufficient memory");
  }

  while(fscanf(fp, "%1024s", line) == 1) {
    ENTRY e, *ep;

    e.key  = line;
    e.data = NULL;
    my_hsearch_r(e, FIND, &ep, &hash);

    if(ep != NULL)  continue;

    e.key  = strdup(line);
    e.data = (void *) ++nlabels;

    if(e.key == NULL || !my_hsearch_r(e, ENTER, &ep, &hash)) {
      Error("Insufficient memory");
    }
  }

  fclose(fp);

  return hash;
}

void ReleaseLabels(Label *labels)
{
  Label *tlptr, *nlptr;

  while(labels) {
    nlptr = labels;
    labels = labels->nextLevel;
    while(nlptr) {
      tlptr = nlptr;
      nlptr=nlptr->next;
      free(tlptr);
    }
  }
}

FILE *OpenInputMLF(const char *in_MLF)
{
  FILE *lfp = NULL;
  char line[1024];
  int i;

  if(in_MLF) {
    if((lfp = my_fopen(in_MLF, "rt", transc_filter)) == NULL) {
      Error("Cannot open input MLF file %s", in_MLF);
    }
/*    if(!strcmp(in_MLF, "-")) {
      lfp = stdin;
    } else if(transc_filter) {
      char *filter = expandFilterCommand(transc_filter, in_MLF);

      if((lfp = popen(filter, "r")) == NULL) {
        Error("Cannot popen input filter %s", filter);
      }
      free(filter);
    } else if((lfp = fopen(in_MLF, "rt")) == NULL) {
      Error("Cannot open input MLF file %s", in_MLF);
    }*/

    if(fgets(line, sizeof(line), lfp) == NULL) {
      Error("Invalid input MLF file %s", in_MLF);
    }

    for(i = strlen(line); i > 0 && isspace(line[i-1]); i--);
    line[i] = '\0';

    if(strcmp(line, "#!MLF!#")
       && strcmp(line, "#!MNF!#") // !!!To be removed, #!MNF!# will be no more supported
       ) {
      Error("Invalid input MLF file %s", in_MLF);
    }
  }
  return lfp;
}

void CloseInputMLF(FILE *lfp)
{
  my_fclose(lfp);
/*  if(lfp != stdout) {
    if(transc_filter) pclose(lfp);
    else                 fclose(lfp);
  }*/
}

FILE *OpenInputLabelFile(
  char *file_name,
  const char *in_lbl_dir,
  const char *in_lbl_ext,
  FILE *MLFfp, const char *in_MLF)
{
  char line[1024];
  FILE *lfp;
  int fn_len;
  int rewinded = 0;

  if(*file_name != '\0') {
    MakeFileName(line, file_name, in_lbl_dir, in_lbl_ext);
    strcpy(file_name, line);
  }

  if(in_MLF) {
    char *chptr, *endptr;

    assert(MLFfp);

    for(;;) { // if file_name doesn't match to name in MLF, try next record
      if(fgets(line, sizeof(line), MLFfp) == NULL) {
        if(feof(MLFfp)) {
          if(*file_name == '\0') return NULL; // No more entries for sequential reading
          if(!rewinded) { // Try again from the begining
            rewinded = 1;
            rewind(MLFfp);
            fgets(line, sizeof(line), MLFfp);
            continue;
          }

          Error("Label/network file %s not found in MLF file %s", file_name, in_MLF);
        } else {
          Error("Cannot read MLF file %s", in_MLF);
        }
      }

      // Trim initial and final spaces and quotas
      for(chptr = line;
          *chptr && isspace(*chptr);
          chptr++);

      for(endptr = line + strlen(line);
          endptr > chptr && isspace(*(endptr-1));
          endptr--);

      if(*chptr == '"' && endptr > chptr+1 && *(endptr-1) == '"') {
        chptr++;
        endptr--;
      }

      *endptr = '\0';

      if(*chptr == '\0') {
        Error("Label/network file name expected in MLF file %s",
              in_MLF);
      }

      if(*file_name == '\0') { // Return next record in MLF
        assert(in_MLF != NULL);
        strcpy(file_name, chptr);
        return MLFfp;
      }

      if(*chptr == '*') {
        if(*++chptr == '/') ++chptr;
        fn_len = strlen(file_name);
        if(endptr - chptr <= fn_len && /* !!! pattern *\foo will accept also /dir/foofoo */
          !strcmp(chptr, file_name + fn_len - (endptr - chptr))) {
          return MLFfp;
        }
      } else if(!strcmp(chptr, file_name)) {
        return MLFfp;
      }

      // file_name doesn't match to name in MLF, try next record
      do { // read everything till line with "."
        if(fgets(line, sizeof(line), MLFfp) == NULL) {
          Error("Invalid MLF file %s (%s)", in_MLF, line);
        }
      } while((chptr = strtok(line, " \n\t")) == NULL || strcmp(chptr, "."));
    }
  }

  // Opening ordinary label file (not record in MLF)
  if((lfp = my_fopen(file_name, "rt", transc_filter)) == NULL) {
    Error("Cannot open input file %s", file_name);
  }

/*  if(!strcmp(file_name, "-")) {
    lfp = stdin;
  } else if(transc_filter) {
    char *filter = expandFilterCommand(transc_filter, file_name);
    if((lfp = popen(filter, "r")) == NULL) {
      Error("Cannot popen input filter %s", filter);
    }
    free(filter);
  } else if((lfp = fopen(file_name, "rt")) == NULL) {
    Error("Cannot open input file %s", file_name);
  }*/

  return lfp;
}

Label *ReadLabels(
  FILE *lfp,
  struct my_hsearch_data *label_hash,
  enum UnknownLabelsAction unknownLabels,
  LabelFormat labelFormat,
  long sampPeriod,
  const char *file_name,
  const char *in_MLF,
  struct LabelStats *stats)
{
  char line[1024];
  Label *labels = NULL;
  ENTRY e, *ep;

  char *chptr, *endptr;
  Label *prev, *current = NULL;
  long long llv;
  double dv;

  if(stats != NULL) {
    stats->endTime      = UNDEF_TIME;
    stats->nLabelsTotal = 0;
    stats->nLabelsRead  = 0;
  }

  for(;;) {
    do {
      if(fgets(line, sizeof(line), lfp) == NULL) {
//        if(in_MLF) Error("Invalid label file %s (%s)", file_name, line);
        return labels;
      }
      chptr = line + strspn(line, " \n\t");
    } while(!*chptr);

    if(chptr[0] == '.' && (chptr[1] == '\0' || isspace(chptr[1]))) {
      return labels;
    }

    prev = current;

    current = (Label *) malloc(sizeof(Label));
    if(!current) Error("Insufficient memory");
    *current = init_label;

    llv = strtoull(chptr, &endptr, 10);
    if(endptr != chptr) {
      long center_shift = labelFormat.CENTRE_TM ? sampPeriod / 2 : 0;

      if(!(labelFormat.TIMES_OFF)) {
        current->start = (llv - center_shift - labelFormat.left_extent) / sampPeriod;
      }
      chptr = endptr;

      llv = strtoull(chptr, &endptr, 10);
      if(endptr != chptr) {
        if(!(labelFormat.TIMES_OFF)) {
          current->stop = (llv + center_shift + labelFormat.right_extent) / sampPeriod;
        }
        chptr = endptr;

        if(stats != NULL) {
          stats->endTime = HIGHER_OF(stats->endTime, current->stop);
        }
      }
    }
    if(getHTKstr(e.key = chptr, &chptr)) {
      Error("%s in label file %s", chptr, file_name);
    }

    my_hsearch_r(e, FIND, &ep, label_hash);

    if(stats != NULL) stats->nLabelsTotal++;
    if(ep == NULL) {
      if(unknownLabels == UL_ERROR) {
        Error("Unknown label '%s' in label file %s", e.key, file_name);
      }
      if(unknownLabels == UL_WARN) {
        Warning("Ignored unknown label '%s' in label file %s", e.key, file_name);
      }
      if(unknownLabels == UL_IGNORE) {
        free(current);
        current = prev;
        continue;
      }

      if(unknownLabels == UL_READ) {
      current->data = current->mpName = NULL;
      } else if(unknownLabels == UL_INSERT) {
  // ??? naco e.key=strdup(e.key) => bordel v pameti
  // current->data = current->mpName = e.data = e.key = strdup(e.key);
  e.data = strdup(e.key);
  current->mpName = (char*) e.data;
  current->data = current->mpName;

        if(e.key == NULL || !my_hsearch_r(e, ENTER, &ep, label_hash)) {
          Error("Insufficient memory");
        }
      } else Error("Fatal: Invalid UnknownLabelsAction value");
    } else {
      current->data = ep->data;
      current->mpName = ep->key;
    }

    if(stats != NULL) stats->nLabelsRead++;
    dv = strtod(chptr, &endptr);
    if(endptr != chptr) {
      current->score = dv;
      chptr = chptr;
    }

    if(prev) {
      prev->next = current;
    } else {
      labels = current;
    }
  }
  return labels;
}

void CloseInputLabelFile(FILE *lfp, const char *out_MLF)
{
  if(out_MLF == NULL) {
    my_fclose(lfp);
/*    if(lfp != stdout) {
      if(transc_filter) pclose(lfp);
      else                 fclose(lfp);
    }*/
  }
}



int   substitution_cost = 10;
int   insertion_cost    =  7;
int   deletion_cost     =  7;
float null_cost         =  0.0;
int   time_med_alg      =  0;

/*
Moves in cost table and associated cost values

 direction of hypotesis
 .
/|\ I    NULL->I
 |  *  *!NULL->correct->0
    | /       !correct->S
    |/
     --* NULL->0
        !NULL->D

   --> direction of (multilevel) transcription

I/D/S   - cost of insertion/deletion/substitution
NULL    - null word in transcription
correct - hypotesis matches to transcription
*/

#define max(x,y) ((x) > (y) ? (x) : (y))

void AlingTranscriptions(Label **aligned_hyps, Label *new_hyp, long *weights)
{
  Label *alg_hyp = *aligned_hyps;
  Label *tlptr, *hlptr, *llptr, *tmplptr, *prevlptr;
  int tlen, hlen, i, j, k = 0;
  unsigned long *cost, ulv;
  char *path, chv = 0;
  
  if(alg_hyp == NULL) { //Hypoteses to which new_hyp is added must not be empty
    alg_hyp = (Label *) malloc(sizeof(Label));
    if(!alg_hyp) Error("Insufficient memory");
    *alg_hyp = init_label;
    alg_hyp->score = null_cost;
  }

  for(tlptr = alg_hyp, tlen = 0; tlptr != NULL; tlptr=tlptr->next) tlen++;
  for(hlptr = new_hyp, hlen = 0; hlptr != NULL; hlptr=hlptr->next) hlen++;

  cost = (unsigned long *) calloc((tlen+1) * (hlen+1), sizeof(unsigned long));
  if(!cost) Error("Insufficient memory");

  path = (char *) malloc((tlen+1) * (hlen+1) * sizeof(char));
  if(!path) Error("Insufficient memory");

  cost[0] = 0;

  ulv = 0;
  for(tlptr = alg_hyp, i=0; tlptr != NULL; tlptr=tlptr->next, i++) {
    for(llptr = tlptr, k=0; llptr !=NULL; llptr = llptr->nextLevel, k++) {
      long w = weights ? weights[k] : 1;
      if(time_med_alg == 0) {
        ulv += (!llptr->mpName ? 0 : deletion_cost) * w;
      } else if(time_med_alg == 1) {
        ulv += (!llptr->mpName ? 1 : llptr->stop-llptr->start) * w;
      } else {
        ulv += (!llptr->mpName ? 1 : max(new_hyp
                                       ? llptr->stop-new_hyp->start
                                       : 0, 0)) * w;
      }
#ifdef KEEP_ORDER
      if(llptr->mpName && new_hyp &&
         llptr->start >= new_hyp->stop) ulv = LONG_MAX;
//      if(llptr->mpName && new_hyp &&
//         llptr->start+llptr->stop >= new_hyp->start+new_hyp->stop) ulv = LONG_MAX;
#endif
    }
    cost[i+1] = ulv;
    path[i+1] = 0;
  }

  ulv = 0;
//  for(i=0, ulv=0; i<hlen; i++) {
  for(hlptr = new_hyp, i=0; hlptr != NULL; hlptr=hlptr->next, i++) {
    for(llptr = alg_hyp, k=0; llptr !=NULL; llptr = llptr->nextLevel, k++) {
//    for(j=0; j < k; j++) { // k = nuber of alg_hyp levels (see cycle above)
      long w = weights ? weights[k] : 1;
      if(time_med_alg == 0) {
        ulv += insertion_cost * w;
      } else if(time_med_alg == 1) {
        ulv += (hlptr->stop-hlptr->start) * w;
      } else {
        ulv += max(hlptr->stop-llptr->start, 0) * w;
      }
#ifdef KEEP_ORDER
      if(llptr->mpName && hlptr->start >= llptr->stop) ulv = LONG_MAX;
//        if(llptr->mpName && hlptr->start+hlptr->stop >= llptr->start+llptr->stop) ulv = LONG_MAX;
#endif
    }
    cost[(i+1) * (tlen+1)] = ulv;
    path[(i+1) * (tlen+1)] = 1;
  }

  for(tlptr = alg_hyp, i=0; tlptr != NULL; tlptr=tlptr->next, i++) {
    for(hlptr =  new_hyp, j=0; hlptr != NULL; hlptr=hlptr->next, j++) {
      unsigned long tcost, hcost, thcost;

      tcost  = cost[(j+1) * (tlen+1) + i];
      hcost  = cost[ j    * (tlen+1) + i + 1];
      thcost = cost[ j    * (tlen+1) + i];

      for(llptr = tlptr, k=0; llptr !=NULL; llptr = llptr->nextLevel, k++) {
        long w = weights ? weights[k] : 1;
        if(time_med_alg == 0) {
          tcost  += (!llptr->mpName ? 0 : deletion_cost)            * w;
          hcost  += insertion_cost                                * w;
          thcost += (!llptr->mpName 
                     ? insertion_cost
                     : llptr->mpName != hlptr->mpName
                       ? substitution_cost : 0) * w;
        } else if(time_med_alg == 1) {
          tcost  += (!llptr->mpName ? 1 : llptr->stop-llptr->start) * w;
          hcost  += (hlptr->stop-hlptr->start)                    * w;
          thcost += (!llptr->mpName
                    ? hlptr->stop-hlptr->start
                    : labs(llptr->start-hlptr->start) +
                      labs(llptr->stop-hlptr->stop)   +
                      (llptr->mpName != hlptr->mpName ? 1 : 0))       * w;
        } else {
          long hroverlap = max(llptr->next && llptr->next->start != UNDEF_TIME
                               ? hlptr->stop - llptr->next->start : 0, 0);
          long troverlap = max(hlptr->next
                               ? llptr->stop-hlptr->next->start : 0, 0);

          long hloverlap = max(llptr->stop != UNDEF_TIME
                               ? llptr->stop-hlptr->start : 0, 0);
          long tloverlap = max(hlptr->stop-llptr->start, 0);

          tcost  += (!llptr->mpName ? 1 : tloverlap + troverlap)    * w;
          hcost  += (hloverlap + hroverlap)                       * w;
          thcost += (!llptr->mpName
                    ? hloverlap + hroverlap
                    : troverlap + hroverlap +
                      (labs(llptr->start - hlptr->start) +
                       labs(llptr->stop  - hlptr->stop)) +
                      (llptr->mpName != hlptr->mpName ? 1 : 0))              * w;
        }
#ifdef KEEP_ORDER
        if(llptr->mpName) {
          if(hlptr->start >= llptr->stop) tcost = LONG_MAX;
          if(llptr->start >= hlptr->stop) hcost = LONG_MAX;
        }
//        if(llptr->mpName) {
//          if(hlptr->start+hlptr->stop >= llptr->start+llptr->stop) tcost = LONG_MAX;
//          if(llptr->start+llptr->stop >= hlptr->start+hlptr->stop) hcost = LONG_MAX;
//        }
#endif
      }

      if(thcost <= hcost && thcost <= tcost) {
        ulv = thcost;
        chv = 2;
      } else if(hcost <= tcost && hcost <= thcost) {
        ulv = hcost;
        chv = 1;
      } else {
        ulv = tcost;
        chv = 0;
      }

      assert(ulv >= 0);

      cost[(j+1) * (tlen+1) + i + 1] = ulv;
      path[(j+1) * (tlen+1) + i + 1] = chv;
    }
  }

  free(cost);

  i = tlen;
  j = hlen;

  for(;;) {
    char tch = path[j * (tlen+1) + i];
    path[j * (tlen+1) + i] = chv;
    chv = tch;

    if(i == 0 && j == 0) break;

    switch(chv) {
      case  0: i--; break;
      case  1: j--; break;
      default: i--; j--;
    }
  }

  switch(path[0]) {
    case  1:
      tlptr = (Label *) malloc(sizeof(Label));
      if(!tlptr) Error("Insufficient memory");
      *tlptr = init_label;
      tlptr->next = alg_hyp;
      if(alg_hyp) tlptr->start = alg_hyp->start;
      tlptr->score = null_cost;
      alg_hyp = prevlptr = tlptr;

      for(llptr = tlptr->next->nextLevel; llptr !=NULL; llptr = llptr->nextLevel) {
        tmplptr = (Label *) malloc(sizeof(Label));
        if(!tmplptr) Error("Insufficient memory");
        *tmplptr = init_label;
        tmplptr->next = llptr;
        tmplptr->start = llptr->start;
        tmplptr->score = null_cost;
        prevlptr->nextLevel = tmplptr;
        prevlptr = tmplptr;
      }

      prevlptr->nextLevel = NULL;
      hlptr = new_hyp;
      i = 0; j = 1;
      break;
    case  0:
      hlptr = (Label *) malloc(sizeof(Label));
      if(!hlptr) Error("Insufficient memory");
      *hlptr = init_label;
      hlptr->next = new_hyp;
      if(new_hyp) hlptr->start = new_hyp->start;
      hlptr->score = null_cost;
      new_hyp = hlptr;
      tlptr = alg_hyp;
      i = 1; j = 0;
      break;
    default:
      tlptr = alg_hyp;
      hlptr =  new_hyp;
      i = j = 1;
  }

  while(i != tlen || j != hlen) {
    chv = path[j * (tlen+1) + i];
    switch(chv) {
      case  1:
        prevlptr = NULL;

        for(llptr = tlptr; llptr !=NULL; llptr = llptr->nextLevel) {
          tmplptr = (Label *) malloc(sizeof(Label));
          if(!tmplptr) Error("Insufficient memory");
          *tmplptr = init_label;
          tmplptr->next = llptr->next;
          llptr->next = tmplptr;
          tmplptr->stop = llptr->stop;
          if(tmplptr->next) tmplptr->start = tmplptr->next->start;
          tmplptr->score = null_cost;
          if(prevlptr) prevlptr->nextLevel = tmplptr;
          prevlptr = tmplptr;
        }

        prevlptr->nextLevel = NULL;
        j++;
        break;
      case  0:
        tmplptr = (Label *) malloc(sizeof(Label));
        if(!tmplptr) Error("Insufficient memory");
        *tmplptr = init_label;
        tmplptr->next = hlptr->next;
        hlptr->next = tmplptr;
        tmplptr->stop = hlptr->stop;
        if(tmplptr->next) tmplptr->start = tmplptr->next->start;
        tmplptr->score = null_cost;
        i++;
        break;
      default:
        i++; j++;
    }
    tlptr = tlptr->next;
    hlptr = hlptr->next;
  }


  hlptr =  new_hyp;
  for(tlptr = alg_hyp; tlptr->nextLevel !=NULL; tlptr = tlptr->nextLevel);
  for(; tlptr !=NULL; tlptr = tlptr->next, hlptr = hlptr->next) {
    tlptr->nextLevel = hlptr;
  }
  
  *aligned_hyps = alg_hyp;
}
