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

#include "common.h"
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <assert.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>



void _Error_(const char *func, const char *file, int line, char *msg, ...) {
  va_list ap;
  va_start(ap, msg);
  fflush(stdout);
  fprintf(stderr, "ERROR (%s:%d) %s: ", file, line, func);
  vfprintf(stderr, msg, ap);
  fputc('\n', stderr);
  va_end(ap);
  exit(1);
}

void _Warning_(const char *func, const char *file, int line, char *msg, ...) {
  va_list ap;
  va_start(ap, msg);
  fflush(stdout);
  fprintf(stderr, "WARNING (%s:%d) %s: ", file, line, func);
  vfprintf(stderr, msg, ap);
  fputc('\n', stderr);
  va_end(ap);
}

void TraceLog(char *msg, ...) {
  va_list ap;
  va_start(ap, msg);
  fflush(stderr);
  vfprintf(stdout, msg, ap);
  fputc('\n', stdout);
  va_end(ap);
}


static char *parmKindNames[] = {
 "WAVEFORM",
 "LPC",
 "LPREFC",
 "LPCEPSTRA",
 "LPDELCEP",
 "IREFC",
 "MFCC",
 "FBANK",
 "MELSPEC",
 "USER",
 "DISCRETE",
 "PLP",
 "ANON"};


int ReadParmKind(const char *str, BOOL checkBrackets)
{
  int  i, parmKind=0;
  int  slen = strlen(str);

  if(checkBrackets) {
    if(str[0] != '<' || str[slen-1] != '>') {
      return -1;
    }
    str++; slen -= 2;
  }
  for(; slen >= 0 && str[slen-2] == '_'; slen -= 2) {
    parmKind |= str[slen-1] == 'E' ? PARAMKIND_E :
                str[slen-1] == 'N' ? PARAMKIND_N :
                str[slen-1] == 'D' ? PARAMKIND_D :
                str[slen-1] == 'A' ? PARAMKIND_A :
                str[slen-1] == 'C' ? PARAMKIND_C :
                str[slen-1] == 'Z' ? PARAMKIND_Z :
                str[slen-1] == 'K' ? PARAMKIND_K :
                str[slen-1] == '0' ? PARAMKIND_0 :
                str[slen-1] == 'V' ? PARAMKIND_V :
                str[slen-1] == 'T' ? PARAMKIND_T : -1;

    if(parmKind == -1) return -1;
  }
  for(i = 0; i < sizeof(parmKindNames) / sizeof(char*); i++) {
    if (!strncmp(str, parmKindNames[i], slen)) {
      return parmKind | i;
    }
  }
  return -1;
}

void MakeFileName(char *outFileName, const char* inFileName,
                  const char *out_dir, const char *out_ext)
{
  const char *base_name, *bname_end = NULL, *chrptr;

//  if(*inFileName == '*' && *++inFileName == '/') ++inFileName;

  base_name = strrchr(inFileName, '/');
  base_name = base_name != NULL ? base_name + 1 : inFileName;
  if(out_ext) bname_end = strrchr(base_name, '.');
  if(!bname_end) bname_end = base_name + strlen(base_name);


  if((chrptr = strstr(inFileName, "/./")) != NULL) {
    // what is in path after /./ serve as base name
    base_name = chrptr + 3;
  }/* else if(*inFileName != '/') {
    // if inFileName isn't absolut path, don't forget directory structure
    base_name = inFileName;
  }*/

  *outFileName = '\0';
  if(out_dir) {
    if(*out_dir) {
      strcat(outFileName, out_dir);
      strcat(outFileName, "/");
    }
    strncat(outFileName, base_name, bname_end-base_name);
  } else {
    strncat(outFileName, inFileName, bname_end-inFileName);
  }

  if(out_ext && *out_ext) {
    strcat(outFileName, ".");
    strcat(outFileName, out_ext);
  }
}

char *strtoupper(char *str)
{
   char *chptr;
   for(chptr = str; *chptr; chptr++) {
     *chptr = toupper(*chptr);
   }
   return str;
}

int qsstrcmp(const void *a, const void *b)
{
  return strcmp(*(char **)a, *(char **)b);
}

int qsptrcmp(const void *a, const void *b)
{
  return *(char **)a - *(char **)b;
}


FLOAT minLogDiff;

void InitLogMath()
{
  minLogDiff = log(DBL_EPSILON);
}

double LogAdd(double x, double y)
{
   double diff;

   if (x < y) {
      diff = x - y;
      x = y;
   } else {
     diff = y - x;
   }

   if(x < LOG_MIN) {
     return LOG_0;
   }

   if(diff < minLogDiff) {
      return  x;
   }

   return x + log(1.0 + exp(diff));
}

double LogSub(double x, double y)
{
   FLOAT diff = y - x;
   assert(diff <= 0.0);

   if(x < LOG_MIN || diff > -EPSILON) { //'diff > -EPSILON' means 'diff == 0.0'
     return LOG_0;
   }

   if(diff < minLogDiff) {
      return  x;
   }

   return x + log(1.0 - exp(diff));
}

FloatInLog FIL_Add(FloatInLog a, FloatInLog b)
{
  FloatInLog ret;
  if(a.negative == b.negative) {
    ret.negative = a.negative;
    ret.logvalue = LogAdd(a.logvalue, b.logvalue);
  } else if(a.logvalue > b.logvalue) {
    ret.negative = a.negative;
    ret.logvalue = LogSub(a.logvalue, b.logvalue);
  } else {
    ret.negative = b.negative;
    ret.logvalue = LogSub(b.logvalue, a.logvalue);
  }
  return ret;
}

FloatInLog FIL_Sub(FloatInLog a, FloatInLog b)
{
  b.negative = ~b.negative;
  return FIL_Add(a, b);
}

FloatInLog FIL_Mul(FloatInLog a, FloatInLog b)
{
  FloatInLog ret;
  ret.negative = a.negative ^ b.negative;
  ret.logvalue = a.logvalue + b.logvalue;
  return ret;
}

FloatInLog FIL_Div(FloatInLog a, FloatInLog b)
{
  b.logvalue = -b.logvalue;
  return FIL_Mul(a, b);
}

void sigmoid_vec(FLOAT *in, FLOAT *out, int size)
{
  while(size--) *out++ = 1.0/(1.0 + _EXP(-*in++));
}

void exp_vec(FLOAT *in, FLOAT *out, int size)
{
  while(size--) *out++ = _EXP(*in++);
}

void log_vec(FLOAT *in, FLOAT *out, int size)
{
  while(size--) *out++ = _LOG(*in++);
}

void sqrt_vec(FLOAT *in, FLOAT *out, int size)
{
  while(size--) *out++ = _SQRT(*in++);
}


void softmax_vec(FLOAT *in, FLOAT *out, int size)
{
  int i;
  FLOAT sum = 0.0;
  for(i = 0; i < size; i++) sum += out[i] = _EXP(in[i]);
  sum = 1.0 / sum;
  for(i = 0; i < size; i++) out[i] *= sum;
}

int ParmKind2Str(int parmKind, char *outstr) {
    if((parmKind & 0x003F) >= sizeof(parmKindNames)/sizeof(parmKindNames[0])) {
      return 0;
    }

    strcpy (outstr, parmKindNames[parmKind & 0x003F]);

    if (parmKind & PARAMKIND_E) strcat (outstr, "_E");
    if (parmKind & PARAMKIND_N) strcat (outstr, "_N");
    if (parmKind & PARAMKIND_D) strcat (outstr, "_D");
    if (parmKind & PARAMKIND_A) strcat (outstr, "_A");
    if (parmKind & PARAMKIND_C) strcat (outstr, "_C");
    if (parmKind & PARAMKIND_Z) strcat (outstr, "_Z");
    if (parmKind & PARAMKIND_K) strcat (outstr, "_K");
    if (parmKind & PARAMKIND_0) strcat (outstr, "_0");
    if (parmKind & PARAMKIND_V) strcat (outstr, "_V");
    if (parmKind & PARAMKIND_T) strcat (outstr, "_T");
    return 1;
}

/*
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))
#define M_SQRT3 1.7320508075688771931766041234368

int cubic_roots(FLOAT  a,  FLOAT  b, FLOAT  c,
                FLOAT *r,  FLOAT *s, FLOAT *t)
{
  FLOAT Q   = (a*a - 3*b) / 9;
  FLOAT R   = (2*a*a*a - 9*a*b + 27*c) / 54;
  FLOAT QQQ = Q*Q*Q;
  FLOAT D   = R*R - QQQ;

  if(D < 0.0) {
    FLOAT sqrtQ = sqrt(Q);
    FLOAT fi = acos(R/sqrt(QQQ));
    *r = -2 * sqrtQ * cos( fi             / 3) - a / 3;
    *s = -2 * sqrtQ * cos((fi + 2 * M_PI) / 3) - a / 3;
    *t = -2 * sqrtQ * cos((fi - 2 * M_PI) / 3) - a / 3;
    return 3;
  } else {
    FLOAT sig = R < 0.0 ? -1.0 : 1.0;
    FLOAT A = -sig * pow((abs(R) + sqrt(D)), 1.0 / 3.0);
    FLOAT B = A == 0 ? 0.0 : Q / A;
    *r =        (A+B) - a / 3;
    *s = -1/2 * (A+B) - a / 3;
    *t = M_SQRT3 / 2 * (A-B);
    return *t == 0 ? 2 : 1;
  }
}
*/

bool isBigEndian()
{
  int a = 1;
  return (bool) ((char *) &a)[0] != 1;
}

int my_hcreate_r(size_t nel, struct my_hsearch_data *tab)
{ 
  memset(&tab->tab, 0, sizeof(tab->tab));
  if((tab->entry = (ENTRY **) malloc(sizeof(ENTRY *)*nel)) == NULL) return 0;
  if(!hcreate_r(nel, &tab->tab)) {
    free(tab->entry);
    return 0;
  }
  tab->nentries = 0;
  tab->tabsize  = nel;
  return 1;
}

int my_hsearch_r(ENTRY item, ACTION action, ENTRY **ret,
                 struct my_hsearch_data *tab)
{
  int i, cc;
  if(action == ENTER && tab->nentries == tab->tabsize) {
    struct hsearch_data newtab;
    ENTRY **epp;

    epp = (ENTRY **) realloc(tab->entry, sizeof(ENTRY *) * tab->tabsize * 2);
    if(epp == NULL) {
      *ret = NULL;
      return 0;
    }

    tab->entry = epp;
    memset(&newtab, 0, sizeof(newtab));
    if(!hcreate_r(tab->tabsize * 2, &newtab)) {
      *ret = NULL;
      return 0;
    }

    tab->tabsize *= 2;

    for(i = 0; i < tab->nentries; i++) {
      cc = hsearch_r(*tab->entry[i], ENTER, ret, &newtab); assert(cc);
      tab->entry[i] = *ret;
    }

    hdestroy_r(&tab->tab);
    tab->tab = newtab;
  }

  cc = hsearch_r(item, action, ret, &tab->tab);
  assert(action == FIND || cc);

  if(action == ENTER && (*ret)->data == item.data) {
    tab->entry[tab->nentries++] = *ret;
  }

  return 1;
}

void my_hdestroy_r(struct my_hsearch_data *tab, int freeKeys)
{
  int i;
  if(freeKeys) {
    for(i = 0; i < tab->nentries; i++) free(tab->entry[i]->key);
  }
  hdestroy_r(&tab->tab);
  free(tab->entry);
}

#include <stdlib.h>
#include <stdio.h>


int fprintHTKstr(FILE *fp, const char *str)
{
  return fprintf(fp, str    == NULL ? "!NULL"   :
                     str[0] == '"'  ? "'%s'"   :
                     str[0] == '\'' ? "\"%s\"" :
                                      "%s", str);
}

int getHTKstr(char *str, char **endPtrOrErrMsg)
{
  char termChar = '\0';
  char *chrptr = str;

  while(isspace(*chrptr)) ++chrptr;

  if(*chrptr == '\'' || *chrptr == '"') {
    termChar = *chrptr;
    chrptr++;
  }

  for(; *chrptr; chrptr++) {
    if(*chrptr == '\'' || *chrptr == '"') {
      if(termChar == *chrptr) {
        termChar = '\0';
        chrptr++;
        break;
      }
    }

    if(isspace(*chrptr) && !termChar) {
      break;
    }

    if(*chrptr == '\\') {
      ++chrptr;
      if(*chrptr == '\0' || (*chrptr    >= '0' && *chrptr <= '7' &&
                            (*++chrptr  <  '0' || *chrptr >  '7' ||
                             *++chrptr  <  '0' || *chrptr >  '7'))) {
        *endPtrOrErrMsg = "Invalid escape sequence";
        return -1;
      }

      if(*chrptr  >= '0' && *chrptr  <= '7') {
        *chrptr = (*chrptr - '0') + (chrptr[-1] - '0') * 8 + (chrptr[-2] - '0') * 64;
      }
    }
    *str++ = *chrptr;
  }

  if(termChar) {
    *endPtrOrErrMsg = "Unterminated quoted string";
     return -2;
  }


  while(isspace(*chrptr)) ++chrptr;
  *endPtrOrErrMsg = chrptr;
  *str = '\0';

  return 0;
}

const char *filter_wldcrd;
char *expandFilterCommand(const char *command, const char *filename)
{
   char *out, *outend;
   const char *chrptr = command;
   int ndollars = 0;
   int fnlen = strlen(filename);

   while(*chrptr++) ndollars += (*chrptr ==  *filter_wldcrd);

   out = (char*) malloc(strlen(command) - ndollars + ndollars * fnlen + 1);
   if(out == NULL) Error("Insufficient memory");

   outend = out;

   for(chrptr = command; *chrptr; chrptr++) {
     if(*chrptr ==  *filter_wldcrd) {
       strcpy(outend, filename);
       outend += fnlen;
     } else {
       *outend++ = *chrptr;
     }
   }
   *outend = '\0';
   return out;
}

#define LINEBUFF_INIT_SIZE 1024
#define LINEBUFF_INCR_SIZE 1024
char *readline(FILE *fp, struct readline_data *data)
{
  char *chrptr;
  int   len;

  if(data->size == 0) {
    data->size = LINEBUFF_INIT_SIZE;
    data->buffer = (char *) malloc(data->size);
    if(data->buffer == NULL) Error("Insufficient memory");
  } else if(data->size == -1) { // EOF reached in previous call
    data->size = 0;
    free(data->buffer);
    return NULL;
  }
  chrptr = data->buffer;
  *chrptr = '\0';

  for(;;) {
    if(fgets(chrptr, data->size-(chrptr-data->buffer), fp) == NULL) {
      if(chrptr != data->buffer) {
        data->size = -1;     // EOF reached but buffer is not empty; next call return NULL
        return data->buffer;
      } else {               // EOF reached and buffer is empty; return NULL imediately
        data->size = 0;
        free(data->buffer);
        return NULL;
      }
    }
    len = strlen(chrptr);
    assert(len >= 0 && len < data->size-(chrptr-data->buffer));

    if(len > 0 && chrptr[len-1] == '\n') {
      chrptr[len-1] = '\0';
      return data->buffer;
    }
    chrptr += len;
    if(chrptr-data->buffer+1 == data->size) {
      data->size += LINEBUFF_INCR_SIZE;
      data->buffer = (char *) realloc(data->buffer, data->size);
      if(data->buffer == NULL) Error("Insufficient memory");
      chrptr = data->buffer+data->size-LINEBUFF_INCR_SIZE-1;
    }
  }
}


static char *HModules[] = {
"HADAPT","HAUDIO","HFB",   "HLABEL","HMAP", "HMEM",  "HMODEL","HNET",  "HPARM",
"HREC",  "HSHELL","HTRAIN","HWAVE", "LCMAP","LGBASE","LMODEL","LPCALC","LWMAP"};

void InsertConfigParam(
  struct my_hsearch_data *config_hash,
  const char *param_name,
  const char *value,
  int optionChar) {
  ENTRY e, *ep;
  int i;
  char *bmod, *emod;


  e.key = (char*) malloc(strlen(param_name)+1);
  for(i=0; *param_name; param_name++) {
    if(*param_name != '-' && *param_name != '_') {
      e.key[i++] = toupper(*param_name);
    }
  }
  e.key[i] = '\0';

  // Find and remove any known module name from param name
  if((emod = strrchr(e.key, ':')) != NULL && emod != e.key) {
    *emod = '\0';
    if((bmod = strrchr(e.key, ':')) == NULL) bmod = e.key;
    else bmod++;

    if(bsearch(&bmod, HModules, sizeof(HModules)/sizeof(*HModules),
              sizeof(*HModules), qsstrcmp)) {
      memmove(bmod, emod+1, strlen(emod+1)+1);
    } else {
      *emod = ':';
    }
  }

  // Replace initial 'H' in application name by 'S'
  if(strrchr(e.key, ':') && e.key[0] == 'H') {
    e.key[0] = 'S';
  }

  my_hsearch_r(e, FIND, &ep, config_hash);

  if((e.data = malloc(strlen(value) + 3)) == NULL) Error("Insufficient memory");
  strcpy(((char *) e.data) + 2, value);

  //1st byte of the value string is read/unread flag
  ((char *) e.data)[0] = '#';
  ((char *) e.data)[1] = optionChar;


  if(ep == NULL) {
    if(e.key == NULL || !my_hsearch_r(e, ENTER, &ep, config_hash)) {
      Error("Insufficient memory");
    }
  } else {
    free(ep->data);
    ep->data = e.data;
    free(e.key);
  }
}

void ReadConfig(const char *file_name, struct my_hsearch_data *config_hash)
{
  char *parptr, *endptr, *chptr, *line;
  struct readline_data rld = {0};
  int line_no = 0;
  FILE *fp;

  if((fp = fopen(file_name, "rb")) == NULL) {
    Error("Cannot open input config file %s", file_name);
  }

  while((line = readline(fp, &rld)) != NULL) {
    line_no++;

    if((chptr = strchr(line, '#')) != NULL) *chptr = '\0';

    for(parptr = line; isspace(*parptr); parptr++);

    if(*parptr == '\0') continue;

    chptr = parptr;
    for(;;) {
      // Replace speces by '_', which is removed in InsertConfigParam
      while(isalnum(*chptr) || *chptr == '_' || *chptr == '-') chptr++;
      while(isspace(*chptr)) *chptr++ = '_';
      if(*chptr != ':') break;
      chptr++;
      while(isspace(*chptr)) *chptr++ = '_';
    }
    if(*chptr != '=') Error("Character '=' expected (%s:%d, char %d)",
                            file_name, line_no, chptr-line+1);
    *chptr = '\0';
    chptr++;
    if(getHTKstr(chptr, &endptr)) {
      Error("%s (%s:%d)", endptr, file_name, line_no);
    }
    if(*endptr) Error("Extra characters at the end of line (%s:%d, char %d)",
                      file_name, line_no, endptr-line+1);

    InsertConfigParam(config_hash, parptr, chptr, 'C');
  }

  if(ferror(fp)) {
    Error("Cannot read input config file", file_name);
  }
  fclose(fp);
}

void CheckCommandLineParamUse(struct my_hsearch_data *config_hash)
{
  int i;
  for(i = 0; i < config_hash->nentries; i++) {
    char *value = (char *) config_hash->entry[i]->data;
    if(value[0] != ' ' && value[1] != 'C') {
      assert(strchr(config_hash->entry[i]->key, ':'));
      Error("Unexpected command line parameter %s",
            strchr(config_hash->entry[i]->key, ':') + 1);
    }
  }
}


void PrintConfig(struct my_hsearch_data *config_hash)
{
  int i;
  char *par, *key, *val;
  printf("\nConfiguration Parameters[%d]\n", config_hash->nentries);
  for(i = 0; i < config_hash->nentries; i++) {
    key = (char *) config_hash->entry[i]->key;
    val = (char *) config_hash->entry[i]->data;
    par = strrchr(key, ':');
    if(par) par++; else par = key;
    printf("%c %-15.*s %-20s = %-30s # -%c\n",
           val[0], par-key, key, par, val+2, val[1]);
  }
  putchar('\n');
}

ENTRY *GetParam(
   struct my_hsearch_data *config_hash,
   const char *param_name)
{
  ENTRY e, *ep;

  e.key = (char *) (param_name - 1);
  do {
    e.key++;
    my_hsearch_r(e, FIND, &ep, config_hash);
  } while(ep == NULL && (e.key = strchr(e.key, ':')) != NULL);
  if(ep != NULL) *(char *) ep->data = ' '; //Mark param as read
  return ep;
}


const char *GetParamStr(
   struct my_hsearch_data *config_hash,
   const char *param_name,
   const char *default_value)
{
  ENTRY *ep = GetParam(config_hash, param_name);
  return ep != NULL ? 2 + (char *) ep->data : default_value;
}

static char *OptOrParStr(ENTRY *ep)
{
  static char str[128];
  int optionChar = ((char *) ep->data)[1];
  if(optionChar == 'C') {
    sprintf(str, "parameter '%s'", ep->key);
  } else if(optionChar == '-') {
    sprintf(str, "option '--%s'", ep->key);
  } else{
    sprintf(str, "option '-%c'", optionChar);
  }
  return str;
}

long GetParamInt(
   struct my_hsearch_data *config_hash,
   const char *param_name,
   long default_value)
{
  char *chrptr;
  ENTRY *ep = GetParam(config_hash, param_name);
  if(ep == NULL) return default_value;

  const char *val = 2 + (char *) ep->data;
  default_value = strtol(val, &chrptr, 0);
  if(!*val || *chrptr) {
    Error("Integer number expected for %s but found '%s'", OptOrParStr(ep),val);
  }
  return default_value;
}

FLOAT GetParamFlt(
   struct my_hsearch_data *config_hash,
   const char *param_name,
   FLOAT default_value)
{
  char *chrptr;
  ENTRY *ep = GetParam(config_hash, param_name);
  if(ep == NULL) return default_value;

  const char *val = 2 + (char *) ep->data;
  default_value = strtod(val, &chrptr);
  if(!*val || *chrptr) {
    Error("Decimal number expected for %s but found '%s'", OptOrParStr(ep),val);
  }
  return default_value;
}

bool GetParamBool(
   struct my_hsearch_data *config_hash,
   const char *param_name,
   bool default_value)
{
  ENTRY *ep = GetParam(config_hash, param_name);
  if(ep == NULL) return default_value;

  const char *val = 2 + (char *) ep->data;
  if(!strcmp(val, "TRUE") || !strcmp(val, "T")) return 1;
  if(strcmp(val, "FALSE") && strcmp(val, "F")) {
    Error("TRUE or FALSE expected for %s but found '%s'", OptOrParStr(ep),val);
  }
  return false;
}

// '...' are pairs: string and corresponding integer value , terminated by NULL
int GetParamEnum(
   struct my_hsearch_data *config_hash,
   const char *param_name,
   int default_value, ...)
{
  ENTRY *ep = GetParam(config_hash, param_name);
  if(ep == NULL) return default_value;

  const char *val = 2 + (char *) ep->data;
  char *s;
  int i, cnt = 0, l = 0;
  va_list ap;

  va_start(ap, default_value);
  while((s = va_arg(ap, char *)) != NULL) {
    l += strlen(s) + 2;
    ++cnt;
    i = va_arg(ap, int);
    if(!strcmp(val, s)) break;
  }
  va_end(ap);
  if(s) return i;

  //To report error, create string listing all possible values
  s = (char*) malloc(l + 1);
  s[0] = '\0';
  va_start(ap, default_value);
  for(i = 0; i < cnt; i++) {
    strcat(s, va_arg(ap, char *));
    va_arg(ap, int);
    if(i < cnt - 2) strcat(s, ", ");
    else if(i == cnt - 2) strcat(s, " or ");
  }
  va_end(ap);
  Error("%s expected for %s but found '%s'", s, OptOrParStr(ep),val);
  return 0;
}

int GetDerivParams(
  struct my_hsearch_data *config_hash,
  int *derivOrder,
  int **derivWinLens,
  int *startFrmExt,
  int *endFrmExt,
  char **CMNPath,
  char **CMNFile,
  const char **CMNMask,
  char **CVNPath,
  char **CVNFile,
  const char **CVNMask,
  const char **CVGFile,
  const char *toolName,
  int pseudoModeule)
{
  const char *str;
  int targetKind;
  char *chrptr, paramName[32];
  const char *CMNDir, *CVNDir;
  strcpy(paramName, toolName);
  strcat(paramName, pseudoModeule == 1 ? "SPARM1:" :
                    pseudoModeule == 2 ? "SPARM2:" : "");
  chrptr = paramName + strlen(paramName);

  strcpy(chrptr, "STARTFRMEXT");
  *startFrmExt = GetParamInt(config_hash, paramName, 0);
  strcpy(chrptr, "ENDFRMEXT");
  *endFrmExt   = GetParamInt(config_hash, paramName, 0);

  *CMNPath = *CVNPath = NULL;
  strcpy(chrptr, "CMEANDIR");
  CMNDir       = GetParamStr(config_hash, paramName, NULL);
  strcpy(chrptr, "CMEANMASK");
  *CMNMask     = GetParamStr(config_hash, paramName, NULL);
  if(*CMNMask != NULL) {
    *CMNPath = (char*) malloc((CMNDir ? strlen(CMNDir) : 0) + npercents(*CMNMask) + 2);
    if(*CMNPath == NULL) Error("Insufficient memory");
    if(CMNDir != NULL) strcat(strcpy(*CMNPath, CMNDir), "/");
    *CMNFile = *CMNPath + strlen(*CMNPath);
  }
  strcpy(chrptr, "VARSCALEDIR");
  CVNDir      = GetParamStr(config_hash, paramName, NULL);
  strcpy(chrptr, "VARSCALEMASK");
  *CVNMask     = GetParamStr(config_hash, paramName, NULL);
  if(*CVNMask != NULL) {
    *CVNPath = (char*) malloc((CVNDir ? strlen(CVNDir) : 0) + npercents(*CVNMask) + 2);
    if(*CVNPath == NULL) Error("Insufficient memory");
    if(CVNDir != NULL) strcat(strcpy(*CVNPath, CVNDir), "/");
    *CVNFile = *CVNPath + strlen(*CVNPath);
  }
  strcpy(chrptr, "VARSCALEFN");
  *CVGFile     = GetParamStr(config_hash, paramName, NULL);
  strcpy(chrptr, "TARGETKIND");
  str = GetParamStr(config_hash, paramName, "ANON");
  targetKind = ReadParmKind(str, FALSE);
  if(targetKind == -1) Error("Invalid TARGETKIND = '%s'", str);

  strcpy(chrptr, "DERIVWINDOWS");
  if((str = GetParamStr(config_hash, paramName, NULL)) != NULL) {
    long lval;
    *derivOrder      = 0;
    *derivWinLens = NULL;
    while((str = strtok((char *) str, " \t_")) != NULL) {
      lval = strtol(str, &chrptr, 0);
      if(!*str || *chrptr) {
        Error("Integers separated by '_' expected for parameter DERIVWINDOWS");
      }
      *derivWinLens = (int *)realloc(*derivWinLens, ++*derivOrder*sizeof(int));
      if(*derivWinLens == NULL) Error("Insufficient memory");
      (*derivWinLens)[*derivOrder-1] = lval;
      str = NULL;
    }
    return targetKind;
  }
  *derivOrder = targetKind & PARAMKIND_T ? 3 :
                targetKind & PARAMKIND_A ? 2 :
                targetKind & PARAMKIND_D ? 1 : 0;

  if(*derivOrder || targetKind != PARAMKIND_ANON) {
   *derivWinLens = (int *) malloc(3 * sizeof(int));
    if(*derivWinLens == NULL) Error("Insufficient memory");

    strcpy(chrptr, "DELTAWINDOW");
    (*derivWinLens)[0] = GetParamInt(config_hash, paramName, 2);
    strcpy(chrptr, "ACCWINDOW");
    (*derivWinLens)[1] = GetParamInt(config_hash, paramName, 2);
    strcpy(chrptr, "THIRDWINDOW");
    (*derivWinLens)[2] = GetParamInt(config_hash, paramName, 2);
    return targetKind;
  }
  *derivWinLens = NULL;
  *derivOrder   = -1;
  return targetKind;
}

FILE *my_fopen(
  const char *file_name,
  const char *type,
  const char *filter)
{
  FILE *fp;

  if(!strcmp(file_name, "-")) {
    fp = *type == 'r' ? stdin : stdout;
  } else if(*file_name == '|') {
    if((fp = popen(++file_name, *type == 'r' ? "r" : "w")) == NULL) {
      Error("Cannot popen %s filter '%s'",
            *type == 'r' ? "input": "output", file_name);
    }
  } else if(filter) {
    char *f = expandFilterCommand(filter, file_name);

    if((fp = popen(f, *type == 'r' ? "r" : "w")) == NULL) {
      Error("Cannot popen %s filter '%s'",
            *type == 'r' ? "input": "output", f);
    }
    free(f);
  } else if((fp = fopen(file_name, type)) == NULL) {
    return NULL;
  }
  return fp;
}

int my_fclose(FILE *fp)
{
  struct stat sb;

  if(fp == stdin || fp == stdout) return 0;

  if(fstat(fileno(fp), &sb)) return EOF;

  if(S_ISFIFO(sb.st_mode)) return pclose(fp);
  else return fclose(fp);
}

const char *script_filter;
const char *parm_filter;
const char *MMF_filter;
const char *parm_ofilter;
const char *hlist_ofilter;
const char *MMF_ofilter;

int ParseOptions(
  int argc,
  char *argv[],
  const char *optionMapping,
  const char *toolName,
  struct my_hsearch_data *cfgHash)
{
  int i, opt = '?', optind;
  BOOL option_must_follow = FALSE;
  char param[1024], *value, *optfmt;
  const char *optarg;
  char *chptr, *bptr, tstr[4] = " -?";
  unsigned long long option_mask = 0;

  #define MARK_OPTION(ch) {if(isalpha(ch)) option_mask |= 1ULL << ((ch) - 'A');}
  #define OPTION_MARK(ch) (isalpha(ch) && ((1ULL << ((ch) - 'A')) & option_mask))
  #define IS_OPTION(str) ((str)[0] == '-' && (isalpha((str)[1]) || (str)[1] == '-'))

  for(optind = 1; optind < argc; optind++) {
    if(!strcmp(argv[optind], "--")) break;
    if(argv[optind][0] != '-' || argv[optind][1] != 'A') continue;
    if(argv[optind][2] != '\0') {
      Error("Unexpected argument '%s' after option '-A'", argv[optind] + 2);
    }
    for(i=0; i < argc; i++) printf("%s ", argv[i]);
    putchar('\n');
    break;
  }
  for(optind = 1; optind < argc; optind++) {
    if(!strcmp(argv[optind], "--")) break;
    if(argv[optind][0] != '-' || argv[optind][1] != 'C') continue;
    if(argv[optind][2] != '\0') {
      ReadConfig(argv[optind] + 2, cfgHash);
    } else if(optind+1 < argc && !IS_OPTION(argv[optind+1])) {
      ReadConfig(argv[++optind], cfgHash);
    } else {
      Error("Config file name expected after option '-C'");
    }
  }
  for(optind = 1; optind < argc; optind++) {
    if(!strcmp(argv[optind], "--")) break;
    if(argv[optind][0] != '-' || argv[optind][1] != '-') continue;
    bptr = (char*) malloc(strlen(toolName) + strlen(argv[optind]+2) + 2);
    if(bptr == NULL) Error("Insufficient memory");
    strcat(strcat(strcpy(bptr, toolName), ":"), argv[optind]+2);
    value = strchr(bptr, '=');
    if(!value) Error("Character '=' expected after option '%s'", argv[optind]);
    *value++ = '\0';
    InsertConfigParam(cfgHash, bptr, value /*? value : "TRUE"*/, '-');
    free(bptr);
  }
  for(optind = 1; optind < argc && IS_OPTION(argv[optind]); optind++) {
    option_must_follow = FALSE;
    tstr[2] = opt = argv[optind][1];
    optarg = argv[optind][2] != '\0' ? argv[optind] + 2 : NULL;

    if(opt == '-' && !optarg) {    // '--' terminates the option list
      return optind+1;
    }
    if(opt == 'C' || opt == '-') { // C, A and long options have been already processed
      if(!optarg) optind++;
      continue;
    }
    if(opt == 'A') continue;

    chptr = strstr(optionMapping, tstr);
    if(chptr == NULL) Error("Invalid command line option '-%c'", opt);

    chptr += 3;
    while(isspace(*chptr)) chptr++;

    if(!chptr || chptr[0] == '-') {// Option without format string will be ignored
      optfmt = " ";
    } else {
      optfmt = chptr;
      while(*chptr && !isspace(*chptr)) chptr++;
      if(!*chptr) Error("Fatal: Unexpected end of optionMap string");
    }
    for(i = 0; !isspace(*optfmt); optfmt++) {
      while(isspace(*chptr)) chptr++;
      value = chptr;
      while(*chptr && !isspace(*chptr)) chptr++;
      assert(chptr-value+1 < sizeof(param));
      strncat(strcat(strcpy(param, toolName), ":"), value, chptr-value);
      param[chptr-value+strlen(toolName)+1] = '\0';
      switch(*optfmt) {
        case 'n': value = strchr(param, '=');
                  if(value) *value = '\0';
                  InsertConfigParam(cfgHash, param,
                                    value ? value + 1: "TRUE", opt);
                  break;
        case 'l':
        case 'o':
        case 'r': i++;
                  if(!optarg && (optind+1==argc || IS_OPTION(argv[optind+1]))) {
                    if(*optfmt == 'r' || *optfmt == 'l') {
                      Error("Argument %d of option '-%c' expected", i, opt);
                    }
                    optfmt = "  "; // Stop reading option arguments
                    break;
                  }
                  if(!optarg) optarg = argv[++optind];
                  if(*optfmt == 'o') {
                    option_must_follow = (BOOL) 1;
                  }
                  bptr = NULL;

                  // For repeated use of option with 'l' (list) format, append
                  // ',' and argument string to existing config parameter value.
                  if(*optfmt == 'l' && OPTION_MARK(opt)) {
                    bptr = strdup(GetParamStr(cfgHash, param, ""));
                    if(bptr == NULL) Error("Insufficient memory");
                    bptr = (char*) realloc(bptr, strlen(bptr) + strlen(optarg) + 2);
                    if(bptr == NULL) Error("Insufficient memory");
                    strcat(strcat(bptr, ","), optarg);
                    optarg = bptr;
                  }
                  MARK_OPTION(opt);
                  InsertConfigParam(cfgHash, param, optarg, opt);
                  free(bptr);
                  optarg = NULL;
                  break;
        default : Error("Fatal: Invalid character '%c' in optionMap after %s",
                        *optfmt, tstr);
      }
    }
    if(optarg) Error("Unexpected argument '%s' after option '-%c'",optarg,opt);
  }

  for(i = optind; i < argc && !IS_OPTION(argv[i]); i++);
  if(i < argc) {
    Error("No option expected after first non-option argument '%s'", argv[optind]);
  }
  if(option_must_follow) {
    Error("Option '-%c' with optional argument must not be the last option",opt);
  }
  return optind;
}

FileListElem **AddFileElem(FileListElem **last, char *fileElem)
{
  char *chrptr;
  *last = (FileListElem *) malloc(sizeof(FileListElem) + strlen(fileElem));
  if(!*last) Error("Insufficient memory");
  chrptr = strcpy((*last)->logical, fileElem);
  for(; *chrptr; chrptr++) if(*chrptr == '\\') *chrptr = '/';
  chrptr = strchr((*last)->logical, '=');
  if(chrptr) *chrptr = '\0';
  (*last)->physical = chrptr ? chrptr+1: (*last)->logical;
  last = &(*last)->next;
  *last = NULL;
  return last;
}

static int memcmpw(const char *nstr, const char *wstr, int len, char **substr)
{
  int i, npercents = 0;
  
  for(i=0; i<len; i++) {
    if((*wstr != *nstr) && (*wstr != '?') && (*wstr != '%')) return 1;
    if(*wstr == '%') (*substr)[npercents++] = *nstr;
    wstr++; nstr++;
  }
  *substr += npercents;
  **substr = '\0';
  return 0;
}

int npercents(const char *str)
{
  int ret = 0;
  while(*str) if(*str++ == '%') ret++;
  return ret;
}

int process_mask(const char *normstr, const char *wildcard, char *substr)
{
  char *hlpptr;
  const char  *endwc = wildcard + strlen(wildcard);
  const char  *endns = normstr + strlen(normstr);

  *substr = '\0';

  if((hlpptr = (char *) memchr(wildcard, '*', endwc-wildcard)) == NULL) {
    return !((endwc-wildcard != endns-normstr) || 
             memcmpw(normstr, wildcard, endns-normstr, &substr));
  }

  if((hlpptr-wildcard > endns-normstr) || 
     memcmpw(normstr, wildcard, hlpptr-wildcard, &substr)) {
    return 0;
  }
  wildcard = hlpptr;
  
  for(;;) {
    while(*wildcard == '*') wildcard++;
    if(!*wildcard) return 1;
    if(!*normstr)  return 0;
    if((hlpptr = (char *) memchr(wildcard, '*', endwc-wildcard)) == NULL) {
      return !((endwc-wildcard > endns-normstr) ||
               memcmpw(endns-(endwc-wildcard), wildcard, endwc-wildcard, &substr));
    }
    
    do {
      if((endns-normstr) < (hlpptr-wildcard)) return 0;
      normstr = (char *) ((*wildcard == '?' || *wildcard == '%')
                ? normstr 
                : memchr(normstr, *wildcard, endns-normstr-(hlpptr-wildcard)+1));
      if(!normstr || !*normstr) return 0;
    } while(memcmpw(normstr++, wildcard, hlpptr-wildcard, &substr));
    
    normstr += hlpptr - wildcard - 1;
    wildcard = hlpptr;
  }
}