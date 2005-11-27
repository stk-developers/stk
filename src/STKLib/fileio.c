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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <ctype.h>
#ifndef WIN32
#include <unistd.h>
#else
#include <io.h>
#include <direct.h>
#endif
#include "fileio.h"

int WriteHTKHeader (FILE * fp_out, HTK_Header header, bool swap)
{
  int cc;

  if(swap) {
    swap4(header.nSamples);
    swap4(header.sampPeriod);
    swap2(header.sampSize);
    swap2(header.sampKind);
  }

  fseek (fp_out, 0L, SEEK_SET);
  cc = fwrite(&header, sizeof(HTK_Header), 1, fp_out);

  if(swap) {
    swap4(header.nSamples);
    swap4(header.sampPeriod);
    swap2(header.sampSize);
    swap2(header.sampKind);
  }

  return cc == 1 ? 0 : -1;
}

int WriteHTKFeature(FILE * fp_out, FLOAT *out, size_t fea_len, bool swap) {
  size_t i, cc = 0;
#if !DOUBLEPRECISION
  if(swap) for(i = 0; i < fea_len; i++) swap4(out[i]);
  cc = fwrite(out, sizeof(FLOAT_32), fea_len, fp_out);
  if(swap) for(i = 0; i < fea_len; i++) swap4(out[i]);
#else
  float f;

  for(i = 0; i < fea_len; i++) {
    f = out[i];
    if(swap) swap4(f);
    cc += fwrite(&f, sizeof(FLOAT_32), 1, fp_out);
  }
#endif

  return cc == fea_len ? 0 : -1;
}

int ReadHTKHeader (FILE * fp_in, HTK_Header * header, bool swap)
{
  if(!fread(&header->nSamples,   sizeof(INT_32),  1, fp_in)) return -1;
  if(!fread(&header->sampPeriod, sizeof(INT_32),  1, fp_in)) return -1;
  if(!fread(&header->sampSize,   sizeof(INT_16),  1, fp_in)) return -1;
  if(!fread(&header->sampKind,   sizeof(UINT_16), 1, fp_in)) return -1;

  if(swap) {
    swap4(header->nSamples);
    swap4(header->sampPeriod);
    swap2(header->sampSize);
    swap2(header->sampKind);
  }

  if(header->sampPeriod < 0 || header->sampPeriod > 100000 ||
     header->nSamples   < 0 || header->sampSize   < 0) {
    return -1;
  }

  return 0;
}


int ReadHTKFeature(FILE * fp_in, FLOAT *in, size_t fea_len, bool swap,
                   bool decompress, FLOAT *scale, FLOAT *bias)
{
  size_t i;
  if(decompress) {
    INT_16 s;
//    FLOAT scale = (xmax - xmin) / (2*32767);
//    FLOAT bias  = (xmax + xmin) / 2;

    for(i = 0; i < fea_len; i++) {
      if(fread(&s, sizeof(INT_16), 1, fp_in) != 1) {
        return -1;
      }
      if(swap) swap2(s);
      in[i] = (s + bias[i]) / scale[i];
    }
    return 0;
  }

#if !DOUBLEPRECISION
  if(fread(in, sizeof(FLOAT_32), fea_len, fp_in) != fea_len) {
    return -1;
  }
  if(swap) for(i = 0; i < fea_len; i++) swap4(in[i]);
#else
  float f;

  for(i = 0; i < fea_len; i++) {
    if(fread(&f, sizeof(FLOAT_32), 1, fp_in) != 1) {
      return -1;
    }
    if(swap) swap4(f);
    in[i] = f;
  }
#endif
  return 0;
}

int Mkdir4File(const char *file_name)
{
  struct stat statbuff;
  char dir_name[260];
  char *chptr;

  if((chptr=strrchr(file_name, '/')) == NULL) {
    return 0;
  }
  strncpy(dir_name, file_name, chptr-file_name);
  dir_name[chptr-file_name]='\0';
  chptr=dir_name;
  if(*chptr == '/') {
    chptr++;
  }
  if(!chptr[0]) {
    return 0;
  }
  for(;;) {
    if((chptr=strchr(chptr, '/')) != NULL) {
      *chptr='\0';
    }
    if((access(dir_name, 0) || stat(dir_name, &statbuff) == -1 ||
       !(S_IFDIR & statbuff.st_mode))) {
#ifndef WIN32
      if(mkdir(dir_name, 0777)) return -1;
#else
      if(mkdir(dir_name)) return -1;
#endif
    }
    if(chptr) {
      *chptr = '/';
    }
    if(!chptr || !chptr[1]) {
      return 0;
    }
    chptr++;
  }
}

/*
FUNCTION

  FLOAT *ReadHTKFeatures(char *file_name, bool swap, int extLeft, int extRight,
                         int targetKind, int derivOrder, int *derivWinLen,
                         HTK_Header *header);

DESCRIPTION

  Function opens HTK feature file, reads its contents, optionally augment
  features with their derivatives (of any order), optionally extend (augmented)
  features by repeating initial or final frames and return pointer to buffer
  with results allocated by malloc.

PARAMETERS

  file_name:
    Name of features file optionally augmented with frame range specification
    in the form 'file_name[start,end]', where 'start' and 'end' are integers.
    If range is specified, only the corresponding segment of feature frames is
    extracted from file. Range [0,0] corresponds to only first frame. The range
    can also extend over feature file boundaries (e.g. start can be negative).
    In such case initial/final frames are repeated for the extending frames, so
    that resulting features have allays end-start+1 frames.

  swap:
    Boolean value specifies whether to swap bytes when reading file or not.

  extLeft:
    Features read from file are extended with extLeft initial frames. Normally,
    these frames are repetitions of the first feature frame in file (with its
    derivative, if derivatives are preset in the file).  However, if segment of
    feature frames is extracted according to range specification, the true
    feature frames from beyond the segment boundary are used, wherever it is
    possible. Note that value of extLeft can be also negative. In such case
    corresponding number of initial frames is discarded.

  extRight:
    The paramerer is complementary to parameter extLeft and has obvious
    meaning. (Controls extensions over the last frame, last frame from file is
    repeated only if necessary).

  targetKind:
    The parameters is used to check whether header->sampKind match to
    requited targetKind and to control suppression of 0'th cepstral or energy
    coefficients accorging to modifiers _E, _0 and _N. Modifiers _D, _A and _T
    are ignored; Computation of derivatioves is controled by parameters
    derivOrder and derivWinLen. Value PARAMKIND_ANON ensures that function do
    not result in targetKind mismatch error and cause no _E or _0 suppression.

  derivOrder:
    Final features will be augmented with their derivatives up to 'derivOrder'
    order. If 'derivOrder' is negative value, no new derivatives are appended
    and derivatives that already present in feature file are preserved.
    Straight features are considered to be of zero order. If some
    derivatives are already present in feature file, these are not computed
    again, only higher order derivatives are appended if required. Note, that
    HTK feature file cannot contain higher order derivatives (e.g. double
    delta) without containing lower ones (e.g. delta). Derivative present
    in feature file that are of higher order than is required are discarded.
    Derivatives are computed in the final stage from (extracted segment of)
    feature frames possibly extended by repeated frames. Derivatives are
    computed using the same formula that is employed also by HTK tools. Lengths
    of windows used for computation of derivatives are passed in parameter
    derivWinLen. To compute derivatives for frames close to boundaries,
    frames before the first and after the last frame (of the extracted segment)
    are considered to be (yet another) repetitions of the first and the
    last frame, respectively. If the segment of frames is extracted according to
    range specification and parameters extLeft and extLeft are set to zero,
    the first and the last frames of the segment are considered to be repeated,
    eventough the true feature frames from beyond the segment boundary can be
    available in the file. Therefore, segment extracted from features that were
    before augmented with derivatives will differ from the same segment augmented
    with derivatives by this function. Difference will be of course only on
    boundaries and only in derivatives. This "incorrect" behavior was chosen
    to fully simulate behavior of HTK tools. To obtain more correct computation
    of derivatives, use parameters extLeft and extRight, which correctly extend
    segment with the true frames (if possible) and in resulting feature matrix
    ignore first extLeft and last extRight frames. For this purpose, both
    extLeft and extRight should be set to sum of all values in the array
    derivWinLen.

  derivWinLen:
    Array of size derivOrder specifying lengths of windows used for computation
    of derivatives. Individual values represents one side context used in
    the computation. The each window length is therefore twice the value
    from array plus one. Value at index zero specify window length for first
    order derivatives (delta), higher indices corresponds to higher order
    derivatives.

  header:
    Structure *header is filled by information read from header of HTK feature
    file. The information is changed according to modifications made to
    features: bits representing added/removed derivatives are set/unset in
    header->sampKind, header->nSamples is set to correspond to the number of
    frames in the resulting feature matrix.

RETURN VALUE

    Function return pointer to buffer (allocated by malloc) with resulting
    feature matrix.
*/

enum CNFileType {CNF_Mean, CNF_Variance, CNF_VarScale};
void ReadCepsNormFile(const char *file, char **last_file, FLOAT **vec_buff,
                      int sampKind, enum CNFileType type, int coefs);

FLOAT *ReadHTKFeatures(
  char *                file_name,
  bool                  swap,
  int                   extLeft,
  int                   extRight,
  int                   targetKind,
  int                   derivOrder,
  int *                 derivWinLen,
  HTK_Header *          header,
  const char *          cmn_file,
  const char *          cvn_file,
  const char *          cvg_file,
  RHFBuffer *           buff)
{
  FLOAT *               feaMx;
  int                   fromFrame;
  int                   toFrame;
  int                   totFrames;
  int                   trg_vec_size;
  int                   src_vec_size;
  int                   i;
  int                   j;
  int                   k;
  int                   e;
  int                   srcDerivOrder;
  int                   loSrcTgzDerOrd;
  int                   coefs;
  int                   trgE;
  int                   trg0;
  int                   trgN;
  int                   srcE;
  int                   src0;
  int                   srcN;
  int                   comp;
  int                   coefSize;
  char *                chptr;

  // remove final spaces from file name
  for(i = strlen(file_name)-1; i >= 0 && isspace(file_name[i]); i--) 
    file_name[i] = '\0';
  
  // read frame range definition if any ( physical_file.fea[s,e] )
  if((chptr = strrchr(file_name, '[')) == NULL ||
     ((i=0), sscanf(chptr, "[%d,%d]%n", &fromFrame, &toFrame, &i), chptr[i] != '\0')) 
  {
    chptr = NULL;
  }
  
  if (chptr != NULL) 
    *chptr = '\0';

  if (strcmp(file_name, "-") && buff->last_file_name &&
      !strcmp(buff->last_file_name, file_name)) 
  {
    *header = buff->last_header;
  } 
  else 
  {
    if (buff->last_file_name) 
    {
      if (buff->fp != stdin) fclose(buff->fp);
      free(buff->last_file_name);
      buff->last_file_name = NULL;
    }
    
    if (!strcmp(file_name, "-")) 
      buff->fp = stdin;
    else 
      buff->fp = fopen(file_name, "rb");
    
    if (buff->fp == NULL) 
      Error("Cannot open feature file: '%s'", file_name);
    
    if (ReadHTKHeader(buff->fp, header, swap)) 
      Error("Invalid HTK header in feature file: '%s'", file_name);
    
    if (header->sampKind & PARAMKIND_C) 
    {
      // File is in compressed form, scale and bias vectors
      // are appended after HTK header.

      int coefs = header->sampSize/sizeof(INT_16);
      buff->A = (float*) realloc(buff->A, coefs * sizeof(FLOAT_32));
      buff->B = (float*) realloc(buff->B, coefs * sizeof(FLOAT_32));
      if(buff->A == NULL || buff->B == NULL) Error("Insufficient memory");

      e  = ReadHTKFeature(buff->fp, buff->A, coefs, swap, 0, 0, 0);
      e |= ReadHTKFeature(buff->fp, buff->B, coefs, swap, 0, 0, 0);
      if(e) Error("Cannot read feature file: '%s'", file_name);
      header->nSamples -= 2 * sizeof(FLOAT_32) / sizeof(INT_16);
    }
    
    if((buff->last_file_name = strdup(file_name)) == NULL) 
    {
      Error("Insufficient memory");
    }
    
    buff->last_header = *header;
  }
  
  if (chptr != NULL) 
    *chptr = '[';

  if (chptr == NULL) 
  { // Range [s,e] was not specified
    fromFrame = 0;
    toFrame   = header->nSamples-1;
  }
  srcDerivOrder = PARAMKIND_T & header->sampKind ? 3 :
                  PARAMKIND_A & header->sampKind ? 2 :
                  PARAMKIND_D & header->sampKind ? 1 : 0;
  srcE = (PARAMKIND_E & header->sampKind) != 0;
  src0 = (PARAMKIND_0 & header->sampKind) != 0;
  srcN =((PARAMKIND_N & header->sampKind) != 0) * (srcE + src0);
  comp =  PARAMKIND_C & header->sampKind;
  header->sampKind &= ~PARAMKIND_C;

  if (targetKind == PARAMKIND_ANON) {
    targetKind = header->sampKind;
  } else if((targetKind & 077) == PARAMKIND_ANON) {
    targetKind &= ~077;
    targetKind |= header->sampKind & 077;
  }
  trgE = (PARAMKIND_E & targetKind) != 0;
  trg0 = (PARAMKIND_0 & targetKind) != 0;
  trgN =((PARAMKIND_N & targetKind) != 0) * (trgE + trg0);


  coefSize = comp ? sizeof(INT_16) : sizeof(FLOAT_32);
  coefs = (header->sampSize/coefSize + srcN) / (srcDerivOrder+1) - srcE - src0;
  src_vec_size = (coefs + srcE + src0) * (srcDerivOrder+1) - srcN;

  //Is coefs dividable by 1 + number of derivatives specified in header
  if(src_vec_size * coefSize != header->sampSize) {
    Error("Invalid HTK header in feature file: '%s'. "
          "sampSize do not match with parmKind", file_name);
  }
  if(derivOrder < 0) derivOrder = srcDerivOrder;


  if((!srcE && trgE) || (!src0 && trg0) || (srcN && !trgN) ||
     (trgN && !trgE && !trg0) || (trgN && !derivOrder) ||
     (srcN && !srcDerivOrder && derivOrder) ||
     ((header->sampKind & 077) != (targetKind & 077) &&
      (header->sampKind & 077) != PARAMKIND_ANON)) 
  {
    char srcParmKind[64], trgParmKind[64];
    ParmKind2Str(header->sampKind, srcParmKind);
    ParmKind2Str(targetKind,       trgParmKind);
    Error("Cannot convert %s to %s", srcParmKind, trgParmKind);
  }

  loSrcTgzDerOrd = LOWER_OF(srcDerivOrder, derivOrder);
  trg_vec_size  = (coefs + trgE + trg0) * (derivOrder+1) - trgN;
  i =  LOWER_OF(fromFrame, extLeft);
  fromFrame -= i;
  extLeft   -= i;

  i =  LOWER_OF(header->nSamples-toFrame-1, extRight);
  toFrame  += i;
  extRight -= i;

  if(fromFrame > toFrame || fromFrame >= header->nSamples || toFrame < 0) {
    Error("Invalid frame range for feature file: '%s'", file_name);
  }
  totFrames = toFrame - fromFrame + 1 + extLeft + extRight;
  feaMx = (FLOAT *) malloc((trg_vec_size * totFrames + 1) * sizeof(FLOAT));
  // + 1 needed for safe reading unwanted _E and _0
  if(feaMx == NULL) Error("Insufficient memory");
  for(i = 0; i <= toFrame - fromFrame; i++) {
    FLOAT *A = buff->A, *B = buff->B;
    FLOAT *mxPtr = feaMx + trg_vec_size * (i+extLeft);
    fseek(buff->fp, sizeof(HTK_Header) + (comp ? src_vec_size * 2 * sizeof(FLOAT_32) : 0) +
              (fromFrame + i) * src_vec_size * coefSize, SEEK_SET);

                        e  = ReadHTKFeature(buff->fp, mxPtr, coefs, swap, comp, A, B);
                        mxPtr += coefs; A += coefs; B+= coefs;
      if(src0 && !srcN) e |= ReadHTKFeature(buff->fp, mxPtr, 1,     swap, comp, A++, B++);
      if(trg0 && !trgN) mxPtr++;
      if(srcE && !srcN) e |= ReadHTKFeature(buff->fp, mxPtr, 1,     swap, comp, A++, B++);
      if(trgE && !trgN) mxPtr++;

    for(j = 0; j < loSrcTgzDerOrd; j++) {
                        e |= ReadHTKFeature(buff->fp, mxPtr, coefs, swap, comp, A, B);
                        mxPtr += coefs; A += coefs; B+= coefs;
      if(src0)          e |= ReadHTKFeature(buff->fp, mxPtr, 1,     swap, comp, A++, B++);
      if(trg0)          mxPtr++;
      if(srcE)          e |= ReadHTKFeature(buff->fp, mxPtr, 1,     swap, comp, A++, B++);
      if(trgE)          mxPtr++;
    }
    if(e) Error("Cannot read feature file: '%s' frame %d/%d", file_name, i, toFrame - fromFrame + 1);
  }

  coefs += trg0 + trgE; // From now, coefs includes also trg0 + trgE !

  for(i = 0; i < extLeft; i++) {
    memcpy(feaMx + trg_vec_size * i,
           feaMx + trg_vec_size * extLeft,
           (coefs * (1+loSrcTgzDerOrd) - trgN) * sizeof(FLOAT));
  }
  for(i = totFrames - extRight; i < totFrames; i++) {
    memcpy(feaMx + trg_vec_size * i,
           feaMx + trg_vec_size * (totFrames - extRight - 1),
           (coefs * (1+loSrcTgzDerOrd) - trgN) * sizeof(FLOAT));
  }
  // Compute missing derivatives
  for(; srcDerivOrder < derivOrder; srcDerivOrder++) {
    int winLen = derivWinLen[srcDerivOrder];
    FLOAT norm = 0.0;
    for(k = 1; k <= winLen; k++) {
      norm += 2 * k * k;
    }
    for(i=0; i < totFrames; i++) {        // for each frame
      for(j=0; j < coefs; j++) {          // for each coefficient
        FLOAT *src = feaMx + i*trg_vec_size + srcDerivOrder*coefs - trgN + j;
        *(src + coefs) = 0.0;
        if(i < winLen || i >= totFrames-winLen) {
          for(k = 1; k <= winLen; k++) {  // boundaries need special treatment
            *(src+coefs) += k*(src[ LOWER_OF(totFrames-1-i,k)*trg_vec_size]
                              -src[-LOWER_OF(i,            k)*trg_vec_size]);
          }
        } else {
          for(k = 1; k <= winLen; k++) {  // otherwice use more efficient code
            *(src+coefs) += k*(src[ k * trg_vec_size]
                              -src[-k * trg_vec_size]);
          }
        }
        *(src + coefs) /= norm;
      }
    }
  }
  header->nSamples = totFrames;
  header->sampSize = trg_vec_size * sizeof(FLOAT_32);
  header->sampKind = targetKind & ~(PARAMKIND_D | PARAMKIND_A | PARAMKIND_T);

  /////////////// Cepstral mean and variance normalization ///////////////////

  if(cmn_file != NULL) {
    ReadCepsNormFile(cmn_file, &buff->last_cmn_file, &buff->cmn,
                     header->sampKind & ~PARAMKIND_Z, CNF_Mean, coefs);
    for(i=0; i < totFrames; i++) {
      for(j=trgN; j < coefs; j++) {
        feaMx[i*trg_vec_size + j - trgN] -= buff->cmn[j];
      }
    }
  }

  header->sampKind |= derivOrder==3 ? PARAMKIND_D | PARAMKIND_A | PARAMKIND_T :
                      derivOrder==2 ? PARAMKIND_D | PARAMKIND_A :
                      derivOrder==1 ? PARAMKIND_D : 0;

  if(cvn_file != NULL) {
    ReadCepsNormFile(cvn_file, &buff->last_cvn_file, &buff->cvn,
                     header->sampKind, CNF_Variance, trg_vec_size);
    for(i=0; i < totFrames; i++) {
      for(j=trgN; j < trg_vec_size; j++) {
        feaMx[i*trg_vec_size + j - trgN] *= buff->cvn[j];
      }
    }
  }
  if(cvg_file != NULL) {
    ReadCepsNormFile(cvg_file, &buff->last_cvg_file, &buff->cvg,
                     -1, CNF_VarScale, trg_vec_size);
    for(i=0; i < totFrames; i++) {
      for(j=trgN; j < trg_vec_size; j++) {
        feaMx[i*trg_vec_size + j - trgN] *= buff->cvg[j];
      }
    }
  }
  return feaMx;
}

void ReadCepsNormFile(const char *file, char **last_file, FLOAT **vec_buff,
                      int sampKind, enum CNFileType type, int coefs)
{
  FILE *fp;
  int i;
  char s1[9], s2[64];
  char *typeStr = (char*) (type == CNF_Mean     ? "MEAN" :
                  type == CNF_Variance ? "VARIANCE" : "VARSCALE");

  char *typeStr2 = (char*) (type == CNF_Mean     ? "CMN" :
                   type == CNF_Variance ? "CVN" : "VarScale");

  if(*last_file != NULL && !strcmp(*last_file, file)) {
    return;
  }
  free(*last_file);
  *last_file=strdup(file);
  *vec_buff = (float*) realloc(*vec_buff, coefs * sizeof(FLOAT));

  if(*last_file == NULL || *vec_buff== NULL) {
    Error("Insufficient memory");
  }
  if((fp = fopen(file, "r")) == NULL) {
    Error("Cannot open %s file: '%s'", typeStr2, file);
  }
  if((type != CNF_VarScale
  && (fscanf(fp, " <%64[^>]> <%64[^>]>", s1, s2) != 2
      || strcmp(strtoupper(s1), "CEPSNORM")
      || ReadParmKind(s2, FALSE) != sampKind))
  || fscanf(fp, " <%64[^>]> %d", s1, &i) != 2
  || strcmp(strtoupper(s1), typeStr)
  || i != coefs) {
    ParmKind2Str(sampKind, s2);
    Error("%s%s%s <%s> %d ... expected in %s file %s",
          type == CNF_VarScale ? "" : "<CEPSNORM> <",
          type == CNF_VarScale ? "" : s2,
          type == CNF_VarScale ? "" : ">",
          typeStr, coefs, typeStr2, file);
  }
  for(i = 0; i < coefs; i++) {
    if(fscanf(fp, " "FLOAT_FMT, *vec_buff+i) != 1) {
      if(fscanf(fp, "%64s", s2) == 1) {
        Error("Decimal number expected but '%s' found in %s file %s",
              s2, typeStr2, file);
      } else if(feof(fp)) {
        Error("Unexpected end of %s file %s", typeStr2, file);
      } else {
        Error("Cannot read %s file %s", typeStr2, file);
      }
    }
    if(type == CNF_Variance)      (*vec_buff)[i] = 1 / sqrt((*vec_buff)[i]);
    else if(type == CNF_VarScale) (*vec_buff)[i] =     sqrt((*vec_buff)[i]);
  }
  if(fscanf(fp, "%64s", s2) == 1) {
    Error("End of file expected but '%s' found in %s file %s",
          s2, typeStr2, file);
  }
}
