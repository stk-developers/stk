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

#ifndef FILEIO_H
#define FILEIO_H

#include "Matrix.h"
#include "common.h"

#include <string>
#include <cstdio>

//using namespace std;
namespace STK
{
  /**
   * Structure for HTK header
   */
  struct HtkHeader
  {
    INT_32    mNSamples;              
    INT_32    mSamplePeriod;
    INT_16    mSampleSize;
    UINT_16   mSampleKind;
  };


  /** 
   * @brief Extension of the HTK header
   */
  struct HtkHeaderExt
  {
    INT_32 mHeaderSize;
    INT_32 mVersion;
    INT_32 mSampSize;
  };  


  struct RHFBuffer
  {
    char*     mpLastFileName;
    char*     mpLastCmnFile;
    char*     mpLastCvnFile;
    char*     mpLastCvgFile;
    FILE*     mpFp;
    FLOAT*    cmn;
    FLOAT*    cvn;
    FLOAT*    cvg;
    HtkHeader last_header;
    FLOAT*    A;
    FLOAT*    B;
  } ;
  
  
  enum CNFileType 
  {
    CNF_Mean, 
    CNF_Variance, 
    CNF_VarScale
  };
  
  //***************************************************************************
  //***************************************************************************
  void ReadCepsNormFile(
    const char*   pFileName, 
    char**        lastFile, 
    FLOAT**       vecBuff,
    int           sampleKind, 
    CNFileType    type, 
    int           coefs);
  
  
  int WriteHTKHeader  (FILE * fp_out, HtkHeader header, bool swap);
  int WriteHTKFeature (FILE * fp_out, FLOAT *out, size_t fea_len, bool swap, bool compress, FLOAT* pScale, FLOAT* pBias);
  int WriteHTKFeatures(FILE * pOutFp, FLOAT * pOut, int nCoeffs, int nSamples, int samplePeriod, int targetKind, bool swap);
  int ReadHTKHeader   (FILE * fp_in, HtkHeader *header, bool swap);
  int ReadHTKFeature  (FILE * fp_in, FLOAT *in, size_t fea_len, bool swap,
                       bool   decompress, FLOAT *A, FLOAT *B);
  
  int Mkdir4File(const char * file_name);
  
  int 
  WriteHTKFeatures(
    FILE *  pOutFp,
    int     samplePeriod,
    int     targetKind,  
    bool    swap,
    Matrix<FLOAT>&        rFeatureMatrix);
  
  FLOAT*
  ReadHTKFeatures(
    char*         file_name,
    bool          swap,
    int           extLeft,
    int           extRight,
    int           targetKind,
    int           derivOrder,
    int*          derivWinLen,
    HtkHeader*    header,
    const char*   cmn_file,
    const char*   cvn_file,
    const char*   cvg_file,
    RHFBuffer*    buff);

  bool 
  ReadHTKFeatures(
    const char*           pFileName,
    bool                  swap,
    int                   extLeft,
    int                   extRight,
    int                   targetKind,
    int                   derivOrder,
    int*                  derivWinLen,
    HtkHeader*            pHeader,
    const char*           pCmnFile,
    const char*           pCvnFile,
    const char*           pCvgFile,
    RHFBuffer*            pBuff,
    Matrix<FLOAT>&        rFeatureMatrix);
      
}; // namespace STK  

#endif // FILEIO_H
