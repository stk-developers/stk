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

//#define VERSION "0.4 "__TIME__" "__DATE__

//#define HAVE_REENTRANT_SEARCH 1
//#define HAVE_SEARCH_H 1

#include "STKLib/Decoder.h"
#include "STKLib/Models.h"
#include "STKLib/fileio.h"
#include "STKLib/labels.h"
#include "STKLib/common.h"
//#include "STKLib/stkstream.h"
#include <stdlib.h>
#include <string.h>
#include <cstdio>
#include <malloc.h>
#include <vector>
#ifndef WIN32
#include <unistd.h>
#else
#include "getopt.h"
#endif

using namespace STK;


void usage(char *progname)
{
  char *tchrptr;
  if ((tchrptr = strrchr(progname, '\\')) != NULL) progname = tchrptr+1;
  if ((tchrptr = strrchr(progname, '/')) != NULL) progname = tchrptr+1;
  fprintf(stderr,
"\nUSAGE: %s [options] DataFiles...\n\n"
" Option                                                     Default\n\n"
" -d s       Dir to find hmm definitions                     Current\n"
" -i s       Output transcriptions to MLF s                  Off\n"
" -l s       Dir to store transcription files                Feature file dir\n"
" -x s       Extension for hmm files                         None\n"
" -y s       Output transcription file extension             rec\n"
" -A         Print command line arguments                    Off\n"
" -C cf      Set config file to cf                           Default\n"
" -D         Display configuration variables                 Off\n"
" -H mmf     Load HMM macro file mmf\n"
" -I mlf     Load master transcription file mlf\n"
" -L dir     Set input transcription dir                     Current\n"
" -S f       Set script file to f                            None\n"
" -T N       Set trace flags to N                            0\n"
" -V         Print version information                       Off\n"
" -X ext     Set input transcription file ext                lab\n"
"\n"
" %s is Copyright (C) 2004-2005 Lukas Burget et al. and\n"
" licensed under the GNU General Public License, version 2.\n"
" Bug reports, feedback, etc, to: burget@fit.vutbr.cz\n"
"\n", progname, progname);
  exit(-1);
}

#define SNAME "SVITE"
char *optionStr =
" -d r   SOURCEMODELDIR"
" -i r   TARGETMLF"
" -l r   TARGETTRANSCDIR"
" -m n   MODELALIGNMENT=TRUE"
" -x r   SOURCEMODELEXT"
" -y r   TARGETTRANSCEXT"
" -D n   PRINTCONFIG=TRUE"
" -H l   SOURCEMMF"
" -I r   SOURCEMLF"
" -L r   SOURCETRANSCDIR"
" -S l   SCRIPT"
" -T r   TRACE"
" -V n   PRINTVERSION=TRUE"
" -X r   SOURCETRANSCEXT";

struct ModelRecord
{ 
  ModelRecord(Hmm* h) : hmm(h), score(0.0) {}
  Hmm*  hmm;
  FLOAT score;
};
  
struct KBest
{
  FLOAT  post;
  size_t index;
};

int main(int argc, char *argv[]) 
{
  HtkHeader                    header;
  ModelSet                      hset;
  FILE *                        sfp;
  FILE *                        lfp = NULL;
  FILE *                        ilfp = NULL;
#ifndef USE_NEW_MATRIX  
  FLOAT*                        obsMx;
#else
  Matrix<FLOAT>                 feature_matrix;
#endif  
  size_t                        i;
  int                           j;
  int                           fcnt = 0;
  Label *                       labels = NULL;
  char                          line[1024];
  char                          label_file[1024];
  int                           nfeature_files = 0;
  
  MyHSearchData        modelHash;
  MyHSearchData        cfgHash;
  
  FileListElem *                feature_files = NULL;
  FileListElem *                file_name = NULL;
  FileListElem **               last_file = &feature_files;
  
//double  word_pruning;
         
  const  char *                 hmm_dir;
  const  char *                 hmm_ext;
  const  char *                 out_lbl_dir;
  const char *                  out_lbl_ext;
  const char *                  in_lbl_dir;
  const char *                  in_lbl_ext;
  const char *                  out_MLF;
  const char *                  in_MLF;
  const char *                  network_file;  
  const char *                  hmm_list;
  char *                        script;
  char *                        mmf;
  const char *                  label_filter;
  const char *                  label_ofilter;
        char *                  cmn_path;
        char *                  cmn_file;
  const char *                  cmn_mask;
        char *                  cvn_path;
        char *                  cvn_file;
  const char *                  cvn_mask;
  const char *                  cvg_file;
  const char *                  mmf_dir;
  const char *                  mmf_mask;

  int                           trace_flag;
  int                           targetKind;
  int                           derivOrder;
  int  *                        derivWinLengths;
  int                           startFrmExt;
  int                           endFrmExt;
  int                           nbest;
  bool                          swap_features;
  
  RHFBuffer                     rhfbuff     = {0};
  STKNetworkOutputFormat        in_lbl_fmt  = {0};
  
  std::vector<ModelRecord> hmmRecords;
  std::string lastLogicalName;
  int nFrames = 0;
  
  
  in_lbl_fmt.mNoTimes = 1;

  if (argc == 1) usage(argv[0]);

  hset.Init();

  if (!my_hcreate_r(100,  &modelHash)
  || !my_hcreate_r(100,  &cfgHash)) {
    Error("Insufficient memory");
  }
  
  j = ParseOptions(argc, argv, optionStr, SNAME, &cfgHash);
  
  for (; j < argc; j++) {
    last_file = AddFileElem(last_file, argv[j]);
    nfeature_files++;
  }
  
  targetKind =   GetDerivParams(&cfgHash, &derivOrder, &derivWinLengths,
                                &startFrmExt, &endFrmExt,
                                &cmn_path, &cmn_file, &cmn_mask,
                                &cvn_path, &cvn_file, &cvn_mask, &cvg_file,
                                SNAME":", 0);
  swap_features=!GetParamBool(&cfgHash,SNAME":NATURALREADORDER",isBigEndian());
  gpFilterWldcrd= GetParamStr(&cfgHash,SNAME":HFILTERWILDCARD", "$");
  gpScriptFilter= GetParamStr(&cfgHash,SNAME":HSCRIPTFILTER",   NULL);
  gpParmFilter  = GetParamStr(&cfgHash,SNAME":HPARMFILTER",     NULL);
  gpHListFilter = GetParamStr(&cfgHash,SNAME":HMMLISTFILTER",   NULL);
  gpMmfFilter   = GetParamStr(&cfgHash,SNAME":HMMDEFFILTER",    NULL);
  label_filter = GetParamStr(&cfgHash, SNAME":HLABELFILTER",    NULL);
  label_ofilter= GetParamStr(&cfgHash, SNAME":HLABELOFILTER",   NULL);
  hmm_list     = GetParamStr(&cfgHash, SNAME":SOURCEHMMLIST",   NULL);
  nbest        = GetParamInt(&cfgHash, SNAME":NBEST",           10);
  network_file = GetParamStr(&cfgHash, SNAME":RECOGNET",        NULL);
  hmm_dir      = GetParamStr(&cfgHash, SNAME":SOURCEMODELDIR",  NULL);
  hmm_ext      = GetParamStr(&cfgHash, SNAME":SOURCEMODELEXT",  NULL);
  in_MLF       = GetParamStr(&cfgHash, SNAME":SOURCEMLF",       NULL);
  in_lbl_dir   = GetParamStr(&cfgHash, SNAME":SOURCETRANSCDIR", NULL);
  in_lbl_ext   = GetParamStr(&cfgHash, SNAME":SOURCETRANSCEXT", NULL);
  out_MLF      = GetParamStr(&cfgHash, SNAME":TARGETMLF",       NULL);
  out_lbl_dir  = GetParamStr(&cfgHash, SNAME":TARGETTRANSCDIR", NULL);
  out_lbl_ext  = GetParamStr(&cfgHash, SNAME":TARGETTRANSCEXT", "rec");
  trace_flag   = GetParamInt(&cfgHash, SNAME":TRACE",           0);
  script =(char*)GetParamStr(&cfgHash, SNAME":SCRIPT",          NULL);
  mmf    =(char*)GetParamStr(&cfgHash, SNAME":SOURCEMMF",       NULL);
  
  mmf_dir      = GetParamStr(&cfgHash, SNAME":MMFDIR",          ".");
  mmf_mask     = GetParamStr(&cfgHash, SNAME":MMFMASK",         NULL);


  if (GetParamBool(&cfgHash, SNAME":PRINTCONFIG", false)) {
    PrintConfig(&cfgHash);
  }
  if (GetParamBool(&cfgHash, SNAME":PRINTVERSION", false)) {
    puts("Version: "VERSION"\n");
  }
  if (!GetParamBool(&cfgHash,SNAME":ACCEPTUNUSEDPARAM", false)) {
    CheckCommandLineParamUse(&cfgHash);
  }

  if (NULL != script)
  {
    for (script=strtok(script, ","); script != NULL; script=strtok(NULL, ",")) {
      if ((sfp = my_fopen(script, "rt", gpScriptFilter)) == NULL) {
        Error("Cannot open script file %s", script);
      }
      while (fscanf(sfp, "%s", line) == 1) {
        last_file = AddFileElem(last_file, line);
        nfeature_files++;
      }
      my_fclose(sfp);
    }
  }
  
  if (NULL != mmf)
  {
    for (mmf=strtok(mmf, ","); mmf != NULL; mmf=strtok(NULL, ",")) {
      hset.ParseMmf(mmf, NULL);
    }
  }
  
  if (hmm_list != NULL) hset.ReadHMMList(hmm_list, hmm_dir, hmm_ext);
  
  transc_filter  = transc_filter != NULL  ? transc_filter  : label_filter;
  transc_ofilter = transc_ofilter != NULL ? transc_ofilter : label_ofilter;

  if (network_file) 
  { 
    IStkStream input_stream;    

    input_stream.open(network_file, std::ios::in, transc_filter ? transc_filter : "");
    
    if (!input_stream.good())
    {
      Error("Cannot open network file: %s", network_file);
    }
    
    ilfp = input_stream.file();
    
    //:TODO:
    // header.mSamplePeriod not initialized yet...has no meaning in this tool anyway
    labels = ReadLabels(ilfp, &modelHash, UL_INSERT, in_lbl_fmt, 
                        header.mSamplePeriod, network_file, NULL, NULL);
  }
  else
  {
    ilfp = OpenInputMLF(in_MLF);
  }
  
  lfp = OpenOutputMLF(out_MLF);
  

  for (file_name = feature_files;; file_name = file_name->mpNext) 
  {

    if (NULL == file_name
    ||  lastLogicalName != file_name->logical) 
    {
      if(!lastLogicalName.empty())
      {
        strcpy(label_file, lastLogicalName.c_str());
        lfp = OpenOutputLabelFile(label_file, out_lbl_dir, out_lbl_ext, lfp, out_MLF);
        fprintf(lfp, "%d\n", nFrames);
        for(i = 0; i < hmmRecords.size(); i++)
        {
          fprintf(lfp, "%s %f\n", hmmRecords[i].hmm->mpMacro->mpName, hmmRecords[i].score / nFrames);
        }
        
        CloseOutputLabelFile(lfp, out_MLF);
      }
      
      if (NULL == file_name)
        break;
    
      lastLogicalName = file_name->logical;
      
      if (!network_file) 
      {

        strcpy(label_file, file_name->logical);
        ilfp = OpenInputLabelFile(label_file, in_lbl_dir,
                                  in_lbl_ext ? in_lbl_ext : "lab",
                                  ilfp, in_MLF);

        labels = ReadLabels(ilfp, &modelHash, UL_INSERT, in_lbl_fmt,
                            header.mSamplePeriod, label_file, in_MLF, NULL);
        CloseInputLabelFile(ilfp, in_MLF);
      }
      
      hmmRecords.clear();
      nFrames = 0;

      for (Label* label = labels; label != NULL; label = label->mpNext) {        
        Macro* macro = FindMacro(&hset.mHmmHash, label->mpName);

        if (NULL == macro) Error("HMM %s not found in HMM set", label->mpName);
        
        hmmRecords.push_back(ModelRecord(reinterpret_cast<Hmm*>(macro->mpData)));
      }
      
      if(hmmRecords.empty()) 
      {
        Error("No models to evaluate is specified for feature file %s", file_name->logical);
      }
      
      if (!network_file) 
        ReleaseLabels(labels);
    }
    

    if (trace_flag & 1) {
      TraceLog("Processing file %d/%d '%s'", ++fcnt,
                 nfeature_files,file_name->mpPhysical);
    }
    
    if (cmn_mask) process_mask(file_name->logical, cmn_mask, cmn_file);
    if (cvn_mask) process_mask(file_name->logical, cvn_mask, cvn_file);
    
    ReadHTKFeatures(file_name->mpPhysical, swap_features,
                    startFrmExt, endFrmExt, targetKind,
                    derivOrder, derivWinLengths, &header,
                    cmn_path, cvn_path, cvg_file, &rhfbuff, feature_matrix);

    if (hset.mInputVectorSize != static_cast<int>(header.mSampleSize / sizeof(float))) 
    {
      Error("Vector size [%d] in '%s' is incompatible with HMM set [%d]",
            header.mSampleSize/sizeof(float), file_name->mpPhysical, hset.mInputVectorSize);
    }
    
    if(mmf_mask != NULL) {
      static std::string lastSpeakerMMF;
      std::string speakerMMF;
      ProcessMask(file_name->logical, mmf_mask, speakerMMF);
        
      if(lastSpeakerMMF != speakerMMF) 
      {
        hset.ParseMmf((std::string(mmf_dir) + "/" + speakerMMF).c_str(), NULL);
        lastSpeakerMMF = speakerMMF;
      }
    }
    
    std::vector<KBest> nBestList(nbest);
    
    int t, k;
    int time = -hset.mTotalDelay;
    hset.ResetXformInstances();

    for(t = 0; t < header.mNSamples; t++) 
    {
      hset.UpdateStacks(feature_matrix[t], ++time, FORWARD);

      if (time <= 0) 
        continue;
      
      State* pState = hmmRecords[0].hmm->mpState[0];
      
      for (i = 0; i < nbest; i++) 
        nBestList[i].post = LOG_0;
        
      for (i = 0; i < pState->mNMixtures; i++)
      {
        Mixture* mix = pState->mpMixture[i].mpEstimates;
        FLOAT* l_obs = XformPass(mix->mpInputXform, feature_matrix[t], time, FORWARD);
        assert(l_obs != NULL);
        FLOAT glike = (mix->mpGaussianSelectionInputXform == NULL ||
                       XformPass(mix->mpGaussianSelectionInputXform, feature_matrix[t], time, FORWARD)[mix->mGaussianSelectionIndex-1] > 0.0)
                      ? Decoder<_DecoderNetwork>::DiagCGaussianDensity(mix, l_obs, NULL) + pState->mpMixture[i].mWeight : LOG_0;

        for(j = 0; nBestList[j].post <= glike && j < nbest; j++) {}
        
        if(j > 0)
        {
          for(k = 1; k < j; k++)
            nBestList[k-1] = nBestList[k];

          nBestList[j-1].post = glike;
          nBestList[j-1].index = i;
        }  
      }
      
      FLOAT totLike = LOG_0;
            
      for(j = nbest-1; j >= 0; j--)
        totLike = LogAdd(totLike, nBestList[j].post);
        
      for(j = nbest-1; j >= 0; j--) {
        FLOAT like = nBestList[j].post;
        nBestList[j].post = exp(like - totLike);
        hmmRecords[0].score += nBestList[j].post * like;
      }
      
      for(i = 1; i < hmmRecords.size(); i++) {
        pState = hmmRecords[i].hmm->mpState[0];
        
        for(j = nbest-1; j >= 0; j--) {
          Mixture* mix = pState->mpMixture[nBestList[j].index].mpEstimates;
          FLOAT* l_obs = XformPass(mix->mpInputXform, feature_matrix[t], time, FORWARD);
          assert(l_obs != NULL);
          
          FLOAT glike = (mix->mpGaussianSelectionInputXform == NULL ||
               XformPass(mix->mpGaussianSelectionInputXform, feature_matrix[t], time, FORWARD)[mix->mGaussianSelectionIndex-1] > 0.0)
              ? Decoder<_DecoderNetwork>::DiagCGaussianDensity(mix, l_obs, NULL) : LOG_0;

          hmmRecords[i].score += nBestList[j].post * (glike + pState->mpMixture[nBestList[j].index].mWeight);
        }        
      }
    }
    
    nFrames += header.mNSamples - hset.mTotalDelay;
    feature_matrix.Destroy();
  }
  if (network_file) 
    ReleaseLabels(labels);

  hset.Release();
  
  my_hdestroy_r(&modelHash,   1);
  
  for (unsigned int i = 0; i < cfgHash.mNEntries; i++) 
    free(cfgHash.mpEntry[i]->data);
  
  my_hdestroy_r(&cfgHash, 1);

  while (feature_files) 
  {
    file_name     = feature_files;
    feature_files = feature_files->mpNext;
    free(file_name);
  }

  return 0;
}
