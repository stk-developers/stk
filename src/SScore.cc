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

#define MODULE_VERSION "0.7 "__TIME__" "__DATE__
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

#include "STKLib/fileio.h"
#include "STKLib/common.h"
#include "STKLib/Models.h"
#include "STKLib/Decoder.h"
#include "STKLib/Score.h"



using namespace STK;


void usage(char *progname)
{
  char *tchrptr;
  if ((tchrptr = strrchr(progname, '\\')) != NULL) progname = tchrptr+1;
  if ((tchrptr = strrchr(progname, '/')) != NULL) progname = tchrptr+1;
  fprintf(stderr,
"\nUSAGE: %s [options] DataFiles...\n\n"
" Option                                                     Default\n\n"
//" -c f   Mixture pruning threshold                         10.0\n"
" -d s       Dir to find hmm definitions                     Current\n"
" -l s       Dir to store output param files                 Feature file dir\n"
" -x s       Extension for hmm files                         None\n"
" -y s       Output param file extension                     Unchanged\n"
" -A         Print command line arguments                    Off\n"
" -C cf      Set config file to cf                           Default\n"
" -D         Display configuration variables                 Off\n"
" -H mmf     Load HMM macro file mmf                         \n"
" -S f       Set script file to f                            None\n"
" -T N       Set trace flags to N                            0\n"
" -V         Print version information                       Off\n"
"\n"
" %s is Copyright (C) 2004-2005 Lukas Burget et al. and\n"
" licensed under the GNU General Public License, version 2.\n"
" Bug reports, feedback, etc, to: burget@fit.vutbr.cz\n"
"\n", progname, progname);
  exit(-1);
}

#define SNAME "SFEACAT"
char *optionStr =
" -d r   SOURCEMODELDIR"
" -l r   TARGETPARAMDIR"
" -x r   SOURCEMODELEXT"
" -y r   TARGETPARAMEXT"
" -D n   PRINTCONFIG=TRUE"
" -H l   SOURCEMMF"
" -S l   SCRIPT"
" -T r   TRACE"
" -V n   PRINTVERSION=TRUE";

int main(int argc, char *argv[]) 
{
  ModelSet                ubm;
  Hmm*                    ubm_hmm;
  Hmm*                    spk_hmm;
  
  FILE*                   sfp;
  FILE*                   ofp = NULL;
  
  Matrix<FLOAT>           feature_matrix;
  
  
  FLOAT*                  obs;
  HtkHeader               header;
  
  GradientSV              grad_sv;
  Matrix<FLOAT>           model_diff;
  OStkStream              o_stream;
  
  
  int                     i;
  int                     fcnt = 0;
  XformInstance*          p_input = NULL;
  char                    p_line[1024];
  char                    p_out_file[1024];
  MyHSearchData  cfg_hash;
  int                     tot_frames = 0;
  int                     vec_size;
  int                     out_size;
  int                     time;
  FileListElem*           feature_files = NULL;
  int                     nfeature_files = 0;
  FileListElem*           file_name = NULL;
  FileListElem**          last_file = &feature_files;


  const char*             src_hmm_list;
  const char*             src_hmm_dir;
  const char*             src_hmm_ext;
        char*             src_mmf;
  const char*             out_dir;
  const char*             out_ext;
  const char*             input_name;
        char*             script;
        char*             cmn_path;
        char*             cmn_file;
  const char*             cmn_mask;
        char*             cvn_path;
        char*             cvn_file;
  const char*             cvn_mask;
  const char*             cvg_file;
  const char*             mmf_dir;
  const char*             mmf_mask;

  int                     trace_flag;
  int                     targetKind;
  int                     derivOrder;
  int*                    derivWinLengths;
  int                     startFrmExt;
  int                     endFrmExt;
  
  bool                    swap_features;
  bool                    swap_fea_out;
  RHFBuffer               rhfbuff = {0};

  const char*             score_type;
        char*             ubm_mmf;
  
  if (argc == 1) 
    usage(argv[0]);

  ubm.Init(MODEL_SET_WITH_ACCUM);

  if (!my_hcreate_r(100,  &cfg_hash))
    Error("Insufficient memory");
  
  i = ParseOptions(argc, argv, optionStr, SNAME, &cfg_hash);
  
  //  htk_compat   = GetParamBool(&cfgHash,SNAME":HTKCOMPAT",       false);
  for (; i < argc; i++) 
  {
    last_file = AddFileElem(last_file, argv[i]);
    nfeature_files++;
  }
  
  targetKind     = GetDerivParams(&cfg_hash, &derivOrder, &derivWinLengths,
                                &startFrmExt, &endFrmExt,
                                &cmn_path, &cmn_file, &cmn_mask,
                                &cvn_path, &cvn_file, &cvn_mask, &cvg_file,
                                SNAME":", 0);

  input_name     = GetParamStr(&cfg_hash, SNAME":SOURCEINPUT",     NULL);
  swap_features  =!GetParamBool(&cfg_hash,SNAME":NATURALREADORDER", isBigEndian());
  swap_fea_out   =!GetParamBool(&cfg_hash,SNAME":NATURALWRITEORDER",isBigEndian());
  gpFilterWldcrd = GetParamStr(&cfg_hash, SNAME":HFILTERWILDCARD", "$");
  gpScriptFilter  = GetParamStr(&cfg_hash, SNAME":HSCRIPTFILTER",   NULL);
  gpParmFilter    = GetParamStr(&cfg_hash, SNAME":HPARMFILTER",     NULL);
  gpHListFilter  = GetParamStr(&cfg_hash, SNAME":HMMLISTFILTER",   NULL);
  gpMmfFilter     = GetParamStr(&cfg_hash, SNAME":HMMDEFFILTER",    NULL);
  gpParmOFilter   = GetParamStr(&cfg_hash, SNAME":HPARMOFILTER",    NULL);
  out_dir        = GetParamStr(&cfg_hash, SNAME":TARGETPARAMDIR",  NULL);
  out_ext        = GetParamStr(&cfg_hash, SNAME":TARGETPARAMEXT",  "gsv");
  trace_flag     = GetParamInt(&cfg_hash, SNAME":TRACE",           0);
  script  = (char*)GetParamStr(&cfg_hash, SNAME":SCRIPT",          NULL);
  src_hmm_list   = GetParamStr(&cfg_hash, SNAME":SOURCEHMMLIST",   NULL);
  src_hmm_dir    = GetParamStr(&cfg_hash, SNAME":SOURCEMODELDIR",  NULL);
  src_hmm_ext    = GetParamStr(&cfg_hash, SNAME":SOURCEMODELEXT",  NULL);
  ubm_mmf = (char*)GetParamStr(&cfg_hash, SNAME":SOURCEMMF",       NULL);

  mmf_dir        = GetParamStr(&cfg_hash, SNAME":MMFDIR",          ".");
  mmf_mask       = GetParamStr(&cfg_hash, SNAME":MMFMASK",         NULL);
  
  
  score_type     = GetParamStr(&cfg_hash, SNAME":SCORETYPE",       NULL);
  

  if (GetParamBool(&cfg_hash, SNAME":PRINTCONFIG", false))
    PrintConfig(&cfg_hash);
  
  if (GetParamBool(&cfg_hash, SNAME":PRINTVERSION", false))
    puts("Version: "MODULE_VERSION"\n");
                                 
  if (!GetParamBool(&cfg_hash,SNAME":ACCEPTUNUSEDPARAM", false))
    CheckCommandLineParamUse(&cfg_hash);

  
  if (NULL != script)
  {
    for (script=strtok(script, ","); script != NULL; script=strtok(NULL, ",")) 
    {
      if ((sfp = my_fopen(script, "rt", gpScriptFilter)) == NULL) 
        Error("Cannot open script file %s", script);
      
      while (fscanf(sfp, "%s", p_line) == 1) 
      {
        last_file = AddFileElem(last_file, p_line);
        nfeature_files++;
      }
      
      my_fclose(sfp);
    }
  }  
  
  // UBM has to be specified
  if (NULL == ubm_mmf)
  {
    Error("UBM not specified");
  }
  
  // UBM
  if (NULL != ubm_mmf)
  {
    for (ubm_mmf=strtok(ubm_mmf, ","); ubm_mmf != NULL; ubm_mmf=strtok(NULL, ",")) 
    {
      ubm.ParseMmf(ubm_mmf, NULL);
    }
  }
  
  if (src_hmm_list)
  {
    ubm.ReadHMMList(src_hmm_list, src_hmm_dir, src_hmm_ext);
  }
  

  if (input_name != NULL) 
  {
    Macro *macro = FindMacro(&ubm.mXformInstanceHash, input_name);
    if (macro == NULL) Error("Undefined source input '%s'", input_name);
    p_input = (XformInstance *) macro->mpData;
  } 
  else if (ubm.mpInputXform) 
  {
    p_input = ubm.mpInputXform;
  }


  
  ubm_hmm = static_cast<Hmm*>(static_cast<Macro*>(ubm.mHmmHash.mpEntry[0]->data)->mpData);


  grad_sv.mMatrix.Init(ubm_hmm->mpState[0]->mNMixtures,
                       ubm.mpInputXform != NULL ? ubm.mpInputXform->OutSize() 
                                                : ubm.mInputVectorSize);
                                           

                      
  for (file_name=feature_files; file_name != NULL; file_name=file_name->mpNext) 
  {
    if (trace_flag & 1) 
    {
      TraceLog("Processing file %d/%d '%s'",
               ++fcnt, nfeature_files, file_name->logical);
    }
                                
    if(cmn_mask) 
      process_mask(file_name->logical, cmn_mask, cmn_file);
      
    if(cvn_mask) 
      process_mask(file_name->logical, cvn_mask, cvn_file);
      
    ReadHTKFeatures(file_name->mpPhysical, swap_features,
                    startFrmExt, endFrmExt, targetKind,
                    derivOrder, derivWinLengths, &header,
                    cmn_path, cvn_path, cvg_file, &rhfbuff, feature_matrix);
                            
    vec_size = header.mSampleSize / sizeof(float);
    out_size = p_input ? p_input->OutSize() : vec_size;

    if (ubm.mInputVectorSize != -1 && ubm.mInputVectorSize != vec_size) 
    {
      Error("Vector size [%d] in '%s' is incompatible with HMM set [%d]",
            header.mSampleSize/sizeof(float), file_name->mpPhysical, ubm.mInputVectorSize);
    }
    
    if(mmf_mask != NULL) 
    {
      static string lastSpeakerMMF;
      string speakerMMF;
      ProcessMask(file_name->logical, mmf_mask, speakerMMF);
        
      if(lastSpeakerMMF != speakerMMF) 
      {
        ubm.ParseMmf((string(mmf_dir) + "/" + speakerMMF).c_str(), NULL);
        lastSpeakerMMF = speakerMMF;
      }
    }
    
    // construct new file name
    MakeFileName(p_out_file, file_name->logical, out_dir, "gsv");
    
    // save old matrix if new needs to be accumulated
    if (grad_sv.mFileName.empty())
    {
      grad_sv.mFileName = p_out_file;
      
    }
    else if (grad_sv.mFileName != p_out_file)
    {
      o_stream.open(grad_sv.mFileName.c_str(), std::ios::out, gpParmOFilter);    
         
      if (!o_stream.good())
      {
        Error("Cannot open output file: '%s'", grad_sv.mFileName.c_str());
      }
      
      o_stream << grad_sv.mMatrix;
      o_stream.close();
      
      grad_sv.mMatrix.Clear();
      grad_sv.mFileName = p_out_file;
    }

    header.mNSamples -= ubm.mTotalDelay;

    time = -ubm.mTotalDelay;
    ubm.ResetXformInstances();
    
    for(i = 0; i < header.mNSamples + ubm.mTotalDelay; i++) 
    {
      ubm.UpdateStacks(feature_matrix[i], ++time, FORWARD);
      if (time <= 0) 
        continue;

      obs = XformPass(p_input, feature_matrix[i], time, FORWARD);

      gradient_supervector_accum(grad_sv.mMatrix, obs, *ubm_hmm);
    }
    
    tot_frames += header.mNSamples;

    if(trace_flag & 1) 
      TraceLog("[%d frames]", header.mNSamples);
      
    feature_matrix.Destroy();
  }
  
  // flush the matrix to a file
  if (!grad_sv.mFileName.empty())
  {
    o_stream.open(grad_sv.mFileName.c_str(), std::ios::out, gpParmOFilter);    
    
    if (!o_stream.good())
    {
      Error("Cannot open output file: '%s'", grad_sv.mFileName.c_str());
    }
    
    o_stream << grad_sv.mMatrix;
    o_stream.close();
  }
  
  
  if (trace_flag & 2) 
    TraceLog("Total number of frames: %d", tot_frames);

  ubm.Release();
  free(derivWinLengths);

  while (feature_files) 
  {
    file_name = feature_files;
    feature_files = feature_files->mpNext;
    free(file_name);
  }
  
  return 0;
}
