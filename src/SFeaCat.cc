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

#define VERSION "0.7 "__TIME__" "__DATE__
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <malloc.h>
#include <assert.h>
#include "STKLib/fileio.h"
#include "STKLib/common.h"
#include "STKLib/Models.h"
#include "STKLib/Viterbi.h"
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
  ModelSet                hset;
  FILE *                  sfp;
  FILE *                  ofp = NULL;
  FLOAT *                 obs_mx;
  FLOAT *                 obs;
  HtkHeader              header;
  int                     i;
  int                     fcnt = 0;
  XformInstance *         p_input = NULL;
  char                    p_line[1024];
  char                    p_out_file[1024];
  MyHSearchData  cfg_hash;
  int                     tot_frames = 0;
  int                     vec_size;
  int                     out_size;
  int                     time;
  FileListElem *          feature_files = NULL;
  int                     nfeature_files = 0;
  FileListElem *          file_name = NULL;
  FileListElem **         last_file = &feature_files;


  const char *            src_hmm_list;
  const char *            src_hmm_dir;
  const char *            src_hmm_ext;
        char *            src_mmf;
  const char *            out_dir;
  const char *            out_ext;
  const char *            input_name;
        char *            script;
        char *            cmn_path;
        char *            cmn_file;
  const char *            cmn_mask;
        char *            cvn_path;
        char *            cvn_file;
  const char *            cvn_mask;
  const char *            cvg_file;
  int                     trace_flag;
  int                     targetKind;
  int                     derivOrder;
  int *                   derivWinLengths;
  int                     startFrmExt;
  int                     endFrmExt;
  
  bool                    swap_features;
  bool                    swap_fea_out;
  RHFBuffer               rhfbuff = {0};

  
  if (argc == 1) 
    usage(argv[0]);

  //InitHMMSet(&hset, 1);
  hset.Init(MODEL_SET_WITH_ACCUM);

  if (!my_hcreate_r(100,  &cfg_hash))
    Error("Insufficient memory");
  
  
  i = ParseOptions(argc, argv, optionStr, SNAME, &cfg_hash);
//  htk_compat   = GetParamBool(&cfgHash,SNAME":HTKCOMPAT",       FALSE);
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
  script_filter  = GetParamStr(&cfg_hash, SNAME":HSCRIPTFILTER",   NULL);
  parm_filter    = GetParamStr(&cfg_hash, SNAME":HPARMFILTER",     NULL);
  gpHListFilter  = GetParamStr(&cfg_hash, SNAME":HMMLISTFILTER",   NULL);
  MMF_filter     = GetParamStr(&cfg_hash, SNAME":HMMDEFFILTER",    NULL);
  parm_ofilter   = GetParamStr(&cfg_hash, SNAME":HPARMOFILTER",    NULL);
  out_dir        = GetParamStr(&cfg_hash, SNAME":TARGETPARAMDIR",  NULL);
  out_ext        = GetParamStr(&cfg_hash, SNAME":TARGETPARAMEXT",  NULL);
  trace_flag     = GetParamInt(&cfg_hash, SNAME":TRACE",           0);
  script  = (char*)GetParamStr(&cfg_hash, SNAME":SCRIPT",          NULL);
  src_hmm_list   = GetParamStr(&cfg_hash, SNAME":SOURCEHMMLIST",   NULL);
  src_hmm_dir    = GetParamStr(&cfg_hash, SNAME":SOURCEMODELDIR",  NULL);
  src_hmm_ext    = GetParamStr(&cfg_hash, SNAME":SOURCEMODELEXT",  NULL);
  src_mmf = (char*)GetParamStr(&cfg_hash, SNAME":SOURCEMMF",       NULL);

  if (GetParamBool(&cfg_hash, SNAME":PRINTCONFIG", FALSE))
    PrintConfig(&cfg_hash);
  
  if (GetParamBool(&cfg_hash, SNAME":PRINTVERSION", FALSE))
    puts("Version: "VERSION"\n");

  if (!GetParamBool(&cfg_hash,SNAME":ACCEPTUNUSEDPARAM", FALSE))
    CheckCommandLineParamUse(&cfg_hash);

  for (script=strtok(script, ","); script != NULL; script=strtok(NULL, ",")) 
  {
    if ((sfp = my_fopen(script, "rt", script_filter)) == NULL) 
      Error("Cannot open script file %s", script);
    
    while (fscanf(sfp, "%s", p_line) == 1) 
    {
      last_file = AddFileElem(last_file, p_line);
      nfeature_files++;
    }
    
    my_fclose(sfp);
  }
  
  for (src_mmf=strtok(src_mmf, ","); src_mmf != NULL; src_mmf=strtok(NULL, ",")) 
  {
    hset.ParseMmf(src_mmf, NULL);
  }
  
  if (src_hmm_list) 
    hset.ReadHMMList(src_hmm_list, src_hmm_dir, src_hmm_ext);

  if (input_name != NULL) 
  {
    Macro *macro = FindMacro(&hset.mXformInstanceHash, input_name);
    if (macro == NULL) Error("Undefined source input '%s'", input_name);
    p_input = (XformInstance *) macro->mpData;
  } 
  else if (hset.mpInputXform) 
  {
    p_input = hset.mpInputXform;
  }

  for (file_name=feature_files; file_name != NULL; file_name=file_name->mpNext) 
  {
    if (trace_flag & 1) TraceLog("Processing file %d/%d '%s'",
                                ++fcnt, nfeature_files, file_name->logical);
                                
    if(cmn_mask) 
      process_mask(file_name->logical, cmn_mask, cmn_file);
      
    if(cvn_mask) 
      process_mask(file_name->logical, cvn_mask, cvn_file);
      
    obs_mx = ReadHTKFeatures(file_name->mpPhysical, swap_features,
                            startFrmExt, endFrmExt, targetKind,
                            derivOrder, derivWinLengths, &header,
                            cmn_path, cvn_path, cvg_file, &rhfbuff);

    vec_size = header.mSampleSize / sizeof(float);
    out_size = p_input ? p_input->mOutSize : vec_size;

    if (hset.mInputVectorSize != -1 && hset.mInputVectorSize != vec_size) 
    {
      Error("Vector size [%d] in '%s' is incompatible with HMM set [%d]",
            header.mSampleSize/sizeof(float), file_name->mpPhysical, hset.mInputVectorSize);
    }
    
    MakeFileName(p_out_file, file_name->logical, out_dir, out_ext);

    if ((ofp = fopen(p_out_file, "wb")) == NULL)
      Error("Cannot open output feature file: '%s'", p_out_file);
    
    header.mSampleKind = p_input ? PARAMKIND_USER : targetKind;
    header.mSampleSize = out_size * sizeof(float);
    header.mNSamples -= hset.mTotalDelay;

    if (WriteHTKHeader(ofp, header, 1))
      Error("Cannot write to output feature file: '%s'", p_out_file);
    
    time = -hset.mTotalDelay;
    hset.ResetXformInstances();

    for(i = 0; i < header.mNSamples + hset.mTotalDelay; i++) 
    {
      hset.UpdateStacks(obs_mx + i * vec_size, ++time, FORWARD);
      
      if (time <= 0) 
        continue;

      obs = XformPass(p_input, obs_mx + i * vec_size, time, FORWARD);

      if (WriteHTKFeature (ofp, obs, out_size, swap_fea_out)) 
        Error("Cannot write to output feature file: '%s'", p_out_file);      
    }
    
    tot_frames += header.mNSamples;

    if(trace_flag & 1) 
      TraceLog("[%d frames]", header.mNSamples);
      
    fclose(ofp);
    free(obs_mx);
  }
  
  if (trace_flag & 2) 
    TraceLog("Total number of frames: %d", tot_frames);

  hset.Release();
  free(derivWinLengths);

  while (feature_files) 
  {
    file_name = feature_files;
    feature_files = feature_files->mpNext;
    free(file_name);
  }
  
  return 0;
}
