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

#define MODULE_VERSION "2.0.7 "__TIME__" "__DATE__" "SVN_ID  


#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <malloc.h>
#include <assert.h>

#include "STKLib/fileio.h"
#include "STKLib/common.h"
#include "STKLib/Models.h"
#include "STKLib/Decoder.h"
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
"\n%s version " MODULE_VERSION "\n"
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
"\n", progname, progname, progname);
  exit(-1);
}

#define SNAME "SFEACAT"
const char *optionStr =
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
  FILE*                   sfp;
  FILE*                   ofp = NULL;
  OStkStream              o_stream;
  
#ifndef USE_NEW_MATRIX  
  FLOAT*                  obs_mx;
#else
  Matrix<FLOAT>           feature_matrix;
#endif
  
  FLOAT*                  obs;
  HtkHeader               header;
  HtkHeaderExt            header_ext;
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
  const char *            mmf_dir;
  const char *            mmf_mask;

  int                     trace_flag;
  int                     targetKind;
  int                     derivOrder;
  int *                   derivWinLengths;
  int                     startFrmExt;
  int                     endFrmExt;
  
  bool                    swap_features;
  bool                    swap_fea_out;
  RHFBuffer               rhfbuff = {0};

  std::string             output_format;
  
  if (argc == 1) 
    usage(argv[0]);

  //InitHMMSet(&hset, 1);
  hset.Init(MODEL_SET_WITH_ACCUM);

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
  hset.mpHListFilter  = GetParamStr(&cfg_hash, SNAME":HMMLISTFILTER",   NULL);
  gpMmfFilter     = GetParamStr(&cfg_hash, SNAME":HMMDEFFILTER",    NULL);
  gpParmOFilter   = GetParamStr(&cfg_hash, SNAME":HPARMOFILTER",    NULL);
  out_dir        = GetParamStr(&cfg_hash, SNAME":TARGETPARAMDIR",  NULL);
  out_ext        = GetParamStr(&cfg_hash, SNAME":TARGETPARAMEXT",  NULL);
  trace_flag     = GetParamInt(&cfg_hash, SNAME":TRACE",           0);
  script  = (char*)GetParamStr(&cfg_hash, SNAME":SCRIPT",          NULL);
  src_hmm_list   = GetParamStr(&cfg_hash, SNAME":SOURCEHMMLIST",   NULL);
  src_hmm_dir    = GetParamStr(&cfg_hash, SNAME":SOURCEMODELDIR",  NULL);
  src_hmm_ext    = GetParamStr(&cfg_hash, SNAME":SOURCEMODELEXT",  NULL);
  src_mmf = (char*)GetParamStr(&cfg_hash, SNAME":SOURCEMMF",       NULL);

  mmf_dir        = GetParamStr(&cfg_hash, SNAME":MMFDIR",          ".");
  mmf_mask       = GetParamStr(&cfg_hash, SNAME":MMFMASK",         NULL);
  output_format  = GetParamStr(&cfg_hash, SNAME":OUTFORMAT",       "HTK");

  bool print_all_options = GetParamBool(&cfg_hash,SNAME":PRINTALLOPTIONS", false);

  if (GetParamBool(&cfg_hash, SNAME":PRINTCONFIG", false))
    PrintConfig(&cfg_hash);
  
  if (GetParamBool(&cfg_hash, SNAME":PRINTVERSION", false))
    puts("Version: "MODULE_VERSION"\n");
                                 
  if (!GetParamBool(&cfg_hash,SNAME":ACCEPTUNUSEDPARAM", false))
    CheckCommandLineParamUse(&cfg_hash);

  if (print_all_options) 
  {
    print_registered_parameters();
  }
  
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
  
  if (NULL != src_mmf)
  {
    for (src_mmf=strtok(src_mmf, ","); src_mmf != NULL; src_mmf=strtok(NULL, ",")) 
    {
      hset.ParseMmf(src_mmf, NULL);
    }
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
                                
    /*
    if(cmn_mask) 
      process_mask(file_name->logical, cmn_mask, cmn_file);
      
    if(cvn_mask) 
      process_mask(file_name->logical, cvn_mask, cvn_file);
    */

    if(cmn_mask) 
      ProcessMask(file_name->logical, cmn_mask, cmn_file);
      
    if(cvn_mask) 
      ProcessMask(file_name->logical, cvn_mask, cvn_file);
      

    ReadHTKFeatures(file_name->mpPhysical, swap_features,
                    startFrmExt, endFrmExt, targetKind,
                    derivOrder, derivWinLengths, &header, &header_ext,
                    cmn_path, cvn_path, cvg_file, &rhfbuff, feature_matrix);
                            

    vec_size = (header.mSampleSize != -1 ? header.mSampleSize / sizeof(float) : header_ext.mSampSize);
    out_size = p_input ? p_input->OutSize() : vec_size;

/*
    printf("vec_size %d\n", vec_size);
    for(int k = 0; k < vec_size; k++)
      printf(" %f", feature_matrix[0][k]);
    printf("\n");
*/

    if (hset.mInputVectorSize != -1 && hset.mInputVectorSize != vec_size) 
    {
      Error("Vector size [%d] in '%s' is incompatible with HMM set [%d]",
            vec_size, file_name->mpPhysical, hset.mInputVectorSize);
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
    
    MakeFileName(p_out_file, file_name->logical, out_dir, out_ext);


    o_stream.open(p_out_file);
    if (!o_stream.good()) {
      Error("Cannot open output feature file: '%s'", p_out_file);
    }

    ofp = o_stream.file();
    // if ((ofp = fopen(p_out_file, "wb")) == NULL)
    //   Error("Cannot open output feature file: '%s'", p_out_file);
    
//    header.mSampleKind = p_input ? PARAMKIND_USER : targetKind;
//    header.mSampleSize = out_size * sizeof(float);
//    header.mNSamples -= hset.mTotalDelay;

    // if (WriteHTKHeader(ofp, header, 0)) {
    //   Error("Cannot write to output feature file: '%s'", p_out_file);
    // }
    // 
    // o_stream.flush();

    time = -hset.mTotalDelay;
    hset.ResetXformInstances();
    
    //FLOAT *out_fea;
    Matrix<FLOAT> out_matrix(header.mNSamples - hset.mTotalDelay, out_size);

    FLOAT *out_mx = (FLOAT*)  malloc((header.mNSamples - hset.mTotalDelay) * out_size * sizeof(FLOAT));

    if (out_mx == NULL) {
      Error("Insufficient memory");
    }

    for(i = 0; i < header.mNSamples; i++) {
      hset.UpdateStacks(feature_matrix[i], ++time, FORWARD);

      if (time <= 0) {
        continue;
      }

      obs = XformPass(p_input, feature_matrix[i], time, FORWARD);
      // memcpy(out_mx + (time - 1) * out_size, obs, sizeof(FLOAT) * out_size);
      memcpy(out_matrix[time - 1], obs, sizeof(FLOAT) * out_size);

//      if (WriteHTKFeature (ofp, obs, out_size, swap_fea_out)) 
//        Error("Cannot write to output feature file: '%s'", p_out_file);      
    }
    
    if (output_format == "HTK") {
      // resolve the _ANON HTK parameter
      int this_target_kind = targetKind;

      if (targetKind == PARAMKIND_ANON) {
        this_target_kind = header.mSampleKind;
      } 
      else if ((targetKind & 077) == PARAMKIND_ANON) {
        this_target_kind &= ~077;
        this_target_kind |= header.mSampleKind & 077;
      }
    
      if (WriteHTKFeatures(
            ofp,
            header.mSamplePeriod,
            p_input ? (PARAMKIND_USER | (this_target_kind & PARAMKIND_C)) : this_target_kind,
            swap_fea_out, 
            out_matrix)) 
      {
        Error("Cannot write to output feature file: '%s'", p_out_file);
      }
    }
    else if (output_format == "RAWASCII") {
      o_stream << out_matrix << std::endl;
    } 
    else {
      Error("Unknown output format: '%s'", output_format.c_str());
    }

    free(out_mx);
    tot_frames += header.mNSamples - hset.mTotalDelay;

    if(trace_flag & 1) 
      TraceLog("[%d frames]", header.mNSamples - hset.mTotalDelay);
      
    o_stream.close();
    feature_matrix.Destroy();
  }
  
  if (trace_flag & 2) {
    TraceLog("Total number of frames: %d", tot_frames);
  }

  hset.Release();
  free(derivWinLengths);

  while (feature_files) {
    file_name = feature_files;
    feature_files = feature_files->mpNext;
    free(file_name);
  }
  
  for (unsigned int i = 0; i < cfg_hash.mNEntries; i++) {
    free(cfg_hash.mpEntry[i]->data);
  }

  my_hdestroy_r(&cfg_hash, 1);

  free(rhfbuff.mpLastFileName);

  return 0;
}
