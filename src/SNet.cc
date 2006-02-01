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
" -o s       Extension for new hmm files                     As src\n"
" -A         Print command line arguments                    Off\n"
" -B         Save HMM macro files as binary                  Off\n"
" -C cf      Set config file to cf                           Default\n"
" -D         Display configuration variables                 Off\n"
" -G fmt     Set source trascription format to fmt           As config\n"
" -H mmf     Load HMM macro file mmf                         \n"
" -I mlf     Load master label file mlf                      \n"
" -L dir     Set input label (or net) dir                    Current\n"
" -M dir     Dir to write HMM macro files                    Current\n"
" -S f       Set script file to f                            None\n"
" -T N       Set trace flags to N                            0\n"
" -V         Print version information                       Off\n"
" -X ext     Set input label (or netwokr) file ext           lab (net)\n"
"\n"
" %s is Copyright (C) 2004-2005 Lukas Burget et al. and\n"
" licensed under the GNU General Public License, version 2.\n"
" Bug reports, feedback, etc, to: burget@fit.vutbr.cz\n"
"\n", progname, progname);
  exit(-1);
}

#define SNAME "SFEACAT"
char *optionStr =
" -o r   TARGETMODELEXT"
" -B n   SAVEBINARY=TRUE"
" -D n   PRINTCONFIG=TRUE" 
" -H l   SOURCEMMF"
" -I r   SOURCEMLF"
" -L r   SOURCETRANSCDIR"
" -M r   TARGETMODELDIR"
" -S l   SCRIPT"
" -T r   TRACE"
" -V n   PRINTVERSION=TRUE"
" -X r   SOURCETRANSCEXT";


int main(int argc, char *argv[]) 
{
  ModelSet hset;
  FILE *sfp, *ilfp;
  FLOAT *obsMx, *obsMx_out, *obs, *obs_out;
  HtkHeader header, header_out;
  Label *labels;
  int i, fcnt = 0;
  XformInstance *input = NULL;
  char line[1024];
  char label_file[1024];
  char *chrptr;
  MyHSearchData labelHash, cfgHash;
  int totFrames = 0;
  int vec_size, out_size, time;
  FileListElem *feature_files = NULL;
  int nfeature_files = 0;
  FileListElem *file_name = NULL;
  FileListElem **last_file = &feature_files;


  const char *src_hmm_dir;
  const char *src_hmm_ext;
        char *src_mmf;
  const char *trg_hmm_dir;
  const char *trg_hmm_ext;
  const char *trg_mmf;
  const char *out_dir;
  const char *out_ext;
  const char *src_lbl_dir;
  const char *src_lbl_ext;
  const char *src_mlf;

  const char *inputName;
  char *script;
  
        char *        cmn_path;
        char *        cmn_file;
  const char *        cmn_mask;
        char *        cvn_path;
        char *        cvn_file;
  const char *        cvn_mask;
  const char *        cvg_file;
        char *        cmn_path_out;
        char *        cmn_file_out;
  const char *        cmn_mask_out;
        char *        cvn_path_out;
        char *        cvn_file_out;
  const char *        cvn_mask_out;
  const char *        cvg_file_out;
  const char *        outlabel_map;
  int  mTraceFlag;
  int  targetKind, targetKind_out;;
  int  derivOrder, derivOrder_out;
  int  *derivWinLengths, *derivWinLengths_out;
  int startFrmExt;
  int endFrmExt;
  int startFrmExt_out;
  int endFrmExt_out;
  bool hmms_binary;
  bool swap_features;
  bool swap_features_out;
  const char *NNet_instance_name = NULL;
  XformInstance *NNet_instance   = NULL;
  XformInstance *NNet_input      = NULL;
  LabelFormat in_lbl_fmt = {0};
  RHFBuffer rhfbuff      = {0};
  RHFBuffer rhfbuff_out  = {0};

  if (argc == 1) usage(argv[0]);

  hset.Init(MODEL_SET_WITH_ACCUM);

  if (!my_hcreate_r(100, &cfgHash)) {
    Error("Insufficient memory");
  }
  i = ParseOptions(argc, argv, optionStr, SNAME, &cfgHash);
//  htk_compat   = GetParamBool(&cfgHash,SNAME":HTKCOMPAT",       FALSE);
  for (; i < argc; i++) {
    last_file = AddFileElem(last_file, argv[i]);
    nfeature_files++;
  }

  outlabel_map = GetParamStr(&cfgHash, SNAME":OUTPUTLABELMAP",  NULL);
  
  targetKind   = GetDerivParams(&cfgHash, &derivOrder, &derivWinLengths,
                                &startFrmExt, &endFrmExt,
                                &cmn_path, &cmn_file, &cmn_mask,
                                &cvn_path, &cvn_file, &cvn_mask, &cvg_file,
                                SNAME":",  !outlabel_map ? 1 : 0);
  if(!outlabel_map)
  {
    targetKind_out = GetDerivParams(&cfgHash, &derivOrder_out, &derivWinLengths_out,
                                    &startFrmExt_out, &endFrmExt_out,
                                    &cmn_path_out, &cmn_file_out, &cmn_mask_out,
                                    &cvn_path_out, &cvn_file_out, &cvn_mask_out,
                                    &cvg_file_out, SNAME":",  2);
  }
  
  in_lbl_fmt.left_extent  = -100 * (long long) (0.5 + 1e5 *
                 GetParamFlt(&cfgHash, SNAME":STARTTIMESHIFT",  0.0));
  in_lbl_fmt.right_extent =  100 * (long long) (0.5 + 1e5 *
                 GetParamFlt(&cfgHash, SNAME":ENDTIMESHIFT",    0.0));
  NNet_instance_name =
                 GetParamStr(&cfgHash, SNAME":SOURCEINPUT",     NULL); // pripadne rename
  swap_features=!GetParamBool(&cfgHash,SNAME":NATURALREADORDER", isBigEndian());
  swap_features_out=swap_features;
//  swap_fea_out =!GetParamBool(&cfgHash,SNAME":NATURALWRITEORDER",isBigEndian());
  gpFilterWldcrd= GetParamStr(&cfgHash, SNAME":HFILTERWILDCARD", "$");
  script_filter= GetParamStr(&cfgHash, SNAME":HSCRIPTFILTER",   NULL);
  parm_filter  = GetParamStr(&cfgHash, SNAME":HPARMFILTER",     NULL);
//  gpHListFilter = GetParamStr(&cfgHash, SNAME":HMMLISTFILTER",   NULL);
  MMF_filter   = GetParamStr(&cfgHash, SNAME":HMMDEFFILTER",    NULL);
//  parm_ofilter = GetParamStr(&cfgHash, SNAME":HPARMOFILTER",    NULL);
  transc_filter= GetParamStr(&cfgHash, SNAME":HLABELFILTER",    NULL);
//  net_filter   = GetParamStr(&cfgHash, SNAME":HNETFILTER",      NULL);
//  dict_filter  = GetParamStr(&cfgHash, SNAME":HDICTFILTER",     NULL);
  MMF_ofilter  = GetParamStr(&cfgHash, SNAME":HMMDEFOFILTER",   NULL);
//  out_dir      = GetParamStr(&cfgHash, SNAME":TARGETPARAMDIR",  NULL);
//  out_ext      = GetParamStr(&cfgHash, SNAME":TARGETPARAMEXT",  NULL);
  src_mlf      = GetParamStr(&cfgHash, SNAME":SOURCEMLF",       NULL);
  src_lbl_dir  = GetParamStr(&cfgHash, SNAME":SOURCETRANSCDIR", NULL);
  src_lbl_ext  = GetParamStr(&cfgHash, SNAME":SOURCETRANSCEXT", "lab");
  mTraceFlag   = GetParamInt(&cfgHash, SNAME":TRACE",           0);
  hmms_binary  = GetParamBool(&cfgHash,SNAME":SAVEBINARY",      FALSE);
  script =(char*)GetParamStr(&cfgHash, SNAME":SCRIPT",          NULL);
//  src_hmm_list = GetParamStr(&cfgHash, SNAME":SOURCEHMMLIST",   NULL);
//  src_hmm_dir  = GetParamStr(&cfgHash, SNAME":SOURCEMODELDIR",  NULL);
//  src_hmm_ext  = GetParamStr(&cfgHash, SNAME":SOURCEMODELEXT",  NULL);
  src_mmf=(char*)GetParamStr(&cfgHash, SNAME":SOURCEMMF",       NULL);
  trg_hmm_dir  = GetParamStr(&cfgHash, SNAME":TARGETMODELDIR",  NULL);
  trg_hmm_ext  = GetParamStr(&cfgHash, SNAME":TARGETMODELEXT",  NULL);
  trg_mmf      = GetParamStr(&cfgHash, SNAME":TARGETMMF",       NULL);
  
//  in_transc_fmt= (TranscriptionFormat) GetParamEnum(&cfgHash,SNAME":SOURCETRANSCFMT",
//                              !network_file && htk_compat ? TF_HTK : TF_STK,
//                              "HTK", TF_HTK, "STK", TF_STK, NULL);


  if (GetParamBool(&cfgHash, SNAME":PRINTCONFIG", FALSE)) {
    PrintConfig(&cfgHash);
  }
  if (GetParamBool(&cfgHash, SNAME":PRINTVERSION", FALSE)) {
    puts("Version: "VERSION"\n");
  }

  if (!GetParamBool(&cfgHash,SNAME":ACCEPTUNUSEDPARAM", FALSE)) {
    CheckCommandLineParamUse(&cfgHash);
  }

  for (script=strtok(script, ","); script != NULL; script=strtok(NULL, ",")) {
    if ((sfp = my_fopen(script, "rt", script_filter)) == NULL) {
      Error("Cannot open script file %s", script);
    }
    while (fscanf(sfp, "%s", line) == 1) {
      last_file = AddFileElem(last_file, line);
      nfeature_files++;
    }
    my_fclose(sfp);
  }
  for (src_mmf=strtok(src_mmf, ","); src_mmf != NULL; src_mmf=strtok(NULL, ",")) {
    hset.ParseMmf(src_mmf, NULL);
  }
//  if (src_hmm_list) ReadHMMList(&hset,     src_hmm_list, src_hmm_dir, src_hmm_ext);

  if (NNet_instance_name != NULL) {
    Macro *macro = FindMacro(&hset.mXformInstanceHash, NNet_instance_name);
    if (macro == NULL) Error("Undefined input '%s'", NNet_instance_name);
    NNet_instance = (XformInstance *) macro->mpData;
  } else if (hset.mpInputXform) {
    NNet_instance = hset.mpInputXform;
  }
  
  //Store pointer to trasformation instance providing the input for NN
  //and remove this transformation from the input of NN_instance, so that
  //NN_instance will not ask this transformation for its input (see XformPass
  //below), which will allow for reading NN input from cache.
  NNet_input = NNet_instance->mpInput;
  NNet_instance->mpInput = NULL;

  if (outlabel_map) {
    labelHash = readLabelList(outlabel_map);
    
    if (NNet_instance->mOutSize != labelHash.mNEntries) {
        Error("Number of entries [%d] in file '%s' does not match with NNet output size [%d]",
              labelHash.mNEntries, outlabel_map, NNet_instance->mOutSize);
      }
    
    //Allocate buffer, to which example of output vector will be created
    //according to labels
    obs_out = (FLOAT *)malloc(NNet_instance->mOutSize * sizeof(FLOAT));
    if (!obs_out) Error("Insufficient memory");
  }
  
  ilfp = OpenInputMLF(src_mlf);
  
    // MAIN FILE LOOP
  for (file_name = feature_files;
      file_name != NULL;
      file_name = outlabel_map ? file_name->mpNext : file_name->mpNext->mpNext) {

    if (mTraceFlag & 1) {
      if (!outlabel_map) {
        TraceLog("Processing file pair %d/%d '%s' <-> %s",  ++fcnt,
        nfeature_files/2, file_name->mpPhysical,file_name->mpNext->logical);
      } else {
        TraceLog("Processing file %d/%d '%s'", ++fcnt, nfeature_files,file_name->logical);
      }
    }
    int nFrames;
    char *phys_fn = (!outlabel_map ? file_name->mpNext : file_name)->mpPhysical;
    char *lgcl_fn = (!outlabel_map ? file_name->mpNext : file_name)->logical;

    // read sentence weight definition if any ( physical_file.fea[s,e]{weight} )
    FLOAT sentWeight = 1.0; // Use this to wight feature frames;
    
    if ((chrptr = strrchr(phys_fn, '{')) != NULL &&
      ((i=0), sscanf(chrptr, "{%f}%n", &sentWeight, &i), chrptr[i] == '\0')) {
      *chrptr = '\0';
    }
    if (cmn_mask) process_mask(lgcl_fn, cmn_mask, cmn_file);
    if (cvn_mask) process_mask(lgcl_fn, cvn_mask, cvn_file);
    obsMx = ReadHTKFeatures(phys_fn, swap_features,
                            startFrmExt, endFrmExt, targetKind,
                            derivOrder, derivWinLengths, &header,
                            cmn_path, cvn_path, cvg_file, &rhfbuff);

    if (hset.mInputVectorSize != header.mSampleSize / sizeof(float)) {
      Error("Vector size [%d] in '%s' is incompatible with source HMM set [%d]",
            header.mSampleSize/sizeof(float), phys_fn, hset.mInputVectorSize);
    }
    nFrames = header.mNSamples - NNet_input->mTotalDelay;
    if (!outlabel_map) { //If output examples are given by features
                         // read the second set of features ...
                         
      if (cmn_mask_out) process_mask(file_name->logical, cmn_mask_out, cmn_file_out);
      if (cvn_mask_out) process_mask(file_name->logical, cvn_mask_out, cvn_file_out);
      obsMx_out = ReadHTKFeatures(file_name->mpPhysical, swap_features_out,
                                  startFrmExt_out, endFrmExt_out, targetKind_out,
                                  derivOrder_out, derivWinLengths_out, &header_out,
                                  cmn_path_out, cvn_path_out, cvg_file_out, &rhfbuff_out);

      if (nFrames != header_out.mNSamples) {
        Error("Mismatch in number of frames in input/output feature file pair: "
              "'%s' <-> '%s'.", file_name->mpNext->mpPhysical, file_name->mpPhysical);
      }
      
    } else {            // ... otherwise, read corresponding label file
    
      strcpy(label_file, file_name->logical);
      ilfp = OpenInputLabelFile(label_file, src_lbl_dir, src_lbl_ext, ilfp, src_mlf);
      labels = ReadLabels(ilfp, &labelHash, UL_WARN, in_lbl_fmt, header.mSamplePeriod,
                          label_file, src_mlf, NULL);      
      CloseInputLabelFile(ilfp, src_mlf);
    }

    // Initialize all transformation instances (including NN_instance)
    hset.ResetXformInstances();
    Label  *lbl_ptr = labels;
    time = 1;
    if (NNet_input) time -= NNet_input->mTotalDelay;
    // Loop over all feature frames.
    for (i = 0; i < header.mNSamples; i++, time++) {
      int j;
      FLOAT *obs = obsMx + i * hset.mInputVectorSize;

      //Get next NN input vector by propagating vector from feature file
      //through NNet_input transformation.
      obs = NNet_input->XformPass(obs, time, FORWARD);
      //Input of NN is not yet read because of NNet_input delay
      if (time <= 0) continue;

      if (outlabel_map) {
        //Create NN output example vector from lables
        for (j = 0; j < NNet_instance->mOutSize; j++) obs_out[j] = 0;
        while (lbl_ptr && lbl_ptr->mStop < time) lbl_ptr = lbl_ptr->mpNext;
        if (lbl_ptr && lbl_ptr->mStart <= time) obs_out[(int) lbl_ptr->mpData - 1] = 1;
      } else {
        //Get NN output example vector from obsMx_out matrix
        obs_out = obsMx_out + (time-1) * NNet_instance->mOutSize;
      }

/////////////////////////////////////////////////////////////////////////////////
//
//    FOR EACH FRAME OF EACH FILE PRINT INPUT TO NN AND ITS DESIRED OUTPUT
//
//             !!!!! TO BE REWRITEN BY SOMETHING MEANINGFUL !!!!!
//
/////////////////////////////////////////////////////////////////////////////////
      for(int j=0; j< NNet_input->mOutSize; j++)    printf("%5.2f ", obs[j]);
      printf("-> ");
      for(int j=0; j< NNet_instance->mOutSize; j++) printf("%5.2f ", obs_out[j]);
      printf("\n");
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
    }

    totFrames  += nFrames;
    //TraceLog("[%d frames]", nFrames);
    free(obsMx);

    if (!outlabel_map) {
      free(obsMx_out);
    } else {
      ReleaseLabels(labels);
    }
  }
  /// END - MAIN FILE LOOP

  if (mTraceFlag & 2) {
    TraceLog("Total number of frames: %d", totFrames);
  }
  NNet_instance->mpInput = NNet_input;

  //writeNN(data, info, prog, matrix, inCache, compCache);
  hset.WriteMmf(trg_mmf, trg_hmm_dir, trg_hmm_ext, hmms_binary);
  hset.Release();

  for (i = 0; i < cfgHash.mNEntries; i++) free(cfgHash.mpEntry[i]->data);
  my_hdestroy_r(&cfgHash, 1);

  free(derivWinLengths);
  if (src_mlf) fclose(ilfp);
  if (outlabel_map) 
  {
    my_hdestroy_r(&labelHash, 1);
    free(obs_out);
  } 
  else 
  {
    free(derivWinLengths_out);
  }
  
  while (feature_files) {
    file_name = feature_files;
    feature_files = feature_files->mpNext;
    free(file_name);
  }
  return 0;
}
