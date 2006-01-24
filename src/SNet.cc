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

/*
int main(int argc, char *argv[]) 
{
  HMMSet hset;
  FILE *sfp, *ofp = NULL;
  FLOAT *obsMx, *obs;
  HtkHeader header;
  int i, fcnt = 0;
  XformInstance *input = NULL;
  char line[1024];
  char outFile[1024];
  MyHSearchData phoneHash, cfgHash;
  int totFrames = 0;
  int vec_size, out_size, time;
  FileListElem *feature_files = NULL;
  int nfeature_files = 0;
  FileListElem *file_name = NULL;
  FileListElem **last_file = &feature_files;


  const char *src_hmm_dir;
  const char *src_hmm_ext;
        char *src_mmf;
  const char *out_dir;
  const char *out_ext;
  const char *inputName;
  char *script;
  int  mTraceFlag;
  int  targetKind;
  int  derivOrder;
  int  *derivWinLengths;
  int startFrmExt;
  int endFrmExt;
  bool swap_features;
  bool swap_fea_out;
//  char *lbl_list_file   = NULL;
  char *NNet_instance_name     = NULL;
  XformInstance *NNet_instance = NULL;
  XformInstance *NNet_input    = NULL;
  LabelFormat in_lbl_fmt = {0};

  if (argc == 1) usage(argv[0]);

  InitHMMSet(&hset, 1);

  if (!my_hcreate_r(100,  &phoneHash)
  || !my_hcreate_r(100,  &cfgHash)) {
    Error("Insufficient memory");
  }
  i = ParseOptions(argc, argv, optionStr, SNAME, &cfgHash);
//  htk_compat   = GetParamBool(&cfgHash,SNAME":HTKCOMPAT",       FALSE);
  for (; i < argc; i++) {
    last_file = AddFileElem(last_file, argv[i]);
    nfeature_files++;
  }
  targetKind   = GetDerivParams(&cfgHash, &derivOrder, &derivWinLengths,
                                &startFrmExt, &endFrmExt, SNAME":",0);
  in_lbl_fmt.left_extent  = -100 * (long long) (0.5 + 1e5 *
                 GetParamFlt(&cfgHash, SNAME":STARTTIMESHIFT",  0.0));
  in_lbl_fmt.right_extent =  100 * (long long) (0.5 + 1e5 *
                 GetParamFlt(&cfgHash, SNAME":ENDTIMESHIFT",    0.0));
  NNet_instance_name =
                 GetParamStr(&cfgHash, SNAME":SOURCEINPUT",     NULL); // pripadne rename
  swap_features=!GetParamBool(&cfgHash,SNAME":NATURALREADORDER", isBigEndian());
//  swap_fea_out =!GetParamBool(&cfgHash,SNAME":NATURALWRITEORDER",isBigEndian());
  gpFilterWldcrd= GetParamStr(&cfgHash, SNAME":HFILTERWILDCARD", "$");
  script_filter= GetParamStr(&cfgHash, SNAME":HSCRIPTFILTER",   NULL);
  parm_filter  = GetParamStr(&cfgHash, SNAME":HPARMFILTER",     NULL);
//  gpHListFilter = GetParamStr(&cfgHash, SNAME":HMMLISTFILTER",   NULL);
  MMF_filter   = GetParamStr(&cfgHash, SNAME":HMMDEFFILTER",    NULL);
//  parm_ofilter = GetParamStr(&cfgHash, SNAME":HPARMOFILTER",    NULL);
  label_filter = GetParamStr(&cfgHash, SNAME":HLABELFILTER",    NULL);
//  net_filter   = GetParamStr(&cfgHash, SNAME":HNETFILTER",      NULL);
//  dict_filter  = GetParamStr(&cfgHash, SNAME":HDICTFILTER",     NULL);
  MMF_ofilter  = GetParamStr(&cfgHash, SNAME":HMMDEFOFILTER",   NULL);
//  out_dir      = GetParamStr(&cfgHash, SNAME":TARGETPARAMDIR",  NULL);
//  out_ext      = GetParamStr(&cfgHash, SNAME":TARGETPARAMEXT",  NULL);
  src_mlf      = GetParamStr(&cfgHash, SNAME":SOURCEMLF",       NULL);
  src_lbl_dir  = GetParamStr(&cfgHash, SNAME":SOURCETRANSCDIR", NULL);
  src_lbl_ext  = GetParamStr(&cfgHash, SNAME":SOURCETRANSCEXT", NULL);
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
    ReadHMMSet(src_mmf, &hset, NULL);
  }
//  if (src_hmm_list) ReadHMMList(&hset,     src_hmm_list, src_hmm_dir, src_hmm_ext);

  if (NNet_instance_name != NULL) {
    Macro *macro = FindMacro(&hset.Xform_instance_hash, NNet_instance_name);
    if (macro == NULL) Error("Undefined input '%s'", NNet_instance_name);
    NNet_instance = (XformInstance *) macro->mpData;
  } else if (hset.inputXform) {
    NNet_instance = hset.inputXform;
  }
  
  //Store pointer to trasformation instance providing the input for NN
  //and remove this transformation from the input of NN_instance, so that
  //NN_instance will not ask this transformation for its input (see XformPass
  //below), which will allow for reading NN input from cache.
  NNet_input = NNet_instance->mpInput;
  NNet_instance->mpInput = NULL;

  if (out_by_labels) {
    //Allocate buffer, to which example of output vector will be created
    //according to labels
    obs_out = (FLOAT *)malloc(NNet_instance->mOutSize * sizeof(FLOAT));
    if (!obs_out) Error("Insufficient memory");
  }
  
  ilfp = OpenInputMLF(src_mlf);
  
    // MAIN FILE LOOP
  for (file_name = feature_files;
      file_name != NULL;
      file_name = out_by_labels ? file_name->mpNext : file_name->mpNext->mpNext) {

    if (mTraceFlag & 1) {
      if (!out_by_labels) {
        TraceLog("Processing file pair %d/%d '%s' <-> %s",  ++fcnt,
        nfeature_files/2, file_name->mpPhysical,file_name->mpNext->logical);
      } else {
        TraceLog("Processing file %d/%d '%s'", ++fcnt, nfeature_files,file_name->logical);
      }
    }
    int nFrames;
    char *phys_fn = (!out_by_labels ? file_name->mpNext : file_name)->mpPhysical;
    char *lgcl_fn = (!out_by_labels ? file_name->mpNext : file_name)->logical;

    // read sentence weight definition if any ( physical_file.fea[s,e]{weight} )
    if ((chrptr = strrchr(phys_fn, '{')) != NULL &&
      ((i=0), sscanf(chrptr, "{%f}%n", &sentWeight, &i), chrptr[i] == '\0')) {
      *chrptr = '\0';
    } else {
      sentWeight = 1.0;
    }
    if (cmn_mask) process_mask(lgcl_fn, cmn_mask, cmn_file);
    if (cvn_mask) process_mask(lgcl_fn, cvn_mask, cvn_file);
    obsMx = ReadHTKFeatures(phys_fn, swap_features,
                            startFrmExt, endFrmExt, targetKind,
                            derivOrder, derivWinLengths, &header,
                            cmn_path, cvn_path, cvg_file, &rhfbuff);

    if (hset.in_vec_size != header.mSampleSize / sizeof(float)) {
      Error("Vector size [%d] in '%s' is incompatible with source HMM set [%d]",
            header.mSampleSize/sizeof(float), phys_fn, hset.in_vec_size);
    }
    nFrames = header.mNSamples - NNet_input->totalDelay;
    if (!out_by_labels) { //If output examples are given by features
                         // read the second set of features ...
                         
      if (cmn_mask_out) process_mask(file_name->logical, cmn_mask_out, cmn_file_out);
      if (cvn_mask_out) process_mask(file_name->logical, cvn_mask_out, cvn_file_out);
      obsMx_out = ReadHTKFeatures(file_name->mpPhysical, swap_features_out,
                                  startFrmExt_out, endFrmExt_out, targetKind_out,
                                  derivOrder_out, derivWinLengths_out, &header_out,
                                  cmn_path_out, cvn_path_out, cvg_file_out, &rhfbuff_out);

      if (NNet_instance->mOutSize != header_out.mSampleSize/sizeof(float)) {
        Error("Vector size [%d] in '%s' is incompatible with NNet output size [%d]",
              header_out.mSampleSize/sizeof(float), file_name->mpPhysical,
              NNet_instance->mOutSize);
      }
      if (nFrames != header_out.mNSamples) {
        Error("Mismatch in number of frames in input/output feature file pair: "
              "'%s' <-> '%s'.", file_name->mpNext->mpPhysical, file_name->mpPhysical);
      }
      
    } else {            // ... otherwise, read corresponding label file
    
      strcpy(label_file, file_name->logical);
      ilfp = OpenInputLabelFile(label_file, src_lbl_dir, src_lbl_ext, ilfp, src_mlf);
      labels = ReadLabels(ilfp, &phoneHash, UL_INSERT, in_lbl_fmt, header.mSamplePeriod,
                          label_file, src_mlf, NULL);      
      CloseInputLabelFile(lfp, src_mlf);
    }

    // Initialize all transformation instances (including NN_instance)
    ResetXformInstances(&hset);
    Label  *lbl_ptr = labels;
    time = 1;
    if (NNet_input) time -= NNet_input->totalDelay;
    // Loop over all feature frames.
    for (i = 0; i < header.mNSamples; i++, time++) {
      int j;
      FLOAT *obs = obsMx + i * hset.in_vec_size;

      //Get next NN input vector by propagating vector from feature file
      //through NNet_input transformation.
      obs = XformPass(NNet_input, obs, time, FORWARD);
      //Input of NN is not yet read because of NNet_input delay
      if (time <= 0) continue;

      if (lbl_list_file) {
        //Create NN output example vector from lables
        for (j = 0; j < NNet_instance->mOutSize; j++) obs_out[j] = 0;
        while (lbl_ptr && lbl_ptr->mStop < time) lbl_ptr = lbl_ptr->mpNext;
        if (lbl_ptr && lbl_ptr->mStart <= time) obs_out[(int) lbl_ptr->mpData - 1] = 1;
      } else {
        //Get NN output example vector from obsMx_out matrix
        obs_out = obsMx_out + (time-1) * NNet_instance->mOutSize;
      }

      // Saving to cache
      memcpy(giveRowPointer(inCache, info->actualCache), obs, inCache->cols*sizeof(FLOAT));

      // Saving output to cache
      memcpy(giveRowPointer(outCache, info->actualCache), obs_out, outCache->cols*sizeof(FLOAT));

      info->actualCache++;

      if (info->actualCache < info->cacheSize)
        continue;

 
    }

    totFrames  += nFrames;
    //TraceLog("[%d frames]", nFrames);
    free(obsMx);

    if (!lbl_list_file) {
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
  WriteHMMSet(out_MMF, out_hmm_dir, out_hmm_ext, hmms_binary, &hset);


  
  
  //ReleaseHMMSet(&hset);
  hset.Release();

//  my_hdestroy_r(&labelHash, 0);
  my_hdestroy_r(&phoneHash, 1);
  for (i = 0; i < cfgHash.mNEntries; i++) free(cfgHash.mpEntry[i]->data);
  my_hdestroy_r(&cfgHash, 1);

  if (network_file) {
    net.Release();
  }
  free(hset.varFloor);
  free(derivWinLengths);
  free(derivWinLengths_out);
  if (src_mlf)  fclose(ilfp);
  while (feature_files) {
    file_name = feature_files;
    feature_files = feature_files->mpNext;
    free(file_name);
  }
  return 0;
}

    
    */
    



#include"SNetLib/nnet.h"
#include"SNetLib/progobj.h"

using namespace SNet;

// ::TODO:: >>> Ask Lukas about Trace!

int main(int argc, char *argv[]){
  ProgObj *p_prog_obj = new ProgObj;
  
  // read command line here
  
  
  
  
  
  p_prog_obj->PrintCommandLine(argc, argv);

  return 0;
}

