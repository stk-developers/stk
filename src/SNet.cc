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
#include "STKLib/hmms.h"
#include "STKLib/viterbi.h"
#ifndef WIN32
#include <unistd.h>
#else
#include "getopt.h"
#endif

void usage(char *progname)
{
  char *tchrptr;
  if((tchrptr = strrchr(progname, '\\')) != NULL) progname = tchrptr+1;
  if((tchrptr = strrchr(progname, '/')) != NULL) progname = tchrptr+1;
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
" -G r   SOURCETRANSCFMT"
" -H l   SOURCEMMF"
" -I r   SOURCEMLF"
" -L r   SOURCETRANSCDIR"
" -M r   TARGETMODELDIR"
" -S l   SCRIPT"
" -T r   TRACE"
" -V n   PRINTVERSION=TRUE"
" -X r   SOURCETRANSCEXT";


int main(int argc, char *argv[]) {
  HMMSet hset;
  FILE *sfp, *ofp = NULL;
  FLOAT *obsMx, *obs;
  HTK_Header header;
  int i, fcnt = 0;
  XformInstance *input = NULL;
  char line[1024];
  char outFile[1024];
  struct my_hsearch_data cfgHash;
  int totFrames = 0;
  int vec_size, out_size, time;
  FileListElem *feature_files = NULL;
  int nfeature_files = 0;
  FileListElem *file_name = NULL;
  FileListElem **last_file = &feature_files;


  const char *src_hmm_list;
  const char *src_hmm_dir;
  const char *src_hmm_ext;
        char *src_mmf;
  const char *out_dir;
  const char *out_ext;
  const char *inputName;
  char *script;
  int  trace_flag;
  int  targetKind;
  int  derivOrder;
  int  *derivWinLengths;
  int startFrmExt;
  int endFrmExt;
  BOOL swap_features;
  BOOL swap_fea_out;

  if(argc == 1) usage(argv[0]);

  InitHMMSet(&hset, 1);

  if(!my_hcreate_r(100,  &cfgHash)) {
    Error("Insufficient memory");
  }
  i = ParseOptions(argc, argv, optionStr, SNAME, &cfgHash);
//  htk_compat   = GetParamBool(&cfgHash,SNAME":HTKCOMPAT",       FALSE);
  for(; i < argc; i++) {
    last_file = AddFileElem(last_file, argv[i]);
    nfeature_files++;
  }
  targetKind   = GetDerivParams(&cfgHash, &derivOrder, &derivWinLengths,
                                &startFrmExt, &endFrmExt, SNAME":",0);
  inputName    = GetParamStr(&cfgHash, SNAME":SOURCEINPUT",     NULL);
  swap_features=!GetParamBool(&cfgHash,SNAME":NATURALREADORDER", isBigEndian());
  swap_fea_out =!GetParamBool(&cfgHash,SNAME":NATURALWRITEORDER",isBigEndian());
  filter_wldcrd= GetParamStr(&cfgHash, SNAME":HFILTERWILDCARD", "$");
  script_filter= GetParamStr(&cfgHash, SNAME":HSCRIPTFILTER",   NULL);
  parm_filter  = GetParamStr(&cfgHash, SNAME":HPARMFILTER",     NULL);
  hlist_filter = GetParamStr(&cfgHash, SNAME":HMMLISTFILTER",   NULL);
  MMF_filter   = GetParamStr(&cfgHash, SNAME":HMMDEFFILTER",    NULL);
  parm_ofilter = GetParamStr(&cfgHash, SNAME":HPARMOFILTER",    NULL);
  out_dir      = GetParamStr(&cfgHash, SNAME":TARGETPARAMDIR",  NULL);
  out_ext      = GetParamStr(&cfgHash, SNAME":TARGETPARAMEXT",  NULL);
  trace_flag   = GetParamInt(&cfgHash, SNAME":TRACE",           0);
  script =(char*)GetParamStr(&cfgHash, SNAME":SCRIPT",          NULL);
  src_hmm_list = GetParamStr(&cfgHash, SNAME":SOURCEHMMLIST",   NULL);
  src_hmm_dir  = GetParamStr(&cfgHash, SNAME":SOURCEMODELDIR",  NULL);
  src_hmm_ext  = GetParamStr(&cfgHash, SNAME":SOURCEMODELEXT",  NULL);
  src_mmf=(char*)GetParamStr(&cfgHash, SNAME":SOURCEMMF",       NULL);

  if(GetParamBool(&cfgHash, SNAME":PRINTCONFIG", FALSE)) {
    PrintConfig(&cfgHash);
  }
  if(GetParamBool(&cfgHash, SNAME":PRINTVERSION", FALSE)) {
    puts("Version: "VERSION"\n");
  }

  if(!GetParamBool(&cfgHash,SNAME":ACCEPTUNUSEDPARAM", FALSE)) {
    CheckCommandLineParamUse(&cfgHash);
  }

  for(script=strtok(script, ","); script != NULL; script=strtok(NULL, ",")) {
    if((sfp = my_fopen(script, "rt", script_filter)) == NULL) {
      Error("Cannot open script file %s", script);
    }
    while(fscanf(sfp, "%s", line) == 1) {
      last_file = AddFileElem(last_file, line);
      nfeature_files++;
    }
    my_fclose(sfp);
  }
  for(src_mmf=strtok(src_mmf, ","); src_mmf != NULL; src_mmf=strtok(NULL, ",")) {
    ReadHMMSet(src_mmf, &hset, NULL);
  }
  if(src_hmm_list) ReadHMMList(&hset,     src_hmm_list, src_hmm_dir, src_hmm_ext);

  if(inputName != NULL) {
    Macro *macro = FindMacro(&hset.Xform_instance_hash, inputName);
    if(macro == NULL) Error("Undefined source input '%s'", inputName);
    input = (XformInstance *) macro->data;
  } else if(hset.inputXform) {
    input = hset.inputXform;
  }

  for(file_name=feature_files; file_name != NULL; file_name=file_name->next) {
    if(trace_flag & 1) TraceLog("Processing file %d/%d '%s'",
                                ++fcnt, nfeature_files, file_name->logical);
    obsMx = ReadHTKFeatures(file_name->physical, swap_features,
                            startFrmExt, endFrmExt, targetKind,
                            derivOrder, derivWinLengths, &header);

    vec_size = header.sampSize / sizeof(float);
    out_size = input ? input->out_size : vec_size;

    if(hset.in_vec_size != -1 && hset.in_vec_size != vec_size) {
      Error("Vector size [%d] in '%s' is incompatible with HMM set [%d]",
            header.sampSize/sizeof(float), file_name->physical, hset.in_vec_size);
    }
    MakeFileName(outFile, file_name->logical, out_dir, out_ext);

    if((ofp = fopen(outFile, "wb")) == NULL) {
      Error("Cannot open output feature file: '%s'", outFile);
    }
    header.sampKind = input ? PARAMKIND_USER : targetKind;
    header.sampSize = out_size * sizeof(float);

    if(WriteHTKHeader(ofp, header, 1)) {
      Error("Cannot write to output feature file: '%s'", outFile);
    }
    time = -hset.totalDelay;
    ResetXformInstances(&hset);

    for(i = 0; i < header.nSamples; i++) {
      UpdateStacks(&hset, obsMx + i * vec_size, ++time, FORWARD);
      if(time <= 0) continue;

      obs = XformPass(input, obsMx + i * vec_size, time, FORWARD);

      if(WriteHTKFeature (ofp, obs, out_size, swap_fea_out)) {
        Error("Cannot write to output feature file: '%s'", outFile);
      }
    }
    totFrames += header.nSamples - hset.totalDelay;

    if(trace_flag & 1) TraceLog("[%d frames]", header.nSamples-hset.totalDelay);
    fclose(ofp);
    free(obsMx);
  }
  if(trace_flag & 2) TraceLog("Total number of frames: %d", totFrames);

  ReleaseHMMSet(&hset);
  free(derivWinLengths);

  while(feature_files) {
    file_name = feature_files;
    feature_files = feature_files->next;
    free(file_name);
  }
  return 0;
}



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
