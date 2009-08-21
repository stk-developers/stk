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

#define SVN_DATE       "$Date: 2008-03-10 16:44:02 +0100 (Mon, 10 Mar 2008) $"
#define SVN_AUTHOR     "$Author: glembek $"
#define SVN_REVISION   "$Revision: 236 $"
#define SVN_ID         "$Id: SNet.cc 236 2008-03-10 15:44:02Z glembek $"


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

#include "SNetLib/progobj.h"

using namespace STK;
using namespace SNet;

void usage(char *progname)
{
  char *tchrptr;
  if ((tchrptr = strrchr(progname, '\\')) != NULL) progname = tchrptr+1;
  if ((tchrptr = strrchr(progname, '/')) != NULL) progname = tchrptr+1;
  fprintf(stderr,
"\n%s version " MODULE_VERSION "\n"
"\nUSAGE: %s [options] DataFiles...\n\n"
" Option                                                     Default\n\n"
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
" -X ext     Set input label (or network) file ext           lab (net)\n"
" -o s       Extension for new hmm files                     As src\n"
"\n"
"\n"
"SNet based parameters:\n"
"--CROSSVALIDATION                                           false\n"
"--CACHESIZE                                                 12000\n"
"--BUNCHSIZE                                                 1000\n"
"--LEARNINGRATE                                              0.008\n"
"--LEARNINGRATELISTMUL (like \"1.0,0.5\")                      NULL\n"
"--CLIENTS                                                   0\n"
"--JOINIP                                                    NULL\n"
"--RANDOMIZE                                                 true\n"
"--SYNCHRONIZE                                               true\n"
"--PORT                                                      2020\n"
"--SEED                                                      0 (means random)\n"
"--NCUMRBU                                                   CLIENTS\n"
"\n"
" %s is Copyright (C) 2004-2005 Stanislav Kontar, Lukas Burget, Ondrej Glembek et al. and\n"
" licensed under the GNU General Public License, version 2.\n"
" Bug reports, feedback, etc, to: burget@fit.vutbr.cz\n"
"\n", progname, progname, progname);
  exit(-1);
}

#define SNAME "SNET"
const char *optionStr =
" -o r   TARGETMODELEXT"
" -B n   SAVEBINARY=TRUE"
" -D n   PRINTCONFIG=TRUE" 
" -H l   SOURCEMMF"
" -I r   SOURCEMLF"
" -L r   SOURCETRANSCDIR"
" -M r   TARGETMODELDIR"
" -S l   SCRIPT"
" -T r   TRACE"
" -X r   SOURCETRANSCEXT";


int main(int argc, char *argv[]) 
{
  ModelSet        hset;
  FILE *          sfp;
  FILE *          ilfp;

  Matrix<FLOAT>   feature_matrix;
  Matrix<FLOAT>   feature_matrix_out;

  FLOAT *         obs_out   = NULL;
  HtkHeader       header, 
                  header_out;
  Label *         labels    = NULL;
  int i, fcnt = 0;
  char line[1024];
  char label_file[1024];
  char *chrptr;
  MyHSearchData labelHash, cfgHash;
  int totFrames = 0;
  int time;
  FileListElem *feature_files = NULL;
  int nfeature_files = 0;
  FileListElem *file_name = NULL;
  FileListElem **last_file = &feature_files;

        char *src_mmf;
  const char *trg_hmm_dir;
  const char *trg_hmm_ext;
  const char *trg_mmf;
  const char *src_lbl_dir;
  const char *src_lbl_ext;
  const char *src_mlf;
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
  int                 trace_flag;
  int                 target_kind = 0, 
                      target_kind_out = 0;
  int  derivOrder, derivOrder_out;
  int  *derivWinLengths, *derivWinLengths_out;
  int startFrmExt;
  int endFrmExt;
  int startFrmExt_out;
  int endFrmExt_out;
  bool hmms_binary;
  bool swap_features;
  bool swap_features_out;
  bool expand_predef_xforms;  
  const char *NNet_instance_name = NULL;
  XformInstance *NNet_instance   = NULL;
  XformInstance *NNet_input      = NULL;
  STKNetworkOutputFormat in_lbl_fmt = {0};
  RHFBuffer rhfbuff      = {0};
  RHFBuffer rhfbuff_out  = {0};
    
  bool cross_validation = false;
  int cache_size;
  int bunch_size;
  float learning_rate;
  char *learning_rate_list;
  int clients;
  char *ip;
  bool randomize;
  bool sync;
  int port;
  int seed;
  int ncumrbu;

  if (argc == 1) usage(argv[0]);

  hset.Init(MODEL_SET_WITH_ACCUM);
  
  // :TODO: Ondra
  // temporary stuff, later to be removed
  // now look in Models.h
  hset.mUseNewMatrix = true;

  if (!my_hcreate_r(100, &cfgHash)) {
    Error("Insufficient memory");
  }
  i = ParseOptions(argc, argv, optionStr, SNAME, &cfgHash);
//  htk_compat   = GetParamBool(&cfgHash,SNAME":HTKCOMPAT",       FALSE);

  for (; i < argc; i++) 
  {
    last_file = AddFileElem(last_file, argv[i]);
    nfeature_files++;
  }

  outlabel_map = GetParamStr(&cfgHash, SNAME":OUTPUTLABELMAP",  NULL);
  
  target_kind   = GetDerivParams(&cfgHash, &derivOrder, &derivWinLengths,
                                &startFrmExt, &endFrmExt,
                                &cmn_path, &cmn_file, &cmn_mask,
                                &cvn_path, &cvn_file, &cvn_mask, &cvg_file,
                                SNAME":",  !outlabel_map ? 1 : 0);
  if(!outlabel_map)
  {
    target_kind_out = GetDerivParams(&cfgHash, &derivOrder_out, &derivWinLengths_out,
                                    &startFrmExt_out, &endFrmExt_out,
                                    &cmn_path_out, &cmn_file_out, &cmn_mask_out,
                                    &cvn_path_out, &cvn_file_out, &cvn_mask_out,
                                    &cvg_file_out, SNAME":",  2);
  }
  
  in_lbl_fmt.mStartTimeShift = 
    GetParamFlt(&cfgHash, SNAME":STARTTIMESHIFT",  0.0);
  in_lbl_fmt.mEndTimeShift = 
    GetParamFlt(&cfgHash, SNAME":ENDTIMESHIFT",    0.0);
  NNet_instance_name = 
    GetParamStr(&cfgHash, SNAME":SOURCEINPUT",     NULL);
  swap_features = 
    !GetParamBool(&cfgHash,SNAME":NATURALREADORDER", isBigEndian());
  swap_features_out=swap_features;
//  swap_fea_out =!GetParamBool(&cfgHash,SNAME":NATURALWRITEORDER",isBigEndian());
  gpFilterWldcrd= GetParamStr(&cfgHash, SNAME":HFILTERWILDCARD", "$");
  gpScriptFilter= GetParamStr(&cfgHash, SNAME":HSCRIPTFILTER",   NULL);
  gpParmFilter  = GetParamStr(&cfgHash, SNAME":HPARMFILTER",     NULL);
//  gpHListFilter = GetParamStr(&cfgHash, SNAME":HMMLISTFILTER",   NULL);
  gpMmfFilter   = GetParamStr(&cfgHash, SNAME":HMMDEFFILTER",    NULL);
//  gpParmOFilter = GetParamStr(&cfgHash, SNAME":HPARMOFILTER",    NULL);
  transc_filter= GetParamStr(&cfgHash, SNAME":HLABELFILTER",    NULL);
//  net_filter   = GetParamStr(&cfgHash, SNAME":HNETFILTER",      NULL);
//  dict_filter  = GetParamStr(&cfgHash, SNAME":HDICTFILTER",     NULL);
  gpMmfOFilter  = GetParamStr(&cfgHash, SNAME":HMMDEFOFILTER",   NULL);
//  out_dir      = GetParamStr(&cfgHash, SNAME":TARGETPARAMDIR",  NULL);
//  out_ext      = GetParamStr(&cfgHash, SNAME":TARGETPARAMEXT",  NULL);
  src_mlf      = GetParamStr(&cfgHash, SNAME":SOURCEMLF",       NULL);
  src_lbl_dir  = GetParamStr(&cfgHash, SNAME":SOURCETRANSCDIR", NULL);
  src_lbl_ext  = GetParamStr(&cfgHash, SNAME":SOURCETRANSCEXT", "lab");
  trace_flag   = GetParamInt(&cfgHash, SNAME":TRACE",           0);
  hmms_binary  = GetParamBool(&cfgHash,SNAME":SAVEBINARY",      false);
  script =(char*)GetParamStr(&cfgHash, SNAME":SCRIPT",          NULL);
//  src_hmm_list = GetParamStr(&cfgHash, SNAME":SOURCEHMMLIST",   NULL);
//  src_hmm_dir  = GetParamStr(&cfgHash, SNAME":SOURCEMODELDIR",  NULL);
//  src_hmm_ext  = GetParamStr(&cfgHash, SNAME":SOURCEMODELEXT",  NULL);
  src_mmf=(char*)GetParamStr(&cfgHash, SNAME":SOURCEMMF",       NULL);
  trg_hmm_dir  = GetParamStr(&cfgHash, SNAME":TARGETMODELDIR",  NULL);
  trg_hmm_ext  = GetParamStr(&cfgHash, SNAME":TARGETMODELEXT",  NULL);
  trg_mmf      = GetParamStr(&cfgHash, SNAME":TARGETMMF",       NULL);
  
  cross_validation  = GetParamBool(&cfgHash,SNAME":CROSSVALIDATION", false);
  cache_size   = GetParamInt(&cfgHash, SNAME":CACHESIZE",            12000);
  bunch_size   = GetParamInt(&cfgHash, SNAME":BUNCHSIZE",            1000);
  learning_rate   = GetParamFlt(&cfgHash, SNAME":LEARNINGRATE",      0.008);
  learning_rate_list =(char*)GetParamStr(&cfgHash, SNAME":LEARNINGRATELISTMUL",          NULL);
  clients   = GetParamInt(&cfgHash, SNAME":CLIENTS",            0);
  ip = (char*)GetParamStr(&cfgHash, SNAME":JOINIP",       NULL);
  randomize =  GetParamBool(&cfgHash,SNAME":RANDOMIZE", true);
  sync = GetParamBool(&cfgHash,SNAME":SYNCHRONIZE", true);
  port = GetParamInt(&cfgHash, SNAME":PORT",            2020);
  seed = GetParamInt(&cfgHash, SNAME":SEED",            0);
  ncumrbu = GetParamInt(&cfgHash, SNAME":NCUMRBU",            clients);
  expand_predef_xforms = GetParamBool(&cfgHash,SNAME":EXPANDPREDEFXFORMS", true);
  // NCUMRBU == number of client update matrixes received before update in async version
  
  if(clients != 0 && script != NULL) Error("Server should not have input data.");
  
//  in_transc_fmt= (TranscriptionFormat) GetParamEnum(&cfgHash,SNAME":SOURCETRANSCFMT",
//                              !network_file && htk_compat ? TF_HTK : TF_STK,
//                              "HTK", TF_HTK, "STK", TF_STK, NULL);

  bool print_all_options = GetParamBool(&cfgHash,SNAME":PRINTALLOPTIONS", false);


  if (GetParamBool(&cfgHash, SNAME":PRINTCONFIG", false)) {
    PrintConfig(&cfgHash);
  }

  if (!GetParamBool(&cfgHash,SNAME":ACCEPTUNUSEDPARAM", false)) {
    CheckCommandLineParamUse(&cfgHash);
  }

  if (print_all_options) 
  {
    print_registered_parameters();
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
  
  if (NULL != src_mmf)
  {
    for (src_mmf=strtok(src_mmf, ","); src_mmf != NULL; src_mmf=strtok(NULL, ",")) {
      hset.ParseMmf(src_mmf, NULL);
    }
    if(expand_predef_xforms)
    {
      hset.ExpandPredefXforms();
    }    
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
    
    if (NNet_instance->OutSize() != labelHash.mNEntries) {
        Error("Number of entries [%d] in file '%s' does not match with NNet output size [%d]",
              labelHash.mNEntries, outlabel_map, NNet_instance->OutSize());
      }

    //Allocate buffer, to which example of output vector will be created
    //according to labels
    obs_out = (FLOAT *)malloc(NNet_instance->OutSize() * sizeof(FLOAT));
    if (!obs_out) Error("Insufficient memory");
  }
  
  ilfp = OpenInputMLF(src_mlf);
  
///***************************************************************************
///Here starts the processing of files
///************************************************************************************************

  int Hits=0;
  // main file loop
  for (file_name = feature_files;
       file_name != NULL;
       file_name = outlabel_map ? file_name->mpNext : file_name->mpNext->mpNext) 
  {

    if (trace_flag & 1) 
    {
      if (!outlabel_map) 
      {
        TraceLog("Processing file pair %d/%d '%s' <-> %s",  ++fcnt,
        nfeature_files/2, file_name->mpPhysical,file_name->mpNext->logical);
      } 
      else 
      {
        TraceLog("Processing file %d/%d '%s'", ++fcnt, nfeature_files,file_name->logical);
      }
    }
    
    int     n_frames;
    char*   phys_fn = (!outlabel_map ? file_name->mpNext : file_name)->mpPhysical;
    char*   lgcl_fn = (!outlabel_map ? file_name->mpNext : file_name)->logical;
    FLOAT   sentWeight = 1.0; // Use this to wight feature frames;
    
    // read sentence weight definition if any ( physical_file.fea[s,e]{weight} )
    if ((chrptr = strrchr(phys_fn, '{')) != NULL 
    && ((i=0), sscanf(chrptr, "{%f}%n", &sentWeight, &i), chrptr[i] == '\0')) 
    {
      *chrptr = '\0';
    }

    if (cmn_mask) process_mask(lgcl_fn, cmn_mask, cmn_file);
    if (cvn_mask) process_mask(lgcl_fn, cvn_mask, cvn_file);

    ReadHTKFeatures(phys_fn, swap_features,
                    startFrmExt, endFrmExt, target_kind,
                    derivOrder, derivWinLengths, &header,
                    cmn_path, cvn_path, cvg_file, &rhfbuff,feature_matrix);

    if ((size_t) hset.mInputVectorSize != header.mSampleSize / sizeof(float)) 
    {
      Error("Vector size [%d] in '%s' is incompatible with source HMM set [%d]",
            header.mSampleSize/sizeof(float), phys_fn, hset.mInputVectorSize);
    }

    n_frames = header.mNSamples - (NNet_input ? NNet_input->mTotalDelay : 0);

    if (!outlabel_map) 
    { // If output examples are given by features
      // read the second set of features ...
      if (cmn_mask_out) process_mask(file_name->logical, cmn_mask_out, cmn_file_out);
      if (cvn_mask_out) process_mask(file_name->logical, cvn_mask_out, cvn_file_out);
      
      ReadHTKFeatures(file_name->mpPhysical, swap_features_out,
                      startFrmExt_out, endFrmExt_out, target_kind_out,
                      derivOrder_out, derivWinLengths_out, &header_out,
                      cmn_path_out, cvn_path_out, cvg_file_out, &rhfbuff_out,
                      feature_matrix_out);

      if (n_frames != header_out.mNSamples) 
      {
        Error("Mismatch in number of frames in input/output feature file pair: "
              "'%s' <-> '%s'.", file_name->mpNext->mpPhysical, file_name->mpPhysical);
      }
      
    } 
    else 
    { // ... otherwise, read corresponding label file           
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
    
    if (NNet_input) 
      time -= (NNet_input ? NNet_input->mTotalDelay : 0);
    
    // Loop over all feature frames.
    for (i = 0; i < header.mNSamples; i++, time++) 
    {
      size_t j;
      FLOAT* obs = feature_matrix[i];

      // Get next NN input vector by propagating vector from feature file
      // through NNet_input transformation.
      obs = XformPass(NNet_input, obs, time, FORWARD);
      
      // Input of NN is not yet read because of NNet_input delay
      if (time <= 0) {
        continue;
      }

      if (outlabel_map) {
        // Create NN output example vector from lables
        for (j = 0; j < NNet_instance->OutSize(); j++) {
          obs_out[j] = 0;
        }

        while (lbl_ptr && lbl_ptr->mStop < time) {
          lbl_ptr = lbl_ptr->mpNext;
        }

        if (lbl_ptr && lbl_ptr->mStart <= time) {
          // TODO: Don't know what this hack is... Find out!!! Caused
          // troubles on 64bit compilation
          obs_out[reinterpret_cast<size_t>(lbl_ptr->mpData) - 1] = 1;
        }
      } 
      else 
      {
        // Get NN output example vector from obsMx_out matrix
        obs_out = feature_matrix_out[time-1];
      }
      
      ///***********************************************************************
      /// For EACH NEW VECTOR - give it to SNet
      // file_name->logical    - filename (logical)
      // i                     - frame number
      // obs                   - xform output
      // obs_out               - wanted output
      // NNet_instance->OutSize() - dimenzionality
      // this returns true, if last vector
      ///*************************************************************************

      // just find max in both vec, add 1 if both on the same place
      // acctually computing accuracy...
      int max_ind_obs=0;
      int max_ind_obs_out=0;
      for (unsigned int kk=0; kk<NNet_input->OutSize(); kk++) {
	if (obs[kk] > obs[max_ind_obs]) {max_ind_obs=kk;}
	if (obs_out[kk] > obs_out[max_ind_obs_out]) {max_ind_obs_out=kk;}
      }
      if (max_ind_obs == max_ind_obs_out) {Hits++;}
      printf ("NN outputs = %i\tmax obs: %i\t label: %i\n", NNet_input->OutSize(),max_ind_obs,max_ind_obs_out);

    } //End of loop over all feature frames


    totFrames  += n_frames;
    //TraceLog("[%d frames]", n_frames);
    
    if (outlabel_map)
    {
      ReleaseLabels(labels);
    }

    feature_matrix_out.Destroy();
    feature_matrix.Destroy();
  }
  // MAIN FILE LOOP END

  ///***************************************************************************     
  /// DELETE SNET
  // delete prog_obj;  
  std::cout << "SNet is Deleted... debug message\n" << std::flush;
  
  if (trace_flag & 2) {
    TraceLog("Total number of frames: %d", totFrames);
  }
  NNet_instance->mpInput = NNet_input;

  if(!cross_validation)  
    hset.WriteMmf(trg_mmf, trg_hmm_dir, trg_hmm_ext, hmms_binary);
    
  // :KLUDGE: Some bug in STK - maybe solved?!
  /// hset.Release();

  for (size_t i = 0; i < cfgHash.mNEntries; i++) free(cfgHash.mpEntry[i]->data);
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
  
  while (feature_files) 
  {
    file_name = feature_files;
    feature_files = feature_files->mpNext;
    free(file_name);
  }
  
  return 0;
}
