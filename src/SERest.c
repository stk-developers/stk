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
#include "STKLib/labels.h"
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
//" -c f   Mixture pruning threshold                         10.0\n"
" -d s       Dir to find hmm definitions                     Current\n"
" -m N       Min examples needed per model                   3\n"
" -o s       Extension for new hmm files                     As src\n"
" -p N       Set parallel mode to N                          Off\n"
" -r         Enable Single Pass Training                     Off\n"
" -s s       Print statistics to file s                      Off\n"
" -t f [i l] Set pruning to f [inc limit]                    Off\n"
" -u tmvwsx  Update t)rans m)eans v)ars w)ghts s)tats x)form tmvwsx\n"
" -v f       Set minimum variance to f                       0.0\n"
" -w f       Set mix weight floor to f*MINMIX                1.0\n"
" -x s       Extension for hmm files                         None\n"
" -A         Print command line arguments                    Off\n"
" -B         Save HMM macro files as binary                  Off\n"
" -C cf      Set config file to cf                           Default\n"
" -D         Display configuration variables                 Off\n"
////" -F fmt   Set source data format to fmt                   as config
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

#define SNAME "SEREST"
char *optionStr =
" -d r   SOURCEMODELDIR"
" -m r   MINEGS"
" -o r   TARGETMODELEXT"
" -p r   PARALLELMODE"
" -r n   SINGLEPASSTRN=TRUE"
" -s r   SAVESTATS"
" -t ror PRUNING PRUNINGINC PRUNINGMAX"
" -u r   UPDATE"
" -v r   MINVAR"
" -w r   MIXWEIGHTFLOOR"
" -x r   SOURCEMODELEXT"
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
  HMMSet hset, *hset_alig = NULL;
  Network net;
  FILE *sfp, *ilfp = NULL;
  FLOAT *obsMx, *obsMx_alig, sentWeight;
  HTK_Header header, header_alig;
  int i, fcnt = 0;
  Label *labels;
  char line[1024];
  char label_file[1024];
  char *chrptr;
//  struct my_hsearch_data labelHash;
  struct my_hsearch_data nonCDphHash, phoneHash, dictHash, cfgHash;

  FLOAT totLogLike      = 0;
  FLOAT totLogPosterior = 0;
  int   totFrames       = 0;

  FileListElem *feature_files = NULL;
  int nfeature_files = 0;
  FileListElem *file_name = NULL;
  FileListElem **last_file = &feature_files;

  UpdateMask update_mask = 0;

  double word_penalty;
  double model_penalty;
  double grammar_scale;
  double transp_scale;
  double outprb_scale;
  double occprb_scale;
  double pronun_scale;
  double state_pruning;
  double stprn_step;
  double stprn_limit;
  double min_variance;
  double min_mix_wght;
  double E_constant;
  double h_constant;
  double I_smoothing;

  const char *cchrptr;
  const char *src_hmm_list;
  const char *src_hmm_dir;
  const char *src_hmm_ext;
        char *src_mmf;
  const char *alg_hmm_list;
  const char *alg_hmm_dir;
  const char *alg_hmm_ext;
        char *alg_mmf;
  const char *trg_hmm_dir;
  const char *trg_hmm_ext;
  const char *trg_mmf;
  const char *src_lbl_dir;
  const char *src_lbl_ext;
  const char *src_mlf;
  const char *network_file;
  const const char *xformList;
  const char *stat_file;
  char *script;
  const char *dictionary;
  const char *label_filter;
  const char *net_filter;
        char *cmn_path;
        char *cmn_file;
  const char *cmn_mask;
        char *cvn_path;
        char *cvn_file;
  const char *cvn_mask;
  const char *cvg_file;
        char *cmn_path_alig;
        char *cmn_file_alig;
  const char *cmn_mask_alig;
        char *cvn_path_alig;
        char *cvn_file_alig;
  const char *cvn_mask_alig;
  const char *cvg_file_alig;
  int  trace_flag;
  int  min_examples;
  int  parallel_mode;
  int  targetKind;
  int  targetKind_alig = PARAMKIND_ANON;
  int  derivOrder;
  int  derivOrder_alig;
  int  *derivWinLengths;
  int  *derivWinLengths_alig = NULL;
  int  startFrmExt;
  int  endFrmExt;
  int  startFrmExt_alig;
  int  endFrmExt_alig;
  BOOL viterbiTrain;
  BOOL xfStatsBin;
  BOOL hmms_binary;
  BOOL alg_mixtures;
  BOOL one_pass_reest;
  BOOL swap_features;
  BOOL swap_features_alig;
  BOOL htk_compat;
  enum {AT_ML=0, AT_MPE, AT_MFE} accum_type;
  enum {UT_ML=0, UT_MMI, UT_MPE} update_type;
  enum {UM_UPDATE=1, UM_DUMP=2, UM_BOTH=UM_UPDATE|UM_DUMP} update_mode;
  enum TranscriptionFormat {TF_HTK, TF_STK, TF_ERR} in_transc_fmt;
  enum NotInDictAction notInDictAction = 0;
  RHFBuffer rhfbuff                  = {0};
  RHFBuffer rhfbuff_alig             = {0};
  ExpansionOptions expOptions        = {0};
  STKNetworkOutputFormat out_net_fmt = {0};
  LabelFormat in_lbl_fmt = {0};
  in_lbl_fmt.TIMES_OFF = 1;

  if(argc == 1) usage(argv[0]);

  InitHMMSet(&hset, 1);

  if(!my_hcreate_r(100,  &dictHash)
  || !my_hcreate_r(100,  &phoneHash)
  || !my_hcreate_r(100,  &cfgHash)) {
    Error("Insufficient memory");
  }
  i = ParseOptions(argc, argv, optionStr, SNAME, &cfgHash);
  htk_compat   = GetParamBool(&cfgHash,SNAME":HTKCOMPAT",       FALSE);
  if(htk_compat) {
    if(argc == i) Error("HMM list file name expected");
    InsertConfigParam(&cfgHash,        SNAME":SOURCEHMMLIST", argv[i++], '-');
  }
  for(; i < argc; i++) {
    last_file = AddFileElem(last_file, argv[i]);
    nfeature_files++;
  }
  one_pass_reest=GetParamBool(&cfgHash,SNAME":SINGLEPASSTRN",   FALSE);
  targetKind   = GetDerivParams(&cfgHash, &derivOrder, &derivWinLengths,
                                &startFrmExt, &endFrmExt,
                                &cmn_path, &cmn_file, &cmn_mask,
                                &cvn_path, &cvn_file, &cvn_mask, &cvg_file,
                                SNAME":", one_pass_reest ? 2 : 0);

  if(one_pass_reest) {
    targetKind_alig = GetDerivParams(&cfgHash, &derivOrder_alig, &derivWinLengths_alig,
                                     &startFrmExt_alig, &endFrmExt_alig,
                                     &cmn_path_alig, &cmn_file_alig, &cmn_mask_alig,
                                     &cvn_path_alig, &cvn_file_alig, &cvn_mask_alig,
                                     &cvg_file_alig, SNAME":", 1);
  }
  xfStatsBin   = GetParamBool(&cfgHash,SNAME":XFORMSTATSBINARY",FALSE);
  xformList    = GetParamStr(&cfgHash, SNAME":XFORMLIST",       NULL);
  viterbiTrain = GetParamBool(&cfgHash,SNAME":VITERBITRAIN",    FALSE);
  expOptions.CD_phone_expansion =
                 GetParamBool(&cfgHash,SNAME":ALLOWXWRDEXP",    FALSE);
  expOptions.respect_pronun_var
               = GetParamBool(&cfgHash,SNAME":RESPECTPRONVARS", FALSE);
  expOptions.strict_timing
               = GetParamBool(&cfgHash,SNAME":EXACTTIMEMERGE",  FALSE);
  expOptions.no_optimization
               =!GetParamBool(&cfgHash,SNAME":MINIMIZENET",     FALSE);
  expOptions.remove_words_nodes
               = GetParamBool(&cfgHash,SNAME":REMEXPWRDNODES",  FALSE);
  in_lbl_fmt.TIMES_OFF =
                !GetParamBool(&cfgHash,SNAME":TIMEPRUNING",    FALSE);
  in_lbl_fmt.left_extent  = -100 * (long long) (0.5 + 1e5 *
                 GetParamFlt(&cfgHash, SNAME":STARTTIMESHIFT",  0.0));
  in_lbl_fmt.right_extent =  100 * (long long) (0.5 + 1e5 *
                 GetParamFlt(&cfgHash, SNAME":ENDTIMESHIFT",    0.0));
  swap_features_alig =
  swap_features=!GetParamBool(&cfgHash,SNAME":NATURALREADORDER",isBigEndian());
  filter_wldcrd= GetParamStr(&cfgHash, SNAME":HFILTERWILDCARD", "$");
  script_filter= GetParamStr(&cfgHash, SNAME":HSCRIPTFILTER",   NULL);
  hlist_filter = GetParamStr(&cfgHash, SNAME":HMMLISTFILTER",   NULL);
  MMF_filter   = GetParamStr(&cfgHash, SNAME":HMMDEFFILTER",    NULL);
  label_filter = GetParamStr(&cfgHash, SNAME":HLABELFILTER",    NULL);
  net_filter   = GetParamStr(&cfgHash, SNAME":HNETFILTER",      NULL);
  dict_filter  = GetParamStr(&cfgHash, SNAME":HDICTFILTER",     NULL);
  MMF_ofilter  = GetParamStr(&cfgHash, SNAME":HMMDEFOFILTER",   NULL);
  dictionary   = GetParamStr(&cfgHash, SNAME":SOURCEDICT",      NULL);
  grammar_scale= GetParamFlt(&cfgHash, SNAME":LMSCALE",         1.0);
  outprb_scale = GetParamFlt(&cfgHash, SNAME":OUTPSCALE",       1.0);
  transp_scale = GetParamFlt(&cfgHash, SNAME":TRANSPSCALE",     1.0);
  pronun_scale = GetParamFlt(&cfgHash, SNAME":PRONUNSCALE",     1.0);
  occprb_scale = GetParamFlt(&cfgHash, SNAME":OCCUPPSCALE",     1.0);
  word_penalty = GetParamFlt(&cfgHash, SNAME":WORDPENALTY",     0.0);
  model_penalty= GetParamFlt(&cfgHash, SNAME":MODELPENALTY",    0.0);
  network_file = GetParamStr(&cfgHash, SNAME":RECOGNET",        NULL);
  src_mlf      = GetParamStr(&cfgHash, SNAME":SOURCEMLF",       NULL);
  src_lbl_dir  = GetParamStr(&cfgHash, SNAME":SOURCETRANSCDIR", NULL);
  src_lbl_ext  = GetParamStr(&cfgHash, SNAME":SOURCETRANSCEXT", NULL);
  state_pruning= GetParamFlt(&cfgHash, SNAME":PRUNING",         0.0);
  stprn_step   = GetParamFlt(&cfgHash, SNAME":PRUNINGINC",      0.0);
  stprn_limit  = GetParamFlt(&cfgHash, SNAME":PRUNINGMAX",      0.0);
  trace_flag   = GetParamInt(&cfgHash, SNAME":TRACE",           0);
  stat_file    = GetParamStr(&cfgHash, SNAME":SAVESTATS",       NULL);
  min_examples = GetParamInt(&cfgHash, SNAME":MINEGS",          3);
  parallel_mode= GetParamInt(&cfgHash, SNAME":PARALLELMODE",   -1);
  min_variance = GetParamFlt(&cfgHash, SNAME":MINVAR",          0.0);
  min_mix_wght = GetParamFlt(&cfgHash, SNAME":MIXWEIGHTFLOOR",  1.0);
  hmms_binary  = GetParamBool(&cfgHash,SNAME":SAVEBINARY",      FALSE);
  script =(char*)GetParamStr(&cfgHash, SNAME":SCRIPT",          NULL);
  src_hmm_list = GetParamStr(&cfgHash, SNAME":SOURCEHMMLIST",   NULL);
  src_hmm_dir  = GetParamStr(&cfgHash, SNAME":SOURCEMODELDIR",  NULL);
  src_hmm_ext  = GetParamStr(&cfgHash, SNAME":SOURCEMODELEXT",  NULL);
  src_mmf=(char*)GetParamStr(&cfgHash, SNAME":SOURCEMMF",       NULL);
  alg_hmm_list = GetParamStr(&cfgHash, SNAME":ALIGNHMMLIST",    NULL);
  alg_hmm_dir  = GetParamStr(&cfgHash, SNAME":ALIGNMODELDIR",   NULL);
  alg_hmm_ext  = GetParamStr(&cfgHash, SNAME":ALIGNMODELEXT",   NULL);
  alg_mmf=(char*)GetParamStr(&cfgHash, SNAME":ALIGNMMF",        NULL);
  alg_mixtures = GetParamBool(&cfgHash,SNAME":ALIGNMIXTURES",   FALSE);
  trg_hmm_dir  = GetParamStr(&cfgHash, SNAME":TARGETMODELDIR",  NULL);
  trg_hmm_ext  = GetParamStr(&cfgHash, SNAME":TARGETMODELEXT",  NULL);
  trg_mmf      = GetParamStr(&cfgHash, SNAME":TARGETMMF",       NULL);

  update_mode  = GetParamEnum(&cfgHash,SNAME":UPDATEMODE",     UM_UPDATE,
                   "UPDATE",UM_UPDATE,"DUMP",UM_DUMP,"BOTH",UM_BOTH, NULL);

  accum_type   = GetParamEnum(&cfgHash,SNAME":ACCUMULATORTYPE", AT_ML,
                              "ML",AT_ML,"MPE",AT_MPE,"MFE",AT_MFE, NULL);

  update_type   = GetParamEnum(&cfgHash,SNAME":UPDATETYPE",     UT_ML,
                              "ML",UT_ML,"MPE",UT_MPE,"MMI",UT_MMI, NULL);

  E_constant   = GetParamFlt(&cfgHash, SNAME":EBWCONSTANTE",    2.0);
  h_constant   = GetParamFlt(&cfgHash, SNAME":EBWCONSTANTH",    2.0);
  I_smoothing  = GetParamFlt(&cfgHash, SNAME":ISMOOTHING",  update_type==UT_MPE
                                                            ?   50 : 200);

  in_transc_fmt= GetParamEnum(&cfgHash,SNAME":SOURCETRANSCFMT",
                              !network_file && htk_compat ? TF_HTK : TF_STK,
                              "HTK", TF_HTK, "STK", TF_STK, NULL);

  cchrptr      = GetParamStr(&cfgHash, SNAME":UPDATE",  "");
  while(*cchrptr) {
    switch(*cchrptr++) {
      case 't': update_mask |= UM_TRANSITION; break;
      case 'm': update_mask |= UM_MEAN;       break;
      case 'v': update_mask |= UM_VARIANCE;   break;
      case 'w': update_mask |= UM_WEIGHT;     break;
      case 's': update_mask |= UM_XFSTATS;    break;
      case 'x': update_mask |= UM_XFORM;      break;
      case 'o': update_mask |= UM_OLDMEANVAR; break;
      default:  Error("Unknown update flag '%c' (tmvwsx)", *cchrptr);
    }
  }
  if(GetParamBool(&cfgHash, SNAME":PRINTCONFIG", FALSE)) {
    PrintConfig(&cfgHash);
  }
  if(GetParamBool(&cfgHash, SNAME":PRINTVERSION", FALSE)) {
    puts("Version: "VERSION"\n");
  }

  GetParamBool(&cfgHash, SNAME":NFRAMEOUTPNORM", FALSE);


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
  if(alg_hmm_list != NULL || alg_mmf != NULL) {
    hset_alig = (HMMSet *) malloc(sizeof(HMMSet));
    if(!hset_alig) Error("Insufficient memory");

    InitHMMSet(hset_alig, 0);

    for(alg_mmf=strtok(alg_mmf, ","); alg_mmf != NULL; alg_mmf=strtok(NULL, ",")) {
      ReadHMMSet(alg_mmf, hset_alig, NULL);
    }
  } else {
    hset_alig = &hset;
  }
  if(parallel_mode == 0) one_pass_reest = 0;

  if((nfeature_files & 1) && one_pass_reest) {
    Error("Single pass re-estimation requires even number (two sets) of feature files");
  }

  char *ut_str = update_type == UT_MMI ? "MMI" :
                 update_type == UT_MPE ? "MPE" : "Chosen type of ";

  if(parallel_mode == -1 &&  update_type != UT_ML) {
    Error("%s update is not possible without using parallel mode", ut_str);
  }
  if((nfeature_files & 1) && update_type != UT_ML && parallel_mode == 0) {
    Error("%s update requires even number (two sets) of accumulator files", ut_str);
  }
  if(src_hmm_list) ReadHMMList(&hset,     src_hmm_list, src_hmm_dir, src_hmm_ext);
  if(alg_hmm_list) ReadHMMList(hset_alig, alg_hmm_list, alg_hmm_dir, alg_hmm_ext);
  nonCDphHash = MakeCIPhoneHash(&hset);

  if(dictionary != NULL) {
    ReadDictionary(dictionary, &dictHash, &phoneHash);
    notInDictAction  = WORD_NOT_IN_DIC_WARN;
    if(expOptions.respect_pronun_var) {
      notInDictAction |= PRON_NOT_IN_DIC_ERROR;
    }
  }
  if(dictHash.nentries == 0) expOptions.no_word_expansion = 1;

  transc_filter = transc_filter != NULL    ? transc_filter :
                  in_transc_fmt == TF_STK  ? net_filter    :
                                             label_filter;

  if(xformList != NULL) ReadXformList(&hset, xformList);

  AllocateAccumulatorsForXformStats(&hset);

  if(network_file) { // Unsupervised training
    //!!! Currently, it is not possible to get here !!!

    ilfp = fopen(network_file, "rt");
    if(ilfp  == NULL) Error("Cannot open network file: %s", network_file);

    Node *node = ReadSTKNetwork(ilfp, &dictHash, &phoneHash, 0, in_lbl_fmt,
                                header.sampPeriod, network_file, NULL);
    NetworkExpansionsAndOptimizations(node, expOptions, out_net_fmt, &dictHash,
                                      &nonCDphHash, &phoneHash);
    InitNetwork(&net, node, hset_alig, &hset);
    fclose(ilfp);
    min_examples = 0;
  } else {
    ilfp = OpenInputMLF(src_mlf);
  }

  hset.MMIUpdate           = update_type;

  ResetAccumsForHMMSet(&hset);

  hset.MMI_E               = E_constant;
  hset.MMI_h               = h_constant;
  hset.MMI_tauI            = I_smoothing;
  hset.gaussLvl2ModelReest = alg_mixtures;
  hset.minOccurances       = min_examples;
  hset.minMixWeight        = min_mix_wght * MIN_WEGIHT;
  hset.updateMask          = update_mask ? update_mask :
                             UM_TRANSITION | UM_MEAN    | UM_VARIANCE |
                             UM_WEIGHT     | UM_XFSTATS | UM_XFORM;

  if((hset.updateMask & (UM_MEAN | UM_VARIANCE)) &&
     !hset.allMixuresUpdatableFromStatAccums) {
    Warning("Statistic are estimated for Xform not being "
            "a single linear Xform on the input of a mixture. "
            "Means and variances will not be updated");
    hset.updateMask &= ~(UM_MEAN | UM_VARIANCE);
  }
  for(file_name = feature_files;
      file_name != NULL;
      file_name = one_pass_reest || (parallel_mode == 0 && update_type != UT_ML)
                  ? file_name->next->next : file_name->next) {
    if(trace_flag & 1) {
      if(one_pass_reest || (parallel_mode == 0 && update_type != UT_ML)) {
        TraceLog("Processing file pair %d/%d '%s' <-> %s",  ++fcnt,
                 nfeature_files / (one_pass_reest ? 2 : 1),
                 file_name->logical,file_name->next->logical);
      } else {
        TraceLog("Processing file %d/%d '%s'", ++fcnt,
                 nfeature_files,file_name->logical);
      }
    }
    if(parallel_mode == 0) {
      FLOAT P, P2;
      long S;

      ReadAccums(file_name->physical, 1.0, &hset, &S, &P, UT_ML);
      totFrames  += S;
      totLogLike += P;

      if(trace_flag & 1) {
        TraceLog("[%d frames] %f", S, P/S);
      }

      if(update_type != UT_ML) {
        // Second set of accums is for compeating models
        ReadAccums(file_name->next->physical, 1.0, &hset, &S, &P2, update_type);
        totLogPosterior += update_type == UT_MPE ? P2 : P - P2;

        if(trace_flag & 1) {
          TraceLog("[%d frames] %f", S, P2/S);
        }
      }
    } else {
      int nFrames;
      char *phys_fn = (one_pass_reest ? file_name->next : file_name)->physical;

      // read sentence weight definition if any ( physical_file.fea[s,e]{weight} )
      if((chrptr = strrchr(phys_fn, '{')) != NULL &&
        ((i=0), sscanf(chrptr, "{%f}%n", &sentWeight, &i), chrptr[i] == '\0')) {
        *chrptr = '\0';
      } else {
        sentWeight = 1.0;
      }

      if(cmn_mask) process_mask(file_name->logical, cmn_mask, cmn_file);
      if(cvn_mask) process_mask(file_name->logical, cvn_mask, cvn_file);
      obsMx = ReadHTKFeatures(phys_fn, swap_features,
                              startFrmExt, endFrmExt, targetKind,
                              derivOrder, derivWinLengths, &header,
                              cmn_path, cvn_path, cvg_file, &rhfbuff);

      if(hset.in_vec_size != header.sampSize / sizeof(float)) {
        Error("Vector size [%d] in '%s' is incompatible with source HMM set [%d]",
              header.sampSize/sizeof(float), phys_fn, hset.in_vec_size);
      }
      nFrames = header.nSamples - hset.totalDelay;

//      for(i = 0; i < nFrames; i++) {
//        int j;
//        for(j = 0; j < hset.in_vec_size; j++) {
//          printf("%5.2f ", obsMx[i * hset.in_vec_size +j]);
//        }
//        puts("");
//      }

      if(one_pass_reest) {
        if(cmn_mask_alig) process_mask(file_name->logical, cmn_mask_alig, cmn_file_alig);
        if(cvn_mask_alig) process_mask(file_name->logical, cvn_mask_alig, cvn_file_alig);
        obsMx_alig = ReadHTKFeatures(file_name->physical, swap_features_alig,
                                     startFrmExt_alig, endFrmExt_alig, targetKind_alig,
                                     derivOrder_alig, derivWinLengths_alig, &header_alig,
                                     cmn_path_alig, cvn_path_alig, cvg_file_alig, &rhfbuff_alig);

        if(hset_alig->in_vec_size != header_alig.sampSize/sizeof(float)) {
          Error("Vector size [%d] in '%s' is incompatible with alignment HMM set [%d]",
                header_alig.sampSize/sizeof(float), file_name->physical, hset_alig->in_vec_size);
        }
        if(nFrames != header_alig.nSamples - hset_alig->totalDelay) {
          if(hset_alig->totalDelay != hset.totalDelay) {
            printf("HPARM1 frames: %d+%ld+%d, alignment model delay: %d\n"
                   "HPARM2 frames: %d+%ld+%d, source model delay: %d\n",
                   startFrmExt_alig, header_alig.nSamples-startFrmExt_alig-endFrmExt_alig,
                   endFrmExt_alig, hset_alig->totalDelay,
                   startFrmExt, header.nSamples-startFrmExt-endFrmExt,
                   endFrmExt, hset.totalDelay);
          }
          Error("Mismatch in number of frames in single pass re-estimation "
                "feature file pair: '%s' <-> '%s'%s",
                file_name->next->physical, file_name->physical,
                hset_alig->totalDelay != hset.totalDelay ?
                ". Consider different delays of alignment/source HMM set!":"");
        }
      } else {
        header_alig = header;
        obsMx_alig  = obsMx;
      }
      if(!network_file) {
        Node *node = NULL;
        strcpy(label_file, file_name->logical);
        ilfp = OpenInputLabelFile(label_file, src_lbl_dir,
                                  src_lbl_ext ? src_lbl_ext :
                                  in_transc_fmt == TF_STK ? "net" : "lab",
                                  ilfp, src_mlf);

        if(in_transc_fmt == TF_HTK) {
          labels = ReadLabels(ilfp, dictionary ? &dictHash : &phoneHash,
                                    dictionary ? UL_ERROR : UL_INSERT, in_lbl_fmt,
                                    header.sampPeriod, label_file, src_mlf, NULL);
          node = MakeNetworkFromLabels(labels, dictionary ? NT_Word : NT_Phone);
          ReleaseLabels(labels);
        } else if(in_transc_fmt == TF_STK) {
          node = ReadSTKNetwork(ilfp, &dictHash, &phoneHash, 0, in_lbl_fmt,
                              header.sampPeriod, label_file, src_mlf);
        } else Error("Too bad. What did you do ?!?");

        NetworkExpansionsAndOptimizations(node, expOptions, out_net_fmt, &dictHash,
                                          &nonCDphHash, &phoneHash);
        InitNetwork(&net, node, hset_alig, &hset);
        CloseInputLabelFile(ilfp, src_mlf);
      }
      net.wPenalty     = word_penalty;
      net.mPenalty     = model_penalty;
      net.lmScale      = grammar_scale;
      net.pronScale    = pronun_scale;
      net.tranScale    = transp_scale;
      net.outpScale    = outprb_scale;
      net.ocpScale     = occprb_scale;
      net.pruningThresh= state_pruning > 0.0 ? state_pruning : -LOG_0;
      net.accumType    = accum_type;

      if(GetParamBool(&cfgHash, SNAME":NFRAMEOUTPNORM", FALSE)) {
        outprb_scale = 1.0 / nFrames;
      }

      FLOAT P;
      for(;;) {
        if(nFrames < 1) {
          Warning("Number of frames smaller than model delay, skipping file %s",
                  file_name->physical);
          P = LOG_MIN;
          break;
        }
printf("Using weight: %f", sentWeight);
        P = !viterbiTrain
          ? BaumWelchReest(&net, obsMx_alig, obsMx, nFrames, sentWeight)
          : ViterbiReest(  &net, obsMx_alig, obsMx, nFrames, sentWeight);

        if(P > LOG_MIN) break;

        if(net.pruningThresh <= LOG_MIN ||
          stprn_step <= 0.0 ||
          (net.pruningThresh += stprn_step) > stprn_limit ) {
          Warning("Overpruning or bad data, skipping file %s",
                  file_name->physical);
          break;
        }

        Warning("Overpruning or bad data in file %s, "
                "trying pruning threshold: %.2f",
                file_name->physical, net.pruningThresh);
      }
      if(P > LOG_MIN) {
       totFrames  += nFrames;
       totLogLike += P;
        if(trace_flag & 1) {
          TraceLog("[%d frames] %f", nFrames, P/nFrames);
        }
      }
      free(obsMx);
      if(one_pass_reest) free(obsMx_alig);
      if(!network_file)  ReleaseNetwork(&net);
    }
  }
  if(trace_flag & 2) {
    TraceLog("Total number of frames: %d\nTotal log likelihood: %e",
              totFrames, totLogLike);
    if(parallel_mode == 0 && update_type != UT_ML) {
      TraceLog("Total log posterior: %e", totLogPosterior);
    }
  }
  if(stat_file) WriteHMMStats(stat_file, &hset);

  if(parallel_mode != 0) DistributeMacroOccurances(&hset);

  if(parallel_mode > 0 || update_mode & UM_DUMP) {
    char accfn[32];
    snprintf(accfn, sizeof(accfn)-1, "SER%d.acc", HIGHER_OF(parallel_mode, 0));
    WriteAccums(accfn, trg_hmm_dir, &hset, totFrames, totLogLike);
  }
  if(parallel_mode <= 0 && update_mode & UM_UPDATE) {
    Macro *macro = FindMacro(&hset.variance_hash, "varFloor1");

    if(macro != NULL || (float) min_variance > 0.0) {
      Variance *tmpvar = macro ? (Variance *) macro->data : NULL;
//      assert(!tmpvar || hset.in_vec_size == tmpvar->vec_size);

      hset.varFloor = (Variance *) malloc(sizeof(Variance)+((tmpvar ? tmpvar->vec_size : hset.in_vec_size)-1)*sizeof(FLOAT));
      if(hset.varFloor == NULL) Error("Insufficient memory");

      hset.varFloor->vec_size = tmpvar ? tmpvar->vec_size : hset.in_vec_size;

      for(i = 0; i < hset.varFloor->vec_size; i++) {
        if(macro) {
          hset.varFloor->vector[i] =
            tmpvar && ((float) min_variance <= 0.0 ||
                        tmpvar->vector[i] < 1/min_variance)
            ? tmpvar->vector[i] : 1 / min_variance;
        }
      }
    }
    // Required by WriteXformStatsAndRunCommands and UpdateHMMSetFromAccums
    ScanHMMSet(&hset, mtm_mean|mtm_variance, NULL, NormalizeStatsForXform, 0);

    if(hset.updateMask & UM_XFSTATS) {
      WriteXformStatsAndRunCommands(trg_hmm_dir, xfStatsBin, &hset);
    }
    if(hset.updateMask != UM_XFSTATS) {
      if(hset.updateMask & UM_XFORM) {
        ReadXformStats(trg_hmm_dir, xfStatsBin, &hset);
      }
      UpdateHMMSetFromAccums(trg_hmm_dir, &hset);
      WriteHMMSet(trg_mmf, trg_hmm_dir, trg_hmm_ext, hmms_binary, &hset);
    }
  }
  if(hset_alig != &hset) {
    ReleaseHMMSet(hset_alig);
    free(hset_alig);
  }
  ReleaseHMMSet(&hset);

//  my_hdestroy_r(&labelHash, 0);
  my_hdestroy_r(&phoneHash, 1);
  my_hdestroy_r(&nonCDphHash, 0);
  FreeDictionary(&dictHash);
  for(i = 0; i < cfgHash.nentries; i++) free(cfgHash.entry[i]->data);
  my_hdestroy_r(&cfgHash, 1);

  if(network_file) {
    ReleaseNetwork(&net);
  }
  free(hset.varFloor);
  free(derivWinLengths);
  free(derivWinLengths_alig);
  if(src_mlf)  fclose(ilfp);
  while(feature_files) {
    file_name = feature_files;
    feature_files = feature_files->next;
    free(file_name);
  }
  return 0;
}

/* Accepted parameters
ACCEPTUNUSEDPARAM
ACCUMULATORTYPE
ACCWINDOW
ALIGNHMMLIST
ALIGNMIXTURES
ALIGNMMF
ALIGNMODELDIR
ALIGNMODELEXT
ALLOWXWRDEXP
DELTAWINDOW
DERIVWINDOWS
ENDTIMESHIFT
EXACTTIMEMERGE
HDICTFILTER
HLABELFILTER
HMMDEFFILTER
HMMDEFOFILTER
HMMLISTFILTER
HNETFILTER
HPARMFILTER
HTKCOMPAT
LABELFORMATING
LMSCALE
MINEGS
MINIMIZENET
MINVAR
MIXWEIGHTFLOOR
MODELPENALTY
NATURALREADORDER
NETFORMATING
OCCUPPSCALE
OUTPSCALE
PARALLELMODE
PRINTCONFIG
PRINTVERSION
PRONUNSCALE
PRUNING
PRUNINGINC
PRUNINGMAX
RECOGNET
REMEXPWRDNODES
RESPECTPRONVARS
SAVEBINARY
SAVESTATS
SCRIPT
SINGLEPASSTRN
SOURCEDICT
SOURCEHMMLIST
SOURCEHMMLIST
SOURCEMLF
SOURCEMMF
SOURCEMODELDIR
SOURCEMODELEXT
SOURCETRANSCDIR
SOURCETRANSCEXT
SOURCETRANSCFMT
STARTTIMESHIFT
TARGETKIND
TARGETMMF
TARGETMODELDIR
TARGETMODELEXT
THIRDWINDOW
TIMEPRUNNING
TRACE
TRANSPSCALE
UPDATE
VITERBITRAIN
WORDPENALTY
XFORMLIST
XFORMSTATSBINARY*/
