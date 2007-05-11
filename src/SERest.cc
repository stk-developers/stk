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

#define MODULE_VERSION "0.7 "__TIME__" "__DATE__" "SVN_ID  


#include "STKLib/fileio.h"
#include "STKLib/common.h"
#include "STKLib/Models.h"
#include "STKLib/Decoder.h"
#include "STKLib/labels.h"
#include "STKLib/stkstream.h"
#include "STKLib/Features.h"


#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <malloc.h>
#include <assert.h>

// Use internal version of getopt function if not available
#ifndef HAVE_UNISTD_H
#  include <unistd.h>
#else
#  include <getopt.h>
#endif

#include <iostream>
#include <string>
#include <sstream>

using namespace std;
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

int main(int argc, char* argv[]) 
{
  ModelSet                        hset;
  ModelSet*                       hset_alig  = NULL;
  ModelSet*                       hset_prior = NULL;
  Decoder<DecoderNetwork>         decoder;

  FeatureRepository               feature_repo;
  FeatureRepository               feature_repo_alig;

  FILE*                           sfp;
  FILE*                           ilfp = NULL;

  Label*                          labels;
  
  int                             i;
  int                             fcnt = 0;
  
  char                            line[1024];
  char                            label_file[1024];
  char*                           chrptr;
  
//  MyHSearchData labelHash;
  MyHSearchData                   nonCDphHash;
  MyHSearchData                   phoneHash;
  MyHSearchData                   dictHash;
  MyHSearchData                   cfgHash;

  FLOAT                           totLogLike      = 0;
  FLOAT                           totLogPosterior = 0;
  int                             totFrames       = 0;

//  FileListElem*                   feature_files = NULL;
//  int                             nfeature_files = 0;
//  FileListElem*                   file_name = NULL;
//  FileListElem**                  last_file = &feature_files;

  int                             update_mask = 0;

  double                          word_penalty;
  double                          model_penalty;
  double                          grammar_scale;
  double                          transp_scale;
  double                          outprb_scale;
  double                          occprb_scale;
  double                          pronun_scale;
  double                          state_pruning;
  double                          stprn_step;
  double                          stprn_limit;
  double                          min_variance;
  double                          min_mix_wght;
  double                          sig_slope;
  double                          E_constant;
  double                          h_constant;
  double                          I_smoothing;
  double                          MAP_tau;
  double                          min_occup;
        
  const char*                     cchrptr;
  const char*                     src_hmm_list;
  const char*                     src_hmm_dir;
  const char*                     src_hmm_ext;
        char*                     src_mmf;
  const char*                     alg_hmm_list;
  const char*                     alg_hmm_dir;
  const char*                     alg_hmm_ext;
        char*                     alg_mmf;
  const char*                     pri_hmm_list;
  const char*                     pri_hmm_dir;
  const char*                     pri_hmm_ext;
        char*                     pri_mmf;
  const char*                     trg_hmm_dir;
  const char*                     trg_hmm_ext;
  const char*                     trg_mmf;
  const char*                     src_lbl_dir;
  const char*                     src_lbl_ext;
  const char*                     src_mlf;
  const char*                     network_file;
  const char*                     xformList;
  const char*                     stat_file;
  char*                           script;
  const char*                     dictionary;
  const char*                     label_filter;
  const char*                     net_filter;
        char*                     cmn_path;
        char*                     cmn_file;
  const char*                     cmn_mask;
        char*                     cvn_path;
        char*                     cvn_file;
  const char*                     cvn_mask;
  const char*                     cvg_file;
        char*                     cmn_path_alig;
        char*                     cmn_file_alig;
  const char*                     cmn_mask_alig;
        char*                     cvn_path_alig;
        char*                     cvn_file_alig;
  const char*                     cvn_mask_alig;
  const char*                     cvg_file_alig;
  const char*                     mmf_dir;
  const char*                     mmf_mask;
  const char*                     mmf_dir_alig;
  const char*                     mmf_mask_alig;
  const char*                     mmf_karkulka;
  bool                            karkulka;
  bool                            cmllr_stats;  
  int                             trace_flag;
  int                             min_examples;
  int                             parallel_mode;
  int                             targetKind;
  int                             targetKind_alig = PARAMKIND_ANON;
  int                             deriv_order;
  int                             deriv_order_alig;
  int*                            deriv_win_lengths;
  int*                            deriv_win_lengths_alig = NULL;
  int                             start_frm_ext;
  int                             end_frm_ext;
  int                             start_frm_ext_alig;
  int                             end_frm_ext_alig;
  int                             max_active;
  int                             min_active;
  bool                            viterbiTrain;
  bool                            xfStatsBin;
  bool                            hmms_binary;
  bool                            save_glob_opt;
  bool                            alg_mixtures;
  bool                            one_pass_reest;
  bool                            swap_features;
  bool                            swap_features_alig;
  bool                            htk_compat;
  bool                            J_smoothing;
  double                          I_max_occup;;
  
  AccumType                       accum_type;
  
  Matrix<FLOAT>                   feature_matrix;
  Matrix<FLOAT>*                  feature_matrix_alig = NULL;
  
  
  enum UpdateType update_type;
  enum Update_Mode {UM_UPDATE=1, UM_DUMP=2, UM_BOTH=UM_UPDATE|UM_DUMP} update_mode;
  enum TranscriptionFormat {TF_HTK, TF_STK, TF_ERR} in_transc_fmt;
  int notInDictAction = (NotInDictActionType) WORD_NOT_IN_DIC_UNSET;
  
  RHFBuffer               rhfbuff           = {0};
  RHFBuffer               rhfbuff_alig      = {0};
  
  rhfbuff.mpLastFileName = NULL;
  rhfbuff_alig.mpLastFileName = NULL;
  
  ExpansionOptions        expOptions        = {0};
  STKNetworkOutputFormat  in_net_fmt        = {0};
//  LabelFormat             in_lbl_fmt        = {0};
//  in_lbl_fmt.TIMES_OFF = 1;


  // show help if no arguments specified
  if (argc == 1) 
    usage(argv[0]);


  
  if (!my_hcreate_r(100,  &dictHash)
   || !my_hcreate_r(100,  &phoneHash)
   || !my_hcreate_r(100,  &cfgHash)) 
  {
    Error("Insufficient memory");
  }
  
  i = ParseOptions(argc, argv, optionStr, SNAME, &cfgHash);
  if (htk_compat) {
    if (argc == i) Error("HMM list file name expected");
    InsertConfigParam(&cfgHash,        SNAME":SOURCEHMMLIST", argv[i++], '-');
  }
  
  htk_compat   = GetParamBool(&cfgHash,SNAME":HTKCOMPAT",       false);  
  script =(char*)GetParamStr(&cfgHash, SNAME":SCRIPT",          NULL);
  parallel_mode= GetParamInt(&cfgHash, SNAME":PARALLELMODE",   -1);
  one_pass_reest=GetParamBool(&cfgHash,SNAME":SINGLEPASSTRN",   false);
  if (parallel_mode == 0) {
    one_pass_reest = false;
  }
  update_type   = (UpdateType) GetParamEnum(&cfgHash,SNAME":UPDATETYPE",   UT_ML,
                              "ML",UT_ML,"MPE",UT_EBW,"MMI",UT_TwoAccumSetEBW, NULL);

  targetKind   = GetDerivParams(&cfgHash, &deriv_order, &deriv_win_lengths,
                                &start_frm_ext, &end_frm_ext,
                                &cmn_path, &cmn_file, &cmn_mask,
                                &cvn_path, &cvn_file, &cvn_mask, &cvg_file,
                                SNAME":", one_pass_reest ? 2 : 0);

  feature_repo.Init(swap_features, start_frm_ext, end_frm_ext,
                    targetKind, deriv_order, deriv_win_lengths,
		    cmn_path, cmn_mask, cvn_path, cvn_mask, cvg_file);
		    
  std::queue<FeatureRepository *> featureRepositoryList;
  featureRepositoryList.push(&feature_repo);

  if (one_pass_reest || (0 == parallel_mode && UT_TwoAccumSetEBW == update_type)) {
    targetKind_alig = GetDerivParams(&cfgHash, &deriv_order_alig, &deriv_win_lengths_alig,
                                     &start_frm_ext_alig, &end_frm_ext_alig,
                                     &cmn_path_alig, &cmn_file_alig, &cmn_mask_alig,
                                     &cvn_path_alig, &cvn_file_alig, &cvn_mask_alig,
                                     &cvg_file_alig, SNAME":", 1);

    feature_repo_alig.Init(swap_features_alig, start_frm_ext_alig, end_frm_ext_alig,
                           targetKind_alig, deriv_order_alig, deriv_win_lengths_alig,
			   cmn_path_alig, cmn_mask_alig, cvn_path_alig, cvn_mask_alig,
                           cvg_file_alig);
			   
    featureRepositoryList.push(&feature_repo_alig);
  }
  // the rest of the parameters are the feature files
  for (; i < argc; i++) {
    featureRepositoryList.front()->AddFile(argv[i]);
    swap(featureRepositoryList.front(), featureRepositoryList.back());
  }
  if (NULL != script) {
    AddFileListToFeatureRepositories(script, gpScriptFilter, featureRepositoryList);
  }
  
  xfStatsBin   = GetParamBool(&cfgHash,SNAME":XFORMSTATSBINARY",false);
  xformList    = GetParamStr(&cfgHash, SNAME":XFORMLIST",       NULL);
  viterbiTrain = GetParamBool(&cfgHash,SNAME":VITERBITRAIN",    false);
  expOptions.mCDPhoneExpansion =
                 GetParamBool(&cfgHash,SNAME":ALLOWXWRDEXP",    false);
  expOptions.mRespectPronunVar
               = GetParamBool(&cfgHash,SNAME":RESPECTPRONVARS", false);
  expOptions.mStrictTiming
               = GetParamBool(&cfgHash,SNAME":EXACTTIMEMERGE",  false);
  expOptions.mNoOptimization
               =!GetParamBool(&cfgHash,SNAME":MINIMIZENET",     false);
  expOptions.mRemoveWordsNodes
               = GetParamBool(&cfgHash,SNAME":REMEXPWRDNODES",  false);
//  in_lbl_fmt.TIMES_OFF =
  in_net_fmt.mNoTimes =
                !GetParamBool(&cfgHash,SNAME":TIMEPRUNING",     false);
//  in_lbl_fmt.left_extent  = -100 * (long long) (0.5 + 1e5 *
  in_net_fmt.mStartTimeShift =
                 GetParamFlt(&cfgHash, SNAME":STARTTIMESHIFT",  0.0);
//  in_lbl_fmt.right_extent =  100 * (long long) (0.5 + 1e5 *
  in_net_fmt.mEndTimeShift =
                 GetParamFlt(&cfgHash, SNAME":ENDTIMESHIFT",    0.0);
  swap_features_alig =
  swap_features=!GetParamBool(&cfgHash,SNAME":NATURALREADORDER",isBigEndian());
  gpFilterWldcrd= GetParamStr(&cfgHash, SNAME":HFILTERWILDCARD","$");
  gpScriptFilter= GetParamStr(&cfgHash, SNAME":HSCRIPTFILTER",  NULL);
  gpHListFilter = GetParamStr(&cfgHash, SNAME":HMMLISTFILTER",  NULL);
  gpMmfFilter   = GetParamStr(&cfgHash, SNAME":HMMDEFFILTER",   NULL);
  label_filter = GetParamStr(&cfgHash, SNAME":HLABELFILTER",    NULL);
  net_filter   = GetParamStr(&cfgHash, SNAME":HNETFILTER",      NULL);
  dict_filter  = GetParamStr(&cfgHash, SNAME":HDICTFILTER",     NULL);
  gpMmfOFilter = GetParamStr(&cfgHash, SNAME":HMMDEFOFILTER",   NULL);
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
  max_active   = GetParamInt(&cfgHash, SNAME":MAXACTIVEMODELS", 0);
  min_active   = GetParamInt(&cfgHash, SNAME":MINACTIVEMODELS", 0);
  trace_flag   = GetParamInt(&cfgHash, SNAME":TRACE",           0);
  stat_file    = GetParamStr(&cfgHash, SNAME":SAVESTATS",       NULL);
  min_examples = GetParamInt(&cfgHash, SNAME":MINEGS",          3);
  min_variance = GetParamFlt(&cfgHash, SNAME":MINVAR",          0.0);
  min_mix_wght = GetParamFlt(&cfgHash, SNAME":MIXWEIGHTFLOOR",  1.0);
  hmms_binary  = GetParamBool(&cfgHash,SNAME":SAVEBINARY",      false);
  save_glob_opt= GetParamBool(&cfgHash,SNAME":SAVEGLOBOPTS",    true);
  cmllr_stats  = GetParamBool(&cfgHash,SNAME":CMLLRSTATS",      false);
  src_hmm_list = GetParamStr(&cfgHash, SNAME":SOURCEHMMLIST",   NULL);
  src_hmm_dir  = GetParamStr(&cfgHash, SNAME":SOURCEMODELDIR",  NULL);
  src_hmm_ext  = GetParamStr(&cfgHash, SNAME":SOURCEMODELEXT",  NULL);
  src_mmf=(char*)GetParamStr(&cfgHash, SNAME":SOURCEMMF",       NULL);
  alg_hmm_list = GetParamStr(&cfgHash, SNAME":ALIGNHMMLIST",    NULL);
  alg_hmm_dir  = GetParamStr(&cfgHash, SNAME":ALIGNMODELDIR",   NULL);
  alg_hmm_ext  = GetParamStr(&cfgHash, SNAME":ALIGNMODELEXT",   NULL);
  alg_mmf=(char*)GetParamStr(&cfgHash, SNAME":ALIGNMMF",        NULL);
  alg_mixtures = GetParamBool(&cfgHash,SNAME":ALIGNMIXTURES",   false);
  pri_hmm_list = GetParamStr(&cfgHash, SNAME":PRIORHMMLIST",    NULL);
  pri_hmm_dir  = GetParamStr(&cfgHash, SNAME":PRIORMODELDIR",   NULL);
  pri_hmm_ext  = GetParamStr(&cfgHash, SNAME":PRIORMODELEXT",   NULL);
  pri_mmf=(char*)GetParamStr(&cfgHash, SNAME":PRIORMMF",        NULL);
  trg_hmm_dir  = GetParamStr(&cfgHash, SNAME":TARGETMODELDIR",  NULL);
  trg_hmm_ext  = GetParamStr(&cfgHash, SNAME":TARGETMODELEXT",  NULL);
  trg_mmf      = GetParamStr(&cfgHash, SNAME":TARGETMMF",       NULL);
  
  mmf_dir      = GetParamStr(&cfgHash, SNAME":MMFDIR",          ".");
  mmf_mask     = GetParamStr(&cfgHash, SNAME":MMFMASK",         NULL);
  mmf_dir_alig = GetParamStr(&cfgHash, SNAME":ALIGNMMFDIR",          ".");
  mmf_mask_alig= GetParamStr(&cfgHash, SNAME":ALIGNMMFMASK",         NULL);
  
  mmf_karkulka = GetParamStr(&cfgHash, SNAME":CWOUTDIR",        ".");
  karkulka     = GetParamBool(&cfgHash,SNAME":CWUPDATE",        false);
  
  
  update_mode  = (Update_Mode) GetParamEnum(&cfgHash,SNAME":UPDATEMODE",    UM_UPDATE,
                   "UPDATE",UM_UPDATE,"DUMP",UM_DUMP,"BOTH",UM_BOTH, NULL);

  accum_type   = (AccumType) GetParamEnum(&cfgHash,SNAME":ACCUMULATORTYPE", AT_ML,
                              "ML",AT_ML,"MPE", AT_MPE, "MMI", AT_MMI,"MCE",AT_MCE,"MFE",AT_MFE, NULL);


  sig_slope    = GetParamFlt(&cfgHash, SNAME":MCESIGSLOPE",    -1.0); // -1.0 ~ off
  E_constant   = GetParamFlt(&cfgHash, SNAME":EBWCONSTANTE",    2.0);
  h_constant   = GetParamFlt(&cfgHash, SNAME":EBWCONSTANTH",    2.0);
  I_smoothing  = GetParamFlt(&cfgHash, SNAME":ISMOOTHING",      200.0);
  J_smoothing  = GetParamBool(&cfgHash,SNAME":JSMOOTHING",      false);
  MAP_tau      = GetParamFlt(&cfgHash, SNAME":MAPTAU",          update_type == UT_ML ? 10 : 0);
  min_occup    = GetParamFlt(&cfgHash, SNAME":MINOCC",          0.0);
  I_max_occup  = GetParamFlt(&cfgHash, SNAME":IMAXOCCUP",       -1.0);

  in_transc_fmt= (TranscriptionFormat) GetParamEnum(&cfgHash,SNAME":SOURCETRANSCFMT",
                              !network_file && htk_compat ? TF_HTK : TF_STK,
                              "HTK", TF_HTK, "STK", TF_STK, NULL);

  cchrptr      = GetParamStr(&cfgHash, SNAME":UPDATE",  "");
  
  for (; *cchrptr; cchrptr++)
  {
    switch (*cchrptr) 
    {
      case 't': update_mask |= UM_TRANSITION; break;
      case 'm': update_mask |= UM_MEAN;       break;
      case 'v': update_mask |= UM_VARIANCE;   break;
      case 'w': update_mask |= UM_WEIGHT;     break;
      case 's': update_mask |= UM_XFSTATS;    break;
      case 'x': update_mask |= UM_XFORM;      break;
      case 'o': update_mask |= UM_OLDMEANVAR; break;
      case 'p': update_mask |= UM_MAP;        break;
      case 'c': update_mask |= UM_CWEIGHTS;   break;
      default:  Error("Unknown update flag '%c' (tmvwsx)", *cchrptr);
    }
  }
  
  if(update_mask == 0)
  {
    update_mask = UM_TRANSITION | UM_MEAN    | UM_VARIANCE |
                  UM_WEIGHT     | UM_XFSTATS | UM_XFORM ;
  }

  // initialize basic ModelSet
  hset.Init(update_type == UT_TwoAccumSetEBW ? MODEL_SET_WITH_TWO_ACCUM_SET : MODEL_SET_WITH_ACCUM);
  hset.mUpdateMask          = update_mask;

  // ***************************************************************************
  // Cluster adaptive training update
  if (update_mask & UM_CWEIGHTS)
  {
    hset.mClusterWeightsOutPath = mmf_karkulka;
  }

  bool print_all_options = GetParamBool(&cfgHash,SNAME":PRINTALLOPTIONS", false);
  
  if (GetParamBool(&cfgHash, SNAME":PRINTCONFIG", false)) 
    PrintConfig(&cfgHash);
  
  if (GetParamBool(&cfgHash, SNAME":PRINTVERSION", false)) 
    puts("Version: "MODULE_VERSION"\n");
  

  GetParamBool(&cfgHash, SNAME":NFRAMEOUTPNORM", false);


  if (!GetParamBool(&cfgHash,SNAME":ACCEPTUNUSEDPARAM", false)) 
  {
    CheckCommandLineParamUse(&cfgHash);
  }
  
  if (print_all_options) 
  {
    print_registered_parameters();
  }


    

  if (NULL != src_mmf)
  {
    for (src_mmf=strtok(src_mmf, ","); src_mmf != NULL; src_mmf=strtok(NULL, ",")) 
    {
      hset.ParseMmf(src_mmf, NULL);
    }
  }
  
  if (alg_hmm_list != NULL || alg_mmf != NULL) 
  {
    hset_alig = new ModelSet;
    hset_alig->Init();

    if (NULL != alg_mmf)
    {
      for (alg_mmf = strtok(alg_mmf, ","); alg_mmf != NULL; alg_mmf=strtok(NULL, ",")) 
      {
        hset_alig->ParseMmf(alg_mmf, NULL);
      }
    }
  } 
  else 
  {
    hset_alig = &hset;
  }
  
  if (pri_hmm_list != NULL || pri_mmf != NULL) 
  {
    hset_prior = new ModelSet;
    hset_prior->Init();

    for (pri_mmf = strtok(pri_mmf, ","); pri_mmf != NULL; pri_mmf=strtok(NULL, ",")) 
    {
      hset_prior->ParseMmf(pri_mmf, NULL);
    }
  } 
  else 
  {
    hset_prior = &hset;
  }
  
  hset.AttachPriors(hset_prior);

  if (one_pass_reest && feature_repo.QueueSize() != feature_repo_alig.QueueSize())
    Error("Single pass re-estimation requires even number (two sets) of feature files");

  if (parallel_mode == -1 &&  update_type == UT_TwoAccumSetEBW)
    Error("Two accumulator set MMI update is not possible without using -p 0");
    
  if (UT_TwoAccumSetEBW == update_type && 0 == parallel_mode
  && feature_repo.QueueSize() != feature_repo_alig.QueueSize())
    Error("Two accumulator set MMI update requires even number of accumulator files");
    
  if (src_hmm_list) 
    hset.ReadHMMList(src_hmm_list, src_hmm_dir ? src_hmm_dir : "", src_hmm_ext ? src_hmm_ext : "");
    
  if (alg_hmm_list) 
    hset_alig->ReadHMMList(alg_hmm_list, alg_hmm_dir, alg_hmm_ext);
    
  if (pri_hmm_list) 
    hset_prior->ReadHMMList(pri_hmm_list, pri_hmm_dir, pri_hmm_ext);
  
  nonCDphHash = hset.MakeCIPhoneHash();

  if (dictionary != NULL) 
  {
    ReadDictionary(dictionary, &dictHash, &phoneHash);
    notInDictAction  = WORD_NOT_IN_DIC_WARN;
    
    if (expOptions.mRespectPronunVar)
      notInDictAction |= PRON_NOT_IN_DIC_ERROR;
  }
  
  if (dictHash.mNEntries == 0) 
    expOptions.mNoWordExpansion = 1;

  transc_filter = transc_filter != NULL    ? transc_filter :
                  in_transc_fmt == TF_STK  ? net_filter    :
                                             label_filter;

  if (xformList != NULL) 
    hset.ReadXformList(xformList);

  hset.mCmllrStats = cmllr_stats;
  hset.AllocateAccumulatorsForXformStats();

  if (network_file) 
  { // Unsupervised training
//    Decoder<DecoderNetwork>::NetworkType::Node*   p_node = NULL;
//    Decoder<DecoderNetwork>::NetworkType              my_net(p_node);
    IStkStream                        input_stream;    
    
    input_stream.open(network_file, ios::in, transc_filter ? transc_filter : "");
    
    if (!input_stream.good())
    {
      Error("Cannot open network file: %s", network_file);
    }
    
    ilfp = input_stream.file();
    
    if (in_transc_fmt == TF_HTK) 
    {
      labels = ReadLabels(
          ilfp, 
          dictionary ? &dictHash : &phoneHash,
          dictionary ? UL_ERROR : UL_INSERT, 
          in_net_fmt,
          feature_repo.CurrentHeader().mSamplePeriod, // :TODO: mSamplePeriod should be probably taken from feature_repo_alig
          network_file, 
          NULL, 
          NULL);
          
//      my_net.BuildFromLabels( labels, dictionary ? NT_WORD : NT_PHONE);
      decoder.rNetwork().BuildFromLabels(labels, dictionary ? NT_WORD : NT_PHONE);

          
      ReleaseLabels(labels);
    }
    else if (in_transc_fmt == TF_STK) 
    {
      ReadSTKNetwork(
        ilfp, 
        &dictHash,
        &phoneHash, 
        notInDictAction, 
        in_net_fmt,
        feature_repo.CurrentHeader().mSamplePeriod, // :TODO: mSamplePeriod should be probably taken from feature_repo_alig
        network_file, 
        NULL,
        false,
        decoder.rNetwork());
    }
    else 
    {
      Error("Too bad. What did you do ?!?");
    }
                                
    decoder.rNetwork().ExpansionsAndOptimizations(
      expOptions, 
      in_net_fmt, 
      &dictHash,
      &nonCDphHash, 
      &phoneHash);
                                      
//    p_node = my_net.pFirst();
//    decoder.Init(p_node, hset_alig, &hset);
    decoder.Init(hset_alig, &hset/*, compactNetworkRepresentation*/);
    
    min_examples = 0;
  } 
  else 
  {
    ilfp = OpenInputMLF(src_mlf);
  }
  
  //
    
  hset.mUpdateType          = update_type;
  hset.mMinVariance         = min_variance;    ///< global minimum variance floor
  hset.MMI_E                = E_constant;
  hset.MMI_h                = h_constant;
  hset.MMI_tauI             = I_smoothing;
  hset.JSmoothing           = J_smoothing;
  hset.mISmoothingMaxOccup  = I_max_occup;
  hset.mMinOccupation       = min_occup;
  hset.mMapTau              = MAP_tau;
  hset.mGaussLvl2ModelReest = alg_mixtures;
  hset.mMinOccurances       = min_examples;
  hset.mMinMixWeight        = min_mix_wght * MIN_WEGIHT;
  hset.mUpdateMask          = update_mask;
  hset.mSaveGlobOpts        = save_glob_opt;
  hset.ResetAccums();  


  if ((hset.mUpdateMask & (UM_MEAN | UM_VARIANCE)) &&
     !hset.mAllMixuresUpdatableFromStatAccums) 
  {
    Warning("Statistic are estimated for Xform not being "
            "a single linear Xform on the input of a mixture. "
            "Means and variances will not be updated");
    hset.mUpdateMask &= ~(UM_MEAN | UM_VARIANCE);
  }
  

  for (feature_repo.Rewind(); !feature_repo.EndOfList(); feature_repo.MoveNext(), feature_repo_alig.MoveNext())
  {
//    FeatureRepository::ListIterator aux_alig_name(feature_repo.pCurrentRecord());

    if (trace_flag & 1) {
      if (one_pass_reest || (0 == parallel_mode && UT_TwoAccumSetEBW == update_type)) {
// get the name of the following feature file
//        aux_alig_name++;

        TraceLog("Processing file pair %d/%d '%s' <-> %s",  ++fcnt,
                 feature_repo.QueueSize(),
                 feature_repo_alig.Current().Logical().c_str(),
		 feature_repo.Current().Logical().c_str()
//                 aux_alig_name->Physical().c_str()
		 );
      }
      else {
        TraceLog("Processing file %d/%d '%s'", ++fcnt,
                 feature_repo.QueueSize(),
                 feature_repo.Current().Logical().c_str());
      }
    }
    
    if (0 == parallel_mode) 
    {
      FLOAT P;
      FLOAT P2;
      long  S;

      hset.ReadAccums(feature_repo_alig.Current(), &S, &P, UT_ML);

      totFrames  += S;
      totLogLike += P;

      if (trace_flag & 1)
        TraceLog("[%d frames] %f", S, P/S);

      if (UT_TwoAccumSetEBW == update_type) 
      {
        // Second set of accums is for compeating models
        hset.ReadAccums(feature_repo.Current(), &S, &P2, update_type);

        totLogPosterior += P - P2;

        if (trace_flag & 1)
          TraceLog("[%d frames] %f", S, P2/S);
      }
    } 
    else 
    {
      int     n_frames;
      int     n_samples_alig;

      // alignment features are read first .....................................
      if (one_pass_reest)
      {
        feature_matrix_alig = new Matrix<FLOAT>;

        feature_repo_alig.ReadFullMatrix(*feature_matrix_alig);
        
        n_samples_alig = feature_repo_alig.CurrentHeader().mNSamples;

        if (hset_alig->mInputVectorSize != static_cast<int>
            (feature_repo_alig.CurrentHeader().mSampleSize / sizeof(float))) 
        {
          Error("Vector size [%d] in '%s' is incompatible with alignment HMM set [%d]",
                feature_repo_alig.CurrentHeader().mSampleSize/sizeof(float),
                feature_repo_alig.Current().Physical().c_str(),
                hset_alig->mInputVectorSize);
        }
      }

      // now reest features are read ...........................................
      feature_repo.ReadFullMatrix(feature_matrix);
      
      if (hset.mInputVectorSize != static_cast<int>
          (feature_repo.CurrentHeader().mSampleSize / sizeof(float))) 
      {
        Error("Vector size [%d] in '%s' is incompatible with source HMM set [%d]",
              feature_repo.CurrentHeader().mSampleSize/sizeof(float), 
              feature_repo.Current().Physical().c_str(), 
              hset.mInputVectorSize);
      }
      
      n_frames = feature_repo.CurrentHeader().mNSamples - hset.mTotalDelay;
        
      if (one_pass_reest 
      && (n_frames != n_samples_alig - hset_alig->mTotalDelay))
      {
        if ((hset_alig->mTotalDelay != hset.mTotalDelay)) 
        {
          printf("HPARM1 frames: %d+%d+%d, alignment model delay: %d\n"
                 "HPARM2 frames: %d+%d+%d, source model delay: %d\n",
                 start_frm_ext_alig,
                 static_cast<int>(n_samples_alig - start_frm_ext_alig
                   - end_frm_ext_alig),
                 end_frm_ext_alig, 
                 hset_alig->mTotalDelay,
                 start_frm_ext, 
                 static_cast<int>(feature_repo.CurrentHeader().mNSamples 
                   - start_frm_ext - end_frm_ext),
                 end_frm_ext, 
                 hset.mTotalDelay);
        }
        
        Error("Mismatch in number of frames in single pass re-estimation "
              "feature file pair: '%s' <-> '%s'%s",
              feature_repo_alig.Current().Physical().c_str(),
              feature_repo.Current().Physical().c_str(),
//              aux_alig_name->Physical().c_str(),
              hset_alig->mTotalDelay != hset.mTotalDelay ?
              ". Consider different delays of alignment/source HMM set!":"");
      }

      // use feature_matrix for alignment
      if (!one_pass_reest) 
      {
        feature_matrix_alig = &feature_matrix;
      }
      
      if(mmf_mask_alig != NULL) 
      {
        static string lastSpeakerMMF;
        string speakerMMF;
        ProcessMask(feature_repo_alig.Current().Logical().c_str(), mmf_mask_alig, speakerMMF);
        
        if(lastSpeakerMMF != speakerMMF) 
        {
          hset_alig->ParseMmf((string(mmf_dir_alig) + "/" + speakerMMF).c_str(), NULL);
          lastSpeakerMMF = speakerMMF;
        }
      }
      
      if(mmf_mask != NULL) 
      {
        static string lastSpeakerMMF;
        string speakerMMF;
        ProcessMask(feature_repo.Current().Logical().c_str(), mmf_mask, speakerMMF);
        
        if(lastSpeakerMMF != speakerMMF) 
        {
          hset.ParseMmf((string(mmf_dir) + "/" + speakerMMF).c_str(), NULL);
          lastSpeakerMMF = speakerMMF;
        }
      }

// //////////////////////////////////////////////////////////////////
      
      
      if (!network_file) 
      {
//        Decoder<DecoderNetwork>::NetworkType::Node* p_node = NULL;
//        Decoder<DecoderNetwork>::NetworkType            my_net(p_node);

        strcpy(label_file, feature_repo.Current().Logical().c_str());
        
        ilfp = OpenInputLabelFile(
            label_file, 
            src_lbl_dir,
            src_lbl_ext ? src_lbl_ext : in_transc_fmt == TF_STK ? "net" : "lab",
            ilfp, 
            src_mlf);

        if (in_transc_fmt == TF_HTK) 
        {
          labels = ReadLabels(
              ilfp, 
              dictionary ? &dictHash : &phoneHash,
              dictionary ? UL_ERROR : UL_INSERT, 
              in_net_fmt,
              feature_repo.CurrentHeader().mSamplePeriod, // :TODO: mSamplePeriod should be probably taken from feature_repo_alig
              label_file, 
              src_mlf, 
              NULL);
              
//          my_net.BuildFromLabels(labels, dictionary ? NT_WORD : NT_PHONE);
          decoder.rNetwork().BuildFromLabels(labels, dictionary ? NT_WORD : NT_PHONE);

          ReleaseLabels(labels);
        } 
        else if (in_transc_fmt == TF_STK) 
        {
          ReadSTKNetwork(
              ilfp, 
              &dictHash, 
              &phoneHash, 
              notInDictAction, 
              in_net_fmt,
              feature_repo.CurrentHeader().mSamplePeriod, // :TODO: mSamplePeriod should be probably taken from feature_repo_alig 
              label_file, 
              src_mlf, false, decoder.rNetwork());
        } 
        else 
        {
          Error("Too bad. What did you do ?!?");
        }

        decoder.rNetwork().ExpansionsAndOptimizations(
          expOptions,
          in_net_fmt,
          &dictHash,
          &nonCDphHash,
          &phoneHash);
            
//        my_net.ExpansionsAndOptimizations(expOptions, in_net_fmt, 
//            &dictHash, &nonCDphHash, &phoneHash);

//        p_node = my_net.pFirst();
//        decoder.Init(p_node, hset_alig, &hset);

        decoder.Init(hset_alig, &hset/*, compactNetworkRepresentation*/);

        CloseInputLabelFile(ilfp, src_mlf);
      }
      
      decoder.mWPenalty        = word_penalty;
      decoder.mMPenalty        = model_penalty;
      decoder.mLmScale         = grammar_scale;
      decoder.mPronScale       = pronun_scale;
      decoder.mTranScale       = transp_scale;
      decoder.mOutpScale       = outprb_scale;
      decoder.mOcpScale        = occprb_scale;
      decoder.mPruningThresh   = state_pruning > 0.0 ? state_pruning : -LOG_0;
      decoder.mMaxActiveModels = max_active;
      decoder.mMinActiveModels = min_active;
      decoder.mAccumType       = accum_type;
      
      double prn_step   = stprn_step;
      double prn_limit  = stprn_limit;

      if (GetParamBool(&cfgHash, SNAME":NFRAMEOUTPNORM", false)) 
      {
        decoder.mOutpScale  = outprb_scale / n_frames;
        decoder.mPruningThresh /= n_frames;
        prn_step       /= n_frames;
        prn_limit      /= n_frames;
      }

      FLOAT P;
      for (;;) 
      {
        if (n_frames < 1) 
        {
          Warning("Number of frames is smaller than model delay, skipping file %s",
                  feature_repo.Current().Logical().c_str());
          P = LOG_MIN;
          break;
        }
        
        if (accum_type == AT_MCE || accum_type == AT_MMI) 
        {
          P = decoder.MCEReest(*feature_matrix_alig, feature_matrix, n_frames, 
              feature_repo.Current().Weight(), sig_slope);
        }
        else 
        {
          P = !viterbiTrain
            ? decoder.BaumWelchReest(*feature_matrix_alig, feature_matrix, 
                n_frames, feature_repo.Current().Weight())
            : decoder.ViterbiReest  (*feature_matrix_alig, feature_matrix, 
                n_frames, feature_repo.Current().Weight());
        }
        
        if (P > LOG_MIN)
          break;

        if (decoder.mPruningThresh <= LOG_MIN ||
          prn_step <= 0.0 ||
          (decoder.mPruningThresh += prn_step) > prn_limit) 
        {
          Warning("Overpruning or bad data, skipping file %s",
                  feature_repo.Current().Logical().c_str());
          break;
        }

        Warning("Overpruning or bad data in file %s, "
                "trying pruning threshold: %.2f",
                feature_repo.Current().Logical().c_str(), decoder.mPruningThresh);
      }
      
      ClusterWeightAccumUserData cwa = {hset.mNClusterWeightVectors, hset.mpGw, 
        hset.mpKw};

      // here
      if (hset.mUpdateMask  & UM_CWEIGHTS)
      {
        hset.Scan(MTM_MIXTURE, NULL, ComputeClusterWeightVectorAccums, &cwa);
      }
      
      if (P > LOG_MIN) 
      {
        totFrames  += n_frames;
        totLogLike += P;
        if (trace_flag & 1) 
          TraceLog("[%d frames] %f (Act: %f/%f)", n_frames, P/n_frames,
                   (float) decoder.mNActiveModelsForUtterance / n_frames, 
                   (float) decoder.mNActiveTokensForUtterance / n_frames);
      }
      
      feature_matrix.Destroy();
      if (one_pass_reest)
        feature_matrix_alig->Destroy();
        
      if (!network_file)  
        decoder.Clear();
    }
  }

  // save unsaved cluster mean weights
  if (hset.mUpdateMask & UM_CWEIGHTS)
  {
    for (int i = 0; i < hset.mNClusterWeightVectors; i++)
    {
      hset.ComputeClusterWeightsVector(i);
      hset.WriteClusterWeightsVector(i);
    }    
    hset.mClusterWeightsStream.close();
  }
  
  
  if (trace_flag & 2) 
  {
    TraceLog("Total number of frames: %d\nTotal log likelihood: %e",
              totFrames, totLogLike);
    
    if (0 == parallel_mode && UT_TwoAccumSetEBW == update_type)
      TraceLog("Total log posterior: %e", totLogPosterior);
  }
  
  if (stat_file) 
    hset.WriteHMMStats(stat_file);
  
  if (parallel_mode != 0) 
    hset.DistributeMacroOccurances();

  if (parallel_mode > 0 || update_mode & UM_DUMP) 
  {  
    char accfn[32];
    snprintf(accfn, sizeof(accfn)-1, "SER%d.acc", HIGHER_OF(parallel_mode, 0));
    hset.WriteAccums(accfn, trg_hmm_dir, totFrames, totLogLike);
  }
  
  if (parallel_mode <= 0 && update_mode & UM_UPDATE) 
  {
    Macro*        macro;
    Variance**    vector_to_update;
    size_t        tmp_var_floor_size;
    string        suffix = "";
    string        macro_name;
    int           m;
    
    
    // we'll go throuch modelset and all XformInstances and assign a
    // varfloor if defined... i=-1 indicates that we're inspecting the global
    // varFloor
    for (m = -1; m < static_cast<int>(hset.mXformInstanceHash.mNEntries); m++) 
    {
      if (m == -1)
      {
        suffix = "";
        tmp_var_floor_size = hset.mInputVectorSize;
        vector_to_update = &(hset.mpVarFloor);
      }
      else
      {
        macro  = reinterpret_cast<Macro*>(hset.mXformInstanceHash.mpEntry[m]->data);
        suffix = macro->mpName;
        XformInstance* xfi =  reinterpret_cast <XformInstance*>(macro->mpData);
        tmp_var_floor_size =  xfi->OutSize();  
        vector_to_update   = &xfi->mpVarFloor;
      }
        
      macro_name = "varFloor1" + suffix;
      macro = FindMacro(&hset.mVarianceHash, macro_name.c_str());
      if(macro) {
        *vector_to_update = (Variance *) macro->mpData;
        
        if((*vector_to_update)->VectorSize() != tmp_var_floor_size) 
        {
          Error("Ivalid size of variance floor vector '%s'", macro_name.c_str());
        }
      }
    }
    
    // Required by WriteXformStatsAndRunCommands and UpdateHMMSetFromAccums
    if(!hset.mCmllrStats) {
      hset.Scan(MTM_MEAN|MTM_VARIANCE, NULL, NormalizeStatsForXform, 0);
    }

    if (hset.mUpdateMask & UM_XFSTATS) 
      hset.WriteXformStatsAndRunCommands(trg_hmm_dir, xfStatsBin);
    
    if (hset.mUpdateMask != UM_XFSTATS) 
    {
      if (hset.mUpdateMask & UM_XFORM) 
      {
        hset.ReadXformStats(trg_hmm_dir, xfStatsBin);
      }
      
      //UpdateHMMSetFromAccums(trg_hmm_dir, &hset);
      hset.UpdateFromAccums(trg_hmm_dir);
      hset.WriteMmf(trg_mmf, trg_hmm_dir, trg_hmm_ext, hmms_binary);
    }
  }
  
  if (hset_alig != &hset) 
  {
    hset_alig->Release();
    delete hset_alig;
  }
  
  if (hset_prior != &hset) 
  {
    hset_prior->Release();
    delete hset_prior;
  }

  hset.Release();

    
//  my_hdestroy_r(&labelHash, 0);
  my_hdestroy_r(&phoneHash, 0);
  my_hdestroy_r(&nonCDphHash, 0);
  FreeDictionary(&dictHash);
  
  for (size_t i = 0; i < cfgHash.mNEntries; i++) 
    free(cfgHash.mpEntry[i]->data);
    
  my_hdestroy_r(&cfgHash, 1);

  if (network_file) 
    decoder.Clear();
  
//  delete hset.mpVarFloor;
  free(deriv_win_lengths);
  free(deriv_win_lengths_alig);
  
  if (src_mlf) 
    fclose(ilfp);
                                                              
  if (NULL != rhfbuff.mpFp)
    fclose(rhfbuff.mpFp);
  
  if (NULL != rhfbuff_alig.mpFp)
    fclose(rhfbuff_alig.mpFp);
    
    
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
