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

#define VERSION "0.6 "__TIME__" "__DATE__

#include "STKLib/Viterbi.h"
#include "STKLib/Models.h"
#include "STKLib/fileio.h"
#include "STKLib/labels.h"
#include "STKLib/common.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <malloc.h>
#ifndef WIN32
#include <unistd.h>
#else
#include "getopt.h"
#endif

using namespace STK;


typedef struct _LRTrace LRTrace;
struct _LRTrace {
  Node *mpWordEnd;
  FLOAT lastLR;
  FLOAT candidateLR;
  long long candidateStartTime;
  long long candidateEndTime;
};

LRTrace *MakeLRTraceTable(Network *net, int *nWords, Token **filler_end)
{
  LRTrace *lrt;
  Node *node;
  int i;

  *nWords = 0;

  for (node=net->mpFirst; node != NULL; node = node->mpNext) {
    if (node->mType == (NT_WORD | NT_STICKY) && node->mpPronun == NULL) break;
  }
  if (node == NULL) {
    Error("Network contains no null node with flag=F (reference model end)");
  }
  *filler_end = node->mpExitToken;

  for (node=net->mpFirst; node != NULL; node = node->mpNext) {
    if (node->mType == (NT_WORD | NT_STICKY) && node->mpPronun != NULL) ++*nWords;
  }
  if (*nWords == 0) {
    Error("Network contains no word nodes with flag=K, (keyword model end)");
  }
  if ((lrt = (LRTrace *) malloc(*nWords * sizeof(LRTrace))) == NULL) {
    Error("Insufficient memory");
  }
  for (i = 0, node=net->mpFirst; node != NULL; node = node->mpNext) {
    if (node->mType == (NT_WORD | NT_STICKY) && node->mpPronun != NULL) {
      lrt[i].mpWordEnd = node;
      i++;
    }
  }
  return lrt;
}

void PutCandidateToLabels(LRTrace *lrt, FLOAT scoreThreshold, FILE *lfp,
                          const char *label_file, const char *out_MLF,
                          LabelFormat out_lbl_fmt, long sampPeriod)
{
  if (lrt->candidateEndTime != 0
     && (scoreThreshold < LOG_MIN || lrt->candidateLR > scoreThreshold)) {
     Label label  = init_label;
     label.mStart = lrt->candidateStartTime;
     label.mStop  = lrt->candidateEndTime;
     label.mpName = lrt->mpWordEnd->mpPronun->mpWord->mpName;
     label.mScore = lrt->candidateLR;

     WriteLabels(lfp, &label, out_lbl_fmt, sampPeriod, label_file, out_MLF);
  }
}



char *longopts =
" --accept-unused-param"
" --acc-window         "
" --allow-xwrd-exp     \n"
" --delta-window       "
" --end-frm-ext        "
" --end-time-shift     \n"
" --evaluation         "
" --exact-time-merge   "
" --hdict-filter       \n"
" --hfilter-wildcard   "
" --hlabel-ofilter     "
" --hmmdef-filter      \n"
" --hmmlist-filter     "
" --hnet-filter        "
" --hparm-filter       \n"
" --keyword-threshold  "
" --label-formating    "
" --lm-scale           \n"
" --minimize-net       "
" --model-penalty      "
" --natural-read-order \n"
" --outp-scale         "
" --print-config       "
" --print-version      \n"
" --pronun-scale       "
" --pruning            "
" --recog-net          \n"
" --rem-exp-wrd-nodes  "
" --respect-pron-vars  "
" --script             \n"
" --source-dict        "
" --source-hmmlist     "
" --source-mlf         \n"
" --source-mmf         "
" --source-model-dir   "
" --source-model-ext   \n"
" --sourcetranscdir    "
" --source-transc-ext  "
" --start-frm-ext      \n"
" --start-time-shift   "
" --target-kind        "
" --target-mlf         \n"
" --target-transc-dir  "
" --target-transc-ext  "
" --third-window       \n"
" --time-pruning       "
" --trace              "
" --transp-scale       \n"
" --word-penalty       \n";


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
" -o s       Output label formating NCST                     None\n"
" -t f       Set pruning to f [inc limit]                    Off\n"
" -w f       Load single reognition network f                Off\n"
" -x s       Extension for hmm files                         None\n"
" -y s       Output transcription file extension             rec\n"
" -A         Print command line arguments                    Off\n"
" -C cf      Set config file to cf                           Default\n"
" -D         Display configuration variables                 Off\n"
//" -G fmt     Set source trascription format to fmt           As config\n"
" -H mmf     Load HMM macro file mmf\n"
" -I mlf     Load master transcription file mlf\n"
" -L dir     Set input transcription dir                     Current\n"
//" -P fmt     Set target transcription format to fmt          As config\n"
" -S f       Set script file to f                            None\n"
" -T N       Set trace flags to N                            0\n"
" -V         Print version information                       Off\n"
" -X ext     Set input transcription file ext                lab\n"
"\n"
"Long options (configuration parameters):\n\n"
"%s\n"
" %s is Copyright (C) 2004-2005 Lukas Burget et al. and\n"
" licensed under the GNU General Public License, version 2.\n"
" Bug reports, feedback, etc, to: burget@fit.vutbr.cz\n"
"\n", progname, longopts, progname);
  exit(-1);
}

#define SNAME "SLRATIO"
char *optionStr =
" -d r   SOURCEMODELDIR"
" -i r   TARGETMLF"
" -l r   TARGETTRANSCDIR"
" -o r   LABELFORMATING"
" -t r   PRUNING"
" -w r   RECOGNET"
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

//" -G r   SOURCETRANSCFMT"
//" -P r   TARGETTRANSCFMT"


int main(int argc, char *argv[]) {
  HtkHeader header;
  ModelSet hset;
  Network net;
  FILE *sfp, *lfp = NULL, *ilfp = NULL;
  FLOAT  *obsMx;
  int i, fcnt = 0;
  char line[1024];
  char label_file[1024];
  const char *cchrptr;
  MyHSearchData nonCDphHash, phoneHash, dictHash, cfgHash;

  FileListElem *feature_files = NULL;
  int nfeature_files = 0;
  FileListElem *file_name = NULL;
  FileListElem **last_file = &feature_files;

  AlignmentType alignment = WORD_ALIGNMENT;

  double score_thresh;
  double word_penalty;
  double model_penalty;
  double grammar_scale;
  double transp_scale;
  double outprb_scale;
  double pronun_scale;
  double occprb_scale;
  double state_pruning;
//double stprn_step;
//double stprn_limit;
//double word_pruning;
  const char *hmm_dir;
  const char *hmm_ext;
  const char *out_lbl_dir;
  const char *out_lbl_ext;
  const char *in_lbl_dir;
  const char *in_lbl_ext;
  const char *out_MLF;
  const char *in_MLF;
  const char *network_file;
  const char *hmm_list;
  const char *dictionary;
  char *script;
  char *mmf;
//  const char *label_filter;
//  const char *net_filter;
//  const char *label_ofilter;
//  const char *net_ofilter;
        char *cmn_path;
        char *cmn_file;
  const char *cmn_mask;
        char *cvn_path;
        char *cvn_file;
  const char *cvn_mask;
  const char *cvg_file;
  int  trace_flag;
  int  targetKind;
  int  derivOrder;
  int  *derivWinLengths;
  int startFrmExt;
  int endFrmExt;
  bool fulleval;
  bool swap_features;
//  enum TranscriptionFormat {TF_HTK, TF_STK} in_transc_fmt, out_transc_fmt;
  int notInDictAction = 0;
  RHFBuffer rhfbuff                    = {0};
  ExpansionOptions expOptions          = {0};
  STKNetworkOutputFormat in_net_fmt    = {0};
//  STKNetworkOutputFormat out_net_fmt   = {0};
  LabelFormat out_lbl_fmt              = {0};
  LabelFormat in_lbl_fmt               = {0};
  in_lbl_fmt.TIMES_OFF = 1;

  int nWords, j;
  LRTrace *lrt = NULL;
  Token *filler_end;

  if (argc == 1) usage(argv[0]);

  //InitHMMSet(&hset, 0);
  hset.Init();

  if (!my_hcreate_r(100,  &dictHash)
  || !my_hcreate_r(100,  &phoneHash)
  || !my_hcreate_r(100,  &cfgHash)) {
    Error("Insufficient memory");
  }
  i = ParseOptions(argc, argv, optionStr, SNAME, &cfgHash);

  for (; i < argc; i++) {
    last_file = AddFileElem(last_file, argv[i]);
    nfeature_files++;
  }
  targetKind =   GetDerivParams(&cfgHash, &derivOrder, &derivWinLengths,
                                &startFrmExt, &endFrmExt,
                                &cmn_path, &cmn_file, &cmn_mask,
                                &cvn_path, &cvn_file, &cvn_mask, &cvg_file,
                                SNAME":", 0);
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
  in_lbl_fmt.TIMES_OFF =
                !GetParamBool(&cfgHash,SNAME":TIMEPRUNING",    false);
  in_lbl_fmt.left_extent  = -100 * (long long) (0.5 + 1e5 *
                 GetParamFlt(&cfgHash, SNAME":STARTTIMESHIFT",  0.0));
  in_lbl_fmt.right_extent =  100 * (long long) (0.5 + 1e5 *
                 GetParamFlt(&cfgHash, SNAME":ENDTIMESHIFT",    0.0));
  fulleval     = GetParamBool(&cfgHash,SNAME":EVALUATION",      false);
  swap_features=!GetParamBool(&cfgHash,SNAME":NATURALREADORDER",isBigEndian());
  gpFilterWldcrd= GetParamStr(&cfgHash, SNAME":HFILTERWILDCARD", "$");
  gpScriptFilter= GetParamStr(&cfgHash, SNAME":HSCRIPTFILTER",   NULL);  
  gpParmFilter  = GetParamStr(&cfgHash, SNAME":HPARMFILTER",     NULL);
  gpHListFilter = GetParamStr(&cfgHash, SNAME":HMMLISTFILTER",   NULL);
  gpMmfFilter   = GetParamStr(&cfgHash, SNAME":HMMDEFFILTER",    NULL);
//  label_filter = GetParamStr(&cfgHash, SNAME":HLABELFILTER",    NULL);
  transc_filter= GetParamStr(&cfgHash, SNAME":HNETFILTER",      NULL);
  dict_filter  = GetParamStr(&cfgHash, SNAME":HDICTFILTER",     NULL);
  transc_ofilter=GetParamStr(&cfgHash, SNAME":HLABELOFILTER",   NULL);
//  net_ofilter  = GetParamStr(&cfgHash, SNAME":HNETOFILTER",     NULL);
  dictionary   = GetParamStr(&cfgHash, SNAME":SOURCEDICT",      NULL);
  hmm_list     = GetParamStr(&cfgHash, SNAME":SOURCEHMMLIST",   NULL);
  grammar_scale= GetParamFlt(&cfgHash, SNAME":LMSCALE",         1.0);
  outprb_scale = GetParamFlt(&cfgHash, SNAME":OUTPSCALE",       1.0);
  transp_scale = GetParamFlt(&cfgHash, SNAME":TRANSPSCALE",     1.0);
  pronun_scale = GetParamFlt(&cfgHash, SNAME":PRONUNSCALE",     1.0);
  occprb_scale = GetParamFlt(&cfgHash, SNAME":OCCUPPSCALE",     1.0);
  word_penalty = GetParamFlt(&cfgHash, SNAME":WORDPENALTY",     0.0);
  model_penalty= GetParamFlt(&cfgHash, SNAME":MODELPENALTY",    0.0);
  score_thresh = GetParamFlt(&cfgHash, SNAME":KEYWORDTHRESHOLD",LOG_0);
  network_file = GetParamStr(&cfgHash, SNAME":RECOGNET",        NULL);
  hmm_dir      = GetParamStr(&cfgHash, SNAME":SOURCEMODELDIR",  NULL);
  hmm_ext      = GetParamStr(&cfgHash, SNAME":SOURCEMODELEXT",  NULL);
  in_MLF       = GetParamStr(&cfgHash, SNAME":SOURCEMLF",       NULL);
  in_lbl_dir   = GetParamStr(&cfgHash, SNAME":SOURCETRANSCDIR", NULL);
  in_lbl_ext   = GetParamStr(&cfgHash, SNAME":SOURCETRANSCEXT", NULL);
  out_MLF      = GetParamStr(&cfgHash, SNAME":TARGETMLF",       NULL);
  out_lbl_dir  = GetParamStr(&cfgHash, SNAME":TARGETTRANSCDIR", NULL);
  out_lbl_ext  = GetParamStr(&cfgHash, SNAME":TARGETTRANSCEXT", "rec");
  state_pruning= GetParamFlt(&cfgHash, SNAME":PRUNING",         0.0);
//  stprn_step   = GetParamFlt(&cfgHash, SNAME":PRUNINGINC",      0.0);
//  stprn_limit  = GetParamFlt(&cfgHash, SNAME":PRUNINGMAX",      0.0);
  trace_flag   = GetParamInt(&cfgHash, SNAME":TRACE",           0);
  script =(char*)GetParamStr(&cfgHash, SNAME":SCRIPT",          NULL);
  mmf    =(char*)GetParamStr(&cfgHash, SNAME":SOURCEMMF",       NULL);

  cchrptr      = GetParamStr(&cfgHash, SNAME":LABELFORMATING",  "");
  while (*cchrptr) {
    switch (*cchrptr++) {
      case 'N': out_lbl_fmt.SCORE_NRM = 1; break;
      case 'S': out_lbl_fmt.SCORE_OFF = 1; break;
      case 'C': out_lbl_fmt.CENTRE_TM = 1; break;
      case 'T': out_lbl_fmt.TIMES_OFF = 1; break;
//      case 'W': out_lbl_fmt.WORDS_OFF = 1; break;
//      case 'M': out_lbl_fmt.MODEL_OFF = 1; break;
//      case 'F': out_lbl_fmt.FRAME_SCR = 1; break;
//      case 'X': out_lbl_fmt.STRIP_TRI = 1; break;
      default:
        Warning("Unknown label formating flag '%c' ignored (NCST)", *cchrptr);
    }
  }
  if (GetParamBool(&cfgHash, SNAME":PRINTCONFIG", false)) {
    PrintConfig(&cfgHash);
  }
  if (GetParamBool(&cfgHash, SNAME":PRINTVERSION", false)) {
    puts("Version: "VERSION"\n");
  }
  if (!GetParamBool(&cfgHash,SNAME":ACCEPTUNUSEDPARAM", false)) {
    CheckCommandLineParamUse(&cfgHash);
  }

  for (script=strtok(script, ","); script != NULL; script=strtok(NULL, ",")) {
    if ((sfp = my_fopen(script, "rt", gpScriptFilter)) == NULL) {
      Error("Cannot open script file %s", optarg);
    }
    while (fscanf(sfp, "%s", line) == 1) {
      last_file = AddFileElem(last_file, line);
      nfeature_files++;
    }
    my_fclose(sfp);
  }
  for (mmf=strtok(mmf, ","); mmf != NULL; mmf=strtok(NULL, ",")) {
    hset.ParseMmf(mmf, NULL);
  }
  if (hmm_list != NULL) hset.ReadHMMList(hmm_list, hmm_dir, hmm_ext);
  nonCDphHash = hset.MakeCIPhoneHash();

  if (dictionary != NULL) {
    ReadDictionary(dictionary, &dictHash, &phoneHash);
    notInDictAction  = WORD_NOT_IN_DIC_WARN;
    if (expOptions.mRespectPronunVar) {
      notInDictAction |= (int) PRON_NOT_IN_DIC_ERROR;
    }
  }
  if (dictHash.mNEntries == 0) expOptions.mNoWordExpansion = 1;

  if (!network_file) {
    ilfp = OpenInputMLF(in_MLF);
  } else {
    ilfp = fopen(network_file, "rt");
    if (ilfp  == NULL) Error("Cannot open network file: %s", network_file);

    Node *node = ReadSTKNetwork(ilfp, &dictHash, &phoneHash,
                                notInDictAction, in_lbl_fmt,
                                header.mSamplePeriod, network_file, NULL);
    NetworkExpansionsAndOptimizations(node, expOptions, in_net_fmt, &dictHash,
                                      &nonCDphHash, &phoneHash);
    //InitNetwork(&net, node, &hset, NULL);
    net.Init(node, &hset, NULL);
    fclose(ilfp);
    lrt = MakeLRTraceTable(&net, &nWords, &filler_end);
  }
  lfp = OpenOutputMLF(out_MLF);

  for (file_name = feature_files; file_name; file_name = file_name->mpNext) {

    if (trace_flag & 1) {
      TraceLog("Processing file %d/%d '%s'", ++fcnt,
                 nfeature_files,file_name->mpPhysical);
    }
    if (cmn_mask) process_mask(file_name->logical, cmn_mask, cmn_file);
    if (cvn_mask) process_mask(file_name->logical, cvn_mask, cvn_file);
    obsMx = ReadHTKFeatures(file_name->mpPhysical, swap_features,
                            startFrmExt, endFrmExt, targetKind,
                            derivOrder, derivWinLengths, &header,
                            cmn_path, cvn_path, cvg_file, &rhfbuff);

    if (hset.mInputVectorSize != static_cast<int>(header.mSampleSize / sizeof(float))) {
      Error("Vector size [%d] in '%s' is incompatible with HMM set [%d]",
            header.mSampleSize/sizeof(float), file_name->mpPhysical, hset.mInputVectorSize);
    }
    if (!network_file) {
      Node *node = NULL;
      strcpy(label_file, file_name->logical);
      ilfp = OpenInputLabelFile(label_file, in_lbl_dir,
                              in_lbl_ext ? in_lbl_ext : "net",
                              ilfp, in_MLF);

      node = ReadSTKNetwork(ilfp, &dictHash, &phoneHash, notInDictAction,
                            in_lbl_fmt, header.mSamplePeriod, label_file, in_MLF);

      NetworkExpansionsAndOptimizations(node, expOptions, in_net_fmt, &dictHash,
                                        &nonCDphHash, &phoneHash);

      //InitNetwork(&net, node, &hset, NULL);
      net.Init(node, &hset, NULL);
      CloseInputLabelFile(ilfp, in_MLF);
      lrt = MakeLRTraceTable(&net, &nWords, &filler_end);
    }
    net.mWPenalty     = word_penalty;
    net.mMPenalty     = model_penalty;
    net.mLmScale      = grammar_scale;
    net.mPronScale    = pronun_scale;
    net.mTranScale    = transp_scale;
    net.mOutpScale    = outprb_scale;
    net.mOcpScale     = occprb_scale;
    net.mAlignment     = alignment;
    net.mPruningThresh = state_pruning > 0.0 ? state_pruning : -LOG_0;

    strcpy(label_file, file_name->logical);
    lfp = OpenOutputLabelFile(label_file,out_lbl_dir,out_lbl_ext,lfp,out_MLF);

    for (i = 0; i < nWords; i++) {
      lrt[i].lastLR = -FLT_MAX;
      lrt[i].candidateLR = -FLT_MAX;
      lrt[i].candidateStartTime = 0;
      lrt[i].candidateEndTime = 0;
    }
    //ViterbiInit(&net);
    net.ViterbiInit();
    net.PassTokenInNetwork = fulleval ? &PassTokenSum : &PassTokenMax;
    net.PassTokenInModel   = fulleval ? &PassTokenSum : &PassTokenMax;

    KillToken(filler_end);
    
    for (j = 0; j < nWords; j++) 
      KillToken(lrt[j].mpWordEnd->mpExitToken);

    if (trace_flag & 2) 
    {
      printf("Node# %13s", "Filler");
      for (j = 0; j < nWords; j++) {
        printf(" %13s", lrt[j].mpWordEnd->mpPronun->mpWord->mpName);
      }
      puts("");
    }
    for (i = 0; i < header.mNSamples; i++) {
      //ViterbiStep(&net, obsMx + i * hset.mInputVectorSize);
      net.ViterbiStep(obsMx + i * hset.mInputVectorSize);

      if (trace_flag & 2) {
        printf("      %13e", filler_end->mLike);
        for (j = 0; j < nWords; j++) {
          printf(" %13e", lrt[j].mpWordEnd->mpExitToken->mLike);
        }
        puts("");
      }
      
      for (j = 0; j < nWords; j++) 
      {
        long long   wordStartTime;
        float       lhRatio;
        Token *     word_end = lrt[j].mpWordEnd->mpExitToken;

        if (!word_end->IsActive() || !filler_end->IsActive()) 
        {
          lrt[j].lastLR = -FLT_MAX;
          KillToken(word_end);
          continue;
        }
        
        lhRatio = word_end->mLike - filler_end->mLike;

        if (lrt[j].lastLR > lhRatio) { // LR is not growing, cannot be max.
          lrt[j].lastLR = lhRatio;
          KillToken(word_end);
          continue;
        }
        lrt[j].lastLR = lhRatio;
        wordStartTime = (long long) (word_end->mpWlr && word_end->mpWlr->mpNext
                        ? word_end->mpWlr->mpNext->mTime : 0);

        KillToken(word_end);
        if (lrt[j].candidateLR > lhRatio && lrt[j].candidateEndTime > wordStartTime) {
          continue;
        }
        if (lrt[j].candidateEndTime <= wordStartTime) {
          PutCandidateToLabels(&lrt[j], score_thresh, lfp, label_file,
                               out_MLF, out_lbl_fmt, header.mSamplePeriod);
        }
        lrt[j].candidateStartTime = wordStartTime;
        lrt[j].candidateEndTime = (long long) net.mTime;
        lrt[j].candidateLR = lhRatio;
      }
      KillToken(filler_end);
    }
    for (j = 0; j < nWords; j++) {
      PutCandidateToLabels(&lrt[j], score_thresh, lfp, label_file,
                           out_MLF, out_lbl_fmt, header.mSamplePeriod);
    }
    //ViterbiDone(&net, NULL);
    net.ViterbiDone(NULL);
    
    free(obsMx);
    CloseOutputLabelFile(lfp, out_MLF);

    if (trace_flag & 1) {
      TraceLog("[%d frames]", header.mNSamples - hset.mTotalDelay);
    }
    if (!network_file) {
      net.Release();
      free(lrt);
    }
  }
  
  if (network_file) {
    net.Release();
    free(lrt);
  }
  
  //ReleaseHMMSet(&hset);
  hset.Release();
  
  my_hdestroy_r(&phoneHash,   1);
  my_hdestroy_r(&nonCDphHash, 0);
  FreeDictionary(&dictHash);

  while (feature_files) {
    file_name     = feature_files;
    feature_files = feature_files->mpNext;
    free(file_name);
  }
  //my_hdestroy_r(&cfgHash, 0);

  return 0;
}

