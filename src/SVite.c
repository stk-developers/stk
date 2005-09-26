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

#define VERSION "0.4 "__TIME__" "__DATE__

#include "STKLib/viterbi.h"
#include "STKLib/hmms.h"
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

void usage(char *progname)
{
  char *tchrptr;
  if((tchrptr = strrchr(progname, '\\')) != NULL) progname = tchrptr+1;
  if((tchrptr = strrchr(progname, '/')) != NULL) progname = tchrptr+1;
  fprintf(stderr,
"\nUSAGE: %s [options] DataFiles...\n\n"
" Option                                                     Default\n\n"
" -a         Align from label files                          On\n"
////"  -b s    def s as utterance boundary word                none\n"
" -d s       Dir to find hmm definitions                     Current\n"
" -f         Output full state alignment                     Off\n"
" -i s       Output transcriptions to MLF s                  Off\n"
" -l s       Dir to store transcription files                Feature file dir\n"
" -m         Output model alignment                          Off\n"
////"  -n i [N] N-best recognition (using i tokens)             off\n"
" -o s       Output label formating NCSTWMXF                 None\n"
" -p f       Inter word trans penalty (log)                  0.0\n"
" -q s       Output network formating JMRVWXalpstv           tvl\n"
////"  -q s    output lattice formating ABtvaldmn              tvaldmn\n"
" -r f       Pronunciation prob scale factor                 1.0\n"
" -s f       Grammar scale factor                            1.0\n"
" -t f [i l] Set pruning to f [inc limit]                    Off\n"
////"  -u i    set pruning max active                          0\n"
//"  -v f    Set word end pruning threshold                  0.0\n"
" -w [f]     Recognise from network                          Off\n"
" -x s       Extension for hmm files                         None\n"
" -y s       Output transcription file extension             rec\n"
////"  -z s    generate lattices with extension s              off\n"
" -A         Print command line arguments                    Off\n"
" -C cf      Set config file to cf                           Default\n"
" -D         Display configuration variables                 Off\n"
" -G fmt     Set source trascription format to fmt           As config\n"
" -H mmf     Load HMM macro file mmf\n"
" -I mlf     Load master transcription file mlf\n"
" -L dir     Set input transcription dir                     Current\n"
" -P fmt     Set target transcription format to fmt          As config\n"
" -S f       Set script file to f                            None\n"
" -T N       Set trace flags to N                            0\n"
" -V         Print version information                       Off\n"
" -X ext     Set input transcription file ext                lab\n"
"\n"
" %s is Copyright (C) 2004-2005 Lukas Burget et al. and\n"
" licensed under the GNU General Public License, version 2.\n"
" Bug reports, feedback, etc, to: burget@fit.vutbr.cz\n"
"\n", progname, progname);
  exit(-1);
}

#define SNAME "SVITE"
char *optionStr =
" -a"
" -d r   SOURCEMODELDIR"
" -f n   STATEALIGNMENT=TRUE"
" -i r   TARGETMLF"
" -l r   TARGETTRANSCDIR"
" -m n   MODELALIGNMENT=TRUE"
" -o r   LABELFORMATING"
" -p r   WORDPENALTY"
" -q r   NETFORMATING"
" -r r   PRONUNSCALE"
" -s r   LMSCALE"
" -t ror PRUNING PRUNINGINC PRUNINGMAX"
" -w o   RECOGNET"
" -x r   SOURCEMODELEXT"
" -y r   TARGETTRANSCEXT"
" -D n   PRINTCONFIG=TRUE"
" -G r   SOURCETRANSCFMT"
" -H l   SOURCEMMF"
" -I r   SOURCEMLF"
" -L r   SOURCETRANSCDIR"
" -P r   TARGETTRANSCFMT"
" -S l   SCRIPT"
" -T r   TRACE"
" -V n   PRINTVERSION=TRUE"
" -X r   SOURCETRANSCEXT";

int main(int argc, char *argv[]) {
  HTK_Header header;
  HMMSet hset;
  Network net;
  FILE *sfp, *lfp = NULL, *ilfp = NULL;
  FLOAT  *obsMx;
  FLOAT like;
  int i, fcnt = 0;
  Label *labels;
  char line[1024];
  char label_file[1024];
  const char *cchrptr;
  struct my_hsearch_data nonCDphHash, phoneHash, dictHash, cfgHash;

  FileListElem *feature_files = NULL;
  int nfeature_files = 0;
  FileListElem *file_name = NULL;
  FileListElem **last_file = &feature_files;

  Alignment alignment = WORD_ALIGNMENT;

  double word_penalty;
  double model_penalty;
  double grammar_scale;
  double transp_scale;
  double outprb_scale;
  double pronun_scale;
  double occprb_scale;
  double state_pruning;
  double stprn_step;
  double stprn_limit;
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
  const char *label_filter;
  const char *net_filter;
  const char *label_ofilter;
  const char *net_ofilter;
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
  BOOL baum_welch;
  BOOL swap_features;
  BOOL htk_compat;
  enum TranscriptionFormat {TF_HTK, TF_STK} in_transc_fmt, out_transc_fmt;
  enum NotInDictAction notInDictAction = 0;
  RHFBuffer rhfbuff                    = {0};
  ExpansionOptions expOptions          = {0};
  STKNetworkOutputFormat in_net_fmt    = {0};
  STKNetworkOutputFormat out_net_fmt   = {0};
  LabelFormat out_lbl_fmt              = {0};
  LabelFormat in_lbl_fmt               = {0};
  in_lbl_fmt.TIMES_OFF = 1;

  if(argc == 1) usage(argv[0]);

  InitHMMSet(&hset, 0);

  if(!my_hcreate_r(100,  &dictHash)
  || !my_hcreate_r(100,  &phoneHash)
  || !my_hcreate_r(100,  &cfgHash)) {
    Error("Insufficient memory");
  }
  i = ParseOptions(argc, argv, optionStr, SNAME, &cfgHash);
  htk_compat = GetParamBool(&cfgHash, SNAME":HTKCOMPAT", FALSE);
  if(htk_compat) {
    if(argc == i) Error("Dictionary file name expected");
    InsertConfigParam(&cfgHash, SNAME":SOURCEDICT", argv[i++], '-');
    if(argc == i) Error("HMM list file name expected");
    InsertConfigParam(&cfgHash, SNAME":SOURCEHMMLIST",    argv[i++], '-');
  }
  for(; i < argc; i++) {
    last_file = AddFileElem(last_file, argv[i]);
    nfeature_files++;
  }
  targetKind =   GetDerivParams(&cfgHash, &derivOrder, &derivWinLengths,
                                &startFrmExt, &endFrmExt,
                                &cmn_path, &cmn_file, &cmn_mask,
                                &cvn_path, &cvn_file, &cvn_mask, &cvg_file,
                                SNAME":", 0);
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
  baum_welch   = GetParamBool(&cfgHash,SNAME":EVALUATION",      FALSE);
  swap_features=!GetParamBool(&cfgHash,SNAME":NATURALREADORDER",isBigEndian());
  filter_wldcrd= GetParamStr(&cfgHash, SNAME":HFILTERWILDCARD", "$");
  script_filter= GetParamStr(&cfgHash, SNAME":HSCRIPTFILTER",   NULL);
  parm_filter  = GetParamStr(&cfgHash, SNAME":HPARMFILTER",     NULL);
  hlist_filter = GetParamStr(&cfgHash, SNAME":HMMLISTFILTER",   NULL);
  MMF_filter   = GetParamStr(&cfgHash, SNAME":HMMDEFFILTER",    NULL);
  label_filter = GetParamStr(&cfgHash, SNAME":HLABELFILTER",    NULL);
  net_filter   = GetParamStr(&cfgHash, SNAME":HNETFILTER",      NULL);
  dict_filter  = GetParamStr(&cfgHash, SNAME":HDICTFILTER",     NULL);
  label_ofilter= GetParamStr(&cfgHash, SNAME":HLABELOFILTER",   NULL);
  net_ofilter  = GetParamStr(&cfgHash, SNAME":HNETOFILTER",     NULL);
  dictionary   = GetParamStr(&cfgHash, SNAME":SOURCEDICT",      NULL);
  hmm_list     = GetParamStr(&cfgHash, SNAME":SOURCEHMMLIST",   NULL);
  grammar_scale= GetParamFlt(&cfgHash, SNAME":LMSCALE",         1.0);
  outprb_scale = GetParamFlt(&cfgHash, SNAME":OUTPSCALE",       1.0);
  transp_scale = GetParamFlt(&cfgHash, SNAME":TRANSPSCALE",     1.0);
  pronun_scale = GetParamFlt(&cfgHash, SNAME":PRONUNSCALE",     1.0);
  occprb_scale = GetParamFlt(&cfgHash, SNAME":OCCUPPSCALE",     1.0);
  word_penalty = GetParamFlt(&cfgHash, SNAME":WORDPENALTY",     0.0);
  model_penalty= GetParamFlt(&cfgHash, SNAME":MODELPENALTY",    0.0);
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
  stprn_step   = GetParamFlt(&cfgHash, SNAME":PRUNINGINC",      0.0);
  stprn_limit  = GetParamFlt(&cfgHash, SNAME":PRUNINGMAX",      0.0);
  trace_flag   = GetParamInt(&cfgHash, SNAME":TRACE",           0);
  script =(char*)GetParamStr(&cfgHash, SNAME":SCRIPT",          NULL);
  mmf    =(char*)GetParamStr(&cfgHash, SNAME":SOURCEMMF",       NULL);

  cchrptr      = GetParamStr(&cfgHash, SNAME":LABELFORMATING",  "");
  while(*cchrptr) {
    switch(*cchrptr++) {
      case 'N': out_lbl_fmt.SCORE_NRM = 1; break;
      case 'S': out_lbl_fmt.SCORE_OFF = 1; break;
      case 'C': out_lbl_fmt.CENTRE_TM = 1; break;
      case 'T': out_lbl_fmt.TIMES_OFF = 1; break;
      case 'W': out_lbl_fmt.WORDS_OFF = 1; break;
      case 'M': out_lbl_fmt.MODEL_OFF = 1; break;
      case 'F': out_lbl_fmt.FRAME_SCR = 1; break;
//      case 'X': out_lbl_fmt.STRIP_TRI = 1; break;
      default:
        Warning("Unknown label formating flag '%c' ignored (NCSTWMF)", *cchrptr);
    }
  }
  cchrptr      = GetParamStr(&cfgHash, SNAME":NETFORMATING",  "");
  if(*cchrptr) {
    out_net_fmt.no_LM_likes    = 1;
    out_net_fmt.no_times       = 1;
    out_net_fmt.no_pronun_vars = 1;
    out_net_fmt.no_acc_likes   = 1;
  }
  while(*cchrptr) {
    switch(*cchrptr++) {
      case 'R': out_net_fmt.base62_labels  = 1; // reticent
                out_net_fmt.lin_node_seqs  = 1;
                out_net_fmt.no_defaults    = 1; break;
      case 'V': out_net_fmt.arc_defs_with_J= 1;
                out_net_fmt.all_field_names= 1; break;
      case 'J': out_net_fmt.arc_defs_to_end= 1; break;
      case 'W': out_net_fmt.no_word_nodes  = 1; break;
      case 'M': out_net_fmt.no_model_nodes = 1; break;
      case 'X': out_net_fmt.strip_triphones= 1; break;
      case 't': out_net_fmt.no_times       = 0; break;
      case 's': out_net_fmt.start_times    = 1; break;
      case 'v': out_net_fmt.no_pronun_vars = 0; break;
      case 'a': out_net_fmt.no_acc_likes   = 0; break;
      case 'l': out_net_fmt.no_LM_likes    = 0; break;
      case 'p': out_net_fmt.aprox_accuracy = 1; break;
      default:
        Warning("Unknown net formating flag '%c' ignored (JMRVWXalpstv)", *cchrptr);
    }
  }
  in_transc_fmt= GetParamEnum(&cfgHash,SNAME":SOURCETRANSCFMT",
                              !network_file && htk_compat ? TF_HTK : TF_STK,
                              "HTK", TF_HTK, "STK", TF_STK, NULL);

  out_transc_fmt=GetParamEnum(&cfgHash,SNAME":TARGETTRANSCFMT",
                              htk_compat ? TF_HTK : TF_STK,
                              "HTK", TF_HTK, "STK", TF_STK, NULL);

  if(GetParamBool(&cfgHash, SNAME":STATEALIGNMENT", FALSE)) {
    alignment |= STATE_ALIGNMENT;
  }
  if(GetParamBool(&cfgHash, SNAME":MODELALIGNMENT", FALSE)) {
    alignment |= MODEL_ALIGNMENT;
  }
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
  for(mmf=strtok(mmf, ","); mmf != NULL; mmf=strtok(NULL, ",")) {
    ReadHMMSet(mmf, &hset, NULL);
  }
  if(hmm_list != NULL) ReadHMMList(&hset, hmm_list, hmm_dir, hmm_ext);
  nonCDphHash = MakeCIPhoneHash(&hset);

  if(dictionary != NULL) {
    ReadDictionary(dictionary, &dictHash, &phoneHash);
    notInDictAction  = WORD_NOT_IN_DIC_WARN;
    if(expOptions.respect_pronun_var) {
      notInDictAction |= PRON_NOT_IN_DIC_ERROR;
    }
  }
  if(dictHash.nentries == 0) expOptions.no_word_expansion = 1;

  transc_filter = transc_filter != NULL   ? transc_filter :
                  in_transc_fmt == TF_STK ? net_filter    : label_filter;

  transc_ofilter = transc_ofilter != NULL ? transc_ofilter :
                   out_transc_fmt == TF_STK  ? net_ofilter : label_ofilter;

  if(dictionary == NULL && in_transc_fmt == TF_HTK) {
    // Word alignment is inpossible in this case.
    if(alignment & WORD_ALIGNMENT) {
      alignment &= ~WORD_ALIGNMENT;
      alignment |= MODEL_ALIGNMENT;
    }
    out_lbl_fmt.WORDS_OFF = 1;
  }
  if(!network_file) {
    ilfp = OpenInputMLF(in_MLF);
  } else {
    ilfp = fopen(network_file, "rt");
    if(ilfp  == NULL) Error("Cannot open network file: %s", network_file);

    Node *node = ReadSTKNetwork(ilfp, &dictHash, &phoneHash,
                                notInDictAction, in_lbl_fmt,
                                header.sampPeriod, network_file, NULL);
    NetworkExpansionsAndOptimizations(node, expOptions, in_net_fmt, &dictHash,
                                      &nonCDphHash, &phoneHash);
    InitNetwork(&net, node, &hset, NULL);
    fclose(ilfp);
  }
  lfp = OpenOutputMLF(out_MLF);

  for(file_name = feature_files; file_name; file_name = file_name->next) {

    if(trace_flag & 1) {
      TraceLog("Processing file %d/%d '%s'", ++fcnt,
                 nfeature_files,file_name->physical);
    }
    if(cmn_mask) process_mask(file_name->logical, cmn_mask, cmn_file);
    if(cvn_mask) process_mask(file_name->logical, cvn_mask, cvn_file);
    obsMx = ReadHTKFeatures(file_name->physical, swap_features,
                            startFrmExt, endFrmExt, targetKind,
                            derivOrder, derivWinLengths, &header,
                            cmn_path, cvn_path, cvg_file, &rhfbuff);

/*  lfp = fopen("xxx.fea", "w");
    header.sampKind = 9;
    WriteHTKHeader(lfp, header, 1);
    WriteHTKFeature(lfp, obsMx, header.nSamples * header.sampSize / sizeof(float), 1);
    fclose(lfp);
    exit(0); */

    if(hset.in_vec_size != header.sampSize / sizeof(float)) {
      Error("Vector size [%d] in '%s' is incompatible with HMM set [%d]",
            header.sampSize/sizeof(float), file_name->physical, hset.in_vec_size);
    }
    if(!network_file) {
      Node *node = NULL;
      strcpy(label_file, file_name->logical);
      ilfp = OpenInputLabelFile(label_file, in_lbl_dir,
                              in_lbl_ext ? in_lbl_ext :
                              in_transc_fmt == TF_STK ? "net" : "lab",
                              ilfp, in_MLF);

      if(in_transc_fmt == TF_HTK) {
        labels = ReadLabels(ilfp, dictionary ? &dictHash : &phoneHash,
                                  dictionary ? UL_ERROR : UL_INSERT, in_lbl_fmt,
                                  header.sampPeriod, label_file, in_MLF, NULL);
        node = MakeNetworkFromLabels(labels, dictionary ? NT_Word : NT_Phone);
        ReleaseLabels(labels);
      } else if(in_transc_fmt == TF_STK) {
        node = ReadSTKNetwork(ilfp, &dictHash, &phoneHash, notInDictAction,
                              in_lbl_fmt, header.sampPeriod, label_file, in_MLF);
      } else Error("Too bad. What did you do ?!?");

      NetworkExpansionsAndOptimizations(node, expOptions, in_net_fmt, &dictHash,
                                        &nonCDphHash, &phoneHash);
      InitNetwork(&net, node, &hset, NULL);


      CloseInputLabelFile(ilfp, in_MLF);
    }

    net.wPenalty     = word_penalty;
    net.mPenalty     = model_penalty;
    net.lmScale      = grammar_scale;
    net.pronScale    = pronun_scale;
    net.tranScale    = transp_scale;
    net.outpScale    = outprb_scale;
    net.ocpScale     = occprb_scale;
    net.alignment    = alignment;
    net.pruningThresh = state_pruning > 0.0 ? state_pruning : -LOG_0;

    if(alignment & STATE_ALIGNMENT && out_lbl_fmt.MODEL_OFF) net.alignment &= ~MODEL_ALIGNMENT;
    if(alignment & MODEL_ALIGNMENT && out_lbl_fmt.WORDS_OFF) net.alignment &= ~WORD_ALIGNMENT;
    if(alignment & STATE_ALIGNMENT && out_lbl_fmt.FRAME_SCR) net.alignment |=  FRAME_ALIGNMENT;

    for(;;) {
      ViterbiInit(&net);
      net.PassTokenInNetwork = baum_welch ? &PassTokenSum : &PassTokenMax;
      net.PassTokenInModel   = baum_welch ? &PassTokenSum : &PassTokenMax;

      for(i = 0; i < header.nSamples; i++) {
        ViterbiStep(&net, obsMx + i * hset.in_vec_size);
      }
      like = ViterbiDone(&net, &labels);

      if(labels) {
        break;
      }
      if(net.pruningThresh <= LOG_MIN ||
         stprn_step <= 0.0 ||
         (net.pruningThresh += stprn_step) > stprn_limit ) {
        Warning("No tokens survived");
        break;
      }
      Warning("No tokens survived, trying pruning threshold: %.2f", net.pruningThresh);
    }
    if(trace_flag & 1 && labels) {
      Label *label;
      int nFrames = header.nSamples - hset.totalDelay;
      for(label = labels; label->nextLevel != NULL; label = label->nextLevel);
      for(; label != NULL; label = label->next) {
        fprintf(stdout, "%s ", label->name);
      }
      TraceLog(" ==  [%d frames] %f", nFrames, like / nFrames);
    }
    free(obsMx);
    strcpy(label_file, file_name->logical);
    lfp = OpenOutputLabelFile(label_file, out_lbl_dir, out_lbl_ext, lfp, out_MLF);

    if(out_transc_fmt == TF_HTK) {
      WriteLabels(lfp, labels, out_lbl_fmt, header.sampPeriod, label_file, out_MLF);
    } else {
      Node *node = MakeNetworkFromLabels(labels,
                                         alignment & (MODEL_ALIGNMENT|STATE_ALIGNMENT)
                                         ? NT_Model : NT_Word);
      WriteSTKNetwork(lfp, node, out_net_fmt, header.sampPeriod, label_file, out_MLF);
      FreeNetwork(node);
    }
    CloseOutputLabelFile(lfp, out_MLF);
    ReleaseLabels(labels);

    if(!network_file) {
      ReleaseNetwork(&net);
    }
  }
  if(network_file) {
    ReleaseNetwork(&net);
  }
  ReleaseHMMSet(&hset);
//  my_hdestroy_r(&labelHash,   0);
  my_hdestroy_r(&phoneHash,   1);
  my_hdestroy_r(&nonCDphHash, 0);
  FreeDictionary(&dictHash);
  for(i = 0; i < cfgHash.nentries; i++) free(cfgHash.entry[i]->data);
  my_hdestroy_r(&cfgHash, 1);

  while(feature_files) {
    file_name     = feature_files;
    feature_files = feature_files->next;
    free(file_name);
  }
  //my_hdestroy_r(&cfgHash, 0);

  return 0;
}

//HVite -T 05 -H models -w wdnet dict words4 MAL_4379315A.fea > htk.log
