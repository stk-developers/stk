/***************************************************************************
 *   copyright           : (C) 2004-2005 by Lukas Burget,UPGM,FIT,VUT,Brno *
 *   email               : burget@fit.vutbr.cz                             *
 ***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#define VERSION "0.2 "__TIME__" "__DATE__
#include "STKLib/net.h"
#include "STKLib/labels.h"
#include "STKLib/common.h"

#ifndef WIN32
#include <unistd.h>
#else
#include "getopt.h"
#endif
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <malloc.h>
#include <assert.h>

#define SIGNIFICANT_PROB_DIFFERENCE (0.01)
using namespace STK;

void usage(char *progname)
{
  char *tchrptr;
  if ((tchrptr = strrchr(progname, '\\')) != NULL) progname = tchrptr+1;
  if ((tchrptr = strrchr(progname, '/')) != NULL) progname = tchrptr+1;
  fprintf(stderr,
"\nUSAGE: %s [options] NetworkFiles...\n\n"
" Option                                                     Default\n\n"
//" -e         Expand monophones to triphones                  Off\n"
" -i s       Output network to MNF s                         Off\n"
" -l s       Dir to store network files                      Label file dir\n"
//"*-m         Expansion respects pronuncioation variants      Off\n"
//" -n f       Output list of all CD model to f                Off\n"
//" -p s       Keep phonem s as context independent            \n"
" -q s       Output network formating JMRVWXalpstv           tvl\n"
//"*-s         Strict time optimization                        Off\n"
//" -t s       Treat s as tee (sp) model                       \n"
//"*-u         Turn off network optimization                   Optimize\n"
//" -w         Remove expanded word nodes                      Off\n"
" -y s       Output network file extension                   net\n"
" -A         Print command line arguments                    Off\n"
" -C cf      Set config file to cf                           Default\n"
" -D         Display configuration variables                 off\n"
" -G fmt     Set source trascription format to fmt           As config\n"
" -I mlf     Load master label file mlf                      \n"
" -L dir     Set input label or network dir                  Current\n"
" -P         Set target network format to fmt                STK\n"
" -S f       Set script file to f                            None\n"
" -T N       Set trace flags to N                            0\n"
" -V         Print version information                       Off\n"
" -X ext     Set input label file ext                        lab\n"
"\n"
" %s is Copyright (C) 2004-2005 Lukas Burget et al. and\n"
" licensed under the GNU General Public License, version 2.\n"
" Bug reports, feedback, etc, to: burget@fit.vutbr.cz\n"
"\n", progname, progname);
  exit(-1);
}

#define SNAME "SEXPAND"
char *optionStr =
" -i r   TARGETMLF"
" -l r   TARGETTRANSCDIR"
//" -m n   RESPECTPRONVARS=TRUE"
//" -n r   TARGETHMMLIST"
//" -p r   CIMODELS"
" -q r   NETFORMATING"
//" -s n   EXACTTIMEMERGE=TRUE"
//" -t r   TEEMODELS"
//" -u n   MINIMIZENET=FALSE"
//" -w n   REMEXPWRDNODES=TRUE"
" -y r   TARGETTRANSCEXT"
" -D n   PRINTCONFIG=TRUE"
" -G r   SOURCETRANSCFMT"
" -I r   SOURCEMLF"
" -L r   SOURCETRANSCDIR"
" -P r   TARGETTRANSCFMT"
" -S l   SCRIPT"
" -T r   TRACE"
" -V n   PRINTVERSION=TRUE"
" -X r   SOURCETRANSCEXT";


int main(int argc, char *argv[]) {
  FILE *sfp = NULL;
  ENTRY e, *ep;
  int  i;
  Label *labels;
  char line[1024];
  char label_file[1024];
  struct my_hsearch_data phoneHash, dictHash = {0};
  struct my_hsearch_data triphHash, nonCDphHash, cfgHash ;

  FileListElem *feature_files=NULL;
  int nfeature_files = 0;
  FileListElem *file_name   = NULL;
  FileListElem **last_file  = &feature_files;

  FILE *out_MLF_fp          = NULL;
  FILE *in_MLF_fp           = NULL;
  const char *cchrptr;
  const char *in_MLF_fn;
  const char *in_lbl_dir;
  const char *in_lbl_ext;
  const char *out_MLF_fn;
  const char *out_lbl_dir;
  const char *out_lbl_ext;
  const char *cd_list_file;
  const char *dictionary;
  const char *label_filter;
  const char *net_filter;
  char *script;
  char *ci_phn;
  char *tee_phn;
  int trace_flag       =  0;
  int fcnt             = 0;
  enum TranscriptionFormat {
    TF_MLF, TF_MNF, TF_HTK, TF_STK, TF_ERR,
    TF_NOF, TF_MOF, //obsolote formats
  } in_transc_fmt = TF_STK, out_transc_fmt = TF_STK;
  int notInDictAction = 0;
  ExpansionOptions expOptions = {0};
  STKNetworkOutputFormat out_net_fmt = {0};
  LabelFormat out_lbl_fmt = {0};
  LabelFormat in_lbl_fmt = {0};

  if (argc == 1) usage(argv[0]);

  if (!my_hcreate_r(10,   &nonCDphHash)
  || !my_hcreate_r(100,  &dictHash)
  || !my_hcreate_r(100,  &phoneHash)
  || !my_hcreate_r(1000, &triphHash)
  || !my_hcreate_r(100,  &cfgHash)) {

  }
  i = ParseOptions(argc, argv, optionStr, SNAME, &cfgHash);
//  htk_compat = GetParamBool(&cfgHash, SNAME":HTKCOMPAT", FALSE);
  for (; i < argc; i++) {
    last_file = AddFileElem(last_file, argv[i]);
    nfeature_files++;
  }
  expOptions.CD_phone_expansion =
                 GetParamBool(&cfgHash,SNAME":ALLOWXWRDEXP",    FALSE);
  expOptions.respect_pronun_var
               = GetParamBool(&cfgHash,SNAME":RESPECTPRONVARS", FALSE);
  expOptions.strict_timing
               = GetParamBool(&cfgHash,SNAME":EXACTTIMEMERGE",  FALSE);
  expOptions.no_optimization
               =!GetParamBool(&cfgHash,SNAME":MINIMIZENET",     TRUE);
  expOptions.remove_words_nodes
               = GetParamBool(&cfgHash,SNAME":REMEXPWRDNODES",  FALSE);
  in_lbl_fmt.left_extent  = -100 * (long long) (0.5 + 1e5 *
                 GetParamFlt(&cfgHash, SNAME":STARTTIMESHIFT",  0.0));
  in_lbl_fmt.right_extent =  100 * (long long) (0.5 + 1e5 *
                 GetParamFlt(&cfgHash, SNAME":ENDTIMESHIFT",    0.0));
  label_filter = GetParamStr(&cfgHash, SNAME":HLABELFILTER",    NULL);
  gpFilterWldcrd= GetParamStr(&cfgHash, SNAME":HFILTERWILDCARD", "$");
  script_filter= GetParamStr(&cfgHash, SNAME":HSCRIPTFILTER",   NULL);
  net_filter   = GetParamStr(&cfgHash, SNAME":HNETFILTER",      NULL);
  dict_filter  = GetParamStr(&cfgHash, SNAME":HDICTFILTER",     NULL);
//label_ofilter= GetParamStr(&cfgHash, SNAME":HLABELOFILTER",   NULL);
  transc_ofilter=GetParamStr(&cfgHash, SNAME":HNETOFILTER",     NULL);
  dictionary   = GetParamStr(&cfgHash, SNAME":SOURCEDICT",      NULL);
  in_MLF_fn    = GetParamStr(&cfgHash, SNAME":SOURCEMLF",       NULL);
  in_lbl_dir   = GetParamStr(&cfgHash, SNAME":SOURCETRANSCDIR", NULL);
  in_lbl_ext   = GetParamStr(&cfgHash, SNAME":SOURCETRANSCEXT", NULL);
  out_MLF_fn   = GetParamStr(&cfgHash, SNAME":TARGETMLF",       NULL);
  out_lbl_dir  = GetParamStr(&cfgHash, SNAME":TARGETTRANSCDIR", NULL);
  out_lbl_ext  = GetParamStr(&cfgHash, SNAME":TARGETTRANSCEXT", "net");
  trace_flag   = GetParamInt(&cfgHash, SNAME":TRACE",           0);
  cd_list_file = GetParamStr(&cfgHash, SNAME":TARGETHMMLIST",   NULL);
  ci_phn =(char*)GetParamStr(&cfgHash, SNAME":CIMODEL",        NULL);
  tee_phn=(char*)GetParamStr(&cfgHash, SNAME":TEEMODEL",       NULL);
  script =(char*)GetParamStr(&cfgHash, SNAME":SCRIPT",          NULL);

  cchrptr      = GetParamStr(&cfgHash, SNAME":NETFORMATING",  "");
  if (*cchrptr) {
    out_net_fmt.no_LM_likes    = 1;
    out_net_fmt.no_times       = 1;
    out_net_fmt.no_pronun_vars = 1;
    out_net_fmt.no_acc_likes   = 1;
  }
  while (*cchrptr) {
    switch (*cchrptr++) {
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
  in_transc_fmt= (TranscriptionFormat) GetParamEnum(&cfgHash,SNAME":SOURCETRANSCFMT", TF_STK,
                              "HTK", TF_HTK, "STK", TF_STK, "net", TF_NOF, NULL);

  out_transc_fmt= (TranscriptionFormat) GetParamEnum(&cfgHash,SNAME":TARGETTRANSCFMT", TF_STK,
                              "STK", TF_STK, "net", TF_NOF, NULL);

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
      Error("Cannot open script file %s", optarg);
    }
    while (fscanf(sfp, "%s", line) == 1) {
      last_file = AddFileElem(last_file, line);
      nfeature_files++;
    }
    my_fclose(sfp);
  }
  for ( ci_phn=strtok(ci_phn, ",");  ci_phn != NULL; ci_phn=strtok(NULL, ",")) {
    e.key = ci_phn;
    my_hsearch_r(e, FIND, &ep, &nonCDphHash);
    if (ep != NULL) continue;
    if ((e.key = strdup(ci_phn)) == NULL) Error("Insufficient memory");
    e.data = (void *) 0;
    my_hsearch_r(e, ENTER, &ep, &nonCDphHash);
  }
  for (tee_phn=strtok(tee_phn, ",");tee_phn != NULL;tee_phn=strtok(NULL, ",")) {
    e.key = tee_phn;
    my_hsearch_r(e, FIND, &ep, &nonCDphHash);
    if (ep != NULL) {
      ep->data = (void *) 1;
      continue;
    }
    if ((e.key = strdup(tee_phn)) == NULL) Error("Insufficient memory");
    e.data = (void *) 1;
    my_hsearch_r(e, ENTER, &ep, &nonCDphHash);
  }
  if (dictionary != NULL) {
    ReadDictionary(dictionary, &dictHash, &phoneHash);
    notInDictAction  = WORD_NOT_IN_DIC_WARN;
    if (expOptions.respect_pronun_var) {
      notInDictAction |= (int) PRON_NOT_IN_DIC_ERROR;
    }
  }
  if (dictHash.nentries == 0) expOptions.no_word_expansion = 1;

  transc_filter = transc_filter != NULL ? transc_filter :
                  in_transc_fmt == TF_STK      ? net_filter    :
                                          label_filter;

  in_lbl_fmt.TIMES_OFF = out_net_fmt.no_times;

  in_MLF_fp  = OpenInputMLF(in_MLF_fn);
  out_MLF_fp = OpenOutputMLF(out_MLF_fn);

  for (file_name=feature_files; file_name != NULL; file_name=file_name->mpNext) {
    if (trace_flag & 1) {
      TraceLog("Processing file %d/%d '%s'", ++fcnt, nfeature_files,file_name->logical);
    }

    strcpy(label_file, file_name->logical);
    in_MLF_fp = OpenInputLabelFile(label_file, in_lbl_dir,
                                   in_lbl_ext, in_MLF_fp, in_MLF_fn);

    if (in_MLF_fn == NULL && IsMLF(in_MLF_fp)) {
      in_transc_fmt = in_transc_fmt == TF_HTK ? TF_MLF :
      in_transc_fmt == TF_STK ? TF_MNF :
      in_transc_fmt == TF_NOF ? TF_MOF : TF_ERR;

      in_MLF_fn = label_file;
      assert(in_transc_fmt != TF_ERR);
    }
    for (;;) { //in cases of MLF or MNF, we must process all records
      Node *node = NULL;

      if (in_transc_fmt == TF_MLF
      || in_transc_fmt == TF_MNF
      || in_transc_fmt == TF_MOF) {
        label_file[0]='\0'; // Ref. MLF is read sequentially record by record
        if (!OpenInputLabelFile(label_file, NULL, NULL, in_MLF_fp, in_MLF_fn)) {
          break; // Whole MLF or MNF processed
        }
        if (trace_flag & 1) TraceLog("Processing file %d '%s'", ++fcnt, label_file);
      }
      if (in_transc_fmt == TF_HTK || in_transc_fmt == TF_MLF) {
        labels = ReadLabels(in_MLF_fp, &dictHash, UL_ERROR, in_lbl_fmt,
                            /*sampleRate*/ 1, label_file, in_MLF_fn, NULL);
        node = MakeNetworkFromLabels(labels, NT);
        ReleaseLabels(labels);
      } else if (in_transc_fmt == TF_STK || in_transc_fmt == TF_MNF) {
        node = ReadSTKNetwork(in_MLF_fp, &dictHash, &phoneHash, notInDictAction,
                              in_lbl_fmt, /*sampleRate*/ 1, label_file, in_MLF_fn);
      } else if (in_transc_fmt == TF_NOF || in_transc_fmt == TF_MOF) {
        node = ReadSTKNetworkInOldFormat(
                 in_MLF_fp, &dictHash, &phoneHash, in_lbl_fmt,
                 /*sampleRate*/ 1, label_file, in_MLF_fn);
      }
      CloseInputLabelFile(in_MLF_fp, in_MLF_fn);
      NetworkExpansionsAndOptimizations(node, expOptions, out_net_fmt, &dictHash,
                                        &nonCDphHash, &triphHash);

      if (out_net_fmt.aprox_accuracy)
        ComputeAproximatePhoneAccuracy(node, 0);

      out_MLF_fp = OpenOutputLabelFile(label_file, out_lbl_dir, out_lbl_ext,
                                      out_MLF_fp, out_MLF_fn);

      if (out_transc_fmt == TF_NOF) {
        WriteSTKNetworkInOldFormat(out_MLF_fp, node, out_lbl_fmt, 1,
                                 label_file, out_MLF_fn);
      } else {
        WriteSTKNetwork(out_MLF_fp, node, out_net_fmt, 1, label_file, out_MLF_fn);
      }
      CloseOutputLabelFile(out_MLF_fp, out_MLF_fn);
      FreeNetwork(node);

      if (in_transc_fmt == TF_HTK
      || in_transc_fmt == TF_STK
      || in_transc_fmt == TF_NOF) {
        break; // We are dealing with single record in these cases
      }
    }
    if (in_transc_fmt == TF_MLF || in_transc_fmt == TF_MNF) {
      CloseInputMLF(in_MLF_fp);
    }
  }
  if (out_MLF_fn != NULL) fclose(out_MLF_fp);

  if (in_MLF_fn != NULL && (in_transc_fmt == TF_HTK || in_transc_fmt == TF_STK)) {
    CloseInputMLF(in_MLF_fp);
  }
  if (cd_list_file != NULL) {
    FILE *fp = fopen(cd_list_file, "wt");
    if (fp == NULL) Error("Cannot open output file: '%s'", cd_list_file);

    for (i = 0; i < triphHash.nentries; i++) {
      if (fprintf(fp, "%s\n", (char *) triphHash.entry[i]->key) < 0) {
        Error("Cannot write to file: '%s'", cd_list_file);
      }
    }
    fclose(fp);
  }
  my_hdestroy_r(&phoneHash, 1);
  my_hdestroy_r(&nonCDphHash, 0);
  FreeDictionary(&dictHash);
  
  for (i = 0; i < cfgHash.nentries; i++) 
    free(cfgHash.entry[i]->data);
  
  my_hdestroy_r(&cfgHash, 1);
  return 0;
}
