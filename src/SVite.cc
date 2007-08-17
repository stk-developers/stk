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


#define MODULE_VERSION "2.0.4 "__TIME__" "__DATE__" "SVN_ID


////////////////////////////////////////////////////////////////////////////////
// STK HEADERS
////////////////////////////////////////////////////////////////////////////////
#include "STKLib/Features.h"

#include "STKLib/Decoder.h"
#include "STKLib/Models.h"
#include "STKLib/fileio.h"
#include "STKLib/labels.h"
#include "STKLib/common.h"
#include "STKLib/Lattice.h"

////////////////////////////////////////////////////////////////////////////////
// THE STANDARD HEADERS
////////////////////////////////////////////////////////////////////////////////
#include <stdlib.h>
#include <string.h>
#include <cstdio>
#include <malloc.h>

#if HAVE_UNISTD_H
# include <unistd.h>
#else
# include <getopt.h>
#endif


// we will be using the STK namespace ..........................................
using namespace STK;

#include <map>
std::map<char *,float> countMap;

//******************************************************************************
//******************************************************************************
void usage(char *progname)
{
  char *tchrptr;
  if ((tchrptr = strrchr(progname, '\\')) != NULL) progname = tchrptr+1;
  if ((tchrptr = strrchr(progname, '/')) != NULL) progname = tchrptr+1;
  fprintf(stderr,
"\n%s version " MODULE_VERSION "\n"
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
////" -q s       Output lattice formating ABtvaldmn              tvaldmn\n"
" -r f       Pronunciation prob scale factor                 1.0\n"
" -s f       Grammar scale factor                            1.0\n"
" -t f [i l] Set pruning to f [inc limit]                    Off\n"
" -u i       set pruning max active                          0\n"
//"  -v f    Set word end pruning threshold                  0.0\n"
" -w [f]     Recognise from network                          Off\n"
" -x s       Extension for hmm files                         None\n"
" -y s       Output transcription file extension             rec\n"
" -z s       generate lattices with extension s              off\n"
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
"\n", progname, progname, progname);
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
" -u r   MAXACTIVEMODELS"
" -w o   RECOGNET"
" -x r   SOURCEMODELEXT"
" -y r   TARGETTRANSCEXT"
" -z r   LATTICEEXT"
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



  ModelSet                      hset;
  Decoder<DecoderNetwork>        decoder;
  Decoder<_CompactDecoderNetwork> compact_decoder;

  Decoder<DecoderNetwork>*       p_decoder;

  FILE *                        lfp = NULL;
  FILE *                        ilfp = NULL;
  
  Matrix<FLOAT>                 feature_matrix;
  FeatureRepository             feature_repo;

  FLOAT                         like;
  int                           i;
  int                           fcnt = 0;
  Label *                       labels;

  Lattice                       lattice;
  char                          line[1024];
  char                          label_file[1024];
  const char *                  cchrptr;
  
  MyHSearchData                 nonCDphHash;
  MyHSearchData                 phoneHash;
  MyHSearchData                 dictHash;
  MyHSearchData                 cfgHash;
  
  int                           alignment = (AlignmentType) WORD_ALIGNMENT;

  double                        word_penalty;
  double                        model_penalty;
  double                        grammar_scale;
  double                        transp_scale;
  double                        outprb_scale;
  double                        pronun_scale;
  double                        occprb_scale;
  double                        posterior_scale = 1.0;
  double                        state_pruning;
  double                        stprn_step;
  double                        stprn_limit;
  double                        poster_prune;
         
  const char*                   hmm_dir;
  const char*                   hmm_ext;
  const char*                   out_lbl_dir;
  const char*                   out_lbl_ext;
  const char*                   in_lbl_dir;
  const char*                   in_lbl_ext;
  const char*                   out_MLF;
  const char*                   in_MLF;
  const char*                   network_file;
  const char*                   hmm_list;
  const char*                   dictionary;
        char*                   script;
        char*                   mmf;
  const char*                   label_filter;
  const char*                   net_filter;
  const char*                   label_ofilter;
  const char*                   net_ofilter;
        char*                   cmn_path;
        char*                   cmn_file;
  const char*                   cmn_mask;
        char*                   cvn_path;
        char*                   cvn_file;
  const char*                   cvn_mask;
  const char*                   cvg_file;
  const char*                   mmf_dir;
  const char*                   mmf_mask;
  const char*                   lat_ext;
  int                           trace_flag;
  int                           targetKind;
  int                           derivOrder;
  int*                          derivWinLengths;
  int                           startFrmExt;
  int                           endFrmExt;
  int                           max_active;
  int                           min_active;
  bool                          baum_welch;
  bool                          swap_features;
  bool                          htk_compat;
  bool                          time_pruning;
  bool                          compactNetworkRepresentation = false;
  enum TranscriptionFormat {TF_HTK, TF_STK, TF_CSTK} in_transc_fmt, out_transc_fmt;
  int                           notInDictAction = WORD_NOT_IN_DIC_UNSET;
  
  ExpansionOptions              expOptions  = {0};
  ExpansionOptions              emptyExpOpts= {0}; 
                                
  STKNetworkOutputFormat        in_net_fmt  = {0};
  STKNetworkOutputFormat        out_net_fmt = {0};
  STKNetworkOutputFormat        tmp_net_fmt = {0};
//  LabelFormat                   out_lbl_fmt = {0};
//  LabelFormat                   in_lbl_fmt  = {0};

  bool  print_all_options; 

  extern int STK::nbest_lattices;



int main(int argc, char *argv[]) 
{
                                emptyExpOpts.mStrictTiming   = true;
                                emptyExpOpts.mNoOptimization = true;
                                emptyExpOpts.mNoWordExpansion= true;  
//  in_lbl_fmt.TIMES_OFF = 1;
  if (argc == 1) 
    usage(argv[0]);

  // initialize the HMM set
  hset.Init();

  if (!my_hcreate_r(100,  &dictHash)
  ||  !my_hcreate_r(100,  &phoneHash)
  ||  !my_hcreate_r(100,  &cfgHash)) 
  {
    Error("Insufficient memory");
  }

  i = ParseOptions(argc, argv, optionStr, SNAME, &cfgHash);

  // if in HTK compatible mode, we expect the source dictionary and dictionary
  // to come
  htk_compat = GetParamBool(&cfgHash, SNAME":HTKCOMPAT", false);
  if (htk_compat) 
  {
    if (argc == i) Error("Dictionary file name expected");
    InsertConfigParam(&cfgHash, SNAME":SOURCEDICT", argv[i++], '-');
    if (argc == i) Error("HMM list file name expected");
    InsertConfigParam(&cfgHash, SNAME":SOURCEHMMLIST",    argv[i++], '-');
  }
  
  // the rest of the parameters are the feature files
  for (; i < argc; i++) 
  {
    feature_repo.AddFile(argv[i]);
  }
  
  // extract the feature parameters
  targetKind =   GetDerivParams(&cfgHash, &derivOrder, &derivWinLengths,
                                &startFrmExt, &endFrmExt,
                                &cmn_path, &cmn_file, &cmn_mask,
                                &cvn_path, &cvn_file, &cvn_mask, &cvg_file,
                                SNAME":", 0);
  // misc parameter extraction
  expOptions.mCDPhoneExpansion =
                 GetParamBool(&cfgHash,SNAME":ALLOWXWRDEXP",    false);
  expOptions.mRespectPronunVar
               = GetParamBool(&cfgHash,SNAME":RESPECTPRONVARS", false);
  expOptions.mStrictTiming
               = GetParamBool(&cfgHash,SNAME":EXACTTIMEMERGE",  false);
  expOptions.mNoWeightPushing
               =!GetParamBool(&cfgHash,SNAME":WEIGHTPUSHING",   true);
  expOptions.mNoOptimization
               =!GetParamBool(&cfgHash,SNAME":MINIMIZENET",     expOptions.mCDPhoneExpansion
                                                                ? true : false);
  expOptions.mRemoveWordsNodes
               = GetParamBool(&cfgHash,SNAME":REMEXPWRDNODES",  false);
//  in_lbl_fmt.TIMES_OFF =
  time_pruning = GetParamBool(&cfgHash,SNAME":TIMEPRUNING",     false);
  in_net_fmt.mNoTimes = !time_pruning;
//  in_lbl_fmt.left_extent  = -100 * (long long) (0.5 + 1e5 *
  in_net_fmt.mStartTimeShift =
                 GetParamFlt(&cfgHash, SNAME":STARTTIMESHIFT",  0.0);
//  in_lbl_fmt.right_extent =  100 * (long long) (0.5 + 1e5 *
  in_net_fmt.mEndTimeShift =
                 GetParamFlt(&cfgHash, SNAME":ENDTIMESHIFT",    0.0);
  in_net_fmt.mNoAcousticLikes = 
               !!GetParamBool(&cfgHash,SNAME":ADDACSCORES",     true);
  baum_welch   = GetParamBool(&cfgHash,SNAME":EVALUATION",      false);
  swap_features=!GetParamBool(&cfgHash,SNAME":NATURALREADORDER",isBigEndian());
  gpFilterWldcrd=GetParamStr(&cfgHash, SNAME":HFILTERWILDCARD", "$");
  gpScriptFilter=GetParamStr(&cfgHash, SNAME":HSCRIPTFILTER",   NULL);
  gpParmFilter = GetParamStr(&cfgHash, SNAME":HPARMFILTER",     NULL);
  gpHListFilter= GetParamStr(&cfgHash, SNAME":HMMLISTFILTER",   NULL);
  gpMmfFilter  = GetParamStr(&cfgHash, SNAME":HMMDEFFILTER",    NULL);
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
  poster_prune = GetParamFlt(&cfgHash, SNAME":POSTERIORPRUNING",0.0);
  max_active   = GetParamInt(&cfgHash, SNAME":MAXACTIVEMODELS", 0);
  min_active   = GetParamInt(&cfgHash, SNAME":MINACTIVEMODELS", 0);
  trace_flag   = GetParamInt(&cfgHash, SNAME":TRACE",           0);
  script =(char*)GetParamStr(&cfgHash, SNAME":SCRIPT",          NULL);
  mmf    =(char*)GetParamStr(&cfgHash, SNAME":SOURCEMMF",       NULL);
  
  mmf_dir      = GetParamStr(&cfgHash, SNAME":MMFDIR",          ".");
  mmf_mask     = GetParamStr(&cfgHash, SNAME":MMFMASK",         NULL);

  STK::nbest_lattices = GetParamInt(&cfgHash, SNAME":NBEST",    0);


  lat_ext      = GetParamStr(&cfgHash, SNAME":LATTICEEXT",      NULL);
  print_all_options = GetParamBool(&cfgHash,SNAME":PRINTALLOPTIONS", false);

  in_transc_fmt = (TranscriptionFormat) 
    GetParamEnum(&cfgHash, SNAME":SOURCETRANSCFMT",
      !network_file && htk_compat ? TF_HTK : TF_STK,
      "HTK", TF_HTK, "STK", TF_STK, "CSTK", TF_CSTK, NULL);
  
  if (in_transc_fmt == TF_CSTK)
  {
    in_transc_fmt = TF_STK;
    compactNetworkRepresentation = true;
  }

  out_transc_fmt = (TranscriptionFormat) 
    GetParamEnum(&cfgHash,SNAME":TARGETTRANSCFMT", 
        htk_compat ?  : TF_STK ,
        "HTK", TF_HTK, "STK", TF_STK, NULL);
	
  if (lat_ext != NULL) {
    out_transc_fmt = TF_STK;
    out_lbl_ext = lat_ext;
  }

  cchrptr      = GetParamStr(&cfgHash, SNAME":LABELFORMATING",  "");
  if (out_transc_fmt == TF_HTK) {
    while (*cchrptr) 
    {
      switch (*cchrptr++) 
      {
        case 'N': out_net_fmt.mScoreNorm       = 1; break;
        case 'S': out_net_fmt.mNoAcousticLikes = 1; break;
        case 'C': out_net_fmt.mCentreTimes     = 1; break;
        case 'T': out_net_fmt.mNoTimes         = 1; break;
        case 'W': out_net_fmt.mNoWordNodes     = 1; break;
        case 'M': out_net_fmt.mNoModelNodes    = 1; break;
        case 'F': out_net_fmt.mFrameScores     = 1; break;
        case 'X': out_net_fmt.mStripTriphones  = 1; break;
        default:
          Warning("Unknown label formating flag '%c' ignored (NCSTWMF)", 
              *cchrptr);
      }
    }
  }
  
  cchrptr      = GetParamStr(&cfgHash, SNAME":NETFORMATING",  "");
  if (out_transc_fmt == TF_STK) {
    if (*cchrptr) 
    {
      out_net_fmt.mNoLMLikes        = 1;
      out_net_fmt.mNoAcousticLikes  = 1;
      out_net_fmt.mNoTimes          = 1;
      out_net_fmt.mNoPronunVars     = 1;
    }
  
    while (*cchrptr) 
    {
      switch (*cchrptr++) 
      {
        case 'R': out_net_fmt.mBase62Labels             = 1; // reticent
                  out_net_fmt.mLinNodeSeqs              = 1;
                  out_net_fmt.mNoDefaults               = 1; break;
        case 'V': out_net_fmt.mArcDefsWithJ             = 1;
                  out_net_fmt.mAllFieldNames            = 1; break;
        case 'J': out_net_fmt.mArcDefsToEnd             = 1; break;
        case 'W': out_net_fmt.mNoWordNodes              = 1; break;
        case 'M': out_net_fmt.mNoModelNodes             = 1; break;
        case 'X': out_net_fmt.mStripTriphones           = 1; break;
        case 't': out_net_fmt.mNoTimes                  = 0; break;
        case 's': out_net_fmt.mStartTimes               = 1; break;
        case 'v': out_net_fmt.mNoPronunVars             = 0; break;
        case 'a': out_net_fmt.mNoAcousticLikes          = 0; break;
        case 'l': out_net_fmt.mNoLMLikes                = 0; break;
        case 'p': out_net_fmt.mAproxAccuracy            = 1; break;
        default:
          Warning("Unknown net formating flag '%c' ignored (JMRVWXalpstv)", 
              *cchrptr);
      }
    }
  }

  if (GetParamBool(&cfgHash, SNAME":STATEALIGNMENT", false)) 
  {
    alignment |= STATE_ALIGNMENT;
  }

  if (GetParamBool(&cfgHash, SNAME":MODELALIGNMENT", false)) 
  {
    alignment |= MODEL_ALIGNMENT;
  }
  
  if (GetParamBool(&cfgHash, SNAME":PRINTCONFIG", false)) 
  {
    PrintConfig(&cfgHash);
  }

  if (GetParamBool(&cfgHash, SNAME":PRINTVERSION", false)) 
  {
    puts("Version: "MODULE_VERSION"\n");
  }

  if (!GetParamBool(&cfgHash,SNAME":ACCEPTUNUSEDPARAM", false)) 
  {
    CheckCommandLineParamUse(&cfgHash);
  }

  if (print_all_options) 
  {
    print_registered_parameters();
  }
 

  // initialize the feature repository
  feature_repo.Init(swap_features, startFrmExt, endFrmExt, targetKind, 
     derivOrder, derivWinLengths, cmn_path, cmn_mask, cvn_path, cvn_mask,
     cvg_file);


  if (NULL != script) 
    feature_repo.AddFileList(script, gpScriptFilter); 
  
  // parse the given MMF file(s)
  if (NULL != mmf)
  {
    for (mmf=strtok(mmf, ","); mmf != NULL; mmf=strtok(NULL, ",")) 
    {
      hset.ParseMmf(mmf, NULL);
    }
  }

  // parse the HMM list
  if (hmm_list != NULL) 
    hset.ReadHMMList(hmm_list, hmm_dir, hmm_ext);
  
  nonCDphHash = hset.MakeCIPhoneHash();

  
  // a dictionary was specified
  if (dictionary != NULL) 
  {
    ReadDictionary(dictionary, &dictHash, &phoneHash);
    notInDictAction  = WORD_NOT_IN_DIC_WARN;

    if (expOptions.mRespectPronunVar) 
    {
      notInDictAction |= PRON_NOT_IN_DIC_ERROR;
    }
  }
  
  if (dictHash.mNEntries == 0) 
    expOptions.mNoWordExpansion = 1;

  transc_filter = transc_filter != NULL   ? transc_filter :
                  in_transc_fmt == TF_STK ? net_filter    : label_filter;

  transc_ofilter = transc_ofilter != NULL ? transc_ofilter :
                   out_transc_fmt == TF_STK  ? net_ofilter : label_ofilter;

  if (dictionary == NULL && in_transc_fmt == TF_HTK) 
  {
    // Word alignment is inpossible in this case.
    if (alignment & WORD_ALIGNMENT) 
    {
      alignment &= ~WORD_ALIGNMENT;
      alignment |= MODEL_ALIGNMENT;
    }
    out_net_fmt.mNoWordNodes = 1;
  }


  if (!compactNetworkRepresentation)
    p_decoder = &decoder;
  else
    p_decoder = reinterpret_cast<Decoder<DecoderNetwork>* > (&compact_decoder);

  
  /*
  IStkStream                      input_stream;    
  input_stream.open(network_file, ios::in, transc_filter ? transc_filter : "");
  
  if (!input_stream.good())
    Error("Cannot open network file: %s", network_file);
  
  ilfp = input_stream.file();

  ReadSTKNetwork(
      ilfp, 
      &dictHash,
      &phoneHash, 
      notInDictAction, 
      in_net_fmt,
      feature_repo.CurrentHeader().mSamplePeriod, 
      network_file, 
      NULL,
      true,
      compact_decoder.rNetwork());

  input_stream.close();

  OStkStream output_stream;
  output_stream.open("test_output.net");

  if (!output_stream.good())
    Error("Cannot open test output");
   
  WriteSTKNetwork(output_stream.file(), compact_decoder.rNetwork(), out_net_fmt, 
      feature_repo.CurrentHeader().mSamplePeriod, "test_label", "test_mlf.mlf",
      0.0, 0.0, 1.0);
  
  output_stream.close();
  compact_decoder.Clear();

  /*
  if (!compactNetworkRepresentation)
    SViteApp< Decoder<DecoderNetwork> >();
  else
    SViteApp< Decoder<_CompactDecoderNetwork> >();
    */

  
  if (network_file) 
  { // Unsupervised training
    IStkStream                      input_stream;    
    
    input_stream.open(network_file, ios::in, transc_filter ? transc_filter : "");

    
    if (!input_stream.good())
    {
      Error("Cannot open network file: %s", network_file);
    }
    
    ilfp = input_stream.file();
    
    if (in_transc_fmt == TF_HTK) 
    {
      //:TODO:
      // header.mSamplePeriod not initialized yet... 
      labels = ReadLabels(
          ilfp, 
          dictionary ? &dictHash : &phoneHash,
          dictionary ? UL_ERROR : UL_INSERT, 
          in_net_fmt,
          feature_repo.CurrentHeader().mSamplePeriod, 
          network_file, 
          NULL, 
          NULL);
              
      decoder.rNetwork().BuildFromLabels(labels, dictionary ? NT_WORD : NT_PHONE);
              
      ReleaseLabels(labels);

      if (!compactNetworkRepresentation)
        decoder.rNetwork().ExpansionsAndOptimizations(
            expOptions, 
            in_net_fmt, 
            &dictHash,
            &nonCDphHash, 
            &phoneHash,
            word_penalty,
            model_penalty,
            grammar_scale,
            posterior_scale);
    }
    else if (in_transc_fmt == TF_STK) 
    {
      //:TODO:
      // header.mSamplePeriod not initialized yet... 

      if (compactNetworkRepresentation)
      {
        ReadSTKNetwork(
          ilfp, 
          &dictHash,
          &phoneHash, 
          notInDictAction, 
          in_net_fmt,
          feature_repo.CurrentHeader().mSamplePeriod, 
          network_file, 
          NULL,
          compactNetworkRepresentation,
          compact_decoder.rNetwork());
      }
      else
      {
        ReadSTKNetwork(
          ilfp, 
          &dictHash,
          &phoneHash, 
          notInDictAction, 
          in_net_fmt,
          feature_repo.CurrentHeader().mSamplePeriod, 
          network_file, 
          NULL,
          compactNetworkRepresentation,
          decoder.rNetwork());

        decoder.rNetwork().ExpansionsAndOptimizations(
          expOptions, 
          in_net_fmt, 
          &dictHash,
          &nonCDphHash, 
          &phoneHash,
          word_penalty,
          model_penalty,
          grammar_scale,
          posterior_scale);
      }
    }
    else 
    {
      Error("Too bad. What did you do ?!?");
    }

    if (compactNetworkRepresentation)
      compact_decoder.Init(&hset, NULL/*, compactNetworkRepresentation*/);
    else
      decoder.Init(&hset, NULL/*, compactNetworkRepresentation*/);

  } 
  else 
  {
    ilfp = OpenInputMLF(in_MLF);
  }

  lfp = OpenOutputMLF(out_MLF);


  //////////////////////////////////////////////////////////////////////////////
  // read consequently all the feature files 
  for (feature_repo.Rewind(); !feature_repo.EndOfList(); feature_repo.MoveNext())
  {
    if (trace_flag & 1) 
    {
      TraceLog("Processing file %d/%d '%s'", ++fcnt, feature_repo.QueueSize(), 
          feature_repo.Current().Physical().c_str());
    }
    
    // read the feature matrix .................................................
    feature_repo.ReadFullMatrix(feature_matrix);

    if (hset.mInputVectorSize != static_cast<int>(feature_matrix.Cols()))
    {
      Error("Vector size [%d] in '%s' is incompatible with HMM set [%d]",
          feature_matrix.Cols(), feature_repo.Current().Physical().c_str(), 
          hset.mInputVectorSize);
    }

    // parse per-datafile models ...............................................
    if (mmf_mask != NULL) 
    {
      static string    last_speaker_mmf;
      string           speaker_mmf;

      ProcessMask(feature_repo.Current().Logical(), mmf_mask, speaker_mmf);
        
      if (last_speaker_mmf != speaker_mmf) 
      {
        hset.ParseMmf((string(mmf_dir) + "/" + speaker_mmf).c_str(), NULL);
        last_speaker_mmf = speaker_mmf;
      }
    }
    
    // read the network file if given ..........................................
    if (!network_file) 
    {
      // construct the name
      strcpy(label_file, feature_repo.Current().Logical().c_str());

      // open the label
      ilfp = OpenInputLabelFile(label_file, in_lbl_dir, 
          in_lbl_ext ? in_lbl_ext :
          in_transc_fmt == TF_STK ? "net" : "lab",
          ilfp, in_MLF);

      // decide what to do depending on input format
      if (in_transc_fmt == TF_HTK) 
      {
        labels = ReadLabels(ilfp, dictionary ? &dictHash : &phoneHash, 
            dictionary ? UL_ERROR : UL_INSERT, in_net_fmt,
            feature_repo.CurrentHeader().mSamplePeriod, label_file, in_MLF, NULL);

        decoder.rNetwork().BuildFromLabels(labels, dictionary ? NT_WORD : NT_PHONE);

        ReleaseLabels(labels);
      } 
      else if (in_transc_fmt == TF_STK) 
      {
        if (trace_flag & 2)
          TraceLog("Recognition network: Loading \"%s\"", label_file);

        if (!compactNetworkRepresentation)
        {
          ReadSTKNetwork(ilfp, &dictHash, &phoneHash, notInDictAction,
              in_net_fmt, feature_repo.CurrentHeader().mSamplePeriod, label_file,
              in_MLF, compactNetworkRepresentation,  decoder.rNetwork());
        }
        else
        {
          ReadSTKNetwork(ilfp, &dictHash, &phoneHash, notInDictAction,
              in_net_fmt, feature_repo.CurrentHeader().mSamplePeriod, label_file,
              in_MLF, compactNetworkRepresentation,  compact_decoder.rNetwork());
        }
      } 
      else 
      {
        Error("Too bad. What did you do ?!?");
      }

      // we perform optimizations of not in compact representation
      if (!compactNetworkRepresentation)
      {
        if (trace_flag & 2)
          TraceLog("Recognition network: Performing expansions and optimizations");

        decoder.rNetwork().ExpansionsAndOptimizations(expOptions, in_net_fmt, &dictHash,
            &nonCDphHash, &phoneHash,
            word_penalty,
            model_penalty,
            grammar_scale,
            posterior_scale);

        decoder.Init(&hset, NULL/*, false*/);
      }
      else
      {
        compact_decoder.Init(&hset, NULL/*, true*/);
      }

      CloseInputLabelFile(ilfp, in_MLF);
    }

    p_decoder->mTimePruning       = time_pruning;
    p_decoder->mWPenalty          = word_penalty;
    p_decoder->mMPenalty          = model_penalty;
    p_decoder->mLmScale           = grammar_scale;
    p_decoder->mPronScale         = pronun_scale;
    p_decoder->mTranScale         = transp_scale;
    p_decoder->mOutpScale         = outprb_scale;
    p_decoder->mOcpScale          = occprb_scale;
    p_decoder->mAlignment         = alignment;
    p_decoder->mPruningThresh     = state_pruning > 0.0 ? state_pruning : -LOG_0;
    p_decoder->mMaxActiveModels   = max_active;
    p_decoder->mMinActiveModels   = min_active;
    p_decoder->mLatticeGeneration = (lat_ext != NULL);
    

    if(p_decoder->mLatticeGeneration)
    {
      p_decoder->mAlignment = WORD_ALIGNMENT | MODEL_ALIGNMENT;
    }
    else
    {
      if (alignment & STATE_ALIGNMENT && out_net_fmt.mNoModelNodes) p_decoder->mAlignment &= ~MODEL_ALIGNMENT;
      if (alignment & MODEL_ALIGNMENT && out_net_fmt.mNoWordNodes)  p_decoder->mAlignment &= ~WORD_ALIGNMENT;
      if (alignment & STATE_ALIGNMENT && out_net_fmt.mFrameScores)  p_decoder->mAlignment |=  FRAME_ALIGNMENT;
    }

    fflush(stdout);
    if (trace_flag & 2)
      TraceLog("Recognition network: Performing recognition...");

    for (;;) 
    {
      if (!compactNetworkRepresentation)
      {
        decoder.ViterbiInit();

        decoder.PassTokenInNetwork = 
          decoder.mLatticeGeneration ? &Decoder<DecoderNetwork>::PassTokenMaxForLattices :
          baum_welch                 ? &Decoder<DecoderNetwork>::PassTokenSum 
                                     : &Decoder<DecoderNetwork>::PassTokenMax;
                                 
        decoder.PassTokenInModel   = baum_welch ? &Decoder<DecoderNetwork>::PassTokenSum 
                                                : &Decoder<DecoderNetwork>::PassTokenMax;

        for (i = 0; i < feature_matrix.Rows(); i++) 
        {
          decoder.ViterbiStep(feature_matrix[i]);
        }
        
        like = decoder.ViterbiDone(&labels, &lattice);
      }
      else
      {
        compact_decoder.ViterbiInit();

        compact_decoder.PassTokenInNetwork = 
          compact_decoder.mLatticeGeneration ? &Decoder<_CompactDecoderNetwork>::PassTokenMaxForLattices :
          baum_welch                 ? &Decoder<_CompactDecoderNetwork>::PassTokenSum 
                                     : &Decoder<_CompactDecoderNetwork>::PassTokenMax;
                                 
        compact_decoder.PassTokenInModel   = baum_welch ? &Decoder<_CompactDecoderNetwork>::PassTokenSum 
                                                : &Decoder<_CompactDecoderNetwork>::PassTokenMax;

        for (i = 0; i < feature_matrix.Rows(); i++) 
        {
          compact_decoder.ViterbiStep(feature_matrix[i]);
        }
        
        like = compact_decoder.ViterbiDone(&labels, &lattice);
      }


      if (labels) 
        break;

      if (p_decoder->mPruningThresh <= LOG_MIN 
      || (stprn_step <= 0.0) 
      || ((p_decoder->mPruningThresh += stprn_step) > stprn_limit)) 
      {
        Warning("No tokens survived");
        break;
      }

      Warning("No tokens survived, trying pruning threshold: %.2f", 
          p_decoder->mPruningThresh);

      lattice.Clear();
    }

    if (trace_flag & 1 && labels) 
    {
      Label* label;
      int    n_frames = feature_matrix.Rows() - hset.mTotalDelay;

      for (label = labels; label->mpNextLevel != NULL; 
           label = label->mpNextLevel)
      {  }

      for (; label != NULL; label = label->mpNext) 
      {
        if(label->mpName != NULL)
          fprintf(stdout, "%s ", label->mpName);
      }

      TraceLog(" ==  [%d frames] %f", n_frames, like / n_frames);

//      TraceLog(" ==  [%d frames] %f (Act: %f/%f)", n_frames, like / n_frames,
//                   (float) p_decoder->mNActiveModelsForUtterance / n_frames, 
//                   (float) p_decoder->mNActiveTokensForUtterance / n_frames);
    }

    if (!lattice.IsEmpty())
    {
      tmp_net_fmt = out_net_fmt;
      tmp_net_fmt.mNoLMLikes        = 0;
      tmp_net_fmt.mNoAcousticLikes  = 0;
      
      if (trace_flag & 2)
        TraceLog("Lattice: Performing optimizations (first pass)...");
	
      lattice.ExpansionsAndOptimizations(emptyExpOpts, tmp_net_fmt, NULL, NULL, 
          NULL, word_penalty, model_penalty, grammar_scale, posterior_scale);

      if (trace_flag & 2)
        TraceLog("Lattice: Computing posterior probabilities...");
      lattice.ForwardBackward(word_penalty, model_penalty, grammar_scale, 
          posterior_scale, true); // set false for "non-viterbi" posterior pruning
      
      if (trace_flag & 2)
        TraceLog("Lattice: Performing posterior pruning...");
	
      lattice.PosteriorPrune(poster_prune  > 0.0 ? poster_prune  :
                             state_pruning > 0.0 ? state_pruning : -LOG_0,
			     word_penalty, model_penalty, grammar_scale, posterior_scale);
			     
//      lattice.PosteriorExpectedCounts(countMap);
      
      lattice.FreePosteriors();


      if (trace_flag & 2) {
        TraceLog("Lattice: Performing optimizations (second pass)...");
      }
      lattice.ExpansionsAndOptimizations(emptyExpOpts, out_net_fmt, NULL, NULL, 
          NULL, word_penalty, model_penalty, grammar_scale, posterior_scale);
    }

    strcpy(label_file, feature_repo.Current().Logical().c_str());
    lfp = OpenOutputLabelFile(label_file, out_lbl_dir, out_lbl_ext, lfp, out_MLF);

    // write the output ........................................................
    //
    if (out_transc_fmt == TF_STK)
    {
      WriteSTKNetwork(lfp, lattice, out_net_fmt, 
          feature_repo.CurrentHeader().mSamplePeriod, label_file, out_MLF,
          p_decoder->mWPenalty, p_decoder->mMPenalty, p_decoder->mLmScale);
          
      // we are not needing the lattice anymore, so free it from memory
      lattice.Clear();
    } 
    else // if (out_transc_fmt == TF_HTK) 
    {
      WriteLabels(lfp, labels, out_net_fmt, 
          feature_repo.CurrentHeader().mSamplePeriod, label_file, out_MLF);
    }
/*    else 
    {
      // create temporary linear network from labels, just to store it
      DecoderNetwork tmp_net(labels, 
          alignment & (MODEL_ALIGNMENT|STATE_ALIGNMENT) ? NT_MODEL : NT_WORD);
      
      WriteSTKNetwork(lfp, tmp_net, out_net_fmt, 
          feature_repo.CurrentHeader().mSamplePeriod, label_file, out_MLF,
          p_decoder->mWPenalty, p_decoder->mMPenalty, p_decoder->mLmScale);
    }*/

    CloseOutputLabelFile(lfp, out_MLF);
    ReleaseLabels(labels);

    if (!network_file) 
    {
      if (!compactNetworkRepresentation)
        decoder.Clear();
      else
        compact_decoder.Clear();
    }
  } // while (!feature_repo.EndOfList())


  // clean up ..................................................................
  if (network_file) 
  {
    if (!compactNetworkRepresentation)
      decoder.Clear();
    else
      compact_decoder.Clear();
  }
  
  for(std::map<char*,float>::iterator iter = countMap.begin(); iter != countMap.end(); ++iter ) {
    std::cout << iter->first << " " << iter->second << std::endl; 
  }

  hset.Release();  
  // my_hdestroy_r(&labelHash,   0);
  my_hdestroy_r(&phoneHash,   1);
  my_hdestroy_r(&nonCDphHash, 0);
  FreeDictionary(&dictHash);
  
  for (unsigned int i = 0; i < cfgHash.mNEntries; i++) 
    free(cfgHash.mpEntry[i]->data);
  
  my_hdestroy_r(&cfgHash, 1);
  
  if (out_MLF) {
    my_fclose(lfp);
  }

  //my_hdestroy_r(&cfgHash, 0);

  return 0;
}

//HVite -T 05 -H models -w wdnet dict words4 MAL_4379315A.fea > htk.log

