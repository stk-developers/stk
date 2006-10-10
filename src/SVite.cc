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


#define MODULE_VERSION "0.4 "__TIME__" "__DATE__" "SVN_ID


////////////////////////////////////////////////////////////////////////////////
// STK HEADERS
////////////////////////////////////////////////////////////////////////////////
#include "STKLib/Features.h"

#include "STKLib/Viterbi.h"
#include "STKLib/Models.h"
#include "STKLib/fileio.h"
#include "STKLib/labels.h"
#include "STKLib/common.h"


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


//******************************************************************************
//******************************************************************************
void usage(char *progname)
{
  char *tchrptr;
  if ((tchrptr = strrchr(progname, '\\')) != NULL) progname = tchrptr+1;
  if ((tchrptr = strrchr(progname, '/')) != NULL) progname = tchrptr+1;
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
////" -q s       Output lattice formating ABtvaldmn              tvaldmn\n"
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

int main(int argc, char *argv[]) 
{
  ModelSet                      hset;
  Network                       net;
  FILE *                        lfp = NULL;
  FILE *                        ilfp = NULL;
  
  Matrix<FLOAT>                 feature_matrix;
  FeatureRepository             feature_repo;

  FLOAT                         like;
  int                           i;
  int                           fcnt = 0;
  Label *                       labels;
  Node *                        pLattice;
  char                          line[1024];
  char                          label_file[1024];
  const char *                  cchrptr;
  
  MyHSearchData        nonCDphHash;
  MyHSearchData        phoneHash;
  MyHSearchData        dictHash;
  MyHSearchData        cfgHash;
  
  FileListElem *                feature_files = NULL;
  FileListElem *                file_name = NULL;
  FileListElem **               last_file = &feature_files;
  
  int                           alignment = (AlignmentType) WORD_ALIGNMENT;

  double                        word_penalty;
  double                        model_penalty;
  double                        grammar_scale;
  double                        transp_scale;
  double                        outprb_scale;
  double                        pronun_scale;
  double                        occprb_scale;
  double                        state_pruning;
  double                        stprn_step;
  double                        stprn_limit;  
//double  word_pruning;
         
  const  char *                 hmm_dir;
  const  char *                 hmm_ext;
  const  char *                 out_lbl_dir;
  const char *                  out_lbl_ext;
  const char *                  in_lbl_dir;
  const char *                  in_lbl_ext;
  const char *                  out_MLF;
  const char *                  in_MLF;
  const char *                  network_file;
  const char *                  hmm_list;
  const char *                  dictionary;
  char *                        script;
  char *                        mmf;
  const char *                  label_filter;
  const char *                  net_filter;
  const char *                  label_ofilter;
  const char *                  net_ofilter;
        char *                  cmn_path;
        char *                  cmn_file;
  const char *                  cmn_mask;
        char *                  cvn_path;
        char *                  cvn_file;
  const char *                  cvn_mask;
  const char *                  cvg_file;
  const char *                  mmf_dir;
  const char *                  mmf_mask;
  const char *                  lat_ext;

  int                           trace_flag;
  int                           targetKind;
  int                           derivOrder;
  int  *                        derivWinLengths;
  int                           startFrmExt;
  int                           endFrmExt;
  bool                          baum_welch;
  bool                          swap_features;
  bool                          htk_compat;
  bool                          compactNetworkRepresentation = false;
  enum TranscriptionFormat {TF_HTK, TF_STK, TF_CSTK} in_transc_fmt, out_transc_fmt;
  int                           notInDictAction = WORD_NOT_IN_DIC_UNSET;
  
  ExpansionOptions              expOptions  = {0};
  STKNetworkOutputFormat        in_net_fmt  = {0};
  STKNetworkOutputFormat        out_net_fmt = {0};
  LabelFormat                   out_lbl_fmt = {0};
  LabelFormat                   in_lbl_fmt  = {0};
  
  
  in_lbl_fmt.TIMES_OFF = 1;

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
    //last_file = AddFileElem(last_file, argv[i]);
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
  expOptions.mNoOptimization
               =!GetParamBool(&cfgHash,SNAME":MINIMIZENET",     expOptions.mCDPhoneExpansion
                                                                ? true : false);
  expOptions.mRemoveWordsNodes
               = GetParamBool(&cfgHash,SNAME":REMEXPWRDNODES",  false);
  in_lbl_fmt.TIMES_OFF =
                !GetParamBool(&cfgHash,SNAME":TIMEPRUNING",    false);
  in_lbl_fmt.left_extent  = -100 * (long long) (0.5 + 1e5 *
                 GetParamFlt(&cfgHash, SNAME":STARTTIMESHIFT",  0.0));
  in_lbl_fmt.right_extent =  100 * (long long) (0.5 + 1e5 *
                 GetParamFlt(&cfgHash, SNAME":ENDTIMESHIFT",    0.0));
  baum_welch   = GetParamBool(&cfgHash,SNAME":EVALUATION",      false);
  swap_features=!GetParamBool(&cfgHash,SNAME":NATURALREADORDER",isBigEndian());
  gpFilterWldcrd= GetParamStr(&cfgHash,SNAME":HFILTERWILDCARD", "$");
  gpScriptFilter= GetParamStr(&cfgHash, SNAME":HSCRIPTFILTER",   NULL);
  gpParmFilter  = GetParamStr(&cfgHash, SNAME":HPARMFILTER",     NULL);
  gpHListFilter = GetParamStr(&cfgHash,SNAME":HMMLISTFILTER",   NULL);
  gpMmfFilter   = GetParamStr(&cfgHash, SNAME":HMMDEFFILTER",    NULL);
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
  
  mmf_dir      = GetParamStr(&cfgHash, SNAME":MMFDIR",          ".");
  mmf_mask     = GetParamStr(&cfgHash, SNAME":MMFMASK",         NULL);
  lat_ext      = GetParamStr(&cfgHash, SNAME":LATTICEEXT",      NULL);


  cchrptr      = GetParamStr(&cfgHash, SNAME":LABELFORMATING",  "");
  while (*cchrptr) 
  {
    switch (*cchrptr++) 
    {
      case 'N': out_lbl_fmt.SCORE_NRM = 1; break;
      case 'S': out_lbl_fmt.SCORE_OFF = 1; break;
      case 'C': out_lbl_fmt.CENTRE_TM = 1; break;
      case 'T': out_lbl_fmt.TIMES_OFF = 1; break;
      case 'W': out_lbl_fmt.WORDS_OFF = 1; break;
      case 'M': out_lbl_fmt.MODEL_OFF = 1; break;
      case 'F': out_lbl_fmt.FRAME_SCR = 1; break;
//      case 'X': out_lbl_fmt.STRIP_TRI = 1; break;
      default:
        Warning("Unknown label formating flag '%c' ignored (NCSTWMF)", 
            *cchrptr);
    }
  }

  cchrptr      = GetParamStr(&cfgHash, SNAME":NETFORMATING",  "");
  
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
        htk_compat ? TF_HTK : TF_STK ,
        "HTK", TF_HTK, "STK", TF_STK, NULL);

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
    out_lbl_fmt.WORDS_OFF = 1;
  }
  
  if (network_file) 
  { // Unsupervised training
    Node *node = NULL;
    IStkStream input_stream;    
    
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
          in_lbl_fmt,
          feature_repo.CurrentHeader().mSamplePeriod, 
          network_file, 
          NULL, 
          NULL);
              
      node = MakeNetworkFromLabels(labels, 
          dictionary ? NT_WORD : NT_PHONE);
              
      ReleaseLabels(labels);
    }
    else if (in_transc_fmt == TF_STK) 
    {
      //:TODO:
      // header.mSamplePeriod not initialized yet... 
      node = ReadSTKNetwork(
         ilfp, 
         &dictHash,
         &phoneHash, 
         notInDictAction, 
         in_lbl_fmt,
         feature_repo.CurrentHeader().mSamplePeriod, 
         network_file, 
         NULL, compactNetworkRepresentation);
    }
    else 
    {
      Error("Too bad. What did you do ?!?");
    }
                                
    if (!compactNetworkRepresentation)
      NetworkExpansionsAndOptimizations(
          node, 
          expOptions, 
          in_net_fmt, 
          &dictHash,
          &nonCDphHash, 
          &phoneHash);
                                      
    net.Init(node, &hset, NULL, compactNetworkRepresentation);
  } 
  else 
  {
    ilfp = OpenInputMLF(in_MLF);
  }

  lfp = OpenOutputMLF(out_MLF);

  // we are going to read from the feature repository
  feature_repo.Rewind();


  //////////////////////////////////////////////////////////////////////////////
  // read consequently all the feature files 
  while (!feature_repo.EndOfList())
  {
    if (trace_flag & 1) 
    {
      TraceLog("Processing file %d/%d '%s'", ++fcnt, feature_repo.QueueSize(), 
          feature_repo.FollowingPhysical().c_str());
    }
    
    // read the feature matrix .................................................
    feature_repo.ReadFullMatrix(feature_matrix);

    if (hset.mInputVectorSize != static_cast<int>(feature_matrix.Cols()))
    {
      Error("Vector size [%d] in '%s' is incompatible with HMM set [%d]",
          feature_matrix.Cols(), feature_repo.CurrentPhysical().c_str(), 
          hset.mInputVectorSize);
    }

    // parse per-datafile models ...............................................
    if (mmf_mask != NULL) 
    {
      static string    last_speaker_mmf;
      string           speaker_mmf;

      ProcessMask(feature_repo.CurrentLogical(), mmf_mask, speaker_mmf);
        
      if (last_speaker_mmf != speaker_mmf) 
      {
        hset.ParseMmf((string(mmf_dir) + "/" + speaker_mmf).c_str(), NULL);
        last_speaker_mmf = speaker_mmf;
      }
    }
    
    // read the network file if given ..........................................
    if (!network_file) 
    {
      Node* node = NULL;

      strcpy(label_file, feature_repo.CurrentLogical().c_str());

      ilfp = OpenInputLabelFile(label_file, in_lbl_dir, 
          in_lbl_ext ? in_lbl_ext :
          in_transc_fmt == TF_STK ? "net" : "lab",
          ilfp, in_MLF);

      if (in_transc_fmt == TF_HTK) 
      {
        labels = ReadLabels(ilfp, dictionary ? &dictHash : &phoneHash, 
            dictionary ? UL_ERROR : UL_INSERT, in_lbl_fmt,
            feature_repo.CurrentHeader().mSamplePeriod, label_file, in_MLF, NULL);

        node = MakeNetworkFromLabels(labels, dictionary ? NT_WORD : NT_PHONE);
        ReleaseLabels(labels);
      } 
      else if (in_transc_fmt == TF_STK) 
      {
        node = ReadSTKNetwork(ilfp, &dictHash, &phoneHash, notInDictAction,
            in_lbl_fmt, feature_repo.CurrentHeader().mSamplePeriod, label_file,
            in_MLF);
      } 
      else 
      {
        Error("Too bad. What did you do ?!?");
      }

      if (!compactNetworkRepresentation)
        NetworkExpansionsAndOptimizations(node, expOptions, in_net_fmt, &dictHash,
            &nonCDphHash, &phoneHash);

      net.Init(node, &hset, NULL, false);

      CloseInputLabelFile(ilfp, in_MLF);
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
    net.mLatticeGeneration = (lat_ext != NULL);

    if(net.mLatticeGeneration)
    {
      net.mAlignment = WORD_ALIGNMENT | MODEL_ALIGNMENT;
    }
    else
    {
      if (alignment & STATE_ALIGNMENT && out_lbl_fmt.MODEL_OFF) net.mAlignment &= ~MODEL_ALIGNMENT;
      if (alignment & MODEL_ALIGNMENT && out_lbl_fmt.WORDS_OFF) net.mAlignment &= ~WORD_ALIGNMENT;
      if (alignment & STATE_ALIGNMENT && out_lbl_fmt.FRAME_SCR) net.mAlignment |=  FRAME_ALIGNMENT;
    }

    for (;;) 
    {
      net.ViterbiInit();
      net.PassTokenInNetwork = net.mLatticeGeneration ? &PassTokenMaxForLattices :
                               baum_welch             ? &PassTokenSum : &PassTokenMax;
                               
      net.PassTokenInModel   = baum_welch ? &PassTokenSum : &PassTokenMax;

      for (i = 0; i < feature_matrix.Rows(); i++) 
      {
        net.ViterbiStep(feature_matrix[i]);
      }
      
      like = net.ViterbiDone(&labels, &pLattice);

      if (labels) 
        break;

      if (net.mPruningThresh <= LOG_MIN 
      || (stprn_step <= 0.0) 
      || ((net.mPruningThresh += stprn_step) > stprn_limit)) 
      {
        Warning("No tokens survived");
        break;
      }

      Warning("No tokens survived, trying pruning threshold: %.2f", 
          net.mPruningThresh);
    }

    if (trace_flag & 1 && labels) 
    {
      Label* label;
      int    n_frames = feature_matrix.Rows() - hset.mTotalDelay;

      for (label = labels; 
          label->mpNextLevel != NULL;
          label = label->mpNextLevel)
      {}

      for (; label != NULL; label = label->mpNext) 
      {
        if(label->mpName != NULL)
          fprintf(stdout, "%s ", label->mpName);
      }

      TraceLog(" ==  [%d frames] %f", n_frames, like / n_frames);
    }

    if (pLattice)
    {
      //NetworkExpansionsAndOptimizations(pLattice, expOptions, out_net_fmt, 
      //    NULL, NULL, NULL);
    }

    strcpy(label_file, feature_repo.CurrentLogical().c_str());
    lfp = OpenOutputLabelFile(label_file, out_lbl_dir, out_lbl_ext, lfp, out_MLF);

    if (pLattice)
    {
      WriteSTKNetwork(lfp, pLattice, out_net_fmt, 
          feature_repo.CurrentHeader().mSamplePeriod, label_file, out_MLF);
          
      FreeNetwork(pLattice);
    } 
    else if (out_transc_fmt == TF_HTK) 
    {
      WriteLabels(lfp, labels, out_lbl_fmt, 
          feature_repo.CurrentHeader().mSamplePeriod, label_file, out_MLF);
    } 
    else 
    {
      Node* node = MakeNetworkFromLabels(labels, 
          alignment & (MODEL_ALIGNMENT|STATE_ALIGNMENT) ? NT_MODEL : NT_WORD);
      
      WriteSTKNetwork(lfp, node, out_net_fmt, 
          feature_repo.CurrentHeader().mSamplePeriod, label_file, out_MLF);

      FreeNetwork(node);
    }

    CloseOutputLabelFile(lfp, out_MLF);
    ReleaseLabels(labels);

    if (!network_file) 
    {
      net.Release();
    }
  } // while (!feature_repo.EndOfList())


  // clean up ..................................................................
  if (network_file) 
  {
    net.Release();
  }
  
  hset.Release();
  
// my_hdestroy_r(&labelHash,   0);
  my_hdestroy_r(&phoneHash,   1);
  my_hdestroy_r(&nonCDphHash, 0);
  FreeDictionary(&dictHash);
  
  for (unsigned int i = 0; i < cfgHash.mNEntries; i++) 
    free(cfgHash.mpEntry[i]->data);
  
  my_hdestroy_r(&cfgHash, 1);

  while (feature_files) 
  {
    file_name     = feature_files;
    feature_files = feature_files->mpNext;
    free(file_name);
  }
  //my_hdestroy_r(&cfgHash, 0);

  return 0;
}

//HVite -T 05 -H models -w wdnet dict words4 MAL_4379315A.fea > htk.log

