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

#define SVN_DATE       "$Date: $"
#define SVN_AUTHOR     "$Author: $"
#define SVN_REVISION   "$Revision: $"
#define SVN_ID         "$Id: $"

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

using namespace STK;


void usage(char *progname)
{
  char *tchrptr;
  if ((tchrptr = strrchr(progname, '\\')) != NULL) progname = tchrptr+1;
  if ((tchrptr = strrchr(progname, '/')) != NULL) progname = tchrptr+1;
  fprintf(stderr,
"\n%s version " MODULE_VERSION "\n"
"\nUSAGE: %s [options] DataFiles...\n\n"
" Option                                                     Default\n\n"
//" -c f   Mixture pruning threshold                         10.0\n"
" -d s       Dir to find hmm definitions                     Current\n"
" -l s       Dir to store output param files                 Feature file dir\n"
" -x s       Extension for hmm files                         None\n"
" -y s       Output param file extension                     Unchanged\n"
" -A         Print command line arguments                    Off\n"
" -C cf      Set config file to cf                           Default\n"
" -D         Display configuration variables                 Off\n"
" -H mmf     Load HMM macro file mmf                         \n"
" -S f       Set script file to f                            None\n"
" -T N       Set trace flags to N                            0\n"
" -V         Print version information                       Off\n"
"\n"
" %s is Copyright (C) 2004-2005 Lukas Burget et al. and\n"
" licensed under the GNU General Public License, version 2.\n"
" Bug reports, feedback, etc, to: burget@fit.vutbr.cz\n"
"\n", progname, progname, progname);
  exit(-1);
}

#define SNAME "SFEACAT"
const char *optionStr =
" -d r   SOURCEMODELDIR"
" -l r   TARGETPARAMDIR"
" -x r   SOURCEMODELEXT"
" -y r   TARGETPARAMEXT"
" -D n   PRINTCONFIG=TRUE"
" -H l   SOURCEMMF"
" -S l   SCRIPT"
" -T r   TRACE"
" -V n   PRINTVERSION=TRUE";

int main(int argc, char *argv[]) 
{
  FILE*                   sfp;
  FILE*                   ofp = NULL;
  
  Matrix<FLOAT>           feature_matrix1;
  Matrix<FLOAT>           feature_matrix2;
  Matrix<FLOAT>           feature_matrix3;
  
  FLOAT*                  obs;
  HtkHeader               header1;
  HtkHeader               header2;
  HtkHeader               header3;
  int                     i;
  int                     fcnt = 0;
  XformInstance*          p_input = NULL;
  char                    p_line[1024];
  char                    p_out_file[1024];
  MyHSearchData  cfg_hash;
  int                     tot_frames = 0;
  FileListElem *          feature_files = NULL;
  int                     nfeature_files = 0;
  FileListElem *          file_name = NULL;
  FileListElem **         last_file = &feature_files;


  const char *            src_hmm_list;
  const char *            src_hmm_dir;
  const char *            src_hmm_ext;
        char *            src_mmf;
  const char *            out_dir;
  const char *            out_ext;
  const char *            input_name;
        char *            script;
        char *            cmn_path;
        char *            cmn_file;
  const char *            cmn_mask;
        char *            cvn_path;
        char *            cvn_file;
  const char *            cvn_mask;
  const char *            cvg_file;
  const char *            mmf_dir;
  const char *            mmf_mask;

  int                     trace_flag;
  int                     targetKind;
  int                     derivOrder;
  int *                   derivWinLengths;
  int                     startFrmExt;
  int                     endFrmExt;
  
  bool                    swap_features;
  bool                    swap_fea_out;
  RHFBuffer               rhfbuff = {0};

  
  if (argc == 1) 
    usage(argv[0]);

  if (!my_hcreate_r(100,  &cfg_hash))
    Error("Insufficient memory");
  
  
  i = ParseOptions(argc, argv, optionStr, SNAME, &cfg_hash);
//  htk_compat   = GetParamBool(&cfgHash,SNAME":HTKCOMPAT",       false);
  for (; i < argc; i++) 
  {
    last_file = AddFileElem(last_file, argv[i]);
    nfeature_files++;
  }
  
  targetKind     = GetDerivParams(&cfg_hash, &derivOrder, &derivWinLengths,
                                &startFrmExt, &endFrmExt,
                                &cmn_path, &cmn_file, &cmn_mask,
                                &cvn_path, &cvn_file, &cvn_mask, &cvg_file,
                                SNAME":", 0);

  input_name     = GetParamStr(&cfg_hash, SNAME":SOURCEINPUT",     NULL);
  swap_features  =!GetParamBool(&cfg_hash,SNAME":NATURALREADORDER", isBigEndian());
  swap_fea_out   =!GetParamBool(&cfg_hash,SNAME":NATURALWRITEORDER",isBigEndian());
  gpFilterWldcrd = GetParamStr(&cfg_hash, SNAME":HFILTERWILDCARD", "$");
  gpScriptFilter  = GetParamStr(&cfg_hash, SNAME":HSCRIPTFILTER",   NULL);
  gpParmFilter    = GetParamStr(&cfg_hash, SNAME":HPARMFILTER",     NULL);
  gpHListFilter  = GetParamStr(&cfg_hash, SNAME":HMMLISTFILTER",   NULL);
  gpMmfFilter     = GetParamStr(&cfg_hash, SNAME":HMMDEFFILTER",    NULL);
  gpParmOFilter   = GetParamStr(&cfg_hash, SNAME":HPARMOFILTER",    NULL);
  out_dir        = GetParamStr(&cfg_hash, SNAME":TARGETPARAMDIR",  NULL);
  out_ext        = GetParamStr(&cfg_hash, SNAME":TARGETPARAMEXT",  NULL);
  trace_flag     = GetParamInt(&cfg_hash, SNAME":TRACE",           0);
  script  = (char*)GetParamStr(&cfg_hash, SNAME":SCRIPT",          NULL);
  src_hmm_list   = GetParamStr(&cfg_hash, SNAME":SOURCEHMMLIST",   NULL);
  src_hmm_dir    = GetParamStr(&cfg_hash, SNAME":SOURCEMODELDIR",  NULL);
  src_hmm_ext    = GetParamStr(&cfg_hash, SNAME":SOURCEMODELEXT",  NULL);
  src_mmf = (char*)GetParamStr(&cfg_hash, SNAME":SOURCEMMF",       NULL);

  mmf_dir        = GetParamStr(&cfg_hash, SNAME":MMFDIR",          ".");
  mmf_mask       = GetParamStr(&cfg_hash, SNAME":MMFMASK",         NULL);

  bool print_all_options = GetParamBool(&cfg_hash,SNAME":PRINTALLOPTIONS", false);

  if (GetParamBool(&cfg_hash, SNAME":PRINTCONFIG", false))
    PrintConfig(&cfg_hash);
  
  if (GetParamBool(&cfg_hash, SNAME":PRINTVERSION", false))
    puts("Version: "MODULE_VERSION"\n");
                                 
  if (!GetParamBool(&cfg_hash,SNAME":ACCEPTUNUSEDPARAM", false))
    CheckCommandLineParamUse(&cfg_hash);

  if (print_all_options) 
  {
    print_registered_parameters();
  }
  
  if (NULL != script)
  {
    for (script=strtok(script, ","); script != NULL; script=strtok(NULL, ",")) 
    {
      IStkStream script_stream(script, std::ios::in, gpScriptFilter);

      if (!script_stream.good()) {
        Error("Cannot open script file %s", script);
      }
      script_stream >> std::ws;
      
      while (!script_stream.eof()) {
        script_stream.getline(p_line, 1023);
        script_stream >> std::ws;

        last_file = AddFileElem(last_file, p_line);
        nfeature_files++;
      }
      
      script_stream.close();
    }
  }
  
  

  for (file_name=feature_files; file_name != NULL; file_name=file_name->mpNext) 
  {
    if (trace_flag & 1) TraceLog("Processing record %d/%d '%s'",
                                ++fcnt, nfeature_files, file_name->logical);
                                
    std::string input_file1;
    std::string input_file2;
    std::string output_file;

    std::stringstream file_record(file_name->mpPhysical);
    file_record >> input_file1 >> input_file2 >> output_file;


    ReadHTKFeatures(input_file1.c_str(), swap_features,
                    startFrmExt, endFrmExt, targetKind,
                    derivOrder, derivWinLengths, &header1,
                    cmn_path, cvn_path, cvg_file, &rhfbuff, feature_matrix1);
                            
    ReadHTKFeatures(input_file2.c_str(), swap_features,
                    startFrmExt, endFrmExt, targetKind,
                    derivOrder, derivWinLengths, &header2,
                    cmn_path, cvn_path, cvg_file, &rhfbuff, feature_matrix2);

    if (feature_matrix1.Rows() != feature_matrix2.Rows()) {
      Error("Numbers of frames differ [%d] != [%d]",
          feature_matrix1.Rows(), feature_matrix2.Rows());
    }

    // Initialize output matrix
    feature_matrix3.Init(feature_matrix1.Rows(), feature_matrix1.Cols() + feature_matrix2.Cols());
                            
    // Stack two input matrices
    for (i = 0; i < feature_matrix1.Rows(); ++i) {
      int j;
      for (j = 0; j < feature_matrix1.Cols(); j++) {
        feature_matrix3[i][j] = feature_matrix1[i][j];
      }

      int offset = j;

      for (j = 0; j < feature_matrix2.Cols(); j++) {
        feature_matrix3[i][j + offset] = feature_matrix2[i][j];
      }
    }


    if ((ofp = fopen(output_file.c_str(), "wb")) == NULL)
      Error("Cannot open output feature file: '%s'", output_file.c_str());

    if (WriteHTKHeader(ofp, header3, 1))
      Error("Cannot write to output feature file: '%s'", output_file.c_str());
    
    targetKind = header1.mSampleKind;
    if (WriteHTKFeatures(
          ofp,
	  header1.mSamplePeriod,
	  p_input ? (PARAMKIND_USER | targetKind & PARAMKIND_C) : targetKind,
	  swap_fea_out,
          feature_matrix3)) 
    {
      Error("Cannot write to output feature file: '%s'", output_file.c_str());
    }

    tot_frames += header1.mNSamples;
    
    if(trace_flag & 1) 
      TraceLog("[%d frames]", header1.mNSamples );
      
    fclose(ofp);

    feature_matrix1.Destroy();
    feature_matrix2.Destroy();
    feature_matrix3.Destroy();
  }
  
  if (trace_flag & 2) 
    TraceLog("Total number of frames: %d", tot_frames);

  free(derivWinLengths);

  while (feature_files) 
  {
    file_name = feature_files;
    feature_files = feature_files->mpNext;
    free(file_name);
  }
  
  for (unsigned int i = 0; i < cfg_hash.mNEntries; i++)
    free(cfg_hash.mpEntry[i]->data);

  my_hdestroy_r(&cfg_hash, 1);

  free(rhfbuff.mpLastFileName);

  return 0;
}
