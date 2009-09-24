
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

#define SVN_DATE       "$Date: 2009-03-27 09:52:26 +0100 (Fri, 27 Mar 2009) $"
#define SVN_AUTHOR     "$Author: glembek $"
#define SVN_REVISION   "$Revision: 276 $"
#define SVN_ID         "$Id: SERest.cc 276 2009-03-27 08:52:26Z glembek $"

#define MODULE_VERSION "2.0.7 "__TIME__" "__DATE__" "SVN_ID  


#include <STKLib/common.h>
#include <STKLib/stkstream.h>
#include <STKLib/BDTree.h>
#include <STKLib/MlfStream.h>
#include <STKLib/Tokenizer.h>



#include <algorithm>
#include <cctype>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <list>

#if !defined HAVE_BOOST
# warning "STree won't run---Requires BOOST"
int main(int argc, char* argv[]) {
  return 1;
}
#else

#include <boost/program_options.hpp>

// #define HAVE_HDF5
#ifdef HAVE_HDF5
# include <hdf5.h>
#endif

// ............................................................................
// shorten namespace to program options
namespace po = boost::program_options;

using namespace STK;
using namespace std;


// typedefs ...................................................................
typedef std::list<BDTree*> ModelList;


// prototypes .................................................................
void
LoadModels(const std::string& rModelListFile, ModelList& rModelList,
    BDTreeBuildTraits& rTraits);


void
AdaptModel(const BDTree& rOrig, const NGramSubset& rNGrams, 
    const std::string& rFileName, BDTreeBuildTraits& rTraits, 
    std::ostream* pOutVectorStream, BDTreeHeader file_header);

#ifdef HAVE_HDF5
hid_t
H5CreateDataset(hid_t h5_file, const std::string& rFileName,
    hid_t type_id, hid_t space_id, hid_t dcpl_id);
#endif

// A helper function to simplify the main part.
template<class T>
std::ostream& operator <<(std::ostream& os, const std::vector<T>& v)
{
  std::copy(v.begin(), v.end(), std::ostream_iterator<T>(std::cout, " ")); 
    return os;
}



/******************************************************************************/
/******************************************************************************/
int main(int argc, char* argv[])
{
  try 
  {
    // general-purpose variables
    bool                    htk_compatible = false;
    string                  config_file;
    string                  action;
    string                  predictor_vocab_fname;
    bool                    predictor_vocab_extended = false;
    VocabularyTable         predictor_vocabulary;
    string                  target_vocab_fname;
    VocabularyTable*        target_vocabulary = NULL;
    vector< string >        positional_parameters;
    BDTreeBuildTraits       new_tree_traits;
    string                  target_model_name;
    string                  score_output;
    string                  source_model_name;
    string                  source_model_list_name;
    string                  target_ARPA_fname;
    bool                    output_ARPA = false;
    int                     binary_format = 0;
    string                  heldout_data_fname;
    FLOAT                   counts_threshold = 0.0;
    NGramSubsets*           heldout_collection_pntr = NULL;
    //VocabularyTable*        word_vocabulary = NULL;
    //VocabularyTable*        stem_vocabulary = NULL;
    //VocabularyTable*        flex_vocabulary = NULL;
    //VocabularyTable*        POS_vocabulary = NULL;
    //VocabularyTable*        tag_vocabulary = NULL;

    InitLogMath();

    //srandom((unsigned int)time(NULL));
    srand((unsigned)time(0));

    // command line options only
    po::options_description generic_options("Command line options");
    generic_options.add_options()
      ("action,A",      po::value<string>(&action),                           "Action to perform")
      ("adaptmap",                                                            "Use MAP adaptation instead of original algorithm") 
      ("adaptr",        po::value<FLOAT>(&new_tree_traits.mAdaptR),           "Adapt r factor")
      ("cfgfile,C",     po::value<string>(&config_file),                      "Configuration file")
      ("countthreshold", po::value<FLOAT >(&counts_threshold),                "Ignore data whos value is less or equal to this threshold")
      ("help",                                                                "Print more detailed help") 
      ("includecounts",                                                       "In connection with outvectors and invector - augment each leaf with counts") 
      ("indexmlf",      po::value<bool>()->default_value(true),               "If true (default) perform indexing of MLF before reading") 
      ("inputmlf,I",    po::value<string>(),                                  "Source data MLF")
      ("invector",      po::value<string>(),                                  "File with supervector of leaf distribution (see --action=pushleaves)")
      ("maxdepth",      po::value<int>(&new_tree_traits.mMaxDepth),           "Maximum tree depth criterion")
      ("minentr",       po::value<FLOAT >(&new_tree_traits.mMinReduction),    "Minimum entropy reduction criterion")
      ("mindata",       po::value<FLOAT >(&new_tree_traits.mMinInData),       "Minimum input data criterion")
      ("MMIalpha",      po::value<FLOAT >(&new_tree_traits.mMMIAlpha),        "MMI alpha constant")
      ("MMIeta",        po::value<FLOAT >(&new_tree_traits.mMMIEta),          "MMI eta constant")
      ("MMIsplit",                                                            "Toggle MMI splitting criteria") 
      ("MMImininc",     po::value<FLOAT >(&new_tree_traits.mMinMMIIncrease),  "MMI max increase")
      ("scoreout",      po::value<string>(&score_output),                     "Score output file")
      ("order",         po::value<int>(&new_tree_traits.mOrder),              "Number of predictors")
      ("outfileversion", po::value<int>(),                                    "Source data script")
      ("outvectors",    po::value<string>(),                                  "Cluster vector output file")
      ("outvectorsh5",  po::value<string>(),                                  "Cluster vector output hdf5 file")
      ("predvoc",       po::value<string>(&predictor_vocab_fname),            "Predictor vocabulary file")
      ("predvocext",    po::value<bool>(&predictor_vocab_extended),           "Whether to use extended predictor vocab file")
      ("script,S",      po::value<string>(),                                  "Source data script")
      ("smoothr",       po::value<FLOAT >(&new_tree_traits.mSmoothR),         "Smooth r factor")
      ("srcmodel",      po::value<string>(&source_model_name),                "Source model file name")
      ("srcmdllst",     po::value<string>(&source_model_list_name),           "Source model list file name")
      ("tgtmdl",        po::value<string>(&target_model_name),                "Target model file name")
      ("tgtvoc",        po::value<string>(&target_vocab_fname),               "Target vocabulary file [same as predictor]")
      ("LVCSR",         po::value<bool>(&new_tree_traits.mLVCSR),             "LVCSR mode: Sparce matrices are used with exchange algorithm")
      ("randtree",      po::value<bool>(&new_tree_traits.mRandomizeTree),     "Use random question set split and predictor position selection when growing the tree")
      ("randpred",      po::value<bool>(&new_tree_traits.mRandomizePredictors), "Use random predictor selection when growing the tree")
      ("randquest",     po::value<bool>(&new_tree_traits.mRandomizeQuestions),"Use random question set split selection when growing the tree")
      ("outputARPA",    po::value<bool>(&output_ARPA),                        "Output scored N-grams in ARPA format")
      ("outARPAfile",   po::value<string>(&target_ARPA_fname),                "Target ARPA format LM")
      ("heldoutdata",   po::value<string>(&heldout_data_fname),               "Heldout data that can be used as stopping criteria at the training phase")
      //("morphvocscript,S", po::value<string>(),                               "Morphological vocabularies script")
      ("morphpred",     po::value<int>(&new_tree_traits.mMorphologicalPredictors),   "Use morphological predictors (for word models)")
      ("binformat",     po::value<int>(&binary_format),                       "Binary format used to store decision trees"); 

    // config file options only
    po::options_description config_options("Configuration");

    // hidden options only
    po::options_description hidden_options("Hidden");

    // command line options will be combination of generic options, config 
    // file options and hidden options
    po::options_description cmdline_options;
    cmdline_options.add(generic_options).add(config_options).add(hidden_options);

    // config file options will be combination of config options and hidden 
    // options
    po::options_description config_file_options;
    config_file_options.add(config_options);

    // visible options
    po::options_description visible_options;
    visible_options.add(generic_options).add(config_options);

    // we can also process positional parameters
    po::positional_options_description p_options;
    if (htk_compatible) {
      p_options.add("tiedlist", 1);
    }
    p_options.add("inputfile", -1);

    po::variables_map var_map;
    po::store(po::command_line_parser(argc, argv)
        .options(cmdline_options).positional(p_options).run(), var_map);
    po::notify(var_map);

    
    // if no parameters given, dump the generic options description
    if (0 == var_map.size()) {
      std::cout << generic_options << std::endl;
      return 0;
    }
    // if help specified, dump all parameters
    else if (var_map.count("help")) {
      std::cout << visible_options << std::endl;
      return 0;
    }
                 
    // parse the config file
    if (!config_file.empty()) {
      IStkStream ifs(config_file.c_str());
      if (!ifs.good()) {
        // TODO: Error message
      }

      po::store(parse_config_file(ifs, config_file_options), var_map);
      ifs.close();
    }

    // END OF COMMAND LINE PARAMETER PARSING ...................................
    

    // SWITCH AMONG ACTIONS ....................................................
    // make the action flag lower-case
    std::transform(action.begin(), action.end(), action.begin(), 
        (int(*)(int)) std::tolower);

    //..........................................................................
    // TRAIN
    if (action == "train") {
      //BDTree* p_newtree;

      // check the params
      if (target_model_name.empty()) {
        throw runtime_error("Target model (--tgtmdl) not specified");
      }

      // load the predictor vocabulary
      if (var_map.count("predvoc")) {
        predictor_vocabulary.LoadFromFile(predictor_vocab_fname,
            predictor_vocab_extended);
        // target vocab will be the same by default
        target_vocabulary = &predictor_vocabulary;
      }
      else {
        throw runtime_error("Predictor vocabulary file not set");
      }

      // if target vocabulary specified, then read it
      if (var_map.count("tgtvoc")) {
        target_vocabulary = new VocabularyTable();
        target_vocabulary->LoadFromFile(target_vocab_fname,
            predictor_vocab_extended);
      }

      NGramSubsets heldout_collection;
      NGramPool heldout_pool(new_tree_traits.mOrder);
      // if file with heldout data is specified, then read it
      if (var_map.count("heldoutdata")) 
      {
        // Create the stats file
        heldout_pool.setPredictorVocab(&predictor_vocabulary);
        heldout_pool.setTargetVocab(target_vocabulary);
        heldout_pool.AddFromFile(heldout_data_fname, 1, counts_threshold);

        std::cout << "Heldout data mass: " << heldout_pool.Mass() << std::endl;
        std::cout << "Heldout token count: " << heldout_pool.TokenCount() << std::endl;

        heldout_collection.push_back(heldout_pool);
	heldout_collection_pntr = &heldout_collection;
      }

      NGramSubsets data_collection;
      std::map<std::string,NGramPool> data_pools;
      std::map<std::string,NGramPool>::iterator i_data_pools;

      new_tree_traits.mMMItoggle = var_map.count("MMIsplit") ? true : false;

      // load the data .........................................................
      //
      // source streams
      IStkStream input_mlf;
      IMlfStream mlf_data_stream(input_mlf);
      IStkStream reg_data_stream;

      // this will be a general stream to read from
      std::istream* p_in_data_stream;  

      // initialize MLF if necessary
      if (var_map.count("inputmlf")) {
        std::string input_mlf_file = var_map["inputmlf"].as<string >();

        input_mlf.open(input_mlf_file.c_str());
        if (! input_mlf.good()) {
          throw runtime_error(string("Error opening master label file ") +
              input_mlf_file);
        }

        // index the MLF if desired
        if (var_map["indexmlf"].as<bool>()) {
          std::cout << "Indexing MLF " << input_mlf_file << " ... " << std::flush;
          mlf_data_stream.Index();
          std::cout << "DONE" << std::endl;
        }
      }


      // retreive the script file name
      string script_file = var_map["script"].as<string >();

      // parse script file
      IStkStream list_stream(script_file.c_str());
      if (!list_stream.good()) {
        throw runtime_error(string("Error opening source list file ") +
            script_file);
      }
      list_stream >> std::ws;

      // browse the list
      string line_buf;
      while (!list_stream.eof()) {
        // read list record
        std::getline(list_stream, line_buf);
        list_stream >> std::ws;

        // parse the logical and physical part and the weights
        FileListElem new_record(line_buf);

        if ((i_data_pools = data_pools.find(new_record.Logical())) == data_pools.end()) {
          i_data_pools = data_pools.insert(data_pools.end(), make_pair(new_record.Logical(), 
                NGramPool(new_tree_traits.mOrder)));

          // Create the stats file
          i_data_pools->second.setPredictorVocab(&predictor_vocabulary);
          i_data_pools->second.setTargetVocab(target_vocabulary);
        }

        // message
        std::cout << "To logical " << new_record.Logical() << " adding physical "
          << new_record.Physical() << " with weight " << new_record.Weight() << std::endl;


        if (var_map.count("inputmlf")) {
          // we'll be adding from stream if MLF specified
          mlf_data_stream.Open(new_record.Physical());
          if (!mlf_data_stream.good()) {
            throw runtime_error(string("Error opening label file ") +
                new_record.Physical());
          }
          p_in_data_stream = &mlf_data_stream;
        }
        else {
          // we'll be adding from stream if MLF specified
          reg_data_stream.open(new_record.Physical().c_str());
          if (!reg_data_stream.good()) {
            throw runtime_error(string("Error opening list file ") +
                new_record.Physical());
          }
          p_in_data_stream = &reg_data_stream;
        }

        // add the data
        i_data_pools->second.AddFromStream(*p_in_data_stream, new_record.Weight(), counts_threshold);

        if (var_map.count("inputmlf")) {
          // Need to close for proper work
          mlf_data_stream.Close();
        }
        else {
          reg_data_stream.close();
        }

        // write some info
        std::cout << "Data mass: " << i_data_pools->second.Mass() << std::endl;
        std::cout << "Token count: " << i_data_pools->second.TokenCount() << std::endl;
      }

      list_stream.close();


      for (i_data_pools=data_pools.begin(); i_data_pools != data_pools.end(); ++i_data_pools) {
        std::cout << "Creating collection for " << i_data_pools->first << std::endl;
        data_collection.push_back(i_data_pools->second);
      }

      // train new tree
      BDTree new_tree(data_collection, heldout_collection_pntr, new_tree_traits, NULL, "");
      
      OStkStream model_stream(target_model_name.c_str());
      if (!model_stream.good()) {
        throw runtime_error(string("Error opening output file ") + 
            target_model_name);
      }

      // fill the header
      BDTreeHeader file_header;
      file_header.mBinary = true;
      file_header.mFileVersion = binary_format;
      file_header.mOrder  = new_tree_traits.mOrder;
      file_header.mPredictorVocabSize = predictor_vocabulary.Size();
      file_header.mVocabSize = target_vocabulary->Size();

      // output header
      if(file_header.mFileVersion)
        file_header.Write_bin1(model_stream);
      else
        file_header.Write(model_stream);

      // output tree
      new_tree.Write(model_stream, file_header);
      // close stream
      model_stream.close();


      BDTreeInfo tree_info;
      new_tree.GetInfo(tree_info, NULL);

      std::cout << "Total nodes:          " << tree_info.mTotalNodes << std::endl;
      std::cout << "Total leaves:         " << tree_info.mTotalLeaves << std::endl;
      std::cout << "Tree average entropy: " << tree_info.mAverageEntropy << std::endl;
    } // action == train

    //..........................................................................
    // SCORE
    else if (action == "score") {
      ModelList           model_list;
      ModelList::iterator i_model;
      std::string         line_buf;

      // load models
      LoadModels(source_model_list_name, model_list, new_tree_traits);

      // load the predictor vocabulary
      if (var_map.count("predvoc")) {
        predictor_vocabulary.LoadFromFile(predictor_vocab_fname,
            predictor_vocab_extended);
        // target vocab will be the same by default
        target_vocabulary = &predictor_vocabulary;
      }
      else {
        throw runtime_error("Predictor vocabulary file not set");
      }

      // if target vocabulary specified, then read it
      if (var_map.count("tgtvoc")) {
        target_vocabulary = new VocabularyTable();
        target_vocabulary->LoadFromFile(target_vocab_fname,
            predictor_vocab_extended);
      }

      // prepare data pool
      NGramPool data(new_tree_traits.mOrder);
      data.setPredictorVocab(&predictor_vocabulary);
      data.setTargetVocab(target_vocabulary);

      // load the data .........................................................
      //
      // source streams
      IStkStream input_mlf;
      IMlfStream mlf_data_stream(input_mlf);
      IStkStream reg_data_stream;

      // this will be a general stream to read from
      std::istream* p_in_data_stream;  

      // initialize MLF if necessary
      if (var_map.count("inputmlf")) {
        std::string input_mlf_file = var_map["inputmlf"].as<string >();

        input_mlf.open(input_mlf_file.c_str());
        if (! input_mlf.good()) {
          throw runtime_error(string("Error opening master label file ") +
              input_mlf_file);
        }

        // index the MLF if desired
        if (var_map["indexmlf"].as<bool>()) {
          std::cout << "Indexing MLF " << input_mlf_file << "...";
          mlf_data_stream.Index();
          std::cout << " DONE" << std::endl;
        }
      }


      // load the data 
      string script_file = var_map["script"].as<string >();

      // parse script file
      IStkStream list_stream(script_file.c_str());
      if (!list_stream.good()) {
        throw runtime_error(string("Error opening source list file ") +
            script_file.c_str());
      }
      list_stream >> std::ws;

      OStkStream score_stream(score_output.c_str());
      if (!score_stream.good()) {
        throw runtime_error(string("Error opening score output file ") +
            score_output.c_str());
      }

      // browse the list
      while (!list_stream.eof()) {
        std::getline(list_stream, line_buf);
        list_stream >> std::ws;

        FileListElem new_record(line_buf);

        std::cout << "Scoring " << new_record.Physical() << " as " 
          << new_record.Logical() << std::endl;
        
        if (var_map.count("inputmlf")) {
          // we'll be adding from stream if MLF specified
          mlf_data_stream.Open(new_record.Physical());
          if (!mlf_data_stream.good()) {
            throw runtime_error(string("Error opening label file ") +
                new_record.Physical());
          }
          p_in_data_stream = &mlf_data_stream;
        }
        else {
          // we'll be adding from stream if MLF specified
          reg_data_stream.open(new_record.Physical().c_str());
          if (!reg_data_stream.good()) {
            throw runtime_error(string("Error opening list file ") +
                new_record.Physical());
          }
          p_in_data_stream = &reg_data_stream;
        }

        // parse data file
        data.AddFromStream(*p_in_data_stream, new_record.Weight(), counts_threshold);

        if (var_map.count("inputmlf")) {
          // Need to close for proper work
          mlf_data_stream.Close();
        }
        else {
          reg_data_stream.close();
        }

        // write intro
        // file name and number of frames
        score_stream << new_record.Logical() << " " << data.Mass();

        for (i_model=model_list.begin(); i_model!=model_list.end(); ++i_model) {
	  if(output_ARPA)
	  {
	    OStkStream lm_stream(target_ARPA_fname.c_str());
	    if (!lm_stream.good()) 
	    {
	      throw runtime_error(string("Error opening output file ") + 
                target_ARPA_fname);
	    }
            (*i_model)->OutputNGramSubsetARPA(data, lm_stream, new_tree_traits.mOrder);

	    lm_stream.close();
	  }
	  else
	  {
	    double e = (*i_model)->ScoreNGramSubset(data);
            score_stream << " " << e;
	  }
        }
        // end line
        score_stream << std::endl;

        data.Clear();
      }

      score_stream.close();
      list_stream.close();
    } // action == score

    //..........................................................................
    // VIEW
    else if (action == "view") {
      // open stream
      IStkStream model_stream(source_model_name.c_str());
      if (!model_stream.good()) {
        throw runtime_error(string("Error opening model file ") + 
            source_model_name);
      }

      // load header
      BDTreeHeader file_header;
      file_header.Read(model_stream);

      // load tree
      BDTree new_tree(model_stream, file_header);

      // recompute non-leaf distributions
      new_tree.RecomputeDists();

      if (new_tree_traits.mSmoothR > 0.0) {
        new_tree.ComputeBackoffDists(NULL, new_tree_traits.mSmoothR);
      }

      if (var_map.count("outvectors")) {
        std::string output_file = var_map["outvectors"].as<string >();
        OStkStream out_vector_stream(output_file.c_str());

        if (! out_vector_stream.good()) {
          throw runtime_error(string("Error opening output vector file") +
              output_file);
        }

        // this vector will be filled with stacked leaf dists
        BasicVector<FLOAT> super_vector;
        // stack the leaves to the supervector
        new_tree.FillLeafSupervector(super_vector, new_tree_traits.mSmoothR > 0.0, var_map.count("includecounts"));
        // write to stream
        out_vector_stream << super_vector << std::endl;
        out_vector_stream.close();
      }

      // dump tree
      new_tree.Dump(std::cout, " ");

      BDTreeInfo tree_info;
      new_tree.GetInfo(tree_info, NULL);

      std::cout << "Total nodes:          " << tree_info.mTotalNodes << std::endl;
      std::cout << "Total leaves:         " << tree_info.mTotalLeaves << std::endl;
      std::cout << "Tree average entropy: " << tree_info.mAverageEntropy << std::endl;
    }

    //..........................................................................
    // ADAPT
    else if (action == "adapt") {
      // if target not specified, logical script name will be used as output 
      // model name
      bool target_specified = true;

      if (source_model_name.empty()) {
        throw runtime_error("Source model not specified");
      }

      if (target_model_name.empty()) {
        if (!var_map.count("outvectors")) {
          throw runtime_error("Neither target model nor --outvectors specified");
        }
        target_specified = false;
      }

      // check whether to use the MAP adaptation
      new_tree_traits.mMapAdapt = var_map.count("adaptmap") ? true : false;

      // load the predictor vocabulary
      if (var_map.count("predvoc")) {
        predictor_vocabulary.LoadFromFile(predictor_vocab_fname,
            predictor_vocab_extended);
        // target vocab will be the same by default
        target_vocabulary = &predictor_vocabulary;
      }
      else {
        throw runtime_error("Predictor vocabulary file not set");
      }

      // if target vocabulary specified, then read it
      if (var_map.count("tgtvoc")) {
        target_vocabulary = new VocabularyTable();
        target_vocabulary->LoadFromFile(target_vocab_fname,
            predictor_vocab_extended);
      }

      // open stream
      IStkStream model_stream(source_model_name.c_str());
      if (!model_stream.good()) {
        throw runtime_error(string("Error opening model file ") + 
            source_model_name);
      }

      // load header
      BDTreeHeader file_header;
      file_header.Read(model_stream);
      new_tree_traits.mOrder = file_header.mOrder;

      // load tree
      BDTree new_tree(model_stream, file_header);
      model_stream.close();

      // recompute non-leaf distributions
      new_tree.RecomputeDists();

      // compute backoffs
      new_tree.ComputeBackoffDists(NULL, new_tree_traits.mSmoothR);


      // load data
      NGramSubsets          data_collection;
      std::list<NGramPool>  data_pools;

      data_pools.push_back(NGramPool(new_tree_traits.mOrder));
      // Create the stats file
      data_pools.back().setPredictorVocab(&predictor_vocabulary);
      data_pools.back().setTargetVocab(target_vocabulary);
      //data_pools.back().setWordVocab(word_vocabulary);
      //data_pools.back().setStemVocab(stem_vocabulary);
      //data_pools.back().setFlexVocab(flex_vocabulary);
      //data_pools.back().setPOSVocab(POS_vocabulary);
      //data_pools.back().setTagVocab(tag_vocabulary);


      // load the data .........................................................
      //
      // source streams
      IStkStream input_mlf;
      IMlfStream mlf_data_stream(input_mlf);
      IStkStream reg_data_stream;

      // this will be a general stream to read from
      std::istream* p_in_data_stream;  

      // initialize MLF if necessary
      if (var_map.count("inputmlf")) {
        std::string input_mlf_file = var_map["inputmlf"].as<string >();

        input_mlf.open(input_mlf_file.c_str());
        if (! input_mlf.good()) {
          throw runtime_error(string("Error opening master label file ") +
              input_mlf_file);
        }

        // index the MLF if desired
        if (var_map["indexmlf"].as<bool>()) {
          std::cout << "Indexing MLF " << input_mlf_file << "...";
          mlf_data_stream.Index();
          std::cout << " DONE" << std::endl;
        }
      }


      // load the data 
      string script_file = var_map["script"].as<string >();

      // parse script file
      IStkStream list_stream(script_file.c_str());
      if (!list_stream.good()) {
        throw runtime_error(string("Error opening source list file ") +
            script_file);
      }
      list_stream >> std::ws;

      // prepare outvector stream if desired
      OStkStream* p_out_vector_stream = NULL;
      if (var_map.count("outvectors")) {
        std::string output_file = var_map["outvectors"].as<string >();
        p_out_vector_stream = new OStkStream(output_file.c_str());

        if (! p_out_vector_stream->good()) {
          throw runtime_error(string("Error opening output vector file") +
              output_file);
        }
      }

      // browse the list
      string line_buf;
      while (!list_stream.eof()) {
        std::getline(list_stream, line_buf);
        list_stream >> std::ws;

        // parse the logical and physical part and the weights
        FileListElem new_record(line_buf);

        if ((!target_specified) && new_record.Physical() == new_record.Logical()) {
          throw runtime_error(string("Target model has not been specified, but"
               " logical and pysical names are identical (") + script_file + ")");
        }

        std::cout << "To logical " << new_record.Logical() << " adding physical "
          << new_record.Physical() << " with weight " << new_record.Weight() 
          << std::endl;

        if (var_map.count("inputmlf")) {
          // we'll be adding from stream if MLF specified
          mlf_data_stream.Open(new_record.Physical());
          if (!mlf_data_stream.good()) {
            throw runtime_error(string("Error opening label file ") +
                new_record.Physical());
          }
          p_in_data_stream = &mlf_data_stream;
        }
        else {
          // we'll be adding from stream if MLF specified
          reg_data_stream.open(new_record.Physical().c_str());
          if (!reg_data_stream.good()) {
            throw runtime_error(string("Error opening list file ") +
                new_record.Physical());
          }
          p_in_data_stream = &reg_data_stream;
        }

        // parse the file
        data_pools.back().AddFromStream(*p_in_data_stream, new_record.Weight(), counts_threshold);

        if (var_map.count("inputmlf")) {
          mlf_data_stream.Close();
        }
        else {
          reg_data_stream.close();
        }

        // if target was not specified, we immediately adapt
        // TODO: not very good, maybe we prefer to collect the names first, 
        // merge same logical names and then adapt - would be more general
        if (!target_specified) {
          ::AdaptModel(new_tree, data_pools.back(), new_record.Logical(), 
              new_tree_traits, p_out_vector_stream, file_header);

          data_pools.back().Clear();
        }
      }

      list_stream.close();

      if (target_specified) {
        data_collection.push_back(data_pools.back());

        if (var_map.count("outfileversion")) {
          file_header.mFileVersion = binary_format;
        }
        else if (binary_format) {
          file_header.mFileVersion = binary_format;
        }
        else if (file_header.mFileVersion != 0 
              && file_header.mFileVersion != 1) {
          file_header.mFileVersion = 0;
        }

        ::AdaptModel(new_tree, data_collection[0], target_model_name, 
            new_tree_traits, p_out_vector_stream, file_header);
      }

      if (NULL != p_out_vector_stream) {
        p_out_vector_stream->close();
        delete p_out_vector_stream;
      }
    }

    //..........................................................................
    // PUSH new values to the leaves from an ascii vector
    else if (action == "pushleaves") {
      if (source_model_name.empty()) {
        throw runtime_error("Source model not specified");
      }

      if (target_model_name.empty()) {
        throw runtime_error("Target model not specified");
      }

      // open stream
      IStkStream model_stream(source_model_name.c_str());
      if (!model_stream.good()) {
        throw runtime_error(string("Error opening model file ") + 
            source_model_name);
      }

      // load header
      BDTreeHeader file_header;
      file_header.Read(model_stream);
      new_tree_traits.mOrder = file_header.mOrder;

      // load tree
      BDTree new_tree(model_stream, file_header);
      model_stream.close();

      BDTreeInfo tree_info;
      new_tree.GetInfo(tree_info, NULL);

      std::cout << "Vocabulary size:      " << file_header.mVocabSize << std::endl;
      std::cout << "Total nodes:          " << tree_info.mTotalNodes << std::endl;
      std::cout << "Total leaves:         " << tree_info.mTotalLeaves << std::endl;
      std::cout << "Tree average entropy: " << tree_info.mAverageEntropy << std::endl;

      // load leaf supervector .................................................
      std::vector<FLOAT> super_vector(tree_info.mTotalLeaves * file_header.mVocabSize);
      //super_vector.reserve();

      if (!var_map.count("invector")) {
        throw runtime_error("--invector was not specified");
      }

      // open the stream
      std::string super_vector_file = var_map["invector"].as<string >();
      IStkStream super_vector_stream(super_vector_file.c_str());
      if (!super_vector_stream.good()) {
        throw runtime_error(string("Error opening supervector file ") + 
            super_vector_file);
      }
      super_vector_stream >> std::ws;

      int i = 0;
      while (!super_vector_stream.eof()) {
        double aux_value;
        super_vector_stream >> aux_value >> std::ws;
        super_vector[i] = aux_value;
        ++i;
        std::cout << i << " - " << aux_value << std::flush << std::endl;
      }

      super_vector_stream.close();
      std::cout << "Finished reading supervector" << std::flush << std::endl;

      // MAIN THING HERE
      // do the leaf update 
      new_tree.PushLeafSupervector(super_vector, false, 
          var_map.count("includecounts"));

      if (var_map.count("outfileversion")) {
        file_header.mFileVersion = binary_format;
      }
      else if (binary_format) {
        file_header.mFileVersion = binary_format;
      }
      else if (file_header.mFileVersion != 0 
            && file_header.mFileVersion != 1) {
        file_header.mFileVersion = 0;
      }

      // write model 
      OStkStream out_model_stream(target_model_name.c_str());
      if (!out_model_stream.good()) {
        throw runtime_error(string("Error opening target model file ") + 
            target_model_name);
      }

      // output header
      if(file_header.mFileVersion == 1) {
        file_header.Write_bin1(out_model_stream);
      }
      else {
        file_header.Write(out_model_stream);
      }

      // write the model and close stream
      new_tree.Write(out_model_stream, file_header);
      out_model_stream.close();
    }

    //..........................................................................
    // CLUSTER
    else if (action == "cluster") {
      if (source_model_name.empty()) {
        throw runtime_error("Source model not specified");
      }

      // load the predictor vocabulary
      if (var_map.count("predvoc")) {
        predictor_vocabulary.LoadFromFile(predictor_vocab_fname,
            predictor_vocab_extended);
        // target vocab will be the same by default
        target_vocabulary = &predictor_vocabulary;
      }
      else {
        throw runtime_error("Predictor vocabulary file not set");
      }

      if (!var_map.count("outvectors") && !var_map.count("outvectorsh5")) {
        throw runtime_error("--outvectors expected");
      }

      // if target vocabulary specified, then read it
      if (var_map.count("tgtvoc")) {
        target_vocabulary = new VocabularyTable();
        target_vocabulary->LoadFromFile(target_vocab_fname,
            predictor_vocab_extended);
      }


      // load model ........................................
      IStkStream model_stream(source_model_name.c_str());
      if (!model_stream.good()) {
        throw runtime_error(string("Error opening model file ") + 
            source_model_name);
      }

      // load header
      BDTreeHeader file_header;
      file_header.Read(model_stream);
      new_tree_traits.mOrder = file_header.mOrder;

      // load tree
      BDTree new_tree(model_stream, file_header);
      model_stream.close();

      // this vector will be filled with stacked leaf dists
      BasicVector<FLOAT> super_vector;

      string script_file = var_map["script"].as<string >();

      if (var_map.count("outvectorsh5")) {
#ifdef HAVE_HDF5
        hsize_t               vec_size;
        hsize_t               dim_size;
        size_t                n_leaves;
        std::vector<BDTree*>  leaves;
        std::string           output_file = var_map["outvectorsh5"].as<string >();

        // compute sizes
        n_leaves = new_tree.GetLeaves(leaves);

        const VecDistribution* p_aux_distr = 
          dynamic_cast<const VecDistribution*>(leaves.front()->cpDistribution());

        if (NULL == p_aux_distr) {
          throw std::runtime_error("Cannot collect counts for non-vector distribution implementations");
        }

        dim_size    = p_aux_distr->Size();
        vec_size    = n_leaves * dim_size;

        // parse script file
        IStkStream list_stream(script_file.c_str());
        if (!list_stream.good()) {
          throw runtime_error(string("Error opening source list file ") +
              script_file);
        }


        hid_t       h5_file;
        hid_t       h5_dataset;            /* file and dataset handles */
        hid_t       h5_datatype;
        hid_t       h5_dataspace;   /* handles */
        hid_t       h5_plist;
        herr_t      h5_status;                             
        unsigned szip_options_mask;
        unsigned szip_pixels_per_block;

        /* Save old error handler */
        herr_t  (*p_h5_old_error_function)(void*);
        void*   p_h5_old_client_data;
        H5Eget_auto(&p_h5_old_error_function, &p_h5_old_client_data);

        /* Turn off error handling */
        H5Eset_auto(NULL, NULL);

        /*
         * Create a new file using H5F_ACC_TRUNC access,
         * default file creation properties, and default file
         * access properties.
         */
        h5_file = H5Fcreate(output_file.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);


        // browse the list
        string line_buf;
        while (!list_stream.eof()) {
          std::getline(list_stream, line_buf);
          list_stream >> std::ws;

          // parse the logical and physical part and the weights
          FileListElem new_record(line_buf);

          // parse the file
          NGramPool  data_pool(file_header.mOrder);
          data_pool.setPredictorVocab(&predictor_vocabulary);
          data_pool.setTargetVocab(target_vocabulary);
          data_pool.AddFromFile(new_record.Physical(), new_record.Weight(), counts_threshold);

          // clear previous counts 
          new_tree.ResetLeaves();

          // cluster the data
          new_tree.CollectCounts(data_pool);

          // stack the leaves to the supervector
          new_tree.FillLeafSupervector(super_vector, false, false);

          /*
           * Describe the size of the array and create the data space for fixed
           * size dataset. 
           */
          h5_dataspace = H5Screate_simple(1, &vec_size, NULL); 

          /* 
           * Define datatype for the data in the file.
           * We will store little endian INT numbers.
           */
          h5_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
          h5_status  = H5Tset_order(h5_datatype, H5T_ORDER_LE);

          /*
           * Define property list for dataset---we want to compress the files
           */
          h5_plist = H5Pcreate(H5P_DATASET_CREATE);
          H5Pset_chunk(h5_plist, 1, &dim_size);
          H5Pset_deflate(h5_plist, 6);

          /*
           * Create a new dataset within the file using defined dataspace and
           * datatype and default dataset creation properties.
           */
          // h5_dataset = H5Dcreate(h5_file, new_record.Logical().c_str(), 
          //     h5_datatype, h5_dataspace, h5_plist);

          h5_dataset = H5CreateDataset(h5_file, new_record.Logical().c_str(), 
              h5_datatype, h5_dataspace, h5_plist);

          // 
          if (h5_dataset < 0) {
            throw runtime_error("Could not create dataset");
          }
          
          /*
           * Write the data to the dataset using default transfer properties.
           */
          h5_status = H5Dwrite(h5_dataset, h5_datatype, H5S_ALL, H5S_ALL, 
              H5P_DEFAULT, super_vector.cpData());

          if (h5_status < 0) {
            throw runtime_error("Could not write dataset");
          }
          
          /*
           * Close/release resources.
           */
          H5Sclose(h5_dataspace);
          H5Tclose(h5_datatype);
          H5Pclose(h5_plist);
          H5Dclose(h5_dataset);
        }

        /*
         * Close/release resources.
         */
        H5Fclose(h5_file);

        /* Restore previous error handler */
        H5Eset_auto(p_h5_old_error_function, p_h5_old_client_data);
	      
        list_stream.close();
#else
        throw runtime_error("HDF5 support has not been enabled");
#endif
      }
      else {
        string output_file = var_map["outvectors"].as<string >();


        // parse script file
        IStkStream list_stream(script_file.c_str());
        if (!list_stream.good()) {
          throw runtime_error(string("Error opening source list file ") +
              script_file);
        }

        // open output file
        OStkStream vector_stream(output_file.c_str());
        if (!vector_stream.good()) {
          throw runtime_error(string("Error opening output file ") +
              output_file);
        }

        // browse the list
        string line_buf;
        while (!list_stream.eof()) {
          std::getline(list_stream, line_buf);
          list_stream >> std::ws;

          // parse the logical and physical part and the weights
          FileListElem new_record(line_buf);

          // parse the file
          NGramPool  data_pool(file_header.mOrder);
          data_pool.setPredictorVocab(&predictor_vocabulary);
          data_pool.setTargetVocab(target_vocabulary);
          data_pool.AddFromFile(new_record.Physical(), new_record.Weight(), counts_threshold);

          // clear previous counts 
          new_tree.ResetLeaves();

          // cluster the data
          new_tree.CollectCounts(data_pool);

          // stack the leaves to the supervector
          new_tree.FillLeafSupervector(super_vector, false, false);

          // write to stream
          vector_stream << new_record.Logical() << " " << super_vector << std::endl;
        }

        list_stream.close();
        vector_stream.close();
      }
    }

    //..........................................................................
    // CRAP
    else {
      throw std::runtime_error("Wrong action specified");
    }

    // CLEAN UP ................................................................
    if (NULL != target_vocabulary && target_vocabulary != &predictor_vocabulary) {
      delete target_vocabulary;
    }
  }
  catch (std::exception& rExc) {
    std::cerr << "Exception thrown" << std::endl;
    std::cerr << rExc.what() << std::endl;
    return 1;
  }

  return 0;
}



/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void
LoadModels(const std::string& rModelListFile, ModelList& rModelList,
    BDTreeBuildTraits& rTraits)
{
  IStkStream    list_stream;
  std::string   line_buf;
  BDTree*       new_tree;
  int           model_count = 0;

  if (rModelListFile.empty()) {
    throw runtime_error("Source model list not specified");
  }

  list_stream.open(rModelListFile.c_str());
  if (!list_stream.good()) {
    throw runtime_error(string("Error opening model list file ") +
        rModelListFile);
  }

  while (!list_stream.eof()) {
    std::getline(list_stream, line_buf);
    list_stream >> std::ws;

    // open stream
    IStkStream model_stream(line_buf.c_str());
    if (!model_stream.good()) {
      throw runtime_error(string("Error opening model file ") + 
          line_buf);
    }

    // load header
    BDTreeHeader file_header;
    file_header.Read(model_stream);

    // for first model remember N-gram order, then only check
    if (model_count == 0) {
      rTraits.mOrder = file_header.mOrder;
    }
    else if (rTraits.mOrder != file_header.mOrder) {
      throw runtime_error(string("Order does not match previous model(s) ") + 
          rModelListFile);
    }

    // load tree
    new_tree = new BDTree(model_stream, file_header);

    // smooth probs if desired
    if (rTraits.mSmoothR > 0.0) {
      // recompute non-leaf distributions
      new_tree->RecomputeDists();
      new_tree->ComputeBackoffDists(NULL, rTraits.mSmoothR);
    }

    // add to list
    rModelList.push_back(new_tree);

    model_stream.close();
  }

  if (list_stream.fail()) {
    throw runtime_error(string("Error reading model list file ") +
        rModelListFile);
  }

  list_stream.close();
}


/******************************************************************************/
/******************************************************************************/
void
AdaptModel(const BDTree& rOrig, const NGramSubset& rNGrams, 
    const std::string& rFileName, BDTreeBuildTraits& rTraits, 
    std::ostream* pOutVectorStream, BDTreeHeader file_header)
{
  BDTree new_tree(rOrig);

  // this step clears the mpDist distributions. they will be filled with data
  // and adapted using mpBackoffDist, which should be now ready, as 
  // ComputeBackoffDists must have been already called
  new_tree.ResetLeaves();

  // adapt using new data
  new_tree.AdaptLeafDists(rNGrams, rTraits);

  BDTreeInfo tree_info;
  new_tree.GetInfo(tree_info, NULL);

  std::cout << "Total nodes:          " << tree_info.mTotalNodes << std::endl;
  std::cout << "Total leaves:         " << tree_info.mTotalLeaves << std::endl;
  std::cout << "Tree average entropy: " << tree_info.mAverageEntropy << std::endl;

  if (NULL != pOutVectorStream) {
    // this vector will be filled with stacked leaf dists
    BasicVector<FLOAT> super_vector;
    // stack the leaves to the supervector
    new_tree.FillLeafSupervector(super_vector, false, false);
    // write to stream
    *pOutVectorStream << rFileName << " " << super_vector << std::endl;
  }
  else {
    // write model 
    OStkStream out_model_stream(rFileName.c_str());
    if (!out_model_stream.good()) {
      throw runtime_error(string("Error opening target model file ") + 
          rFileName);
    }

    // output header
    if(file_header.mFileVersion == 1)
      file_header.Write_bin1(out_model_stream);
    else
      file_header.Write(out_model_stream);

    new_tree.Write(out_model_stream, file_header);

    out_model_stream.close();
  }
}

#ifdef HAVE_HDF5
hid_t 
H5RCreateDataset(hid_t h5_location, const char* pFileName,
    hid_t type_id, hid_t space_id, hid_t dcpl_id)
{
  hid_t h5_dataset;

  // check wheather slashes exist in the middle somewhere
  const char* p_last_pos = std::strchr(pFileName, '/');

  // make sure that the current char is not slash
  assert(p_last_pos != pFileName);

  // if no slash is found, we can directly create the dataset,
  // otherwise we'll create the parent groups
  if (p_last_pos == NULL) {
    h5_dataset = H5Dcreate(h5_location, pFileName,
      type_id, space_id, dcpl_id);
  }
  else {
    size_t  len = p_last_pos - pFileName;
    char*   p_group_name = new char[len+1];
    hid_t   h5_grp;

    // extract the current group name
    std::strncpy(p_group_name, pFileName, len);
    p_group_name[len] = '\0';

    // try to open the group (within the current location)
    h5_grp = H5Gopen(h5_location, p_group_name);

    // if we didn't manage to open the group, we'll create it and we'll set
    // the current location accordingly
    if (h5_grp < 0) {
      std::cout << "Creating group " << p_group_name << std::endl;
      h5_location = H5Gcreate(h5_location, p_group_name, -1);
    }
    else {
      h5_location = h5_grp;
    }

    delete [] p_group_name;

    h5_dataset = H5RCreateDataset(h5_location, p_last_pos+1, type_id,
        space_id, dcpl_id);
    
    H5Gclose(h5_location);
  }

  return h5_dataset;
}

hid_t
H5CreateDataset(hid_t h5_file, const std::string& rFileName,
    hid_t type_id, hid_t space_id, hid_t dcpl_id)
{
  const char* p_file_name = rFileName.c_str();

  while (*p_file_name == '/' || isspace(*p_file_name)) {
    p_file_name++;
  }

  return H5RCreateDataset(h5_file, p_file_name, type_id, space_id, dcpl_id);
}
#endif

#endif //boost 

