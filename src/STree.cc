#include <STKLib/stkstream.h>
#include <STKLib/BDTree.h>

#include <algorithm>
#include <boost/program_options.hpp>
#include <cctype>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <list>


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
    VocabularyTable         predictor_vocabulary;
    string                  target_vocab_fname;
    VocabularyTable*        target_vocabulary = NULL;
    vector< string >        positional_parameters;
    BDTreeBuildTraits       new_tree_traits;
    string                  target_model_name;
    string                  score_output;
    string                  source_model_name;
    string                  source_model_list_name;


    // command line options only
    po::options_description generic_options("Command line options");
    generic_options.add_options()
      ("action,A",      po::value<string>(&action),                           "Action to perform")
      ("adaptr",        po::value<double>(&new_tree_traits.mAdaptR),          "Adapt r factor")
      ("cfgfile,C",     po::value<string>(&config_file),                      "Configuration file")
      ("help",                                                                "Print more detailed help") 
      ("maxdepth",      po::value<int>(&new_tree_traits.mMaxDepth),           "Maximum tree depth criterion")
      ("minentr",       po::value<double>(&new_tree_traits.mMinReduction),    "Minimum entropy reduction criterion")
      ("mindata",       po::value<double>(&new_tree_traits.mMinInData),       "Minimum input data criterion")
      ("MMIsplit",                                                            "Toggle MMI splitting criteria") 
      ("multiscript",   po::value<string>(),                                  "MMI data collection")
      ("MMImininc",     po::value<double>(&new_tree_traits.mMinMMIIncrease),  "MMI max increase")
      ("multilang",                                                           "Toggle multilang splitting criteria") 
      ("scoreout",      po::value<string>(&score_output),                     "Score output file")
      ("order",         po::value<int>(&new_tree_traits.mOrder),              "Number of predictors")
      ("predvoc",       po::value<string>(&predictor_vocab_fname),            "Predictor vocabulary file")
      ("script,S",      po::value<string>(),                                  "Source data script")
      ("smoothr",       po::value<double>(&new_tree_traits.mSmoothR),         "Smooth r factor")
      ("srcmodel",      po::value<string>(&source_model_name),                "Source model file name")
      ("srcmdllst",     po::value<string>(&source_model_list_name),           "Source model list file name")
      ("tgtmdl",        po::value<string>(&target_model_name),                "Target model file name")
      ("tgtvoc",        po::value<string>(&target_vocab_fname),               "Target vocabulary file [same as predictor]")
      ;

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
      BDTree* p_newtree;

      // check the params
      if (target_model_name.empty()) {
        throw runtime_error("Target model (--targetmodel) not specified");
      }

      // load the predictor vocabulary
      if (var_map.count("predvoc")) {
        predictor_vocabulary.LoadFromFile(predictor_vocab_fname);
        // target vocab will be the same by default
        target_vocabulary = &predictor_vocabulary;
      }
      else {
        throw runtime_error("Predictor vocabulary file not set");
      }

      // if target vocabulary specified, then read it
      if (var_map.count("tgtvoc")) {
        target_vocabulary = new VocabularyTable();
        target_vocabulary->LoadFromFile(target_vocab_fname);
      }

      NGramSubsets data_collection;
      std::list<NGramPool> data_pools;
      new_tree_traits.mMMItoggle = var_map.count("MMIsplit") ? true : false;
      new_tree_traits.mMultiScript = var_map.count("multiscript") ? true : false;

      if (new_tree_traits.mMMItoggle || new_tree_traits.mMultiScript) {
        if (!var_map.count("multiscript")) {
          throw runtime_error(string("--multiscript not specified"));
        }

        std::string multi_script_file = var_map["multiscript"].as<string >();

        // open MMI script list
        IStkStream multi_script_stream(multi_script_file.c_str());
        if (!multi_script_stream.good()) {
          throw runtime_error(string("Error opening multiscript file ") +
              multi_script_file);
        }

        string line_buf;
        while (!multi_script_stream.eof()) {
          std::getline(multi_script_stream, line_buf);
          multi_script_stream >> std::ws;

          // create new pool
          data_pools.push_back(NGramPool(new_tree_traits.mOrder));
          // Create the stats file
          data_pools.back().setPredictorVocab(&predictor_vocabulary);
          data_pools.back().setTargetVocab(target_vocabulary);
          std::cout << "Adding " << line_buf << " to multiscript list" << std::endl;

          // parse script file
          IStkStream list_stream(line_buf.c_str());
          if (!list_stream.good()) {
            throw runtime_error(string("Error opening source list file ") +
                line_buf);
          }
          list_stream >> std::ws;

          // browse the list
          string script_line_buf;
          while (!list_stream.eof()) {
            std::getline(list_stream, script_line_buf);
            list_stream >> std::ws;

            // parse the file
            data_pools.back().AddFromFile(script_line_buf, 1.0);
          }

          list_stream.close();

          data_collection.push_back(data_pools.back());
        }

        multi_script_stream.close();
      }
      else {
        data_pools.push_back(NGramPool(new_tree_traits.mOrder));

        // Create the stats file
        data_pools.back().setPredictorVocab(&predictor_vocabulary);
        data_pools.back().setTargetVocab(target_vocabulary);

        // load the data 
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
          std::getline(list_stream, line_buf);
          list_stream >> std::ws;

          // parse the file
          data_pools.back().AddFromFile(line_buf, 1.0);
        }

        list_stream.close();

        data_collection.push_back(data_pools.back());
      }

      // train new tree
      BDTree new_tree(data_collection, new_tree_traits, NULL, "");
      
      OStkStream model_stream(target_model_name.c_str());
      if (!model_stream.good()) {
        throw runtime_error(string("Error opening output file ") + 
            target_model_name);
      }

      // fill the header
      BDTreeHeader file_header;
      file_header.mBinary = true;
      file_header.mOrder  = new_tree_traits.mOrder;
      file_header.mVocabSize = predictor_vocabulary.Size();

      // output header
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
        predictor_vocabulary.LoadFromFile(predictor_vocab_fname);
        // target vocab will be the same by default
        target_vocabulary = &predictor_vocabulary;
      }
      else {
        throw runtime_error("Predictor vocabulary file not set");
      }

      // if target vocabulary specified, then read it
      if (var_map.count("tgtvoc")) {
        target_vocabulary = new VocabularyTable();
        target_vocabulary->LoadFromFile(target_vocab_fname);
      }

      // prepare data pool
      NGramPool data(new_tree_traits.mOrder);
      data.setPredictorVocab(&predictor_vocabulary);
      data.setTargetVocab(target_vocabulary);

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

        // parse data file
        data.AddFromFile(line_buf, 1.0);

        // write intro
        // file name and number of frames
        score_stream << line_buf << " " << data.Mass();

        for (i_model=model_list.begin(); i_model!=model_list.end(); ++i_model) {
          double e = (*i_model)->ScoreNGramSubset(data);
          score_stream << " " << e;
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
      IStkStream model_stream(target_model_name.c_str());
      if (!model_stream.good()) {
        throw runtime_error(string("Error opening model file ") + 
            source_model_name);
      }

      // load header
      BDTreeHeader file_header;
      file_header.Read(model_stream);

      // load tree
      BDTree new_tree(model_stream, file_header);

      // smooth probs if desired
      if (new_tree_traits.mSmoothR > 0.0) {
        // recompute non-leaf distributions
        new_tree.RecomputeDists();
        new_tree.ComputeBackoffDists(NULL, new_tree_traits.mSmoothR);
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
      if (source_model_name.empty()) {
        throw runtime_error("Source model not specified");
      }

      if (target_model_name.empty()) {
        throw runtime_error("Target model not specified");
      }

      // load the predictor vocabulary
      if (var_map.count("predvoc")) {
        predictor_vocabulary.LoadFromFile(predictor_vocab_fname);
        // target vocab will be the same by default
        target_vocabulary = &predictor_vocabulary;
      }
      else {
        throw runtime_error("Predictor vocabulary file not set");
      }

      // if target vocabulary specified, then read it
      if (var_map.count("tgtvoc")) {
        target_vocabulary = new VocabularyTable();
        target_vocabulary->LoadFromFile(target_vocab_fname);
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

      // load data
      NGramSubsets          data_collection;
      std::list<NGramPool>  data_pools;

      data_pools.push_back(NGramPool(new_tree_traits.mOrder));
      // Create the stats file
      data_pools.back().setPredictorVocab(&predictor_vocabulary);
      data_pools.back().setTargetVocab(target_vocabulary);

      // load the data 
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
        std::getline(list_stream, line_buf);
        list_stream >> std::ws;

        // parse the file
        data_pools.back().AddFromFile(line_buf, 1.0);
      }

      list_stream.close();

      data_collection.push_back(data_pools.back());

      // recompute non-leaf distributions
      new_tree.RecomputeDists();

      // compute backoffs
      new_tree.ComputeBackoffDists(NULL, new_tree_traits.mSmoothR);

      // clear new distributions
      new_tree.ResetLeaves();

      // adapt using new data
      new_tree.AdaptLeafDists(data_collection[0], new_tree_traits.mAdaptR);

      // write model 
      OStkStream out_model_stream(target_model_name.c_str());
      if (!out_model_stream.good()) {
        throw runtime_error(string("Error opening target model file ") + 
            target_model_name);
      }

      file_header.Write(out_model_stream);
      new_tree.Write(out_model_stream, file_header);

      out_model_stream.close();
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
          rModelListFile);
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

