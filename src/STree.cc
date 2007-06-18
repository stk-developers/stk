#include <STKLib/stkstream.h>
#include <STKLib/BDTree.h>

#include <algorithm>
#include <boost/program_options.hpp>
#include <cctype>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>


// shorten namespace to program options
namespace po = boost::program_options;

using namespace STK;
using namespace std;

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
    VocabularyTable*        target_vocabulary;
    vector< string >        positional_parameters;


    // command line options only
    po::options_description generic_options("Command line options");
    generic_options.add_options()
      ("help",          "Print more detailed help") 
      ("action,A",      po::value<string>(&action),         "Action to perform")
      ("configfile,C",  po::value<string>(&config_file),    "Configuration file")
      ("predictorvocab", po::value<string>(&predictor_vocab_fname),    "Predictor vocabulary file")
      ("targetvocab",   po::value<string>(&target_vocab_fname),    "Target vocabulary file [same as predictor]")
      ("order",         po::value<int>(),                        "Number of predictors")
      ("script,S",      po::value<vector<string> >(),  "Script files")
      ;

    // config file options only
    po::options_description config_options("Configuration");
    config_options.add_options()
      ("sourcehmm,H",       po::value<string>(), "Source hmm") 
      ;

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
    
    // load the predictor vocabulary
    if (var_map.count("predictorvocab")) {
      predictor_vocabulary.LoadFromFile(predictor_vocab_fname);
      // target vocab will be the same by default
      target_vocabulary = &predictor_vocabulary;
    }
    else {
      throw runtime_error("Predictor vocabulary file not set");
    }

    // if target vocabulary specified, then read it
    if (var_map.count("targetvocab")) {
      target_vocabulary = new VocabularyTable();
      target_vocabulary->LoadFromFile(target_vocab_fname);
    }


    // SWITCH AMONG ACTIONS ....................................................
    // make the action flag lower-case
    std::transform(action.begin(), action.end(), action.begin(), 
        (int(*)(int)) std::tolower);

    if (action == "train") {
      BDTree* p_newtree;
      // Create the stats file
      NGramPool data(3);
      data.setPredictorVocab(&predictor_vocabulary);
      data.setTargetVocab(target_vocabulary);

      // load the data 
      vector<string>::const_iterator it_file;
      const vector<string>& script = 
        var_map["script"].as<vector<string> >();

      // parse all script files
      for (it_file=script.begin(); it_file != script.end(); ++it_file) {
        std::cout << "Parsing script file: " << *it_file << std::endl;

        IStkStream list_stream(it_file->c_str());
        if (!list_stream.good()) {
          throw runtime_error(string("Error opening script file ") +
              *it_file);
        }
        list_stream >> std::ws;

        // browse the list
        string line_buf;
        while (!list_stream.eof()) {
          std::getline(list_stream, line_buf);
          list_stream >> std::ws;

          // parse the file
          data.AddFromFile(line_buf, 1.0);
        }


        list_stream.close();
      }


    } // action == train

    else if (action == "score") {
    } // action == score

    else {
      throw std::runtime_error("Wrong action specified");
    }

    // CLEAN UP ................................................................
    if (target_vocabulary != &predictor_vocabulary) {
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


