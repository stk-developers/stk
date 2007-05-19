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
  // general-purpose variables
  bool                    htk_compatible = false;
  std::string             config_file;
  std::string             action;
  std::vector< std::string > 
                          positional_parameters;

  try {
    // command line options only
    po::options_description generic_options("Command line options");
    generic_options.add_options()
      ("help", "Print this help message") 
      ("action,A",     po::value<std::string>(&action),         "Action to perform")
      ("configfile,C", po::value<std::string>(&config_file),    "Configuration file")
      ("order",        po::value<int>(),                        "Number of predictors")
      ("script,S",     po::value<std::vector<std::string> >(),  "Script files")
      ;

    // config file options only
    po::options_description config_options("Configuration");
    config_options.add_options()
      ("sourcehmm,H",       po::value<std::string>(), "Source hmm") 
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
    if (var_map.count("help")) {
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
    
    // make the action flag lower-case
    std::transform(action.begin(), action.end(), action.begin(), 
        (int(*)(int)) std::tolower);

    // SWITCH AMONG ACTIONS ....................................................
    if (action == "train") {
      BDTree* p_newtree;

    } // action == train

    else if (action == "score") {
    } // action == score

    else {
      throw std::runtime_error("Wrong action specified");
    }
  }
  catch (std::exception& rExc) {
    std::cerr << rExc.what() << std::endl;
    return 1;
  }

  return 0;
}

