#include <STKLib/common.h>
#include <STKLib/stkstream.h>
#include <STKLib/BDTree.h>
#include <STKLib/MlfStream.h>
#include <STKLib/Tokenizer.h>
#include "STKLib/Features.h"


#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <list>

#if !defined HAVE_BOOST
# error Requires BOOST
#else

#include <boost/program_options.hpp>

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
    bool                    swap_features;
    int                     start_frm_ext;
    int                     end_frm_ext;
    string                  target_kind;

    FeatureRepository       feature_repo;
    string                  out_stat_file;
    string                  script_file;
    string                  script_file_filter;
    Matrix<FLOAT>           feature_matrix;

    double                  zero_stats(0.0);
    BasicVector<FLOAT>      first_stats;
    Matrix<FLOAT>           second_stats;

    string                  config_file;
    vector< string >        positional_parameters;

    InitLogMath();

    //srandom((unsigned int)time(NULL));
    srand((unsigned)time(0));


    // command line options only
    po::options_description generic_options("Command line options");
    generic_options.add_options()
      ("cfgfile,C",     po::value<string>(&config_file),                      "Configuration file")
      ("help",                                                                "Print more detailed help") 

      ("outstatfile",   po::value<string>(),                                  "Output the stats here")
      ("script,S",      po::value<string>(),                                  "Source data script")
      ("scriptfilter",  po::value<string>()->default_value(""),               "Script file filter")
      ("targetkind",    po::value<string>(&target_kind)->default_value(""),   "Target kind code");


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

    po::variables_map var_map;
    po::store(po::command_line_parser(argc, argv)
        .options(cmdline_options).run(), var_map);
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

    // SOME CHECKING ...........................................................
    if (!var_map.count("script")) {
      throw runtime_error("Script file not specified");
    }

    if (!var_map.count("outstatfile")) {
      throw runtime_error("Output file not specified");
    }

    out_stat_file = var_map["outstatfile"].as<string>();
    script_file   = var_map["script"].as<string>();


    // END OF COMMAND LINE PARAMETER PARSING ...................................
    int target_kind_code = ReadParmKind(target_kind.c_str(), false);

    // swap_features, start_frm_ext, end_frm_ext, targetKind, deriv_order, 
    // deriv_win_lengths, cmn_path, cmn_mask, cvn_path, cvn_mask, cvg_file
    feature_repo.Init(!isBigEndian(), 0, 0, target_kind_code,
        0, NULL, NULL, NULL, NULL, NULL, NULL);

    // retreive the script file name
    script_file = var_map["script"].as<string >();

    feature_repo.AddFileList(script_file.c_str(), script_file_filter == "" ?
        NULL : script_file_filter.c_str());

    BasicVector<double> aux_vec;

    for (feature_repo.Rewind(); !feature_repo.EndOfList(); 
        feature_repo.MoveNext()) {
      std::cout << feature_repo.Current().Logical() << " = " <<
        feature_repo.Current().Physical() << std::endl;
      
      // now reest features are read ...........................................
      feature_repo.ReadFullMatrix(feature_matrix);

      if (0 == first_stats.Length()) {
        first_stats.Init(feature_matrix.Cols());
        second_stats.Init(feature_matrix.Cols(), feature_matrix.Cols());
      }
      // check the dimensionality of the current feature file
      else if (first_stats.Length() != feature_matrix.Cols()) {
        throw runtime_error(std::string("Dimensionality of the feature file ")
            + feature_repo.Current().Physical() + " differs from the previous");
      }


      // adding a size_t
      zero_stats += feature_matrix.Rows();
      first_stats.AddColSum(feature_matrix);
      second_stats.AddCMtMMul(1.0, feature_matrix, feature_matrix);

      feature_matrix.Destroy();
    }

    OStkStream out_stream(out_stat_file.c_str());
    if (!out_stream.good()) {
      throw runtime_error(string("Error opening output file ") + 
          out_stat_file);
    }

    out_stream << 3 << " ";
    out_stream << 1 << " ";
    out_stream << first_stats.Length() << " " ;
    out_stream << second_stats.Rows() << " " ;
    out_stream << 0 << std::endl; // write extra zero

    out_stream << zero_stats << std::endl << first_stats << std::endl << second_stats;

    out_stream.close();
  }
  catch (std::exception& rExc) {
    std::cerr << "Exception thrown" << std::endl;
    std::cerr << rExc.what() << std::endl;
    return 1;
  }

  return 0;
}



#endif //boost 

