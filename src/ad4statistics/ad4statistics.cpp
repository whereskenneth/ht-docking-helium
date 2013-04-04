#include <iostream>
#include <vector>
#include <string>



#include <boost/program_options.hpp>
namespace po=boost::program_options;


#include "statistics_functions.h"
#include "write.h"

int main(int argc, char* argv[]) {

  po::options_description desc("Options:");

  desc.add_options()
    ("help, h","Produce this message")
    ("input,i",po::value< std::vector<std::string> >()->multitoken(), "Input dlg files to run calculation statistics")
    ("output,o",po::value< std::string >(), "Output ascii file containing docking convergence statistics")
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,desc), vm);
  po::notify(vm);

  std::vector<std::string> InputFileList;
  std::string OutputFile;


  if (vm.count("input")) {

    InputFileList = vm["input"].as<std::vector<std::string> >();

  }
  else {

    std::cerr << "No input files given! Quitting..." << std::endl;
    return 1;

  }
  
  if (vm.count("output")) {

    OutputFile = vm["output"].as<std::string>();

  }

  else {

    std::cerr << "No output file given. Redirecting to standard out." << std::endl;

  }




  std::vector<float> energy_levels;
  std::vector<int> bins;
  

  run_convergence_statistics(InputFileList, &energy_levels, &bins);

  //write_to_output(OutputFile, energy_levels, bins);






  return 0;












}


  
