#include <iostream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <stdio.h>
#include <time.h>

#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/descriptor.h>
#include <openbabel/generic.h>
#include <boost/program_options.hpp>

#include "WriteToSDF.h"
namespace po=boost::program_options;



int main(int argc, char* argv[]) {


  po::options_description desc("Options:");
	desc.add_options()
		("help,h", "Produce this message")
		("input,i", po::value< std::vector<std::string> >()->multitoken(), "Input dlg files to convert to sdf")
    ("output,o", po::value< std::string >(), "Output filename for sdf file")
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

    std::cerr << "No input file given! Quitting..." << std::endl;
    return 1;
  }

  if (vm.count("output")) {
  
    OutputFile=vm["output"].as<std::string>();

  }

  else {

    std::cerr << "No output file given. Redirecting to standard output" << std::endl;

  }



  
    


  if (writesdf(InputFileList,OutputFile)) {
    return 0;
  }



  else {
  
    return 1;
  }


} 












