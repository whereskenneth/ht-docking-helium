#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>



void write_to_output(std::string OutputFile, std::vector<float> energy_levels, std::vector<int> bins) {



  std::string Line;
  std::vector<std::string> previous_data;

  std::vector<float> old_energies;
  std::vector<int> old_bins;

  float tmp_energy;
  int tmp_bin;


  ifstream ifile(OutputFile.c_str());
  if (!ifile) {

    ofstream ofile(OutputFile.c_str());

  }
  else {

    while (!ifile.eof()) {
      
      ifile >> tmp_bin >> tmp_energy;
      old_bins.push_back(tmp_bin);
      old_energies.push_back(tmp_energy);

    }

    get_old_data(previous_data,&old_energies,&old_bins);

  }    


}

