#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>



#include "statistics_functions.h"

 
void parse_generation(std::vector<std::string> generation_array,
                 std::vector<float>* tmp_energy_levels,
                 std::vector<int>* tmp_bins) {


  int dock_counter=0;
  int generation_counter=0;
  int bin_value;
  int old_bin_value=0;
  int max_bin_value=0;
  float energy_value;
  float min_energy_value;
  std::vector<float> unnormalized_tmp_energy_levels;

  for (int i=0; i < generation_array.size(); i++) {

    //std::cout << generation_array[i].substr(11,generation_array[i].find("Oldest")-11) << " ";
    //std::cout << generation_array[i].substr(generation_array[i].find("Oldest")+16,7) << std::endl;
    
    energy_value=atof(generation_array[i].substr(generation_array[i].find("Oldest")+16,7).c_str());
    bin_value=atoi(generation_array[i].substr(11,generation_array[i].find("Oldest")-11).c_str());
    //std:: cout << unnormalized_tmp_energy_levels[0] << " ";

    if (bin_value > max_bin_value) {
    // Initialize the bins
      (*tmp_bins).push_back(bin_value); 
      (*tmp_energy_levels).push_back(0.0);
      unnormalized_tmp_energy_levels.push_back(energy_value);
      max_bin_value = bin_value;
      if (energy_value < min_energy_value) {
        min_energy_value=energy_value;
      }
      old_bin_value=bin_value;

    }

    else if (bin_value < old_bin_value) {
    // We start over adding stuff to the bins
    // First normalize what we have
      generation_counter=0;     
      for (int j=0; j<unnormalized_tmp_energy_levels.size(); j++) {
        (*tmp_energy_levels)[j]+=(unnormalized_tmp_energy_levels[j]/min_energy_value-(*tmp_energy_levels)[j])/(dock_counter+1);
      }
      //std::cout << unnormalized_tmp_energy_levels[0] << std::endl;
      //std::cout << unnormalized_tmp_energy_levels[0]/min_energy_value << std::endl;
      //std::cout << (*tmp_energy_levels)[0] << std::endl;
      dock_counter++;
      unnormalized_tmp_energy_levels.clear();     
      min_energy_value = energy_value;
      
      unnormalized_tmp_energy_levels.push_back(energy_value);
      //std::cout << unnormalized_tmp_energy_levels[0] << std::endl;
      generation_counter++;

      old_bin_value=bin_value;
      
    }
    
    else if (bin_value > old_bin_value) {


      unnormalized_tmp_energy_levels.push_back(energy_value);
      generation_counter++;
      if (energy_value < min_energy_value) {
        min_energy_value=energy_value;
      }   

      old_bin_value=bin_value;
     

    }


  } // End For loop through generation data


  // Do a final addition of values to tmp energy levels

  for (int j=0; j<unnormalized_tmp_energy_levels.size(); j++) {
      (*tmp_energy_levels)[j]+=(unnormalized_tmp_energy_levels[j]/min_energy_value-(*tmp_energy_levels)[j])/(dock_counter+1);

  }



}



std::vector<std::string> parse_dlg_for_generation(std::string InputFile) {

  std::ifstream t( InputFile.c_str(), std::ifstream::in);
  
  std::string Line;
  
  std::vector<std::string> generation_string;

  while (getline(t,Line)) {

    if ( Line.find("Generation:") != std::string::npos) {
      std::stringstream buffer(Line);
      generation_string.push_back(buffer.str());
      buffer.flush();
    }

  }


  return generation_string;
    
   
}






void run_convergence_statistics( std::vector<std::string> InputFileList, 
                            std::vector<float>* energy_levels, 
                            std::vector<int>* bins) {



 

  for (int i=0; i<InputFileList.size(); i++) {


    std::vector<std::string> generation_array;
    std::vector<float> tmp_energy_levels;
    std::vector<int> tmp_bins;

    generation_array=parse_dlg_for_generation(InputFileList[i]);

    parse_generation(generation_array,&tmp_energy_levels,&tmp_bins);
    
    if (tmp_bins.size() > (*bins).size()) {
      (*bins) = tmp_bins;
    }
    
    for (int j=0; j<tmp_energy_levels.size(); j++) {
  
      if ( j > (*energy_levels).size() || i==0) {
        (*energy_levels).push_back(tmp_energy_levels[j]);  
      }
      else {
        (*energy_levels)[j]+=tmp_energy_levels[j];
      }

    }
  }
}
