



void parse_generation(std::vector<std::string> generation_array,
                      std::vector<float>* tmp_energy_levels,
                      std::vector<int>* tmp_bins);

std::vector<std::string> parse_dlg_for_generation(std::string InputFile);

void run_convergence_statistics( std::vector<std::string> InputFileList, 
                                std::vector<float>* energy_levels, 
                                std::vector<int>* bins);
