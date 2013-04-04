#include <iostream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <stdio.h>
#include <time.h>
#include <fstream>
#include <sstream>

#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/descriptor.h>
#include <openbabel/generic.h>
#include "WriteToSDF.h"



void process_and_convert(std::vector<std::vector<std::string> > pdbqt_file_array, std::vector<std::string> SDTags, std::vector<OpenBabel::OBMol>* mol_array, std::string OutputFilename, std::string receptor_name) {

  
  // SDTags :
  // 1 - Cluster 
  // 2 - conformer 
  // 3- binding energy 
  // 4-Ki


  std::string cluster;
  std::string num_conformers;
  std::string binding_energy;
  std::string Ki;
  float Ki_numeric;
  std::string Ki_units;
  std::string Ki_final;
  std::string mol_title;



  OpenBabel::OBConversion obconversion;
  OpenBabel::OBMol mol;

  obconversion.SetOutFormat("SDF");
  obconversion.SetInFormat("PDBQT");


  //std::cout << pdbqt_file_array.size() << std::endl;
  for (int i=0; i<pdbqt_file_array.size(); i++) {



    std::stringstream pdbqt_stream;


    for (int j=0; j<pdbqt_file_array[i].size(); j++) {
     

      if (pdbqt_file_array[i][j].find("Cluster Rank") != std::string::npos) {
        int cluster_int = atoi(pdbqt_file_array[i][j].substr(23,1).c_str());
        std::stringstream cluster_stream; 
        cluster_stream << cluster_int;
        cluster=cluster_stream.str();
      }

      if (pdbqt_file_array[i][j].find("Number of conformations in this cluster") != std::string::npos) {
        int num_conformers_numeric = atoi(pdbqt_file_array[i][j].substr(49,2).c_str());
        std::stringstream conf_stream;
        conf_stream << num_conformers_numeric;
        num_conformers=conf_stream.str();
        

      }

      if (pdbqt_file_array[i][j].find("Free Energy of Binding") != std::string::npos) {
  
        float binding_energy_fl = atof(pdbqt_file_array[i][j].substr(46,8).c_str());
        std::stringstream binding_energy_ss;
        binding_energy_ss << binding_energy_fl;
        binding_energy=binding_energy_ss.str();
      }

      if (pdbqt_file_array[i][j].find("Estimated Inhibition Constant") != std::string::npos) {

        std::stringstream Ki_final_string;
        Ki=pdbqt_file_array[i][j].substr(46,10);
        Ki_units=Ki.substr(8,2);
        Ki_numeric=atof(Ki.substr(0,7).c_str());

        if (Ki_units.find("mM") != std::string::npos) {
          Ki_numeric*=1000;
        }
        if (Ki_units.find("nM") != std::string::npos) {
          Ki_numeric/=1000; 
        }

        Ki_final_string << Ki_numeric;
        Ki_final = Ki_final_string.str();

      }


      if (pdbqt_file_array[i][j].find("USER    NEWDPF move") != std::string::npos) {
        std::string title;
        title = pdbqt_file_array[i][j].substr(20,50);
        mol_title=title.erase(title.find("."));
        title.clear();
      }

      //std::cout << pdbqt_file_array[i][j] << std::endl;
      if ( pdbqt_file_array[i][j].find("USER") == std::string::npos) {
        pdbqt_stream << pdbqt_file_array[i][j] << std::endl;
      }
      


    }
    //std::cout << pdbqt_stream.str() << std::endl;

    obconversion.SetInStream(&pdbqt_stream); 

    obconversion.Read(&mol);
    std::stringstream output;    
    std::ofstream ofs(OutputFilename.c_str());

    if (OutputFilename.size()) {
      obconversion.SetOutStream(&ofs);
    }  
    else {  
      obconversion.SetOutStream(&output);  
    }


    mol.SetTitle(mol_title);

    OpenBabel::OBPairData *ClusterNum = new OpenBabel::OBPairData;
    ClusterNum->SetAttribute(SDTags[0]);
    ClusterNum->SetValue(cluster);
    mol.SetData(ClusterNum);
  
    OpenBabel::OBPairData *NumConformers = new OpenBabel::OBPairData;
    NumConformers->SetAttribute(SDTags[1]);
    NumConformers->SetValue(num_conformers);
    mol.SetData(NumConformers);

    OpenBabel::OBPairData *BindingEnergy = new OpenBabel::OBPairData;
    BindingEnergy->SetAttribute(SDTags[2] + " (kcal/mol");
    BindingEnergy->SetValue(binding_energy);
    mol.SetData(BindingEnergy);
  
    OpenBabel::OBPairData *Affinity = new OpenBabel::OBPairData;
    Affinity->SetAttribute(SDTags[3]);
    Affinity->SetValue(Ki_final);
    mol.SetData(Affinity);
    
    OpenBabel::OBPairData *Receptor = new OpenBabel::OBPairData;
    Receptor->SetAttribute("Receptor");
    Receptor->SetValue(receptor_name);
    mol.SetData(Receptor);




    mol.DeleteData("MODEL");
    

    obconversion.Write(&mol);  
    std::cout << output.str(); 
   
    std::cout.flush();
    pdbqt_stream.flush();
  }
}









int parse_dlg(std::vector<std::string> DLG_file, std::vector<std::vector<std::string> >* mol_vector_in_pdbqt, std::string* receptor_name){ 

 
  std::vector<std::string> pdbqt_temp_file;
  std::vector<std::string> atom_type_array;

  int counter =0;
  
  int atom_identifier;
  for (int i=0; i < DLG_file.size(); i++) {


    if (DLG_file[i].find("DPF> fld") != std::string::npos) {


      *receptor_name = DLG_file[i].erase(0,DLG_file[i].find("fld"));
      *receptor_name = (*receptor_name).erase((*receptor_name).find("."));
      *receptor_name = (*receptor_name).erase(0,(*receptor_name).find(" ")+1);

      //std::cout << *receptor_name << std::endl;

    }

        
    
    bool found_docked=false;
    if (DLG_file[i].find("DOCKED") != std::string::npos) {
      found_docked=true;
    }
    
    if (DLG_file[i].find("INPUT-LIGAND-PDBQT") != std::string::npos) {
      if (DLG_file[i].find("ATOM") != std::string::npos) {


        found_docked=true;
        //std::cout << DLG_file[i].substr(96,3) << std::endl;
        atom_type_array.push_back(DLG_file[i].substr(96,3));
        
      }
    }


    
    if ( (DLG_file[i].find("MODEL") != std::string::npos || DLG_file[i].find("ATOM") != std::string::npos || DLG_file[i].find("TER\n") != std::string::npos || DLG_file[i].find("ENDMDL") != std::string::npos || DLG_file[i].find("USER") != std::string::npos) && !found_docked) {

      if (DLG_file[i].find("ATOM") != std::string::npos) {

        atom_identifier=atoi(DLG_file[i].substr(8,3).c_str());


        pdbqt_temp_file.push_back(DLG_file[i].substr(0,76) + atom_type_array[atom_identifier-1]);
      }
      else {
        
        pdbqt_temp_file.push_back(DLG_file[i]);
      
      }
         
      //std::cout << pdbqt_temp_file[counter] << std::endl;
      //counter++;
      //std::cout << DLG_file[i] << std::endl;
      if (DLG_file[i].find("ENDMDL") != std::string::npos) {
        (*mol_vector_in_pdbqt).push_back(pdbqt_temp_file);
        pdbqt_temp_file.clear();
      }

    } //end if actual pdbqt file
  
  }// end for

  return 0;


}

  






std::vector<std::string> read_input_dlg(std::string InputFile) {

  std::ifstream t( InputFile.c_str(), std::ifstream::in);
  std::string Line;
  std::vector<std::string> dlg_string_file;

  while (getline(t,Line))
  {
    std::stringstream buffer(Line);
    dlg_string_file.push_back(buffer.str());
    buffer.flush();
  }



  return dlg_string_file;
}







int writesdf(std::vector<std::string> InputFileList, std::string OutputFile) {


  OpenBabel::OBConversion obconversion;
  std::vector<OpenBabel::OBMol> mol_array;
  std::vector<std::vector<std::string> > mol_vector_in_pdbqt;







  

  std::string receptor_name;
  std::vector<std::string> SDTags;
  std::vector<std::string> SDData;
  //Each line from DLG file will be stored as a string in a vector element
  std::vector<std::string> DLG_file;

  SDTags.push_back("ClusterRank");
  SDTags.push_back("Conformers_In_Cluster");
  SDTags.push_back("Binding_Energy");  
  SDTags.push_back("Estimated_Ki (uM)");

  for (int i=0; i<InputFileList.size(); i++) {

    DLG_file=read_input_dlg(InputFileList[i]);
    
    parse_dlg(DLG_file, &mol_vector_in_pdbqt, &receptor_name);
    process_and_convert(mol_vector_in_pdbqt, SDTags, &mol_array,OutputFile,receptor_name);

    mol_vector_in_pdbqt.clear();
  }

  return 1;

 
}
