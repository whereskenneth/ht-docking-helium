#!/bin/bash

if [ ! -d processed_results ]; then
	mkdir processed_results
fi


counter=1;
total=$(ls docking_results | wc -l)
for i in docking_results/*.tar.gz
do
  
	
  if [ ! -d docking_results/tmp ]
   then mkdir docking_results/tmp
   else 
		rm docking_results/tmp/*
		echo "Cleaning up..."
  fi

  
	
  
  echo "Converting $counter of $total docking results..."
  
  echo "Untarring docking results..."  
  tar -xzf $i -C docking_results/tmp

	echo "Converting dlgs to sdfs..."
  dlgtosdf -i docking_results/tmp/*.dlg >> docking_results/final_results.sdf

	mv $i processed_results

  counter=$((counter+1))
done

rm -rf docking_results/tmp

cp receptor/*.pdbqt docking_results
