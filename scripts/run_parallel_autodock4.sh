#!/bin/bash


num_processors_for_docking=108
sge_name_for_script=ad42_mpi
sge_name_for_queue=AS

let overheard_time=0
let docking_time=0
let total_ligands=0
let log_indicator=$RANDOM
let counter=1



if [ ! -d logs ]; then
	mkdir logs
fi


ls ready-for-docking > undocked_list.txt

let total_files=$(cat undocked_list.txt | wc -l)

for i in $(cat undocked_list.txt)
do

	total_ligands=$(($total_ligands + $(zcat ready-for-docking/$i | tar -tf - | grep dpf | wc -l)))
	echo "Processing ligands up to $total_ligands compounds"

	echo "Processing $counter of $total_files"
	counter=$((counter+1))
	T="$(date +%s)"
	basename=$(basename $i .tar.gz)
	if [ ! -d Dockings ]; then mkdir Dockings; fi
	if [ ! -d docking_results ]; then mkdir docking_results; fi
	tar -xzf ready-for-docking/$i -C Dockings

	ls Dockings/	> ligand_list.txt

	echo '#!/bin/bash' > cluster_sub_script.sh
	echo "module load intel" >> cluster_sub_script.sh
	echo "module load openmpi_intel" >> cluster_sub_script.sh
	echo 'mpirun autodock-intel-serial-mpi $PWD/ligand_list.txt $PWD/Dockings/ $PWD/logs/ default_seed reuse_maps' >> cluster_sub_script.sh

	T="$(($(date +%s)-T))"
	overhead_time=$((overhead_time+T))

	T="$(date +%s)"
	qsub -V -N $sge_name_for_script -pe orte $num_processors_for_docking -cwd -j y -o /nfsscratch/autodocklogs-${log_indicator}.log -q $sge_name_for_queue cluster_sub_script.sh

	
	echo "Running jobs on cluster..."
	job_waiter=$(qstat | grep $sge_name_for_script)
	while [ "$job_waiter" ]; do sleep 20; job_waiter=$(qstat | grep $sge_name_for_script); done
	
	echo "Finished!"
	T="$(($(date +%s)-T))"
	docking_time="$((docking_time+T))"


	T="$(date +%s)"
	echo "Packing docking logs..."
	mkdir tmp_dlglinks

	for j in $(find $PWD/Dockings -name "*.dlg")
	do
		ln -s $j tmp_dlglinks
	done

	cd tmp_dlglinks
	tar -czhf ${basename}-results.tar.gz * --remove-files
	mv ${basename}-results.tar.gz ../docking_results
	cd ../
	rm -rf tmp_dlglinks

	echo "Finished!"

	echo "Cleaning up docking directory..."

	rm -rf Dockings/*

	echo "Finished!"

	T="$(($(date +%s)-T))"
	overhead_time="$((overhead_time+T))"

	echo $i >> finished_docking_list.txt
	sed -i '1d' undocked_list.txt

done	

total_time="$((docking_time+overhead_time))"
	
echo "Docked $total_ligands ligands in a total of $((total_time/60)) minutes" >> timings.txt

echo "Overhead time = $((overhead_time/60)) minutes" >> timings.txt
echo "Docking time = $((docking_time/60)) minutes" >> timings.txt

	






