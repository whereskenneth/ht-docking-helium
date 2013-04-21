#!/bin/bash


# Some control numbers for the script, does not need to be changed. Can be if 
# wanted, for optimization purposes

sge_name_for_script=prepdpf4
number_of_processors_anticipated=108
max_ligand_list_size=2000
min_ligand_list_size=1000
limit_num_ligands=100000
ga_num_evals=50000
ga_pop_size=150
ga_runs=1
sge_queue_custom=AS


T="$(date +%s)"
# Need to know how many ligands we will be dealing with

let global_num_ligands=0
let global_total_num_ligands=0

if [ $(($(find ligand-archives -name "*.tar.gz" | wc -l))) -gt 0 ]; then sleep 0
else 
	echo "No ligand archives found. Exiting"
	exit
fi



echo "Getting total ligand count..."
for i in ligand-archives/*.tar.gz
do
  
  global_num_ligands=$(zcat $i | tar -tf - | wc -l)
	global_total_num_ligands=$((global_total_num_ligands+global_num_ligands))
done
echo "Total ligand count: $global_total_num_ligands"




# Function that creates a list of ligand files, divied up by how many processors and compounds there are. 

function prepare_for_dpf4 () {

	if [ $((total_num_ligands/number_of_processors_anticipated)) -gt $max_ligand_list_size ]; then
		let ligand_list_size=$max_ligand_list_size
	elif [ $((total_num_ligands/number_of_processors_anticipated)) -lt $min_ligand_list_size ]; then
		let ligand_list_size=$min_ligand_list_size
	else
		let ligand_list_size=$((total_num_ligands/$number_of_processors_anticipated+1))
	fi

	cd ligands
	let counter_l=0
	let counter_j=1

	if [ -d ../ligand-list-for-dpf4 ]; then rm ../ligand-list-for-dpf4/*
	else
		mkdir ../ligand-list-for-dpf4
	fi


  # put each ligand file's name into a list of ~2K compounds,
  # depending on number of procs and number of ligs
	for i in *.pdbqt
	do
		echo $i >> ../ligand-list-for-dpf4/ligand-list-for-dpf4-${counter_l}
		if [ $(($counter_j % $ligand_list_size)) -eq 0 ]; then
			counter_l=$((counter_l+1))
		fi
		counter_j=$((counter_j+1))
	done
	
	cd ../

}


# This function creates a shell script that will be submitted to SGE to process the ligands into
# folders which are meaningful to the parallel autodock. For each list of ~2K ligands produced in  
# the function prepare_for_dpf4, a shell script is spawned that will iterate through the list. Each 
# of these shell scripts is submitted to whichever queue is specified in the control portion of this script.

function run_dpf4 () {


	
	if [ -d dpf4tmp ]; then rm dpf4tmp/*
	else
		mkdir dpf4tmp
	fi
	
  basename=$PWD

  if [ ! -d ready-for-docking ]; then mkdir ready-for-docking; fi
  if [ ! -d Dockings ]; then mkdir Dockings; fi
  if [ ! -d dpf4s ]; then mkdir dpf4s; fi

	cd ligand-list-for-dpf4
	for i in ligand-list-for-dpf4-*
	do
		touch ../dpf4tmp/$i.sh
		echo "#!/bin/bash" >> ../dpf4tmp/$i.sh
		echo "module load intel" >> ../dpf4tmp/$i.sh
		echo "module load python26" >> ../dpf4tmp/$i.sh
    echo 'random_folder=$RANDOM-$RANDOM' >> ../dpf4tmp/$i.sh
    echo 'mkdir Dockings/docking-${random_folder}' >> ../dpf4tmp/$i.sh
		echo -ne 'for j in $(cat ' >> ../dpf4tmp/$i.sh
		echo -ne "ligand-list-for-dpf4/$i" >> ../dpf4tmp/$i.sh
		echo ');' >> ../dpf4tmp/$i.sh
		echo 'do ' >> ../dpf4tmp/$i.sh
    echo '  ligand_name=$(basename $j .pdbqt)' >> ../dpf4tmp/$i.sh
    echo '  receptor_name=$(basename receptor/*.pdbqt .pdbqt)' >> ../dpf4tmp/$i.sh
		echo -ne '  prepare_dpf42.py -l ligands/${j} -r receptor/1zuw.pdbqt -o dpf4s/${ligand_name}_${receptor_name}.dpf' >> ../dpf4tmp/$i.sh
	  echo " -p ga_num_evals=$ga_num_evals -p ga_pop_size=$ga_pop_size -p ga_run=$ga_runs -p rmstol=2.0" >> ../dpf4tmp/$i.sh
    echo '  mkdir Dockings/docking-${random_folder}/${ligand_name}_${receptor_name}' >> ../dpf4tmp/$i.sh
    echo '  mv dpf4s/${ligand_name}_${receptor_name}.dpf Dockings/docking-${random_folder}/${ligand_name}_${receptor_name}'  >> ../dpf4tmp/$i.sh
    echo -ne "  ln -s $basename/" >> ../dpf4tmp/$i.sh
    echo 'receptor/* Dockings/docking-${random_folder}/${ligand_name}_${receptor_name}'  >> ../dpf4tmp/$i.sh
    echo '  mv ligands/${ligand_name}.pdbqt Dockings/docking-${random_folder}/${ligand_name}_${receptor_name}'  >> ../dpf4tmp/$i.sh
		echo 'done' >> ../dpf4tmp/$i.sh
		echo 'cd Dockings/docking-${random_folder}' >> ../dpf4tmp/$i.sh
    echo 'tar -czf ../../ready-for-docking/docking_folder-${random_folder}.tar.gz * --remove-files' >> ../dpf4tmp/$i.sh
		echo 'cd ../../' >> ../dpf4tmp/$i.sh
    echo 'rm -rf Dockings/docking-${random_folder}' >> ../dpf4tmp/$i.sh
    
	done
	
	cd ../

  
  # Submit each shell script to the cluster as an independent, serial job. When the lustre file system becomes 
  # available again, we will change the above script to output to lustre (or create a symlink from the folder
  # dpf4s to some folder on lustre.
	for i in dpf4tmp/*.sh
	do
	  
		qsub -V -N $sge_name_for_script -cwd -j y -o ~/scratch/$sge_name_for_script.log -q $sge_queue_custom $i
	done 


	# This line queries the queue, looks for the job named $sge_name_for_script, 
  # and looks every 20 seconds until there are no  
  # more jobs with then name $sge_name_for_script. The function then exits.
	echo "Running jobs on cluster..."
	sleep 10
	job_waiter=$(qstat | grep $sge_name_for_script)
	while [ "$job_waiter" ]; do sleep 20; job_waiter=$(qstat | grep $sge_name_for_script); done
	
	echo "Finished!"

}





let number_of_ligands=0
let total_number_of_ligands=0
let total_docking_archives=1
let ligand_counter_for_ouput=0

ligand_archive_list_notifier=$RANDOM


if [ -d ligands ]; then
  rm -rf ligands/*
else
  mkdir ligands
fi


for i in ligand-archives/*.tar.gz
do
	echo $i >> ligand_archive_list-${ligand_archive_list_notifier}.txt
done

# Typically, each archive will house 100K compounds in pdbqt format. Typically, we have around 50 processors 
# free at any given time. This means that the best way to divy up each archive is to split it into subsets
# containing ~2K ligands apiece. The time it takes to process 2K compounds is around 10 minutes. 

exec 3<&0
exec 0<ligand_archive_list-${ligand_archive_list_notifier}.txt

while read line
do


	# We will fill our ligand directory up to 100K ligands, but no more so that file access isn't slowed down
 	if [ $total_number_of_ligands -lt $limit_num_ligands ]; then
		echo "Extracting $line..."
		tar -xzf $line -C ligands
		number_of_ligands=$(zcat $line  | tar -tf - | wc -l)
		ligand_counter_for_output=$((ligand_counter_for_output+number_of_ligands))
		total_number_of_ligands=$((total_number_of_ligands+number_of_ligands))
    #sleep 0
	else 
		echo "Preparing dpf4s for $ligand_counter_for_output out of $global_total_num_ligands ligands"
		prepare_for_dpf4
		run_dpf4
		echo "Extracting $line..."
		tar -xzf $line -C ligands
		number_of_ligands=$(zcat $line  | tar -tf - | wc -l)
		ligand_counter_for_output=$((ligand_counter_for_output+number_of_ligands))
		total_number_of_ligands=$number_of_ligands
    total_docking_archives=$((total_docking_archives+1))
 	fi
	

done


# Finally, perform the last run through of ligands for processing.
prepare_for_dpf4
run_dpf4


exec 0<&3 

#echo "ligands/${LIGAND_FILENAME}"
#echo "dpf4s/$(basename $LIGAND_FILENAME .pdbqt)"


if [ ! -d logs ]; then mkdir logs; fi
rm -rf dpf4tmp
rm -rf dpf4s
rm -rf ligand_archive_list-${ligand_archive_list_notifier}.txt
rm -rf ligand-list-for-dpf4

T="$(($(date +%s)-T))"

echo "Processed $global_total_num_ligands ligands" > timings.txt
echo "Timing in minutes to process: $((T/60))" >> timings.txt


