# PBS submission script for running tasks in parallel
# declare a name for this job to be sample_job
#PBS -N lamodel  
# request 1 node
#PBS -l nodes=1:ppn=1
#PBS -l cput=9999:00:00
#PBS -l walltime=1:00:00
#$ -cwd 

# By default, PBS scripts execute in your home directory, not the 
# directory from which they were submitted. The following line 
# places you in the directory from which the job was submitted.  

#cd $PBS_O_WORKDIR
# run the program


./lamodel $LAPARAMS


