#!/bin/bash

last_job_id=
num_jobs=3

# Parse the flags and arguments

while getopts "n:j:" opt; do
  case $opt in
    n) num_jobs=$OPTARG      ;;  
    j) last_job_id=$OPTARG   ;;  
    *) echo 'error' >&2
       exit 1
  esac
done


# Loop over the number of jobs
for i in $(seq 1 $num_jobs); do
  # Define the job script
  script_file="run.sbatch"
  
# Submit the job to the queue

  if [ -z "$last_job_id" ]; then

    # if a job id is not passed, the first job is queued to start immediately
    if [ $i == 1 ] ; then 
        last_job_id=$(sbatch --parsable $script_file)
        echo $last_job_id
    else
        last_job_id=$(sbatch --parsable --dependency=afterany:$last_job_id $script_file)
        echo $last_job_id
    fi  

  else

    last_job_id=$(sbatch --parsable --dependency=afterany:$last_job_id $script_file)
    echo $last_job_id
  
  fi  
done

