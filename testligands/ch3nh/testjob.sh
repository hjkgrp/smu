#! -S bin/bash
#$ -l h_rt=24:00:00 # specify maximum runtime
#$ -l h_rss=8G	     # specify requested memory
#$ -q cpus	     # specify which queue to submit the job to
#$ -l cpus=1	     # Always set to 1 no matter how many gpus are used.	

sleep 30
echo "hello"
