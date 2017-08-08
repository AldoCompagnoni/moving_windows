#!/bin/bash 

#$  -N  array_test
#$  -S  /bin/bash
#$  -l  h_rt=86400
#$  -l  h_vmem=8G
#$  -cwd

# IMPORTAN: number of parallel environment needs to be SAME as binding linear number.
# first of the following two arguments specifies the number of cores to use.
#$  -pe smp 8
#$  -binding linear:8

#$ -o /work/$USER/$JOB_NAME-$JOB_ID-$TASK_ID.log
#$ -e /work/$USER/$JOB_NAME-$JOB_ID-$TASK_ID.err

module load R/3.3.1-1

Rscript array_test.R \
	"$SGE_TASK_ID" \
	"/data/idiv_knight/sApropos/" \
	"/work/$USER/$JOB_NAME-$JOB_ID-$SGE_TASK_ID"
