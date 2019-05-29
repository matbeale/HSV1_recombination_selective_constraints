#!/bin/bash


#usage : Extract_Pathseek_Stats.sh <samplename>


samplename=$1

# basefolder=/home/ubuntu/analyses/snippy-pipelined/pipeline-0.5_dedup/GoodCov_Genomes/dedup_Core
### basefolder=/home/ubuntu/analyses/snippy-pipelined/pipeline-0.5_dedup/lowcov_genomes
basefolder=/home/ubuntu/analyses/snippy-pipelined/pipeline-0.5_dedup/GoodCov_Genomes
runfolder=$basefolder/$samplename
workingfolder=$runfolder/coverage-stats

# Test if files actually exist
if [ ! -f $runfolder/$samplename.bam ]; then
echo "File $workingfolder/$samplename.depth does not exist"
echo "Usage: qsub VarScan.sh <samplename>"
exit 1
fi

mkdir -v -p $workingfolder/

samtools depth $runfolder/$samplename.bam > $workingfolder/$samplename.depth

# R script that determines coverage stats and also outputs a plot.  This version has modifications specific to the snippy reference pipeline that allow for missing regions due to the repeat regions.  The missing regions are excluded from coverage calculations.
Rscript ~/scripts/coverage-plot-2.R $workingfolder/$samplename.depth


































