#!/bin/bash

# Maps a set of pre-filtered viral reads to reference using snippy, then removes low coverage genomes, Generates a core snp genome, phylogeny, annotation with snpeff, and extracts gene sequences, aligns to HSV2 outgroup and calculates DoS per dataset

#usage : pipelinexxx.sh samplefile

mincov=5
minfrac=0.5
threads=7
basefolder=/home/ubuntu/analyses/snippy-pipelined
runfolder=$basefolder/pipeline-0.5_dedup

mkdir -v -p $runfolder/

echo
echo "Running snippy to call snps with >$mincov reads and >$minfrac proportion for all samples in $1"
echo
while read seqname read1 read2; do snippy --cpus $threads --mincov $mincov --minfrac $minfrac --prefix "$seqname" --force --outdir $runfolder/"$seqname" --ref ~/references/NC_001806.2.fasta --R1 ${read1} --R2 ${read2} ; done < $1
echo 

echo "Annotating all vcfs using snpEff"
echo
while read seqname read1 read2; do cd $runfolder/${seqname}/ ; java -jar ~/programs/snpEff/snpEff.jar NC_001806.2 $runfolder/${seqname}/${seqname}.filt.vcf > $runfolder/${seqname}/${seqname}.filt.ann.vcf ; mv $runfolder/${seqname}/snpEff_genes.txt $runfolder/${seqname}/${seqname}_snpEff_genes.txt ; mv $runfolder/${seqname}/snpEff_summary.html $runfolder/${seqname}/${seqname}_snpEff_summary.html ; done < $1

echo 
echo "Moving known low coverage genomes to separate folder to exclude from downstream analysis"
mkdir -p -v $runfolder/lowcov_genomes/
mv $runfolder/HSV1-nCSF2/ $runfolder/lowcov_genomes/HSV1-nCSF2/
mv $runfolder/HSV1-CSF9/ $runfolder/lowcov_genomes/HSV1-CSF9/
mv $runfolder/HSV1-CSF11/ $runfolder/lowcov_genomes/HSV1-CSF11/
mv $runfolder/HSV1-nCSF12/ $runfolder/lowcov_genomes/HSV1-nCSF12/

echo "Moving good coverage genomes to separate folder"
mkdir -p -v $runfolder/GoodCov_Genomes/
mv $runfolder/HSV1-*/ $runfolder/GoodCov_Genomes/

echo
echo "Now calling core snps from genomes with good coverage"
cd $runfolder/GoodCov_Genomes/
mkdir -p -v $runfolder/GoodCov_Genomes/Core_Snps/
while read seqname read1 read2; do cp $runfolder/GoodCov_Genomes/$seqname/$seqname.tab $runfolder/GoodCov_Genomes/$seqname/snps.tab ; done < $1
snippy-core --mincov=$mincov --prefix=HSV1_Core $runfolder/GoodCov_Genomes/HSV1-*/
mv $runfolder/GoodCov_Genomes/HSV1_Core*  $runfolder/GoodCov_Genomes/Core_Snps/

echo
echo "Now redo Core SNPs for deduplicated samples for PopGen"
mkdir -p -v $runfolder/GoodCov_Genomes/dedup_Core/
while read line ; do mv $runfolder/GoodCov_Genomes/$line/ $runfolder/GoodCov_Genomes/dedup_Core/$line/ ; done < ~/data/HSV1_deduplicated.samplelist
cd $runfolder/GoodCov_Genomes/dedup_Core/
snippy-core --mincov=$mincov --prefix=HSV1_Core-dedup $runfolder/GoodCov_Genomes/dedup_Core/HSV1-*/

echo 
echo "Compile all Core VCFs for group analysis, and run snpEff on deduplicated VCF"
mkdir -p -v $runfolder/GoodCov_Genomes/dedup_Core/VCFs/
find $runfolder/GoodCov_Genomes/dedup_Core/ -name "*.filt.ann.vcf" -exec cp "{}" $runfolder/GoodCov_Genomes/dedup_Core/VCFs/ \;
cp $runfolder/GoodCov_Genomes/dedup_Core/HSV1_Core-dedup.vcf $runfolder/GoodCov_Genomes/dedup_Core/VCFs/
cd $runfolder/GoodCov_Genomes/dedup_Core/VCFs/
java -jar ~/programs/snpEff/snpEff.jar NC_001806.2 $runfolder/GoodCov_Genomes/dedup_Core/VCFs/HSV1_Core-dedup.vcf > $runfolder/GoodCov_Genomes/dedup_Core/VCFs/HSV1_Core-dedup.ann.vcf 

echo
echo "generate Core genomes and annotate core VCFs for CSF and non-CSF datasets"
snippy-core --mincov=$mincov --prefix=HSV1_Core-dedup_CSF $runfolder/GoodCov_Genomes/dedup_Core/HSV1-CSF*/
snippy-core --mincov=$mincov --prefix=HSV1_Core-dedup_nCSF $runfolder/GoodCov_Genomes/dedup_Core/HSV1-nCSF*/

java -jar ~/programs/snpEff/snpEff.jar NC_001806.2 $runfolder/GoodCov_Genomes/dedup_Core/VCFs/HSV1_Core_CSF.vcf > $runfolder/GoodCov_Genomes/dedup_Core/VCFs/HSV1_Core-dedup_CSF.ann.vcf


java -jar ~/programs/snpEff/snpEff.jar NC_001806.2 $runfolder/GoodCov_Genomes/dedup_Core/VCFs/HSV1_Core_nCSF.vcf > $runfolder/GoodCov_Ge    nomes/dedup_Core/VCFs/HSV1_Core-dedup_nCSF.ann.vcf

echo 
echo "Extract gene sequences from low complexity consensus seqs (subs)"
mkdir -p -v $runfolder/GoodCov_Genomes/dedup_Core/fastas/
find $runfolder/GoodCov_Genomes/dedup_Core/ -name "*.consensus.subs.fa" -exec cp "{}" $runfolder/GoodCov_Genomes/dedup_Core/fastas/ \;
module load cufflinks
mkdir -p -v $runfolder/GoodCov_Genomes/dedup_Core/genefastas/
while read seqname read1 read2; do gffread ~/references/NC_001806.2.CDS.gff -O -g $runfolder/GoodCov_Genomes/dedup_Core/fastas/$seqname.consensus.subs.fa -x $runfolder/GoodCov_Genomes/dedup_Core/genefastas/$seqname.cds.fasta; perl -i -pe "s/^\>/\>$seqname/g" $runfolder/GoodCov_Genomes/dedup_Core/genefastas/$seqname.cds.fasta ; perl -i -pe "s/$seqname\w+\s+/$seqname\_/g" $runfolder/GoodCov_Genomes/dedup_Core/genefastas/$seqname.cds.fasta  ; done < $1 # generate CDS gene seqs for each wgs fasta and add in sequence name to each one

echo 
echo "Separating out genes into individual files"
for i in $runfolder/GoodCov_Genomes/dedup_Core/genefastas/*.cds.fasta; do while read line; do ~/scripts/selectSeqs.pl -m "$line$" $i >> $runfolder/GoodCov_Genomes/dedup_Core/genefastas/$line.fasta; done < ~/data/HSV1_nonR-gene-names_subsetted.txt; done
while read line ; do cp $runfolder/GoodCov_Genomes/dedup_Core/genefastas/$line.fasta $runfolder/GoodCov_Genomes/dedup_Core/genefastas/$line.HSV2.fasta ; done < ~/data/HSV1_nonR-gene-names_subsetted.txt


echo 
echo "Adding in HSV2 gene sequences from NC_001798.2.HSV2.cds.fasta"
while read line; do ~/scripts/selectSeqs.pl -m "$line$" ~/references/NC_001798.2.HSV2.cds.fasta >> $runfolder/GoodCov_Genomes/dedup_Core/genefastas/$line.HSV2.fasta; done < ~/data/HSV1_nonR-gene-names_subsetted.txt

echo 
echo "Now align all genes using PRANK"
while read line; do prank -d=${runfolder}/GoodCov_Genomes/dedup_Core/genefastas/$line.HSV2.fasta -o=${runfolder}/GoodCov_Genomes/dedup_Core/genefastas/$line.HSV2.prank -f=fasta -codon -DNA -iterate=10 -F ; done < ~/data/HSV1_nonR-gene-names_subsetted.txt

echo 
echo "Now use core snp alignment (not deduplicated) in RAxML"
~/scripts/Fasta2Phylip.pl $runfolder/GoodCov_Genomes/Core_Snps/HSV1_Core.aln $runfolder/GoodCov_Genomes/Core_Snps/HSV1_Core.phy
mkdir -p -v $runfolder/GoodCov_Genomes/Core_Snps/raxml/
raxmlHPC-PTHREADS -s $runfolder/GoodCov_Genomes/Core_Snps/HSV1_Core.phy -n HSV1_Core$minfrac -T $threads -m GTRCAT -f a -x 12345 -p 12345 -N 1000 -w $runfolder/GoodCov_Genomes/Core_Snps/raxml/



