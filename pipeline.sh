# Tuxedo Pipeline with Modifications and Annotations

### Preliminary steps:
* The genome fasta file and annotation gff file downloaded from NCBI or Ensemble
  * Which database to use depends on downstream tasks, if using STAR Ensemble is usually better after its release from embargo and publication of a recent paper, on 07/06/2016.  
* The Ion Torrent Sequencer outputs the read files as Binary Alignment Map, or “bam” files, and for most analyses, this must be converted to the standard “fastq” read format.  
* Optional: The Ion Torrent Server converted the files to fastq files, which were then downloaded and copied to the genelabs1.lsu.edu “Partek Flow” computer.  
  * For simplicity, Partek Flow was used to remove adapters and low quality bases with the program cutadapt, then following GABP and a technical note from Life Technologies, the RNA fastq read files were mapped with Tophat2 to accommodate splice junctions
* Otherwise, the following pipeline and commands were used:

# prior to running tophat2, use:
gffread gfffile.gff -T -o gfffile.gff.gtf
# then spot check and scan file, generally downstream HTSEQ usage is 
# --type=CDS --idattr=gene_name
# Tophat will find features automatically
/gpfs/gpfs0/tools/biobuilds/biobuilds-2016.04/bin/tophat2 -p 18 -o FASTQFILE.fastq.tophatout -G /gpfs/gpfs0/home/jcaskey/ThuneCatfish/tophat-catfish/Catfish.gtf Catfish.fasta FASTQFILE.fastq

# Followed by split pipeline:
picardjdk SortSam I=FASTQFILE.fastq.tophatout/accepted_hits.bam O=FASTQFILE.fastq.tophataligned.coordsort.bam SO=coordinate
##
bedtools bamtofastq -i FASTQFILE.fastq.tophatout/unmapped.bam -fq FASTQFILE.unaligned.fastq 
#
bowtie2 -x Catfish.fasta -p 18 --local --very-sensitive-local -U FASTQFILE.unaligned.fastq -S FASTQFILE.unaligned.fastq.unalignedBt2.sam
#
picardjdk SortSam I=FASTQFILE.unaligned.fastq.unalignedBt2.sam O=FASTQFILE.unalignedBt2.sam.coordsort.bam SO=coordinate
## 
## ReorderSam required b/c contigs and indices do not match with #
## MergeSamFiles command
## 
picardjdk ReorderSam I=FASTQFILE.fastq.tophataligned.coordsort.bam O=FASTQFILE.fastq.tophataligned.coordsort.bam.reord.bam R=Catfish.fasta 
#
picardjdk ReorderSam I=FASTQFILE.unalignedBt2.sam.coordsort.bam O=FASTQFILE.unalignedBt2.sam.coordsort.reord.bam R=Catfish.fasta
##
## can use for loop as follows:
## 
## j=1
## NM=“”
## for i in *.bam ; do if !(($j%2)) ; then echo “submitting ${i} and ${NM} for merge job ${j}” ; bsub -e picardMerge.${j}.err -J “PicardMerge Job ${j}” -q alignment picardjdk MergeSamFiles I=${i} I=${NM} O=OUTFILE.bam SO=coordinate ; else NM=${i} ; fi ; let j=$j+1 ; done
picardjdk MergeSamFiles I=FASTQFILE.fastq.tophataligned.coordsort.bam.reord.bam I=FASTQFILE.unalignedBt2.sam.coordsort.reord.bam O=FASTQFILE.tophataligned.coordsort.reord.bam.merged.bam AS=false SO=coordinate
##
## CODA1: Cufflinks output ####
################################
## use cuffdiff with the following options 
## in the tophat-catfish directory
## library-type specific to Ion Torrent
## Only RPKM works with Biobuilds on Delta Cluster
## repeat for each, also can use for loop as follows:
## 
## j=1
## for i in *.merged.bam ; do bsub -e cufflnk${j}.err -J “Cufflink Job ${j}” -q alignment … ; let j=$j+1 ; done ;
## 
cufflinks -p 18 --library-type fr-secondstrand -G /gpfs/gpfs0/home/jcaskey/ThuneCatfish/tophat-catfish/Catfish.gtf -o FASTQFILE.tophataligned.coordsort.reord.bam.merged.bam.cufflinkout -b /gpfs/gpfs0/home/jcaskey/ThuneCatfish/tophat-catfish/Catfish.fasta FASTQFILE.tophataligned.coordsort.reord.bam.merged.bam
##
## for cuffmerge, use echo $file >> filename.txt for each
## PWDVAR=$(pwd)
## for i in *.cufflinks_out ; do echo "${PWDVAR}/${i}/transcripts.gtf" >> cuffmergeList.txt ; done
## optionally use for loop 
## then check
cuffmerge -p 18 -g Catfish.gtf -s Catfish.fasta -o OUTFILE_out filename.txt
##
## for cuffdiff, first assign each variable to variable name, example:
## WTBAM1=R_2015_03_10_12_53_03_user_LSU-44-Lidia_catfish_RNAseq_sampleC_wt_79_trimmed.coordsort.bam.reord.bam.merged.bam.bam.coordsort.bam
## then do
bsub -J “CuffDiff Job" -e CJobErr1.err -q alignment /gpfs/gpfs0/tools/biobuilds/biobuilds-2016.04/bin/cuffdiff -p 18 --library-type fr-secondstrand -o updated_wtDesG_cuffdiffout -L WT,DesG -b Catfish.fasta updated_wtDesG_cuffmergeOut/merged.gtf $WTBAM1,$WTBAM2,$WTBAM3 $DesG1BAM,$DesG2BAM,$DesG3BAM
##
## CODA2: DESeq2 Output ####
############################
# use htseq-count 
# with the following options in htseq-catfish directory:
# note a for loop can be employed similar to described above
# note also that union seems to work the best 
# and is recommended by Huber et al.
# bsub -J “htseqjob 1" -e CJup2.err -q alignment -o BAMFILE.bam.union.htseqout.txt /gpfs/gpfs0/tools/biobuilds/biobuilds-2016.04/bin/htseq-count --format=bam --type=CDS --idattr=gene_name BAMFILENAME.coordsort.merged.bam Catfish.gtf 
# The htseqoutput will need to be parsed, 
# because it will contain a header and footer.  
# the head and foot values will vary, and note that 
# the command is GNU-specific (it will not work on BSD)
tail -n +34 BAMFILE.coordsort.merged.bam.union.htseqout.txt | head -n -7 > BAMFILE.coordsort.merged.bam.union.htseqout.parsed.txt
# perl or manually editing the file are alternatives
# then use with DESeq2, as per Dr. Kern, Dunnet’s test-multiple hypothesis test for groups, since goal is to compare each group back to WT (no interaction variable).  In the case of DESeq2, use nBinom with Wald’s test-load each dataset, set design to treatment, where treatment is the mutation, for example WT vs. DesG, treatment is mutant or WT, repeat for uninf vs. wt, wt vs. 65st, etc.  Another possibility is LRT, but probably isn’t necessary.

# For above, cufflinks/cuffdiff assesses expression between experimental points, and creates combinatorial expression curves among transcripts, then classified into increasing, decreasing, flat, or mixed.  The Jensen-Shannon divergence was used between points, and a matrix of the abundance estimate
