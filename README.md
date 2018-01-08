# Tuxedo Pipeline with Modifications and Annotations

### Preliminary steps:
* The genome fasta file and annotation gff file downloaded from NCBI or Ensemble
  * Which database to use depends on downstream tasks, if using STAR Ensemble is usually better after its release from embargo and publication of a recent paper, on 07/06/2016.  
* The Ion Torrent Sequencer outputs the read files as Binary Alignment Map, or “bam” files, and for most analyses, this must be converted to the standard “fastq” read format.  
* Optional: The Ion Torrent Server converted the files to fastq files, which were then downloaded and copied to the genelabs1.lsu.edu “Partek Flow” computer.  
  * For simplicity, Partek Flow was used to remove adapters and low quality bases with the program cutadapt, then following [Generally Accepted](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-91) [Bioinformatics](http://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-4-r36) [Practices](http://bib.oxfordjournals.org/content/14/6/671.long) and a [technical note from Life Technologies](https://ioncommunity.thermofisher.com/docs/DOC-7062), the RNA fastq read files were mapped with Tophat2 to accommodate splice junctions
* Otherwise, [the following pipeline and commands were used](https://github.com/disulfidebond/tuxedo_pipeline_with_modifications/blob/master/pipeline.sh).
