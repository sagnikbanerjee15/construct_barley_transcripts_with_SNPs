# Construct barley transcripts with SNPs

De novo recontruction of barley transcription factors from CI16151 to preserve SNPs and indels

The sequence of 25 barley genes were reconstructed from RNA-Seq data collected from a infection time-course experiment from a single cultivar (**CI16151**). Raw reads, from a total of 18 CI16151 samples, were mapped to the barley genome using a STAR with the `alignIntronMin` and  `alignIntronMax` set to 20 and 10,000 respectively. The genomic locus for each of the 25 genes were extracted from the V2 version of Barley Morex annotation [citation reqd.]. Aligned CI16151 raw reads were then subset for only those extracted loci and converted into a fastq file. Trinity was used to assemble the reawreads *de novo* to preserve the SNPs and indels. Out of 25 genes, 11 genes were found to contain at least one SNP.
