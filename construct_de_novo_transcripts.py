
import sys
import os

def downloadCI16151data():
    srr_id_filename = "/project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/raw_data/list_of_ids"
    cmd  = f" /project/maizegdb/sagnik/FINDER/Finder//utils/downloadAndDumpFastqFromSRA.py "
    cmd += f"--sra {srr_id_filename} "
    cmd += f" --cpu 60 "
    cmd += f" --output /project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/raw_data/ "
    cmd += f" 1> /project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/raw_data/download.output "
    cmd += f" 2> /project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/raw_data/download.error "
    os.system(cmd)

def alignReadsToBarley():
    srr_id_filename = "/project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/raw_data/list_of_ids"
    srr_ids = open(srr_id_filename,"r").read().split("\n")[:-1]
    for srr_id in srr_ids:
        print(srr_id)
        sys.stdout.flush()
        cmd  = "STAR "
        cmd += " --runThreadN 60 "
        cmd += " --genomeDir /project/maizegdb/sagnik/data/finder/Hordeum_vulgare/genome/star_index "
        cmd += " --readFilesIn /project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/raw_data/"+srr_id+".fastq "
        cmd += " --alignIntronMin 20  --alignIntronMax 10000 "
        cmd += " --limitBAMsortRAM 107374182400 "
        cmd += " --outSAMattributes NH HI AS nM NM MD jM jI XS "
        cmd += " --genomeLoad LoadAndKeep "
        cmd += " --outSAMtype BAM SortedByCoordinate "
        cmd += " --outFileNamePrefix /project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/STAR_alignments/"+srr_id+"_ "
        cmd += " 1> /project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/STAR_alignments/"+srr_id+".output "
        cmd += " 2> /project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/STAR_alignments/"+srr_id+".error "
        if os.path.exists("/project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/STAR_alignments/"+srr_id+"_Aligned.sortedByCoord.out.bam")==False:
            os.system(cmd)
        

def constructDeNovoTranscript(gene):
    gff3_filename = "/project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/Barley_Morex_V2_gene_annotation_PGSB.all.gff3"
    cmd = f"cat {gff3_filename} |grep {gene}|grep gene > /project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/{gene}_loc"
    os.system(cmd)
    line = open(f"/project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/{gene}_loc","r").read().strip()
    line = line.strip().split("\t")
    chromosome = line[0]
    start = line[3]
    end = line[4]
    
    # Retrieve region from merged bam file
    cmd  = f"samtools view -@ 60 /project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/STAR_alignments/merged.bam {chromosome}:{start}-{end}"
    cmd += "|awk '{print \"@\"$1\"/1\"\"\\n\"$10\"\\n+\\n\"$11}' "
    cmd += f"> /project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/STAR_alignments/merged_{gene}.fastq "
    os.system(cmd)
    """
    cmd = "spades.py --rna "
    cmd += f" -s /project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/STAR_alignments/CI16151_merged_{gene}.fastq " 
    cmd += " -t 60 "
    cmd += " -k 127 "
    cmd += f" -o /project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/contigs/CI16151_merged_{gene} "
    if os.path.exists(f"/project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/contigs/CI16151_merged_{gene}/transcripts.fasta")==False:
        os.system(cmd)
    """
    # Align the long contigs to barley genome using gmap
    """cmd  = "gmapl "
    cmd += " -D /project/maizegdb/sagnik/data/finder/Hordeum_vulgare/genome/ "
    cmd += " -d gmap_index "
    cmd += " --min-intronlength=20 "
    cmd += " --max-intronlength-middle=10000 "
    cmd += " --max-intronlength-ends=10000 "
    cmd += " --no-chimeras "
    cmd += " -t 60 "
    cmd += " -f samse "
    cmd += f" /project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/contigs/CI16151_merged_{gene}/transripts.fasta "
    cmd += f" > /project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/contigs/CI16151_merged_{gene}_gmap_aligned.sam "
    if os.path.exists(f"/project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/contigs/CI16151_merged_{gene}_gmap_aligned.sam")==False:
        os.system(cmd)
    """
    """
    cmd  = f" Trinity "
    cmd += f" --seqType fq "
    cmd += f" --single /project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/STAR_alignments/CI16151_merged_{gene}.fastq "
    cmd += f" --CPU 40 "
    cmd += f" --KMER_SIZE 32 "
    cmd += f" --max_memory 45G "
    cmd += f" --output /project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/contigs/CI16151_merged_{gene}_trinity "
    cmd += f" --full_cleanup "
    cmd += f"1> /project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/contigs/CI16151_merged_{gene}_trinity.output "
    cmd += f"2> /project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/contigs/CI16151_merged_{gene}_trinity.error "
    os.system(cmd)
    """
    fhw=open(f"/project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/contigs/merged_{gene}_soapdenovo.configuration_file","w")
    fhw.write(f"""#maximal read length
max_rd_len=151
[LIB]
#maximal read length in this lib
rd_len_cutof=151
#average insert size
avg_ins=200
#if sequence needs to be reversed 
reverse_seq=0
#in which part(s) the reads are used
asm_flags=3
#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)
map_len=32
#fastq file for read 1 
q=/project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/STAR_alignments/CI16151_merged_{gene}.fastq""")
    fhw.close()
    
    cmd  = f" SOAPdenovo-Trans-127mer "
    cmd += f" all "
    cmd += f" -o /project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/contigs/merged_{gene}_soapdenovo "
    cmd += f" -K 127 "
    cmd += f" -p 40 "
    cmd += f" -s /project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/contigs/merged_{gene}_soapdenovo.configuration_file "
    cmd += f" 1> /project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/contigs/merged_{gene}_soapdenovo.output "
    cmd += f" 2> /project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/contigs/merged_{gene}_soapdenovo.error "
    os.system(cmd)
    
    cmd  = "STARlong "
    cmd += " --runThreadN 60 "
    cmd += " --genomeDir /project/maizegdb/sagnik/data/finder/Hordeum_vulgare/genome/star_index "
    cmd += f" --readFilesIn /project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/contigs/merged_{gene}/transcripts.fasta "
    cmd += " --alignIntronMin 20  --alignIntronMax 10000 "
    cmd += " --limitBAMsortRAM 107374182400 "
    cmd += " --outSAMattributes NH HI AS nM NM MD jM jI XS "
    cmd += " --genomeLoad LoadAndKeep "
    cmd += " --outSAMtype BAM SortedByCoordinate "
    cmd += f" --outFileNamePrefix /project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/contigs/merged_{gene}_ "
    cmd += f" 1> /project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/contigs/merged_{gene}.output "
    cmd += f" 2> /project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/contigs/merged_{gene}.error " 
    #if os.path.exists(f"/project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/contigs/merged_{gene}_Aligned.sortedByCoord.out.bam")==False:
    os.system(cmd)
    
    #if os.path.exists(f"/project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/contigs/merged_{gene}_Aligned.sortedByCoord.out.bam.csi")==False:
    cmd = f"samtools index -c /project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/contigs/merged_{gene}_Aligned.sortedByCoord.out.bam"
    os.system(cmd)
        

def mergeAllSamples():
    srr_id_filename = "/project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/raw_data/list_of_ids"
    srr_ids = open(srr_id_filename,"r").read().split("\n")[:-1]
    cmd = "samtools merge -@ 60 /project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/STAR_alignments/merged.bam " 
    for srr_id in srr_ids:
        cmd += "/project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/STAR_alignments/"+srr_id+"_Aligned.sortedByCoord.out.bam "
    if os.path.exists("/project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/STAR_alignments/merged.bam")==False:
        os.system(cmd)
    cmd="samtools index -c /project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/STAR_alignments/merged.bam "
    if os.path.exists("/project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/STAR_alignments/merged.bam.csi")==False:
        os.system(cmd)
    

def main():
    horvu_genes = """HORVU.MOREX.r2.4HG0338120
HORVU.MOREX.r2.1HG0006350
HORVU.MOREX.r2.3HG0266360
HORVU.MOREX.r2.5HG0362710
HORVU.MOREX.r2.4HG0338560
HORVU.MOREX.r2.6HG0455190
HORVU.MOREX.r2.2HG0174690
HORVU.MOREX.r2.7HG0616250
HORVU.MOREX.r2.7HG0616260
HORVU.MOREX.r2.7HG0611210
HORVU.MOREX.r2.4HG0276000
HORVU.MOREX.r2.6HG0464980
HORVU.MOREX.r2.3HG0252860
HORVU.MOREX.r2.6HG0464870
HORVU.MOREX.r2.4HG0298980
HORVU.MOREX.r2.4HG0316820
HORVU.MOREX.r2.5HG0385430
HORVU.MOREX.r2.5HG0396940
HORVU.MOREX.r2.7HG0611160
HORVU.MOREX.r2.5HG0437940
HORVU.MOREX.r2.1HG0058670"""
    horvu_genes = horvu_genes.split("\n")
    
    downloadCI16151data()
    
    alignReadsToBarley()
    
    mergeAllSamples()
    
    for gene in horvu_genes:
        constructDeNovoTranscript(gene)
    
    # Merge all alignments
    cmd="samtools merge -f -@ 60 /project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/final.bam "
    for gene in horvu_genes:
        cmd+= f"/project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/contigs/merged_{gene}_Aligned.sortedByCoord.out.bam "
    os.system(cmd)
    """
    cmd = "samtools merge -@ 60 /project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/merged_from_90_samples.bam  "
    for gene in horvu_genes:
        cmd += f" /project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/STAR_alignments/CI16151_merged_{gene}.fastq "
    os.system(cmd)
    cmd="samtools index -c -@ 60 /project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/merged_from_18_samples.bam"
    os.system(cmd)
    """
    cmd=" samtools index -c -@ 60 /project/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/final.bam "
    os.system(cmd)
        
    
if __name__ == "__main__":
    main()