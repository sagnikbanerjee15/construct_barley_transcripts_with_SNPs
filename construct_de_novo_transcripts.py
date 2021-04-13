
import sys
import os

def downloadCI16151data():
    srr_id_filename = "/90daydata/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/raw_data/list_of_ids"
    cmd  = f" /90daydata/maizegdb/sagnik/Finder//utils/downloadAndDumpFastqFromSRA.py "
    cmd += f"-s {srr_id_filename}"
    cmd += f"-n 32"
    cmd += f"-o /90daydata/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/raw_data/ "
    cmd += f"1> /90daydata/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/raw_data/download.output "
    cmd += f"1> /90daydata/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/raw_data/download.error "
    os.system(cmd)

def alignReadsToBarleyBlumeria():
    srr_id_filename = "/work/LAS/rpwise-lab/sagnik/RNA-Seq_90_sample/data/raw_data/list_of_ids"
    srr_ids = open(srr_id_filename,"r").read().split("\n")[:-1]
    for srr_id in srr_ids:
        cmd  = "STAR "
        cmd += " --runThreadN 32 "
        cmd += " --genomeDir /90daydata/maizegdb/sagnik/data/finder/Hordeum_vulgare/genome/star_index "
        cmd += " --readFilesIn /90daydata/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/raw_data/"+srr_id+".fastq "
        cmd += " --alignIntronMin 20  --alignIntronMax 10000 "
        cmd += " --limitBAMsortRAM 107374182400 "
        cmd += " --outSAMattributes NH HI AS nM NM MD jM jI XS "
        cmd += " --genomeLoad LoadAndKeep "
        cmd += " --outFileNamePrefix /90daydata/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/raw_data//STAR_alignments/"+srr_id+"_ "
        cmd += " > /90daydata/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/raw_data/"+srr_id+".output "
        cmd += " > /90daydata/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/raw_data/"+srr_id+".error "
        if os.path.exists("/90daydata/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/raw_data/"+srr_id+"_Aligned.sortedByCoord.bam"):
            os.system(cmd)
        

def constructDeNovoTranscript(gene):
    gff3_filename = "/90daydata/maizegdb/sagnik/construct_barley_transcripts_with_SNPs/Barley_Morex_V2_gene_annotation_PGSB.all.gff3"

def mergeAllSamples():
    srr_id_filename = "/work/LAS/rpwise-lab/sagnik/RNA-Seq_90_sample/data/raw_data/list_of_ids"
    srr_ids = open(srr_id_filename,"r").read().split("\n")[:-1]
    cmd = "samtools merge -@ 60 "
    for srr_id in srr_ids:
        cmd += ""
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
    
    alignReadsToBarleyBlumeria()
    
    mergeAllSamples()
    
    for gene in horvu_genes:
        constructDeNovoTranscript(gene)
    
if __name__ == "__main__":
    main()