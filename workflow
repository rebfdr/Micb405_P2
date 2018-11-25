####################################
# Genome Annotation
####################################

#!/bin/bash

output_dir=/projects/micb405/resources/project_2/2018/SaanichInlet_120m/MetaBAT2_SaanichInlet_120m/
prokka_output_dir=~/prokka_output
for f in $output_dir/MedQPlus_MAGs/*fa
do
    sid=$( basename $f | sed 's/.fa//g' )
    bin_num=$( basename $f | sed 's/.fa//g' | sed 's/SaanichInlet_120m.//g')
    tax=$(grep -w $sid $output_dir/gtdbtk_output/gtdbtk.*.classification_pplacer.tsv | awk '{ print $2 }' | awk -F";" '{ print $1 }' | sed 's/d__//g')
    echo $sid,$tax,$bin_num
    prokka --kingdom $tax --outdir $prokka_output_dir/$bin_num/ --force $f
done

##################################
# Metabolic Reconstruction
##################################
#use Github to download the .faa files onto desktop
scp -r rdremedios_mb18@orca-wg.bcgsc.ca:/home/rdremedios_mb18/prokka_output/* ~/Desktop/micb405_P2
## then upload each .faa file into KAAS using the metagenome default options
## download the resulting .txt file

#create file with only annotated rows
grep '\sK' ~/Desktop/micb405_P2/kaas_output_txt/*.txt > *.cleaned.txt

### or ###

# move .faa files into one .faa file in putty
cat prokka_ouput/*/*.faa > prokka_final.faa

#use Github to download the .faa file onto desktop
scp -r rdremedios_mb18@orca-wg.bcgsc.ca:/home/rdremedios_mb18/prokka_final.faa ~/Desktop/micb405_P2
## then upload each .faa file into KAAS using the metagenome default options
## download the resulting .txt file

#create file with only annotated rows
grep '\sK' ~/Desktop/micb405_P2/final_kass.txt > ~/Desktop/micb405_P2/final_kass.cleaned.txt

##################################
# Transcriptional Activity
##################################
#concatenate files for reference
 cat prokka_output/*/*.ffn > ref.fasta
 
 #index reference
 bwa index ref.fasta
 
 #align reads, generate sam files for each cruise
 for filename in SI042_120m.qtrim.artifact.rRNA.clean.fastq SI048_120m.qtrim.artifact.rRNA.clean.fastq.gz SI072_120m.qtrim.artifact.rRNA.clean.fastq.gz SI073_120m.qtrim.artifact.rRNA.clean.fastq.gz SI075_120m.qtrim.artifact.rRNA.clean.fastq.gz ; do
    echo "Processing alignment for - $filename"
    bwa mem -t 2 ~/project2/ref.fasta /projects/micb405/resources/project_2/2018/Metatranscriptomes/${filename} > ~/project2/${filename}___outfile.sam
    echo "Done alignment - $filename"
done

#generate csv files
/projects/micb405/resources/project_2/2018/rpkm \
-c project2/ref.fasta \
-a project2/SI042_120m.qtrim.artifact.rRNA.clean.fastq.gz___outfile.sam \
-o project2/rpkm_outputs/SI042_120m.csv

#download .csv files
scp -r rdremedios_mb18@orca-wg.bcgsc.ca:/home/rdremedios_mb18/project2/rpkm_outputs ~/Desktop/micb405_P2

####################################
#figure generation
####################################

for f in prokka_output/*/*faa
do
prokka_id=$( head -1 $f | awk -F_ '{ print $1 }' | sed 's/^>//g' )
mag_id=$( echo $f | sed 's/.faa//g')
echo $prokka_id,$mag_id
done >Prokka_MAG_map.csv

#download Prokka_MAG_map.csv file
scp rdremedios_mb18@orca-wg.bcgsc.ca:/home/rdremedios_mb18/Prokka_MAG_map.csv ~/Desktop/micb405_P2

#concatenate RPKM.csv files from each cruise

#download the other two files from gtdbk

#load files in ko
ko <- read.table("C:/Users/Rebecca/Desktop/micb405_P2/final_kass.cleaned.txt") %>% 
  dplyr::rename(orf = V1) %>% 
  dplyr::rename(ko = V2)

rpkm <- read.table("C:/Users/Rebecca/Desktop/micb405_P2/rpkm_ouputs/final_rpkm.csv", sep=',') %>% 
  dplyr::rename(orf = V1) %>% 
  dplyr::rename(rpkm = V2)
  
prokka_mag_map <- read.table("C:/Users/Rebecca/Desktop/micb405_P2/Prokka_MAG_map.csv", header=F, sep=',') %>% 
  dplyr::rename(prokka_id = V1) %>% 
  dplyr::rename(mag = V2)

arc_class <- read.table("C:/Users/Rebecca/Desktop/micb405_P2/gtdbtk.ar122.classification_pplacer.tsv", sep="\t")

bac_class <- read.table("C:/Users/Rebecca/Desktop/micb405_P2/gtdbtk.bac120.classification_pplacer.tsv", sep="\t")
