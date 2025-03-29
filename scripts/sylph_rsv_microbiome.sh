# this is profiling using sylph : https://github.com/bluenote-1577/sylph
# databases are downloaded into ~/Databases/BacterialDatabases/ ../FungalDatabases ../ViralDatabases
# the databases are preconfigured sketches for use in sylph

cd ~/Analysis/TrimmedMSSFiles/

files=$(cat ~/Documents/Code/RSV-microbiome-2025/sample_list.txt)
for f in $files; do

cd ~/Analysis/TrimmedMSSFiles/

sylph sketch -1 trimmed_read1_${f}.fastq.gz -2 trimmed_read2_${f}.fastq.gz -d Sketches/${f}.sylsp

done

files=$(cat ~/Documents/Code/RSV-microbiome-2025/sample_list.txt)
for f in $files; do

cd ~/Analysis/TrimmedMSSFiles/


# Bacterial query and profiling

sylph query  ~/Databases/BacterialDatabases/gtdb-r220-c200-dbv1.syldb Sketches/${f}.sylsp/trimmed_read1_${f}.fastq.gz.paired.sylsp -t 32 > SylphProfiles/${f}_ani_queries.tsv
sylph profile  ~/Databases/BacterialDatabases/gtdb-r220-c200-dbv1.syldb Sketches/${f}.sylsp/trimmed_read1_${f}.fastq.gz.paired.sylsp -t 32 > SylphProfiles/${f}_profiled_bacteria.tsv


# Viral profiling
sylph profile ~/Databases/ViralDatabases/imgvr_c200_v0.3.0.syldb Sketches/${f}.sylsp/trimmed_read1_${f}.fastq.gz.paired.sylsp -t 24 > SylphProfiles/${f}_profiled_viruses.tsv 

# the problem here is I don't know how to match up the contig name with virus
# this is the portal for the viral database and one of them should have features that map to the contigs


# Fungal profiling
sylph profile ~/Databases/FungalDatabases/fungi-refseq-2023nov28-c200-v0.3.syldb Sketches/${f}.sylsp/trimmed_read1_${f}.fastq.gz.paired.sylsp -t 24 > SylphProfiles/${f}_profiled_fungi.tsv


done
