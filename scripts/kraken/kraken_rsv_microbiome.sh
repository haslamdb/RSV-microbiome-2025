# This script runs Kraken2 and Bracken on RSV microbiome samples
# Database is located at ~/Databases/KrakenPlusPF

# Read sample list
files=$(cat ~/Documents/Code/RSV-microbiome-2025/sample_list.txt)

# Run Kraken2 and Bracken on all samples. I can't get new Kraken2 versions to run, so going with the older version
for f in $files; do
    echo "Processing sample ${f} with Kraken2..."
    
    # Run Kraken2 with paired-end reads
    /usr/local/Kraken2_Old/kraken2 --paired \
        --threads 12 \
        --db ~/Databases/KrakenPlusPF \
        --output ~/Documents/Alignments/KrakenAlignments/KrakenNew/${f}.kraken \
        --report ~/Documents/Alignments/KrakenAlignments/KrakenNew/${f}.kreport \
        ~/Analysis/TrimmedMSSFiles/trimmed_read1_${f}.fastq.gz \
        ~/Analysis/TrimmedMSSFiles/trimmed_read2_${f}.fastq.gz
    
    # Run Bracken for abundance estimation at species level
    echo "Running Bracken for sample ${f}..."
    bracken -d ~/Databases/KrakenPlusPF \
        -i ~/Documents/Alignments/KrakenAlignments/KrakenNew/${f}.kreport \
        -o ~/Documents/Alignments/KrakenAlignments/KrakenNew/${f}_S_abundance.txt \
        -w ~/Documents/Alignments/KrakenAlignments/KrakenNew/${f}_bracken.kreport \
        -r 150 -l S -t 10
    
    # Run Bracken for genus level
    bracken -d ~/Databases/KrakenPlusPF \
        -i ~/Documents/Alignments/KrakenAlignments/KrakenNew/${f}.kreport \
        -o ~/Documents/Alignments/KrakenAlignments/KrakenNew/${f}_G_abundance.txt \
        -w ~/Documents/Alignments/KrakenAlignments/KrakenNew/${f}_G_bracken.kreport \
        -r 150 -l G -t 10
    
    # Run Bracken for family level
    bracken -d ~/Databases/KrakenPlusPF \
        -i ~/Documents/Alignments/KrakenAlignments/KrakenNew/${f}.kreport \
        -o ~/Documents/Alignments/KrakenAlignments/KrakenNew/${f}_F_abundance.txt \
        -w ~/Documents/Alignments/KrakenAlignments/KrakenNew/${f}_F_bracken.kreport \
        -r 150 -l F -t 10
    
    echo "Completed processing sample ${f}"
done

echo "All samples have been processed with Kraken2 and Bracken."
