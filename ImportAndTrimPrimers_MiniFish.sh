#!/bin/bash

# Usage: ./script_name.sh plate1 plate2 plate3 ...

# Function to process statistics from primer trimming
process_stats() {
    local plate=$1
    local input_file=$2
    local output_file=$3
    
    grep "Total written (filtered):" "$input_file" | \
    sed 's/.*(\([0-9.]*\)%)/\1/' | \
    awk -v plate="$plate" '
    BEGIN {
        sum = 0
        sum_sq = 0
        n = 0
    }
    {
        sum += $1
        sum_sq += $1 * $1
        n++
    }
    END {
        if (n > 0) {
            mean = sum / n
            variance = (sum_sq / n) - (mean * mean)
            printf "Plate %s:\nMean: %.2f%%\nVariance: %.2f\n", plate, mean, variance
        }
    }' >> "$output_file"
}

# Check if any arguments were provided
if [ $# -eq 0 ]; then
    echo "No plate numbers provided"
    echo "Usage: $0 plate1 plate2 plate3 ..."
    exit 1
fi

source /Users/gc547/miniconda3/bin/activate qiime2-amplicon-2024.10 

echo "Processing $# plate numbers"
date_suffix=$(date +%Y%m%d_%H%M%S)

# Create or clear the statistics files
stats_file1="all_plates_statistics_trim1_${date_suffix}.txt"
stats_file2="all_plates_statistics_trim2_${date_suffix}.txt"
> "$stats_file1"
> "$stats_file2"

# Single loop through all arguments
for K in "$@"
do
    echo "Starting processing of plate: $K at $(date)"

    # Import
    qiime tools import\
    --type 'SampleData[PairedEndSequencesWithQuality]'\
    --input-path /Volumes/Data_SS1/MiFish/Plate$K/reads/ \
    --input-format CasavaOneEightSingleLanePerSampleDirFmt\
    --output-path demux_Plate$K.qza || { echo "Import failed for Plate $K"; continue; }

    # Create demultiplexing summary
    qiime demux summarize \
    --i-data demux_Plate$K.qza \
    --o-visualization demux_Plate$K.qzv || { echo "Demux summarize failed for Plate $K"; continue; }

    # Trim 3' primers
    echo "Trimming 3' primers from plate: $K"
    qiime cutadapt trim-paired \
    --i-demultiplexed-sequences demux_Plate$K.qza \
    --p-adapter-f CAACCACGAAAGTGGCTTTA \
    --p-adapter-r CAACTTGGGCCTCTCGTATAAC \
    --o-trimmed-sequences trimd_Plate$K.qza \
    --verbose > cutadapt_out_Plate$K.txt || { echo "First trimming failed for Plate $K"; continue; }

    # Process statistics for first trim
    process_stats "$K" "cutadapt_out_Plate$K.txt" "$stats_file1"

    # Trim 5' primers
    echo "Trimming 5' primers from plate: $K"
    qiime cutadapt trim-paired \
    --i-demultiplexed-sequences trimd_Plate$K.qza \
    --p-front-f GTTATACGAGAGGCCCAAGTTG \
    --p-front-r TAAAGCCACTTTCGTGGTTG \
    --o-trimmed-sequences trimd2_Plate$K.qza \
    --verbose > cutadapt_out2_Plate$K.txt || { echo "Second trimming failed for Plate $K"; continue; }

    # Process statistics for second trim
    process_stats "$K" "cutadapt_out2_Plate$K.txt" "$stats_file2"

    # Denoise with Dada2
    echo "Denoising plate: $K"
    qiime dada2 denoise-paired \
    --i-demultiplexed-seqs trimd2_Plate$K.qza \
    --p-trunc-len-f 91 \
    --p-trunc-len-r 91 \
    --p-trim-left-f 0 \
    --p-trim-left-r 0 \
    --p-min-overlap 50 \
    --p-n-threads 8 \
    --o-representative-sequences rep-seqs_Plate$K.qza \
    --o-table table_Plate$K.qza \
    --o-denoising-stats denoise_Plate$K.qza \
    || { echo "Denoising failed for Plate $K"; continue; }  # Added error checking

    # Create denoising visualisation files
    echo "Denoising finished, creating qzv files for plate: $K"
    qiime metadata tabulate \
    --m-input-file denoise_Plate$K.qza \
    --o-visualization denoise_Plate$K.qzv \
    || { echo "Denoising stats visualization failed for Plate $K"; continue; }

    qiime feature-table tabulate-seqs \
    --i-data rep-seqs_Plate$K.qza \
    --o-visualization rep-seqs_Plate$K.qzv \
    || { echo "Feature table visualization failed for Plate $K"; continue; }

    # Optional: Add feature table summary
    qiime feature-table summarize \
    --i-table table_Plate$K.qza \
    --o-visualization table_Plate$K.qzv \
    || { echo "Table summary visualization failed for Plate $K"; continue; }

    echo "Completed processing Plate $K at $(date)"
done

# Add summary of all processing steps
echo "All plates have been processed at $(date)"
echo "Results have been written to:"
echo "- Primer trimming statistics: $stats_file1 and $stats_file2"
echo "- DADA2 output files: rep-seqs_Plate*.qza, table_Plate*.qza, denoise_Plate*.qza"
echo "- Visualization files: *.qzv"