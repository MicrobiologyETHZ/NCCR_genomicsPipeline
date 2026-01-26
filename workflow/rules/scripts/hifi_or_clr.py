#!/usr/bin/env python3
"""
PacBio Read Type Detection Script
Determines if reads are HiFi (CCS) or CLR based on quality and length characteristics

Usage: python detect_pacbio_type.py <fastq_file>
"""

import sys
import gzip
from pathlib import Path
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt
import argparse


def analyze_pacbio_reads(fastq_file, sample_size=10000):
    """
    Analyze PacBio reads to determine if they are HiFi or CLR

    HiFi characteristics:
    - Higher mean quality scores (Q20+, typically Q30+)
    - More uniform length distribution
    - Lengths typically 10-25kb
    - Quality scores more consistent across read length

    CLR characteristics:
    - Lower mean quality scores (Q10-15)
    - Broader length distribution
    - Can be very long (>50kb) but lower quality
    - Quality often decreases toward read ends
    """

    print(f"Analyzing {fastq_file}...")

    # Determine if file is gzipped
    opener = gzip.open if str(fastq_file).endswith('.gz') else open
    mode = 'rt' if str(fastq_file).endswith('.gz') else 'r'

    lengths = []
    qualities = []
    mean_qualities = []

    try:
        with opener(fastq_file, mode) as handle:
            for i, record in enumerate(SeqIO.parse(handle, "fastq")):
                if i >= sample_size:
                    break

                length = len(record.seq)
                quality_scores = record.letter_annotations["phred_quality"]
                mean_qual = np.mean(quality_scores)

                lengths.append(length)
                qualities.extend(quality_scores)
                mean_qualities.append(mean_qual)

                if i % 1000 == 0:
                    print(f"Processed {i} reads...")

    except Exception as e:
        print(f"Error reading file: {e}")
        return None

    if not lengths:
        print("No reads found in file!")
        return None

    # Calculate statistics
    stats = {
        'total_reads_analyzed': len(lengths),
        'mean_length': np.mean(lengths),
        'median_length': np.median(lengths),
        'length_std': np.std(lengths),
        'min_length': np.min(lengths),
        'max_length': np.max(lengths),
        'mean_quality': np.mean(mean_qualities),
        'median_quality': np.median(mean_qualities),
        'quality_std': np.std(mean_qualities),
        'min_quality': np.min(mean_qualities),
        'max_quality': np.max(mean_qualities),
        'reads_above_q20': sum(1 for q in mean_qualities if q >= 20),
        'reads_above_q30': sum(1 for q in mean_qualities if q >= 30),
        'reads_10kb_plus': sum(1 for l in lengths if l >= 10000),
        'reads_25kb_plus': sum(1 for l in lengths if l >= 25000)
    }

    # Calculate percentages
    total_reads = stats['total_reads_analyzed']
    stats['percent_q20_plus'] = (stats['reads_above_q20'] / total_reads) * 100
    stats['percent_q30_plus'] = (stats['reads_above_q30'] / total_reads) * 100
    stats['percent_10kb_plus'] = (stats['reads_10kb_plus'] / total_reads) * 100
    stats['percent_25kb_plus'] = (stats['reads_25kb_plus'] / total_reads) * 100

    return stats, lengths, mean_qualities


def classify_read_type(stats):
    """
    Classify reads as HiFi or CLR based on statistical characteristics
    """

    # Decision criteria
    hifi_score = 0
    clr_score = 0

    # Quality-based criteria
    if stats['mean_quality'] >= 25:
        hifi_score += 3
    elif stats['mean_quality'] >= 20:
        hifi_score += 2
    elif stats['mean_quality'] >= 15:
        hifi_score += 1
    else:
        clr_score += 2

    # High quality read percentage
    if stats['percent_q30_plus'] >= 70:
        hifi_score += 3
    elif stats['percent_q20_plus'] >= 80:
        hifi_score += 2
    elif stats['percent_q20_plus'] >= 50:
        hifi_score += 1
    else:
        clr_score += 2

    # Length distribution (HiFi tends to be more uniform)
    length_cv = stats['length_std'] / \
        stats['mean_length']  # Coefficient of variation
    if length_cv < 0.5:
        hifi_score += 2
    elif length_cv < 0.8:
        hifi_score += 1
    else:
        clr_score += 1

    # Typical length ranges
    if 8000 <= stats['median_length'] <= 25000:
        hifi_score += 2
    elif stats['median_length'] > 30000:
        clr_score += 1

    # Make decision
    if hifi_score > clr_score:
        return "HiFi", hifi_score, clr_score
    else:
        return "CLR", hifi_score, clr_score


def generate_plots(lengths, mean_qualities, output_prefix):
    """Generate diagnostic plots"""

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))

    # Length distribution
    ax1.hist(lengths, bins=50, alpha=0.7, edgecolor='black')
    ax1.set_xlabel('Read Length (bp)')
    ax1.set_ylabel('Frequency')
    ax1.set_title('Read Length Distribution')
    ax1.axvline(np.median(lengths), color='red', linestyle='--',
                label=f'Median: {np.median(lengths):.0f} bp')
    ax1.legend()

    # Quality distribution
    ax2.hist(mean_qualities, bins=50, alpha=0.7,
             edgecolor='black', color='orange')
    ax2.set_xlabel('Mean Quality Score')
    ax2.set_ylabel('Frequency')
    ax2.set_title('Quality Score Distribution')
    ax2.axvline(np.median(mean_qualities), color='red', linestyle='--',
                label=f'Median: Q{np.median(mean_qualities):.1f}')
    ax2.axvline(20, color='green', linestyle=':', label='Q20')
    ax2.axvline(30, color='blue', linestyle=':', label='Q30')
    ax2.legend()

    # Length vs Quality scatter
    sample_indices = np.random.choice(
        len(lengths), min(5000, len(lengths)), replace=False)
    sample_lengths = [lengths[i] for i in sample_indices]
    sample_qualities = [mean_qualities[i] for i in sample_indices]

    ax3.scatter(sample_lengths, sample_qualities, alpha=0.5, s=1)
    ax3.set_xlabel('Read Length (bp)')
    ax3.set_ylabel('Mean Quality Score')
    ax3.set_title('Length vs Quality Relationship')
    ax3.axhline(20, color='green', linestyle=':', alpha=0.7, label='Q20')
    ax3.axhline(30, color='blue', linestyle=':', alpha=0.7, label='Q30')
    ax3.legend()

    # Cumulative length distribution
    sorted_lengths = sorted(lengths, reverse=True)
    cumulative_lengths = np.cumsum(sorted_lengths)
    total_bases = cumulative_lengths[-1]

    ax4.plot(range(len(sorted_lengths)),
             cumulative_lengths / total_bases * 100)
    ax4.set_xlabel('Number of Reads')
    ax4.set_ylabel('Cumulative Bases (%)')
    ax4.set_title('Cumulative Length Distribution')
    ax4.grid(True, alpha=0.3)

    # Add N50 line
    n50_index = np.where(cumulative_lengths >= total_bases * 0.5)[0][0]
    n50_length = sorted_lengths[n50_index]
    ax4.axhline(50, color='red', linestyle='--',
                alpha=0.7, label=f'N50: {n50_length:,} bp')
    ax4.legend()

    plt.tight_layout()
    plt.savefig(f"{output_prefix}_read_analysis.png",
                dpi=300, bbox_inches='tight')
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Determine PacBio read type (HiFi vs CLR)')
    parser.add_argument('fastq_file', help='PacBio FASTQ file to analyze')
    parser.add_argument('--sample-size', type=int, default=10000,
                        help='Number of reads to analyze (default: 10000)')
    parser.add_argument('--output-prefix', default='pacbio_analysis',
                        help='Prefix for output files')
    parser.add_argument('--generate-plots', action='store_true',
                        help='Generate diagnostic plots')

    args = parser.parse_args()

    # Check if file exists
    if not Path(args.fastq_file).exists():
        print(f"Error: File {args.fastq_file} not found!")
        sys.exit(1)

    # Analyze reads
    result = analyze_pacbio_reads(args.fastq_file, args.sample_size)
    if result is None:
        sys.exit(1)

    stats, lengths, mean_qualities = result

    # Classify read type
    read_type, hifi_score, clr_score = classify_read_type(stats)

    # Print results
    print("\n" + "="*60)
    print("PacBio Read Type Analysis Results")
    print("="*60)
    print(f"File analyzed: {args.fastq_file}")
    print(f"Reads analyzed: {stats['total_reads_analyzed']:,}")
    print()

    print("Length Statistics:")
    print(f"  Mean length: {stats['mean_length']:,.0f} bp")
    print(f"  Median length: {stats['median_length']:,.0f} bp")
    print(
        f"  Length range: {stats['min_length']:,} - {stats['max_length']:,} bp")
    print(f"  Reads ≥10kb: {stats['percent_10kb_plus']:.1f}%")
    print(f"  Reads ≥25kb: {stats['percent_25kb_plus']:.1f}%")
    print()

    print("Quality Statistics:")
    print(f"  Mean quality: Q{stats['mean_quality']:.1f}")
    print(f"  Median quality: Q{stats['median_quality']:.1f}")
    print(
        f"  Quality range: Q{stats['min_quality']:.1f} - Q{stats['max_quality']:.1f}")
    print(f"  Reads ≥Q20: {stats['percent_q20_plus']:.1f}%")
    print(f"  Reads ≥Q30: {stats['percent_q30_plus']:.1f}%")
    print()

    print("Classification:")
    print(f"  Predicted read type: {read_type}")
    print(f"  Confidence scores - HiFi: {hifi_score}, CLR: {clr_score}")
    print()

    # Interpretation
    if read_type == "HiFi":
        print("Interpretation:")
        print("  ✓ These appear to be HiFi (CCS) reads")
        print("  ✓ High accuracy, suitable for high-quality assembly")
        print("  ✓ Recommend using --pacbio-hifi flag in metaFlye")
        print("  ✓ Minimal error correction needed")
    else:
        print("Interpretation:")
        print("  ⚠ These appear to be CLR (Continuous Long Reads)")
        print("  ⚠ Lower accuracy, will need more error correction")
        print("  ⚠ Recommend using --pacbio-raw flag in metaFlye")
        print("  ⚠ Consider more aggressive polishing")

    print()

    # Generate plots if requested
    if args.generate_plots:
        try:
            generate_plots(lengths, mean_qualities, args.output_prefix)
            print(
                f"Diagnostic plots saved as {args.output_prefix}_read_analysis.png")
        except ImportError:
            print("Warning: matplotlib not available, skipping plots")
        except Exception as e:
            print(f"Error generating plots: {e}")

    # Save detailed results
    with open(f"{args.output_prefix}_results.txt", 'w') as f:
        f.write("PacBio Read Analysis Results\n")
        f.write("="*30 + "\n")
        f.write(f"File: {args.fastq_file}\n")
        f.write(f"Predicted type: {read_type}\n")
        f.write(f"Confidence: HiFi={hifi_score}, CLR={clr_score}\n")
        f.write("\nDetailed Statistics:\n")
        for key, value in stats.items():
            f.write(f"{key}: {value}\n")

    print(f"Detailed results saved to {args.output_prefix}_results.txt")


if __name__ == "__main__":
    main()
