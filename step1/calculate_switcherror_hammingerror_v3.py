#!/usr/bin/env python3
import sys
from collections import defaultdict

# --- Script Description ---
# This script evaluates the accuracy of SNP phasing by comparing a predicted phasing file
# against a ground truth SNP set. It calculates two key metrics for phasing quality:
# 1.  Switch Error Rate: Measures the rate of incorrect long-range connections,
#     indicating how often the predicted haplotype switches from one true haplotype
#     to another along the chromosome.
# 2.  Hamming Error Rate: Measures the rate of local, single-point phasing errors,
#     reflecting the overall purity of a phased block.
#
# The script is designed to work with custom file formats as described in the functions below.

def parse_predicted_phasing(filepath):
    """
    Parses the predicted phasing file (e.g., from whatshap, processed).

    This function reads a file where each line represents a phased heterozygous SNP.
    It extracts the genotype (GT, e.g., "0|1") and a phase set identifier (PS),
    which groups SNPs that are phased together.

    Args:
        filepath (str): Path to the predicted phasing file.
                        Expected format: "CHROM POS REF ALT GT:AD:DP:GQ:PL:PS"

    Returns:
        dict: A dictionary mapping a unique SNP identifier ("CHROM_POS") to its
              phasing information ("GT:PS_ID").
    """
    phasing_db = {}
    print(f"INFO: Parsing predicted phasing file: {filepath}")
    with open(filepath, 'r') as f:
        for line_num, line in enumerate(f, 1):
            parts = line.strip().split()
            if len(parts) < 5:
                # Skip malformed lines
                continue
            
            chrom, pos, _, _, gt_info = parts
            key = f"{chrom}_{pos}"
            
            gt_parts = gt_info.split(':')
            gt = gt_parts[0]
            
            # The PS tag is typically the last field in the FORMAT string from whatshap
            ps_id = gt_parts[-1]
            
            if gt not in ["0|1", "1|0"]:
                # We only evaluate heterozygous SNPs that have been phased.
                continue
            
            phasing_db[key] = f"{gt}:{ps_id}"

    print(f"INFO: Parsed {len(phasing_db)} phased heterozygous variants from prediction file.")
    return phasing_db

def build_phase_blocks(ground_truth_filepath, phasing_db):
    """
    Constructs phase blocks for evaluation by linking the ground truth with predictions.

    This function reads a ground truth SNP file (e.g., from nucmer comparison of
    haplotype assemblies). For each SNP in the ground truth, it finds the
    corresponding prediction from `phasing_db` and groups them into phase blocks
    based on their PS identifier.

    Args:
        ground_truth_filepath (str): Path to the ground truth SNP file.
                                     Expected format: "REF_CHR REF_POS REF_BASE ALT_CHR ALT_POS ALT_BASE"
        phasing_db (dict): The dictionary of predicted phasings from `parse_predicted_phasing`.

    Returns:
        dict: A dictionary where keys are phase block IDs ("CHROM-PS_ID") and
              values are lists of SNP information dictionaries.
    """
    phase_blocks = defaultdict(list)
    print(f"INFO: Building phase blocks for evaluation using ground truth: {ground_truth_filepath}")
    
    with open(ground_truth_filepath, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 6:
                continue
            
            ref_chrom, ref_pos, true_ref_base, _, _, true_alt_base = parts
            key = f"{ref_chrom}_{ref_pos}"

            # Only consider SNPs that are present in the prediction set
            if key not in phasing_db:
                continue
            
            # Retrieve prediction info
            info = phasing_db[key]
            predicted_gt, ps_id = info.split(':')
            
            block_id = f"{ref_chrom}-{ps_id}"
            
            # Store the essential information for error calculation
            phase_blocks[block_id].append({
                "predicted_gt": predicted_gt,
                "true_ref": true_ref_base,
                "true_alt": true_alt_base,
                "position": int(ref_pos) # Store position for sorting
            })

    # It's crucial to sort SNPs within each block by position for correct switch error calculation
    for block_id in phase_blocks:
        phase_blocks[block_id].sort(key=lambda x: x["position"])

    print(f"INFO: Constructed {len(phase_blocks)} phase blocks containing comparable SNPs.")
    return phase_blocks

def standard_error_calculation(phase_blocks):
    """
    Calculates Switch and Hamming errors using standard, position-aware algorithms.

    Args:
        phase_blocks (dict): The dictionary of phase blocks to evaluate.
    """
    total_snps = 0
    total_switch_errors = 0
    total_hamming_errors = 0

    print("\n--- Detailed Block Statistics (Standard Method) ---")
    print(f"{'BLOCK-ID':<30} {'No_of_SNPs':>12} {'Switch_Errors':>15} {'Hamming_Errors':>18}")

    for block_id, variants in sorted(phase_blocks.items()):
        num_snps = len(variants)
        if num_snps == 0:
            continue

        # --- Hamming Error Calculation ---
        # Logic: Determine the majority phasing orientation for the entire block (0|1 vs 1|0).
        # Any SNP phased in the minority orientation is considered a Hamming Error.
        # This measures the overall "purity" or consistency of the block.
        votes_for_01 = sum(1 for v in variants if v['predicted_gt'] == '0|1')
        votes_for_10 = num_snps - votes_for_01
        
        # Determine the majority orientation for this block
        block_majority_is_01 = votes_for_01 >= votes_for_10
        
        num_hamming_errors_in_block = 0
        for snp in variants:
            snp_is_01 = (snp['predicted_gt'] == '0|1')
            if snp_is_01 != block_majority_is_01:
                num_hamming_errors_in_block += 1
        
        # --- Switch Error Calculation ---
        # Logic: Iterate through adjacent pairs of SNPs in the block. A switch error
        # occurs if the phasing orientation flips from one SNP to the next.
        # This measures the fragmentation or long-range continuity of the block.
        num_switch_errors_in_block = 0
        if num_snps > 1:
            for i in range(1, num_snps):
                prev_gt_is_01 = (variants[i-1]['predicted_gt'] == '0|1')
                curr_gt_is_01 = (variants[i]['predicted_gt'] == '0|1')
                if prev_gt_is_01 != curr_gt_is_01:
                    num_switch_errors_in_block += 1
        
        # Aggregate totals
        total_snps += num_snps
        total_switch_errors += num_switch_errors_in_block
        total_hamming_errors += num_hamming_errors_in_block
        
        print(f"{block_id:<30} {num_snps:>12} {num_switch_errors_in_block:>15} {num_hamming_errors_in_block:>18}")

    # --- Final Summary Report ---
    print("\n--- Final Summary (Standard Method) ---")
    if total_snps > 0:
        # Switch Error Rate is calculated over all possible switch locations.
        # For N SNPs in a block, there are N-1 possible locations for a switch.
        possible_switch_locations = total_snps - len(phase_blocks) if total_snps > len(phase_blocks) else 0
        final_switch_rate = (total_switch_errors / possible_switch_locations) * 100 if possible_switch_locations > 0 else 0.0
        
        # Hamming Error Rate is the simple fraction of erroneous sites.
        final_hamming_rate = (total_hamming_errors / total_snps) * 100

        print(f"Total Number of SNPs Compared: {total_snps}")
        print("-" * 50)
        print("Switch Error Metrics (Standard):")
        print(f"  Total Switches (flips between adjacent SNPs): {total_switch_errors}")
        print(f"  Switch Error Rate: {final_switch_rate:.2f}%  (Calculated as: switches / (total SNPs - total blocks))")
        print("-" * 50)
        print("Hamming Error Metrics (Standard):")
        print(f"  Total Hamming Errors (mismatches to block majority): {total_hamming_errors}")
        print(f"  Hamming Error Rate: {final_hamming_rate:.2f}% (Calculated as: errors / total SNPs)")
        print("-" * 50)
    else:
        print("No comparable SNPs found. Please check your input files and ensure chromosome names match.")

def main():
    """Main function to run the phasing evaluation."""
    if len(sys.argv) != 3:
        print("Usage: python evaluate_phasing_standard.py <predicted_phasing_file> <ground_truth_snp_file>")
        print("  - <predicted_phasing_file>: e.g., 'hifi.phase.txt'")
        print("  - <ground_truth_snp_file>:  e.g., 'hap1234.hapVar.snps'")
        sys.exit(1)
        
    predicted_file, truth_file = sys.argv[1], sys.argv[2]
    
    print("!!! IMPORTANT: Ensure chromosome names are consistent between input files !!!\n")
    
    try:
        # Step 1: Parse the predicted phasing data into a searchable format.
        predicted_data = parse_predicted_phasing(predicted_file)
        
        # Step 2: Build phase blocks by matching ground truth SNPs with predictions.
        block_data = build_phase_blocks(truth_file, predicted_data)
        
        # Step 3: Calculate and report errors for the constructed blocks.
        standard_error_calculation(block_data)
        
    except FileNotFoundError as e:
        print(f"ERROR: File not found - {e}")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
