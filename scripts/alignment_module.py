#!/usr/bin/env python3
"""
SARS-CoV-2 Alignment Module
==========================
Aligns sequences against WIV04 reference using MAFFT.
"""

import os
import subprocess
import tempfile
from typing import List, Tuple
from sequence_parser import load_sequences, SequenceInfo

class AlignmentEngine:
    """Handles sequence alignment using MAFFT"""
    
    def __init__(self, reference_path: str):
        self.reference_path = reference_path
        self.verify_mafft()
    
    def verify_mafft(self):
        """Check if MAFFT is installed"""
        try:
            result = subprocess.run(['mafft', '--version'], 
                                  capture_output=True, text=True)
            print(f"MAFFT found: {result.stderr.strip()}")
        except FileNotFoundError:
            raise RuntimeError("MAFFT not found. Please install MAFFT.")
    
    def create_combined_fasta(self, sequences: List[Tuple[SequenceInfo, str]], 
                             output_path: str):
        """Combine reference + new sequences for alignment"""
        with open(output_path, 'w') as f:
            # Write reference first
            with open(self.reference_path, 'r') as ref:
                f.write(ref.read())
                if not ref.read().endswith('\n'):
                    f.write('\n')
            
            # Write new sequences
            for seq_info, sequence in sequences:
                f.write(f">{seq_info.accession}|{seq_info.pangolin_lineage}\n")
                f.write(f"{sequence}\n")
    
    def run_alignment(self, input_fasta: str, output_fasta: str, 
                     algorithm: str = "auto") -> bool:
        """Run MAFFT alignment"""
        
        # Choose MAFFT algorithm based on sequence count
        if algorithm == "auto":
            # For 228 sequences, use --auto (fast and accurate)
            mafft_cmd = ['mafft', '--auto', '--quiet', input_fasta]
        else:
            mafft_cmd = ['mafft', algorithm, '--quiet', input_fasta]
        
        try:
            print(f"Running MAFFT alignment...")
            with open(output_fasta, 'w') as out_file:
                result = subprocess.run(mafft_cmd, stdout=out_file, 
                                      stderr=subprocess.PIPE, text=True)
            
            if result.returncode != 0:
                print(f"MAFFT error: {result.stderr}")
                return False
            
            print(f"Alignment completed: {output_fasta}")
            return True
            
        except Exception as e:
            print(f"Alignment failed: {e}")
            return False
    
    def align_sequences(self, sequences: List[Tuple[SequenceInfo, str]], 
                       output_dir: str) -> str:
        """Main alignment workflow"""
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # File paths
        combined_fasta = os.path.join(output_dir, "combined_sequences.fasta")
        aligned_fasta = os.path.join(output_dir, "aligned_sequences.fasta")
        
        # Step 1: Combine reference + new sequences
        print(f"Combining {len(sequences)} sequences with reference...")
        self.create_combined_fasta(sequences, combined_fasta)
        
        # Step 2: Run alignment
        success = self.run_alignment(combined_fasta, aligned_fasta)
        
        if success:
            self.validate_alignment(aligned_fasta, len(sequences) + 1)  # +1 for reference
            return aligned_fasta
        else:
            raise RuntimeError("Alignment failed")
    
    def validate_alignment(self, aligned_fasta: str, expected_count: int):
        """Validate alignment output"""
        seq_count = 0
        alignment_length = 0
        
        with open(aligned_fasta, 'r') as f:
            current_seq = ""
            for line in f:
                if line.startswith('>'):
                    if current_seq and alignment_length == 0:
                        alignment_length = len(current_seq.replace('-', ''))
                    seq_count += 1
                    current_seq = ""
                else:
                    current_seq += line.strip()
        
        print(f"Alignment validation:")
        print(f"  - Sequences: {seq_count}/{expected_count}")
        print(f"  - Reference length: ~{alignment_length} bp")
        
        if seq_count != expected_count:
            print(f"Warning: Expected {expected_count} sequences, got {seq_count}")

def test_alignment():
    """Test alignment with your data"""
    
    # Load sequences
    sequences = load_sequences('../data/MAIN.fasta')
    print(f"Loaded {len(sequences)} sequences for alignment")
    
    # Initialize alignment engine
    aligner = AlignmentEngine('../temp_ref.fasta')
    
    # Run alignment (this will take a few minutes)
    try:
        aligned_file = aligner.align_sequences(sequences, '../outputs/alignments')
        print(f"SUCCESS: Aligned sequences saved to {aligned_file}")
        
        # Quick preview
        print("\nFirst few lines of alignment:")
        with open(aligned_file, 'r') as f:
            for i, line in enumerate(f):
                if i < 5:
                    print(line.strip())
                else:
                    break
                    
    except Exception as e:
        print(f"FAILED: {e}")

if __name__ == "__main__":
    test_alignment()
