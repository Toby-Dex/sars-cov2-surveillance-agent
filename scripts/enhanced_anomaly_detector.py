#!/usr/bin/env python3
"""
Enhanced SARS-CoV-2 Anomaly Detection Engine
==========================================
Detects sequences with pandemic potential using advanced mutation analysis.
"""

import os
import numpy as np
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from sequence_parser import SequenceInfo, parse_fasta_header
from collections import defaultdict, Counter
import re

@dataclass
class AnomalyScore:
    """Container for anomaly detection results"""
    sequence_id: str
    lineage: str
    total_score: float
    recombinant_score: float
    geographic_score: float
    temporal_score: float
    genetic_score: float
    mutation_score: float  # NEW: specific mutation pattern score
    risk_level: str
    flags: List[str]

class MutationDetector:
    """Detects pandemic-associated mutation patterns"""
    
    def __init__(self):
        # Critical mutation patterns with weights
        self.high_risk_mutations = {
            'PRRA': 0.4,      # Furin cleavage site
            'P681H': 0.3,     # Omicron transmissibility
            'P681R': 0.3,     # Delta transmissibility
            'N501Y': 0.3,     # Enhanced ACE2 binding
            'E484K': 0.4,     # Immune escape (Beta)
            'E484A': 0.3,     # Immune escape (Omicron)
            'L452R': 0.3,     # Delta signature
            'Q493E': 0.5,     # XEC/KP.3 - very recent concern
            'R346T': 0.3,     # JN.1 immune escape
            'F456L': 0.3,     # Enhanced binding
            'K417N': 0.2,     # Beta signature
            'K417T': 0.2,     # Omicron signature
        }
        
        # Lineage-specific expected mutations
        self.lineage_signatures = {
            'B.1.1.7': ['N501Y', 'P681H'],           # Alpha
            'B.1.351': ['N501Y', 'E484K', 'K417N'],  # Beta
            'B.1.617.2': ['L452R', 'P681R'],         # Delta
            'B.1.1.529': ['N501Y', 'E484A', 'K417T'], # Omicron
            'XEC.2': ['Q493E', 'F456L'],             # XEC
            'JN.1': ['R346T', 'F456L'],              # JN.1
        }
    
    def extract_spike_region(self, aligned_sequence: str, reference: str) -> str:
        """Extract spike protein region for mutation analysis"""
        # Approximate spike boundaries (positions 21563-25384 in reference)
        spike_start_ratio = 21563 / 29903  # ~0.72
        spike_end_ratio = 25384 / 29903    # ~0.85
        
        seq_length = len(aligned_sequence)
        start_pos = int(seq_length * spike_start_ratio)
        end_pos = int(seq_length * spike_end_ratio)
        
        return aligned_sequence[start_pos:end_pos]
    
    def detect_mutation_patterns(self, sequence: str, lineage: str) -> Tuple[float, List[str]]:
        """Detect high-risk mutation patterns"""
        flags = []
        score = 0.0
        
        # Check for high-risk mutations in sequence
        mutations_found = []
        for mutation, weight in self.high_risk_mutations.items():
            if self.pattern_in_sequence(mutation, sequence):
                score += weight
                mutations_found.append(mutation)
                flags.append(f"High-risk mutation detected: {mutation}")
        
        # Check lineage consistency
        if lineage in self.lineage_signatures:
            consistency_score, consistency_flags = self.check_lineage_consistency(
                sequence, lineage, mutations_found
            )
            score += consistency_score
            flags.extend(consistency_flags)
        
        # Bonus for multiple concerning mutations
        if len(mutations_found) > 2:
            bonus = 0.2 * (len(mutations_found) - 2)
            score += bonus
            flags.append(f"Multiple high-risk mutations: {len(mutations_found)} detected")
        
        # Special concern for very recent mutations
        recent_mutations = ['Q493E', 'R346T', 'F456L']
        recent_found = [m for m in mutations_found if m in recent_mutations]
        if recent_found:
            score += 0.3
            flags.append(f"Recent high-concern mutations: {', '.join(recent_found)}")
        
        return min(score, 1.0), flags
    
    def pattern_in_sequence(self, pattern: str, sequence: str) -> bool:
        """Check if mutation pattern exists in sequence"""
        # Simplified pattern matching - looks for amino acid patterns
        # In production, you'd do proper nucleotide->amino acid translation
        
        # Remove gaps and check for pattern
        clean_seq = sequence.replace('-', '').upper()
        pattern_upper = pattern.upper()
        
        # Direct pattern match
        if pattern_upper in clean_seq:
            return True
        
        # Check for nucleotide patterns that code for these amino acids
        nucleotide_patterns = {
            'PRRA': ['CCGCGAGCA', 'CCTCGGGCA', 'CCCCGGGCA'],  # Various PRRA codons
            'P681H': ['CCTCAT', 'CCTCAC'],  # P->H at 681
            'P681R': ['CCTCGT', 'CCTAGA'],  # P->R at 681
            'Q493E': ['CAAGAG', 'CAGGAG'],  # Q->E at 493
        }
        
        if pattern in nucleotide_patterns:
            for nt_pattern in nucleotide_patterns[pattern]:
                if nt_pattern in clean_seq:
                    return True
        
        return False
    
    def check_lineage_consistency(self, sequence: str, lineage: str, 
                                mutations_found: List[str]) -> Tuple[float, List[str]]:
        """Check if mutations match expected lineage signatures"""
        flags = []
        score = 0.0
        
        if lineage not in self.lineage_signatures:
            return 0.0, flags
        
        expected = set(self.lineage_signatures[lineage])
        found = set(mutations_found)
        
        missing = expected - found
        unexpected = found - expected
        
        # Flag missing expected mutations
        if missing:
            score += 0.2 * len(missing)
            flags.append(f"Missing {lineage} signatures: {', '.join(missing)}")
        
        # Flag unexpected mutations for this lineage
        if unexpected:
            score += 0.3
            flags.append(f"Unexpected mutations in {lineage}: {', '.join(list(unexpected)[:3])}")
        
        return min(score, 0.5), flags

class EnhancedAnomalyDetector:
    """Enhanced anomaly detector with mutation pattern analysis"""
    
    def __init__(self, aligned_fasta_path: str):
        self.aligned_fasta_path = aligned_fasta_path
        self.sequences = self.load_aligned_sequences()
        self.reference_genome = self.extract_reference()
        self.mutation_detector = MutationDetector()
        
        # Define pandemic-risk patterns
        self.recombinant_lineages = {'XEC', 'XEC.1', 'XEC.2', 'XEC.8', 'XFG'}
        self.high_transmissibility_lineages = {'JN.1', 'KP.2', 'KP.3', 'BA.2.86'}
        self.established_lineages = {'B.1.1.7', 'B.1.617.2', 'B.1.1.529', 'B.1.351'}
        
        print(f"Loaded {len(self.sequences)} aligned sequences for enhanced analysis")
    
    def load_aligned_sequences(self) -> Dict[str, Tuple[SequenceInfo, str]]:
        """Load aligned sequences into memory"""
        sequences = {}
        current_header = None
        current_seq = []
        
        with open(self.aligned_fasta_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Process previous sequence
                    if current_header and not current_header.startswith('>NC_045512'):
                        header_clean = current_header.strip('>')
                        if '|' in header_clean:
                            parts = header_clean.split('|')
                            accession = parts[0]
                            lineage = parts[1] if len(parts) > 1 else "Unknown"
                            
                            # Create SequenceInfo with simplified data
                            seq_info = SequenceInfo(
                                accession=accession,
                                title=f"SARS-CoV-2 {lineage}",
                                length=len(''.join(current_seq)),
                                collection_date="2024",
                                organism="SARS-CoV-2",
                                pangolin_lineage=lineage
                            )
                            sequence = ''.join(current_seq)
                            sequences[seq_info.accession] = (seq_info, sequence)
                    
                    # Start new sequence
                    current_header = line
                    current_seq = []
                else:
                    current_seq.append(line)
            
            # Process last sequence
            if current_header and not current_header.startswith('>NC_045512'):
                header_clean = current_header.strip('>')
                if '|' in header_clean:
                    parts = header_clean.split('|')
                    accession = parts[0]
                    lineage = parts[1] if len(parts) > 1 else "Unknown"
                    
                    seq_info = SequenceInfo(
                        accession=accession,
                        title=f"SARS-CoV-2 {lineage}",
                        length=len(''.join(current_seq)),
                        collection_date="2024",
                        organism="SARS-CoV-2",
                        pangolin_lineage=lineage
                    )
                    sequence = ''.join(current_seq)
                    sequences[seq_info.accession] = (seq_info, sequence)
        
        return sequences
    
    def extract_reference(self) -> str:
        """Extract reference genome from aligned file"""
        with open(self.aligned_fasta_path, 'r') as f:
            in_reference = False
            reference_seq = []
            
            for line in f:
                line = line.strip()
                if line.startswith('>NC_045512.2'):
                    in_reference = True
                    continue
                elif line.startswith('>') and in_reference:
                    break
                elif in_reference:
                    reference_seq.append(line)
            
            return ''.join(reference_seq)
    
    def calculate_genetic_distance(self, seq1: str, seq2: str) -> float:
        """Calculate genetic distance between two aligned sequences"""
        if len(seq1) != len(seq2):
            return 1.0
        
        differences = sum(1 for a, b in zip(seq1, seq2) if a != b and a != '-' and b != '-')
        valid_positions = sum(1 for a, b in zip(seq1, seq2) if a != '-' and b != '-')
        
        return differences / valid_positions if valid_positions > 0 else 1.0
    
    def detect_recombinant_signatures(self, seq_info: SequenceInfo, sequence: str) -> Tuple[float, List[str]]:
        """Detect recombination signatures"""
        flags = []
        score = 0.0
        
        # Check if lineage is known recombinant
        if any(rec in seq_info.pangolin_lineage for rec in self.recombinant_lineages):
            score += 0.8
            flags.append(f"Known recombinant lineage: {seq_info.pangolin_lineage}")
        
        # Check for unusual gap patterns
        leading_gaps = len(sequence) - len(sequence.lstrip('-'))
        if leading_gaps > 50:
            score += 0.3
            flags.append(f"Unusual insertion pattern: {leading_gaps} leading gaps")
        
        return min(score, 1.0), flags
    
    def calculate_composite_score(self, seq_info: SequenceInfo, sequence: str) -> AnomalyScore:
        """Calculate enhanced anomaly score with mutation analysis"""
        
        # Individual component scores
        recombinant_score, recomb_flags = self.detect_recombinant_signatures(seq_info, sequence)
        
        # NEW: Enhanced mutation pattern detection
        mutation_score, mutation_flags = self.mutation_detector.detect_mutation_patterns(
            sequence, seq_info.pangolin_lineage
        )
        
        # Genetic distance from reference
        genetic_score = self.calculate_genetic_distance(sequence, self.reference_genome)
        
        # Geographic and temporal (simplified for now)
        geographic_score = 0.0
        temporal_score = 0.0
        
        # Enhanced weighted composite score
        weights = {
            'recombinant': 0.3,   # Reduced to make room for mutation analysis
            'mutation': 0.4,      # NEW: Highest weight for specific mutations
            'genetic': 0.2,
            'geographic': 0.05,
            'temporal': 0.05
        }
        
        total_score = (
            recombinant_score * weights['recombinant'] +
            mutation_score * weights['mutation'] +
            genetic_score * weights['genetic'] +
            geographic_score * weights['geographic'] +
            temporal_score * weights['temporal']
        )
        
        # Enhanced risk level determination
        if total_score >= 0.8:
            risk_level = "CRITICAL"
        elif total_score >= 0.6:
            risk_level = "HIGH"
        elif total_score >= 0.4:
            risk_level = "MEDIUM"
        else:
            risk_level = "LOW"
        
        all_flags = recomb_flags + mutation_flags
        
        return AnomalyScore(
            sequence_id=seq_info.accession,
            lineage=seq_info.pangolin_lineage,
            total_score=total_score,
            recombinant_score=recombinant_score,
            geographic_score=geographic_score,
            temporal_score=temporal_score,
            genetic_score=genetic_score,
            mutation_score=mutation_score,  # NEW
            risk_level=risk_level,
            flags=all_flags
        )
    
    def run_analysis(self) -> List[AnomalyScore]:
        """Run enhanced anomaly detection analysis"""
        results = []
        
        print("Running enhanced anomaly detection with mutation analysis...")
        
        for accession, (seq_info, sequence) in self.sequences.items():
            score = self.calculate_composite_score(seq_info, sequence)
            results.append(score)
        
        # Sort by total score (highest risk first)
        results.sort(key=lambda x: x.total_score, reverse=True)
        
        return results
    
    def generate_enhanced_report(self, results: List[AnomalyScore], output_path: str):
        """Generate enhanced surveillance report"""
        
        critical_risk = [r for r in results if r.risk_level == "CRITICAL"]
        high_risk = [r for r in results if r.risk_level == "HIGH"]
        medium_risk = [r for r in results if r.risk_level == "MEDIUM"]
        
        with open(output_path, 'w') as f:
            f.write("Enhanced SARS-CoV-2 Variant Surveillance Report\n")
            f.write("="*55 + "\n\n")
            
            f.write(f"Analysis Summary:\n")
            f.write(f"- Total sequences analyzed: {len(results)}\n")
            f.write(f"- CRITICAL risk sequences: {len(critical_risk)}\n")
            f.write(f"- HIGH risk sequences: {len(high_risk)}\n")
            f.write(f"- MEDIUM risk sequences: {len(medium_risk)}\n")
            f.write(f"- LOW risk sequences: {len(results) - len(critical_risk) - len(high_risk) - len(medium_risk)}\n\n")
            
            # Critical and High risk sequences
            for risk_category, risk_list in [("CRITICAL", critical_risk), ("HIGH", high_risk)]:
                if risk_list:
                    f.write(f"{risk_category} RISK SEQUENCES:\n")
                    f.write("-" * (len(risk_category) + 16) + "\n")
                    for result in risk_list[:10]:  # Top 10 per category
                        f.write(f"\nSequence: {result.sequence_id}\n")
                        f.write(f"Lineage: {result.lineage}\n")
                        f.write(f"Total Score: {result.total_score:.3f}\n")
                        f.write(f"Components: R={result.recombinant_score:.2f}, ")
                        f.write(f"M={result.mutation_score:.2f}, ")
                        f.write(f"G={result.genetic_score:.2f}\n")
                        if result.flags:
                            f.write(f"Flags: {'; '.join(result.flags)}\n")
                    f.write("\n")
            
            # Top 15 overall
            f.write(f"TOP 15 HIGHEST SCORING SEQUENCES:\n")
            f.write("-" * 35 + "\n")
            for i, result in enumerate(results[:15]):
                f.write(f"{i+1:2d}. {result.sequence_id} ({result.lineage}) - ")
                f.write(f"Score: {result.total_score:.3f} - {result.risk_level}\n")
        
        print(f"Enhanced report saved to: {output_path}")

def test_enhanced_detection():
    """Test enhanced anomaly detection"""
    
    # Initialize enhanced detector
    detector = EnhancedAnomalyDetector('../outputs/alignments/aligned_sequences.fasta')
    
    # Run enhanced analysis
    results = detector.run_analysis()
    
    # Generate enhanced report
    os.makedirs('../outputs/reports', exist_ok=True)
    detector.generate_enhanced_report(results, '../outputs/reports/enhanced_surveillance_report.txt')
    
    # Print summary
    print(f"\nEnhanced Analysis Complete!")
    print(f"Total sequences: {len(results)}")
    print(f"Critical risk: {len([r for r in results if r.risk_level == 'CRITICAL'])}")
    print(f"High risk: {len([r for r in results if r.risk_level == 'HIGH'])}")
    print(f"Medium risk: {len([r for r in results if r.risk_level == 'MEDIUM'])}")
    
    print("\nTop 5 highest risk sequences:")
    for i, result in enumerate(results[:5]):
        print(f"{i+1}. {result.sequence_id} ({result.lineage})")
        print(f"   Score: {result.total_score:.3f} (Mutation: {result.mutation_score:.3f})")
        if result.flags:
            print(f"   Key flags: {result.flags[0] if result.flags else 'None'}")

if __name__ == "__main__":
    test_enhanced_detection()
