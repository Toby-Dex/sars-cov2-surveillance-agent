#!/usr/bin/env python3
"""
SARS-CoV-2 Anomaly Detection Engine
==================================
Detects sequences with pandemic potential based on multiple criteria.
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
    risk_level: str
    flags: List[str]

class AnomalyDetector:
    """Detects sequences with pandemic potential"""
    
    def __init__(self, aligned_fasta_path: str):
        self.aligned_fasta_path = aligned_fasta_path
        self.sequences = self.load_aligned_sequences()
        self.reference_genome = self.extract_reference()
        
        # Define pandemic-risk patterns
        self.recombinant_lineages = {'XEC', 'XEC.1', 'XEC.2', 'XEC.8', 'XFG'}
        self.high_transmissibility_lineages = {'JN.1', 'KP.2', 'KP.3', 'BA.2.86'}
        self.established_lineages = {'B.1.1.7', 'B.1.617.2', 'B.1.1.529', 'B.1.351'}
        
        print(f"Loaded {len(self.sequences)} aligned sequences for analysis")
    
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
            return 1.0  # Maximum distance for length mismatch
        
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
        
        # Check for unusual gap patterns (like XEC.2 leading gaps we saw)
        leading_gaps = len(sequence) - len(sequence.lstrip('-'))
        if leading_gaps > 50:  # Unusual insertion pattern
            score += 0.3
            flags.append(f"Unusual insertion pattern: {leading_gaps} leading gaps")
        
        # Check for internal gap clusters (potential recombination breakpoints)
        gap_clusters = self.find_gap_clusters(sequence)
        if len(gap_clusters) > 3:
            score += 0.2 * len(gap_clusters)
            flags.append(f"Multiple gap clusters: {len(gap_clusters)} regions")
        
        return min(score, 1.0), flags
    
    def find_gap_clusters(self, sequence: str) -> List[Tuple[int, int]]:
        """Find clusters of gaps in sequence"""
        clusters = []
        in_cluster = False
        cluster_start = 0
        
        for i, char in enumerate(sequence):
            if char == '-':
                if not in_cluster:
                    cluster_start = i
                    in_cluster = True
            else:
                if in_cluster and i - cluster_start > 5:  # Minimum cluster size
                    clusters.append((cluster_start, i))
                in_cluster = False
        
        return clusters
    
    def analyze_geographic_anomaly(self, seq_info: SequenceInfo) -> Tuple[float, List[str]]:
        """Detect geographic clustering anomalies"""
        flags = []
        score = 0.0
        
        if not seq_info.country:
            return 0.0, flags
        
        # Analyze lineage-country combinations
        lineage_countries = defaultdict(set)
        for _, (other_info, _) in self.sequences.items():
            if other_info.country and other_info.pangolin_lineage:
                lineage_countries[other_info.pangolin_lineage].add(other_info.country)
        
        # Check if this lineage-country combo is unusual
        same_lineage_countries = lineage_countries.get(seq_info.pangolin_lineage, set())
        
        if len(same_lineage_countries) > 0 and seq_info.country not in same_lineage_countries:
            # This lineage appearing in new geographic region
            score += 0.4
            flags.append(f"Lineage {seq_info.pangolin_lineage} in new region: {seq_info.country}")
        
        # Check for high-risk lineages appearing globally
        if (seq_info.pangolin_lineage in self.recombinant_lineages and 
            len(same_lineage_countries) > 5):
            score += 0.3
            flags.append(f"Recombinant lineage with global spread")
        
        return min(score, 1.0), flags
    
    def analyze_temporal_anomaly(self, seq_info: SequenceInfo, sequence: str) -> Tuple[float, List[str]]:
        """Detect temporal evolution inconsistencies"""
        flags = []
        score = 0.0
        
        # Calculate genetic distance from reference
        genetic_dist = self.calculate_genetic_distance(sequence, self.reference_genome)
        
        # Parse collection date
        try:
            if '-' in seq_info.collection_date:
                year_month = seq_info.collection_date.split('-')
                year = int(year_month[0])
                month = int(year_month[1]) if len(year_month) > 1 else 6
            else:
                year = int(seq_info.collection_date)
                month = 6
        except:
            return 0.0, flags
        
        # Expected evolutionary rate (roughly 2e-3 substitutions per site per year)
        years_since_2020 = year + (month / 12.0) - 2020
        expected_distance = years_since_2020 * 0.002
        
        # Flag if genetic distance is much higher than expected
        if genetic_dist > expected_distance * 2:
            score += 0.5
            flags.append(f"Accelerated evolution: {genetic_dist:.4f} vs expected {expected_distance:.4f}")
        
        # Flag if very recent but highly divergent
        if year >= 2024 and genetic_dist > 0.01:
            score += 0.3
            flags.append(f"Recent highly divergent sequence: {year}")
        
        return min(score, 1.0), flags
    
    def calculate_composite_score(self, seq_info: SequenceInfo, sequence: str) -> AnomalyScore:
        """Calculate overall anomaly score"""
        
        # Individual component scores
        recombinant_score, recomb_flags = self.detect_recombinant_signatures(seq_info, sequence)
        geographic_score, geo_flags = self.analyze_geographic_anomaly(seq_info)
        temporal_score, temp_flags = self.analyze_temporal_anomaly(seq_info, sequence)
        
        # Genetic distance from reference
        genetic_score = self.calculate_genetic_distance(sequence, self.reference_genome)
        
        # Weighted composite score
        weights = {
            'recombinant': 0.4,  # Highest weight - recombinants are key
            'geographic': 0.2,   
            'temporal': 0.2,
            'genetic': 0.2
        }
        
        total_score = (
            recombinant_score * weights['recombinant'] +
            geographic_score * weights['geographic'] +
            temporal_score * weights['temporal'] +
            genetic_score * weights['genetic']
        )
        
        # Determine risk level
        if total_score >= 0.7:
            risk_level = "HIGH"
        elif total_score >= 0.4:
            risk_level = "MEDIUM"
        else:
            risk_level = "LOW"
        
        all_flags = recomb_flags + geo_flags + temp_flags
        
        return AnomalyScore(
            sequence_id=seq_info.accession,
            lineage=seq_info.pangolin_lineage,
            total_score=total_score,
            recombinant_score=recombinant_score,
            geographic_score=geographic_score,
            temporal_score=temporal_score,
            genetic_score=genetic_score,
            risk_level=risk_level,
            flags=all_flags
        )
    
    def run_analysis(self) -> List[AnomalyScore]:
        """Run complete anomaly detection analysis"""
        results = []
        
        print("Running anomaly detection analysis...")
        
        for accession, (seq_info, sequence) in self.sequences.items():
            score = self.calculate_composite_score(seq_info, sequence)
            results.append(score)
        
        # Sort by total score (highest risk first)
        results.sort(key=lambda x: x.total_score, reverse=True)
        
        return results
    
    def generate_report(self, results: List[AnomalyScore], output_path: str):
        """Generate surveillance report"""
        
        high_risk = [r for r in results if r.risk_level == "HIGH"]
        medium_risk = [r for r in results if r.risk_level == "MEDIUM"]
        
        with open(output_path, 'w') as f:
            f.write("SARS-CoV-2 Variant Surveillance Report\n")
            f.write("="*50 + "\n\n")
            
            f.write(f"Analysis Summary:\n")
            f.write(f"- Total sequences analyzed: {len(results)}\n")
            f.write(f"- HIGH risk sequences: {len(high_risk)}\n")
            f.write(f"- MEDIUM risk sequences: {len(medium_risk)}\n")
            f.write(f"- LOW risk sequences: {len(results) - len(high_risk) - len(medium_risk)}\n\n")
            
            # High risk sequences
            if high_risk:
                f.write("HIGH RISK SEQUENCES:\n")
                f.write("-" * 30 + "\n")
                for result in high_risk:
                    f.write(f"\nSequence: {result.sequence_id}\n")
                    f.write(f"Lineage: {result.lineage}\n")
                    f.write(f"Total Score: {result.total_score:.3f}\n")
                    f.write(f"Components: R={result.recombinant_score:.2f}, ")
                    f.write(f"G={result.geographic_score:.2f}, ")
                    f.write(f"T={result.temporal_score:.2f}, ")
                    f.write(f"D={result.genetic_score:.2f}\n")
                    if result.flags:
                        f.write(f"Flags: {'; '.join(result.flags)}\n")
            
            # Top 10 overall
            f.write(f"\n\nTOP 10 HIGHEST SCORING SEQUENCES:\n")
            f.write("-" * 35 + "\n")
            for i, result in enumerate(results[:10]):
                f.write(f"{i+1:2d}. {result.sequence_id} ({result.lineage}) - ")
                f.write(f"Score: {result.total_score:.3f} - {result.risk_level}\n")
        
        print(f"Report saved to: {output_path}")

def test_anomaly_detection():
    """Test anomaly detection on aligned sequences"""
    
    # Initialize detector
    detector = AnomalyDetector('../outputs/alignments/aligned_sequences.fasta')
    
    # Run analysis
    results = detector.run_analysis()
    
    # Generate report
    os.makedirs('../outputs/reports', exist_ok=True)
    detector.generate_report(results, '../outputs/reports/surveillance_report.txt')
    
    # Print summary
    print(f"\nAnalysis complete!")
    print(f"Total sequences: {len(results)}")
    print(f"High risk: {len([r for r in results if r.risk_level == 'HIGH'])}")
    print(f"Medium risk: {len([r for r in results if r.risk_level == 'MEDIUM'])}")
    
    print("\nTop 5 highest risk sequences:")
    for i, result in enumerate(results[:5]):
        print(f"{i+1}. {result.sequence_id} ({result.lineage}) - Score: {result.total_score:.3f}")

if __name__ == "__main__":
    test_anomaly_detection()
