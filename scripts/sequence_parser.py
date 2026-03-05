#!/usr/bin/env python3
"""
SARS-CoV-2 Sequence Parser
=========================
Extracts metadata from FASTA headers for surveillance analysis.
Enhanced with sequence gap handling and quality filtering.
"""
import re
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass

def clean_sequence(sequence: str) -> Tuple[str, Dict[str, float]]:
    """
    Clean sequence and calculate quality metrics
    Returns: (cleaned_sequence, quality_stats)
    """
    original_length = len(sequence)
    if original_length == 0:
        return "", {"completeness": 0.0, "gap_percentage": 100.0, "ambiguous_percentage": 100.0}
    
    # Convert to uppercase and remove whitespace
    cleaned = re.sub(r'\s+', '', sequence.upper())
    
    # Count different types of characters
    valid_bases = len(re.findall(r'[ATGC]', cleaned))
    gaps = len(re.findall(r'[N-]', cleaned))
    ambiguous = len(re.findall(r'[RYWSMKHBVD]', cleaned))
    
    # Calculate quality metrics
    quality_stats = {
        "completeness": valid_bases / original_length,
        "gap_percentage": gaps / original_length * 100,
        "ambiguous_percentage": ambiguous / original_length * 100,
        "original_length": original_length,
        "cleaned_length": len(cleaned)
    }
    
    # Replace ambiguous bases with N for consistency
    cleaned = re.sub(r'[RYWSMKHBVD]', 'N', cleaned)
    
    return cleaned, quality_stats

def is_sequence_usable(quality_stats: Dict[str, float], min_completeness: float = 0.80) -> bool:
    """
    Determine if sequence meets quality thresholds for surveillance analysis
    """
    return (
        quality_stats["completeness"] >= min_completeness and
        quality_stats["original_length"] >= 25000  # Minimum reasonable SARS-CoV-2 length
    )

@dataclass
class SequenceInfo:
    """Container for parsed sequence metadata"""
    accession: str
    title: str
    length: int
    collection_date: str
    organism: str
    pangolin_lineage: str
    country: Optional[str] = None
    
    # Enhanced with quality tracking
    quality_stats: Optional[Dict[str, float]] = None
    is_usable: bool = True
    
    def extract_country_from_title(self) -> Optional[str]:
        """Extract country code from virus title"""
        # Pattern: SARS-CoV-2/Country/Identifier/Year
        match = re.search(r'SARS-CoV-2/[^/]*/([^/]+)/', self.title)
        return match.group(1) if match else None

def parse_fasta_header(header: str) -> SequenceInfo:
    """Parse your custom NCBI FASTA header format"""
    parts = header.strip('>').split('|')
    
    if len(parts) < 6:
        raise ValueError(f"Invalid header format: {header}")
    
    seq_info = SequenceInfo(
        accession=parts[0].strip(),
        title=parts[1].strip(),
        length=int(parts[2]) if parts[2].strip().isdigit() else 0,
        collection_date=parts[3].strip(),
        organism=parts[4].strip(),
        pangolin_lineage=parts[5].strip()
    )
    
    # Extract country from title
    seq_info.country = seq_info.extract_country_from_title()
    return seq_info

def load_sequences(fasta_path: str, filter_low_quality: bool = True) -> List[Tuple[SequenceInfo, str]]:
    """Load and parse FASTA file with enhanced quality filtering"""
    sequences = []
    current_header = None
    current_seq = []
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Process previous sequence
                if current_header:
                    seq_info = parse_fasta_header(current_header)
                    raw_sequence = ''.join(current_seq)
                    
                    # Clean sequence and get quality stats
                    cleaned_sequence, quality_stats = clean_sequence(raw_sequence)
                    seq_info.quality_stats = quality_stats
                    seq_info.is_usable = is_sequence_usable(quality_stats)
                    
                    # Only add if usable (or if filtering is disabled)
                    if not filter_low_quality or seq_info.is_usable:
                        sequences.append((seq_info, cleaned_sequence))
                
                # Start new sequence
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
        
        # Process last sequence
        if current_header:
            seq_info = parse_fasta_header(current_header)
            raw_sequence = ''.join(current_seq)
            
            # Clean sequence and get quality stats
            cleaned_sequence, quality_stats = clean_sequence(raw_sequence)
            seq_info.quality_stats = quality_stats
            seq_info.is_usable = is_sequence_usable(quality_stats)
            
            # Only add if usable (or if filtering is disabled)
            if not filter_low_quality or seq_info.is_usable:
                sequences.append((seq_info, cleaned_sequence))
    
    return sequences

def get_sequence_summary(sequences: List[Tuple[SequenceInfo, str]]) -> Dict[str, float]:
    """Generate summary statistics for loaded sequences"""
    if not sequences:
        return {}
    
    total_completeness = sum(seq_info.quality_stats["completeness"] for seq_info, _ in sequences)
    total_gap_percentage = sum(seq_info.quality_stats["gap_percentage"] for seq_info, _ in sequences)
    
    return {
        "total_sequences": len(sequences),
        "avg_completeness": total_completeness / len(sequences),
        "avg_gap_percentage": total_gap_percentage / len(sequences),
        "usable_sequences": sum(1 for seq_info, _ in sequences if seq_info.is_usable)
    }

# Test function
def test_parser():
    """Test the parser with your data"""
    print("🧬 Testing enhanced SARS-CoV-2 sequence parser...")
    
    # Load with quality filtering
    sequences = load_sequences('../data/MAIN.fasta')
    summary = get_sequence_summary(sequences)
    
    print(f"✅ Loaded {summary['total_sequences']} usable sequences")
    print(f"📊 Average completeness: {summary['avg_completeness']:.2%}")
    print(f"🔍 Average gap percentage: {summary['avg_gap_percentage']:.1f}%")
    
    print("\n🔬 First 3 sequences:")
    
    for i, (seq_info, sequence) in enumerate(sequences[:3]):
        stats = seq_info.quality_stats
        print(f"\n{i+1}. {seq_info.accession}")
        print(f"   Lineage: {seq_info.pangolin_lineage}")
        print(f"   Country: {seq_info.country}")
        print(f"   Date: {seq_info.collection_date}")
        print(f"   Completeness: {stats['completeness']:.2%}")
        print(f"   Gap %: {stats['gap_percentage']:.1f}%")
        print(f"   Length: {stats['cleaned_length']} bp")
        print(f"   Sequence preview: {sequence[:50]}...")
        
        # Show gap handling example
        if 'N' in sequence[:100]:
            gap_example = sequence[:100]
            gap_positions = [i for i, base in enumerate(gap_example) if base == 'N']
            print(f"   Gap positions in first 100bp: {gap_positions[:10]}")

if __name__ == "__main__":
    test_parser()
