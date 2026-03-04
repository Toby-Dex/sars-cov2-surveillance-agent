#!/usr/bin/env python3
"""
SARS-CoV-2 Sequence Parser
=========================
Extracts metadata from FASTA headers for surveillance analysis.
"""

import re
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass

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

def load_sequences(fasta_path: str) -> List[Tuple[SequenceInfo, str]]:
    """Load and parse FASTA file"""
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
                    sequence = ''.join(current_seq)
                    sequences.append((seq_info, sequence))
                
                # Start new sequence
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
        
        # Process last sequence
        if current_header:
            seq_info = parse_fasta_header(current_header)
            sequence = ''.join(current_seq)
            sequences.append((seq_info, sequence))
    
    return sequences

# Test function
def test_parser():
    """Test the parser with your data"""
    sequences = load_sequences('../data/MAIN.fasta')
    
    print(f"Loaded {len(sequences)} sequences")
    print("\nFirst 3 sequences:")
    
    for i, (seq_info, sequence) in enumerate(sequences[:3]):
        print(f"\n{i+1}. {seq_info.accession}")
        print(f"   Lineage: {seq_info.pangolin_lineage}")
        print(f"   Country: {seq_info.country}")
        print(f"   Date: {seq_info.collection_date}")
        print(f"   Length: {seq_info.length}")
        print(f"   Sequence preview: {sequence[:50]}...")

if __name__ == "__main__":
    test_parser()
