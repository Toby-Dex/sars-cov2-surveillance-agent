#!/usr/bin/env python3
"""
SARS-CoV-2 Sequence Parser (Enhanced)
=====================================
Extracts metadata from FASTA headers for surveillance analysis.
Enhanced with sequence gap handling, quality filtering, and flexible header parsing.
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
    """
    Parse FASTA header with flexible format support
    
    Supported formats:
    1. Standard 6-field: >AccessionID|Title|Length|Date|Organism|Lineage
    2. Simple 2-field: >SequenceName|Lineage  
    3. GISAID format: >EPI_ISL_123456|hCoV-19/Country/Sample/2024
    4. Basic format: >SequenceName
    """
    parts = header.strip('>').split('|')
    
    # Handle different header formats
    if len(parts) >= 6:
        # Standard 6-field format (your original format)
        seq_info = SequenceInfo(
            accession=parts[0].strip(),
            title=parts[1].strip(),
            length=int(parts[2]) if parts[2].strip().isdigit() else 0,
            collection_date=parts[3].strip(),
            organism=parts[4].strip(),
            pangolin_lineage=parts[5].strip()
        )
    
    elif len(parts) == 2:
        # Simple format: >SequenceName|Lineage
        seq_info = SequenceInfo(
            accession=parts[0].strip(),
            title="",
            length=0,
            collection_date="Unknown",
            organism="SARS-CoV-2",
            pangolin_lineage=parts[1].strip()
        )
    
    elif len(parts) >= 3:
        # Try to parse as partial format
        accession = parts[0].strip()
        title = parts[1].strip() if len(parts) > 1 else ""
        
        # Try to extract lineage from title or use "Unknown"
        lineage = extract_lineage_from_text(title) or "Unknown"
        
        seq_info = SequenceInfo(
            accession=accession,
            title=title,
            length=int(parts[2]) if len(parts) > 2 and parts[2].strip().isdigit() else 0,
            collection_date=parts[3].strip() if len(parts) > 3 else "Unknown",
            organism=parts[4].strip() if len(parts) > 4 else "SARS-CoV-2",
            pangolin_lineage=lineage
        )
    
    else:
        # Single field format: >SequenceName
        accession = parts[0].strip()
        lineage = extract_lineage_from_text(accession) or "Unknown"
        
        seq_info = SequenceInfo(
            accession=accession,
            title="",
            length=0,
            collection_date="Unknown",
            organism="SARS-CoV-2",
            pangolin_lineage=lineage
        )
    
    # Extract country from title if available
    if seq_info.title:
        seq_info.country = seq_info.extract_country_from_title()
    
    return seq_info

def extract_lineage_from_text(text: str) -> Optional[str]:
    """Extract SARS-CoV-2 lineage from text using common patterns"""
    if not text:
        return None
    
    # Common lineage patterns
    patterns = [
        r'\b([A-Z]+(?:\.[0-9]+)+)\b',  # B.1.1.7, XBB.1.5, etc.
        r'\b([A-Z]{2,3}(?:\.[0-9]+)*)\b',  # XEC, JN.1, etc.
        r'\b(Alpha|Beta|Gamma|Delta|Omicron)\b'  # Greek variant names
    ]
    
    for pattern in patterns:
        match = re.search(pattern, text, re.IGNORECASE)
        if match:
            candidate = match.group(1)
            # Validate that this looks like a real lineage
            if is_valid_lineage(candidate):
                return candidate
    
    return None

def is_valid_lineage(lineage: str) -> bool:
    """Check if string looks like a valid SARS-CoV-2 lineage"""
    if not lineage:
        return False
    
    # Known valid patterns
    valid_patterns = [
        r'^[A-Z]+(\.[0-9]+)+$',  # B.1.1.7 style
        r'^[A-Z]{2,3}(\.[0-9]+)*$',  # XBB, XEC style
        r'^(Alpha|Beta|Gamma|Delta|Omicron)$'  # Greek names
    ]
    
    return any(re.match(pattern, lineage, re.IGNORECASE) for pattern in valid_patterns)

def load_sequences(fasta_path: str, filter_low_quality: bool = True, min_completeness: float = 0.80) -> List[Tuple[SequenceInfo, str]]:
    """Load and parse FASTA file with enhanced quality filtering and flexible header parsing"""
    sequences = []
    current_header = None
    current_seq = []
    
    try:
        with open(fasta_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if line.startswith('>'):
                    # Process previous sequence
                    if current_header:
                        try:
                            seq_info = parse_fasta_header(current_header)
                            raw_sequence = ''.join(current_seq)
                            
                            # Clean sequence and get quality stats
                            cleaned_sequence, quality_stats = clean_sequence(raw_sequence)
                            seq_info.quality_stats = quality_stats
                            seq_info.is_usable = is_sequence_usable(quality_stats, min_completeness)
                            
                            # Only add if usable (or if filtering is disabled)
                            if not filter_low_quality or seq_info.is_usable:
                                sequences.append((seq_info, cleaned_sequence))
                        
                        except Exception as e:
                            print(f"Warning: Skipped sequence at line {line_num}: {e}")
                            continue
                    
                    # Start new sequence
                    current_header = line
                    current_seq = []
                else:
                    current_seq.append(line)
            
            # Process last sequence
            if current_header:
                try:
                    seq_info = parse_fasta_header(current_header)
                    raw_sequence = ''.join(current_seq)
                    
                    # Clean sequence and get quality stats
                    cleaned_sequence, quality_stats = clean_sequence(raw_sequence)
                    seq_info.quality_stats = quality_stats
                    seq_info.is_usable = is_sequence_usable(quality_stats, min_completeness)
                    
                    # Only add if usable (or if filtering is disabled)
                    if not filter_low_quality or seq_info.is_usable:
                        sequences.append((seq_info, cleaned_sequence))
                
                except Exception as e:
                    print(f"Warning: Skipped last sequence: {e}")
    
    except FileNotFoundError:
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
    except Exception as e:
        raise Exception(f"Error reading FASTA file: {e}")
    
    return sequences

def get_sequence_summary(sequences: List[Tuple[SequenceInfo, str]]) -> Dict[str, float]:
    """Generate summary statistics for loaded sequences"""
    if not sequences:
        return {
            "total_sequences": 0,
            "avg_completeness": 0.0,
            "avg_gap_percentage": 0.0,
            "usable_sequences": 0
        }
    
    total_completeness = sum(seq_info.quality_stats["completeness"] for seq_info, _ in sequences if seq_info.quality_stats)
    total_gap_percentage = sum(seq_info.quality_stats["gap_percentage"] for seq_info, _ in sequences if seq_info.quality_stats)
    
    valid_stats = [seq_info for seq_info, _ in sequences if seq_info.quality_stats]
    
    return {
        "total_sequences": len(sequences),
        "avg_completeness": total_completeness / len(valid_stats) if valid_stats else 0.0,
        "avg_gap_percentage": total_gap_percentage / len(valid_stats) if valid_stats else 0.0,
        "usable_sequences": sum(1 for seq_info, _ in sequences if seq_info.is_usable)
    }

def test_header_parsing():
    """Test the flexible header parsing with various formats"""
    test_headers = [
        # Your original 6-field format
        ">OQ123456|SARS-CoV-2/USA/Example/2024|29903|2024-01-15|SARS-CoV-2|XEC.2",
        
        # Simple 2-field format (your problematic case)
        ">Test_Sequence|B.1.1.7",
        
        # Other formats
        ">EPI_ISL_123456|hCoV-19/USA/Sample/2024",
        ">MW123456.1",
        ">Sample_XBB15|XBB.1.5"
    ]
    
    print("🧪 Testing flexible header parsing:")
    print("=" * 50)
    
    for header in test_headers:
        try:
            seq_info = parse_fasta_header(header)
            print(f"\nHeader: {header}")
            print(f"  ✅ Accession: {seq_info.accession}")
            print(f"  ✅ Lineage: {seq_info.pangolin_lineage}")
            print(f"  ✅ Title: {seq_info.title}")
            print(f"  ✅ Date: {seq_info.collection_date}")
            print(f"  ✅ Organism: {seq_info.organism}")
        except Exception as e:
            print(f"\nHeader: {header}")
            print(f"  ❌ Error: {e}")

# Test function
def test_parser():
    """Test the parser with your data"""
    print("🧬 Testing enhanced SARS-CoV-2 sequence parser...")
    
    # First test header parsing
    test_header_parsing()
    
    # Then test with actual file if it exists
    try:
        sequences = load_sequences('../data/MAIN.fasta')
        summary = get_sequence_summary(sequences)
        
        print(f"\n📊 File Analysis Results:")
        print(f"✅ Loaded {summary['total_sequences']} sequences")
        print(f"📊 Average completeness: {summary['avg_completeness']:.2%}")
        print(f"🔍 Average gap percentage: {summary['avg_gap_percentage']:.1f}%")
        print(f"✅ Usable sequences: {summary['usable_sequences']}")
        
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
    
    except FileNotFoundError:
        print(f"\n⚠️  MAIN.fasta not found - testing with sample headers only")
    except Exception as e:
        print(f"\n❌ Error testing with file: {e}")

if __name__ == "__main__":
    test_parser()
