#!/usr/bin/env python3
"""
Simple Anomaly Detector Test
"""

def load_simple_sequences(aligned_fasta_path):
    """Load sequences with simplified header parsing"""
    sequences = {}
    current_header = None
    current_seq = []
    
    with open(aligned_fasta_path, 'r') as f:
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
                        sequence = ''.join(current_seq)
                        sequences[accession] = {'lineage': lineage, 'sequence': sequence}
                
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
                sequence = ''.join(current_seq)
                sequences[accession] = {'lineage': lineage, 'sequence': sequence}
    
    return sequences

def simple_anomaly_test():
    """Quick test of sequence loading and basic analysis"""
    
    # Load sequences
    sequences = load_simple_sequences('../outputs/alignments/aligned_sequences.fasta')
    
    print(f"Loaded {len(sequences)} sequences")
    
    # Count by lineage
    lineage_counts = {}
    for acc, data in sequences.items():
        lineage = data['lineage']
        lineage_counts[lineage] = lineage_counts.get(lineage, 0) + 1
    
    print("\nLineage distribution:")
    for lineage, count in sorted(lineage_counts.items()):
        print(f"  {lineage}: {count}")
    
    # Flag XEC.2 sequences as high risk
    xec_sequences = [acc for acc, data in sequences.items() if 'XEC' in data['lineage']]
    
    print(f"\nHigh-risk XEC sequences found: {len(xec_sequences)}")
    print("First 5 XEC sequences:")
    for i, acc in enumerate(xec_sequences[:5]):
        lineage = sequences[acc]['lineage']
        print(f"  {i+1}. {acc} ({lineage})")
    
    return sequences

if __name__ == "__main__":
    sequences = simple_anomaly_test()
