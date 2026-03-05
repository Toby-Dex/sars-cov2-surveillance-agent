#!/usr/bin/env python3
"""
SARS-CoV-2 Surveillance Agent
=======================================================
Integrated surveillance platform with sequence analysis and anomaly detection
"""

import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import sys
import os
import tempfile
import re
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass

# Add scripts directory to path
sys.path.append(str(Path(__file__).parent / "scripts"))

# Quick fix for flexible header parsing - embedded directly in app
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
    quality_stats: Optional[Dict[str, float]] = None
    is_usable: bool = True
    
    def extract_country_from_title(self) -> Optional[str]:
        """Extract country code from virus title"""
        if not self.title:
            return None
        
        # Multiple patterns for country extraction
        patterns = [
            r'SARS-CoV-2/([^/]+)/',  # SARS-CoV-2/Country/...
            r'/([A-Z]{2,3})/',       # /USA/ or /GBR/
            r'hCoV-19/([^/]+)/',     # GISAID format
        ]
        
        for pattern in patterns:
            match = re.search(pattern, self.title)
            if match:
                country = match.group(1).strip()
                # Clean up common country codes
                country_mapping = {
                    'USA': 'United States',
                    'GBR': 'United Kingdom', 
                    'UK': 'United Kingdom',
                    'CHN': 'China',
                    'IND': 'India',
                    'DEU': 'Germany',
                    'FRA': 'France',
                    'CAN': 'Canada',
                    'AUS': 'Australia',
                    'JPN': 'Japan'
                }
                return country_mapping.get(country, country)
        
        return None

def parse_fasta_header(header: str) -> SequenceInfo:
    """
    Parse FASTA header with flexible format support
    
    Supported formats:
    1. Standard 6-field: >AccessionID|Title|Length|Date|Organism|Lineage
    2. Simple 2-field: >SequenceName|Lineage  
    3. Basic format: >SequenceName
    """
    parts = header.strip('>').split('|')
    
    # Handle different header formats
    if len(parts) >= 6:
        # Standard 6-field format (original format)
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
        # Partial format - do your best
        seq_info = SequenceInfo(
            accession=parts[0].strip(),
            title=parts[1].strip() if len(parts) > 1 else "",
            length=int(parts[2]) if len(parts) > 2 and parts[2].strip().isdigit() else 0,
            collection_date=parts[3].strip() if len(parts) > 3 else "Unknown",
            organism=parts[4].strip() if len(parts) > 4 else "SARS-CoV-2",
            pangolin_lineage=parts[5].strip() if len(parts) > 5 else "Unknown"
        )
    
    else:
        # Single field format: >SequenceName
        seq_info = SequenceInfo(
            accession=parts[0].strip(),
            title="",
            length=0,
            collection_date="Unknown",
            organism="SARS-CoV-2", 
            pangolin_lineage="Unknown"
        )
    
    # Extract country from title if available
    if seq_info.title:
        seq_info.country = seq_info.extract_country_from_title()
    
    return seq_info

def clean_sequence(sequence: str) -> Tuple[str, Dict[str, float]]:
    """Clean sequence and calculate quality metrics"""
    original_length = len(sequence)
    if original_length == 0:
        return "", {"completeness": 0.0, "gap_percentage": 100.0, "ambiguous_percentage": 100.0}
    
    cleaned = re.sub(r'\s+', '', sequence.upper())
    
    valid_bases = len(re.findall(r'[ATGC]', cleaned))
    gaps = len(re.findall(r'[N-]', cleaned))
    ambiguous = len(re.findall(r'[RYWSMKHBVD]', cleaned))
    
    quality_stats = {
        "completeness": valid_bases / original_length,
        "gap_percentage": gaps / original_length * 100,
        "ambiguous_percentage": ambiguous / original_length * 100,
        "original_length": original_length,
        "cleaned_length": len(cleaned)
    }
    
    cleaned = re.sub(r'[RYWSMKHBVD]', 'N', cleaned)
    return cleaned, quality_stats

def is_sequence_usable(quality_stats: Dict[str, float], min_completeness: float = 0.80) -> bool:
    """Determine if sequence meets quality thresholds"""
    return (
        quality_stats["completeness"] >= min_completeness and
        quality_stats["original_length"] >= 25000
    )

def load_sequences(fasta_path: str, filter_low_quality: bool = True, min_completeness: float = 0.80) -> List[Tuple[SequenceInfo, str]]:
    """Load and parse FASTA file with flexible header parsing"""
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
                            
                            cleaned_sequence, quality_stats = clean_sequence(raw_sequence)
                            seq_info.quality_stats = quality_stats
                            seq_info.is_usable = is_sequence_usable(quality_stats, min_completeness)
                            
                            if not filter_low_quality or seq_info.is_usable:
                                sequences.append((seq_info, cleaned_sequence))
                        
                        except Exception as e:
                            print(f"Warning: Skipped sequence at line {line_num}: {e}")
                            continue
                    
                    current_header = line
                    current_seq = []
                else:
                    current_seq.append(line)
            
            # Process last sequence
            if current_header:
                try:
                    seq_info = parse_fasta_header(current_header)
                    raw_sequence = ''.join(current_seq)
                    
                    cleaned_sequence, quality_stats = clean_sequence(raw_sequence)
                    seq_info.quality_stats = quality_stats
                    seq_info.is_usable = is_sequence_usable(quality_stats, min_completeness)
                    
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
    
    valid_stats = [seq_info for seq_info, _ in sequences if seq_info.quality_stats]
    
    if valid_stats:
        total_completeness = sum(seq_info.quality_stats["completeness"] for seq_info in valid_stats)
        total_gap_percentage = sum(seq_info.quality_stats["gap_percentage"] for seq_info in valid_stats)
        
        return {
            "total_sequences": len(sequences),
            "avg_completeness": total_completeness / len(valid_stats),
            "avg_gap_percentage": total_gap_percentage / len(valid_stats),
            "usable_sequences": sum(1 for seq_info, _ in sequences if seq_info.is_usable)
        }
    else:
        return {
            "total_sequences": len(sequences),
            "avg_completeness": 0.0,
            "avg_gap_percentage": 0.0,
            "usable_sequences": 0
        }

# Try to import additional modules
try:
    from scripts.enhanced_anomaly_detector import EnhancedAnomalyDetector, AnomalyScore
    from scripts.alignment_module import AlignmentEngine
    USING_ORIGINAL_MODULES = True
except ImportError:
    class EnhancedAnomalyDetector:
        def __init__(self, *args): pass
        def run_analysis(self): return []
    
    class AlignmentEngine:
        def __init__(self, *args): pass
    
    class AnomalyScore:
        def __init__(self):
            self.sequence_id = "Unknown"
            self.lineage = "Unknown"
            self.risk_level = "LOW"
            self.total_score = 0.0
            self.mutation_score = 0.0
            self.recombinant_score = 0.0
            self.flags = []
    
    USING_ORIGINAL_MODULES = False

# Page configuration
st.set_page_config(
    page_title="SARS-CoV-2 Surveillance Agent",
    page_icon="🦠",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS with light sky blue background
st.markdown("""
<style>
    .stApp {
        background-color: #87CEEB;
    }
    
    .main-header {
        background: linear-gradient(90deg, #1f77b4, #ff7f0e);
        color: white;
        padding: 1rem;
        border-radius: 10px;
        text-align: center;
        margin-bottom: 2rem;
    }
    .metric-card {
        background: #f8f9fa;
        padding: 1rem;
        border-radius: 8px;
        border-left: 4px solid #1f77b4;
    }
    .risk-critical { border-left-color: #dc3545; background: #fff5f5; }
    .risk-high { border-left-color: #fd7e14; background: #fff8f1; }
    .risk-medium { border-left-color: #ffc107; background: #fffbf0; }
    .risk-low { border-left-color: #28a745; background: #f8fff9; }
    .parser-info {
        background: #e3f2fd;
        padding: 0.5rem;
        border-radius: 5px;
        border-left: 3px solid #1976d2;
        margin-bottom: 1rem;
    }
    
    /* Override Streamlit's default white backgrounds */
    .block-container {
        background-color: rgba(255, 255, 255, 0.9);
        border-radius: 10px;
        padding: 2rem;
        margin: 1rem;
    }
</style>
""", unsafe_allow_html=True)

def main():
    """Main application"""
    
    # Header
    st.markdown("""
    <div class="main-header">
        <h1>SARS-CoV-2 Surveillance Agent</h1>
        <p>Real-time phylogenetic analysis and anomaly detection for pandemic preparedness</p>
        <small>Developed by Tobi Lawal | Emory University</small>
    </div>
    """, unsafe_allow_html=True)
    
    # Show parser status
    if not USING_ORIGINAL_MODULES:
        st.markdown("""
        <div class="parser-info">
            <strong>Enhanced Mode:</strong> Using flexible sequence parser with support for multiple header formats
            including simple formats like <code>>SequenceName|Lineage</code>
        </div>
        """, unsafe_allow_html=True)
    
    # Sidebar navigation
    with st.sidebar:
        st.markdown("## Navigation")
        page = st.selectbox("Select Analysis Module", [
            "Dashboard Overview",
            "Sequence Upload & Analysis", 
            "Anomaly Detection",
            "Surveillance Reports",
            "About"
        ])
        
        st.markdown("---")
        st.markdown("### System Status")
        
        # Check system components
        components_status = check_system_components()
        for component, status in components_status.items():
            if status:
                st.success(f"{component}")
            else:
                st.error(f"{component}")
    
    # Route to selected page
    if page == "Dashboard Overview":
        dashboard_page()
    elif page == "Sequence Upload & Analysis":
        upload_analysis_page()
    elif page == "Anomaly Detection":
        anomaly_detection_page()
    elif page == "Surveillance Reports":
        reports_page()
    elif page == "About":
        about_page()

def check_system_components() -> Dict[str, bool]:
    """Check if system components are available"""
    components = {}
    
    # Check data files
    components["MAIN.fasta Data"] = Path("data/MAIN.fasta").exists()
    components["Reference Genome"] = Path("temp_ref.fasta").exists()
    
    # Check modules
    if USING_ORIGINAL_MODULES:
        components["Original Sequence Parser"] = True
        try:
            from enhanced_anomaly_detector import EnhancedAnomalyDetector
            components["Anomaly Detector"] = True
        except:
            components["Anomaly Detector"] = False
        
        try:
            from alignment_module import AlignmentEngine
            components["Alignment Engine"] = True
        except:
            components["Alignment Engine"] = False
    else:
        components["Enhanced Sequence Parser"] = True
        components["Anomaly Detector"] = False
        components["Alignment Engine"] = False
    
    return components

def dashboard_page():
    """Main dashboard with surveillance overview"""
    st.header("Surveillance Dashboard")
    
    try:
        # Load existing surveillance data
        if Path("data/MAIN.fasta").exists():
            with st.spinner("Loading surveillance data..."):
                sequences = load_sequences("data/MAIN.fasta")
                summary = get_sequence_summary(sequences)
            
            # Key metrics
            col1, col2, col3, col4 = st.columns(4)
            with col1:
                st.metric("Total Sequences", f"{summary['total_sequences']:,}")
            with col2:
                st.metric("Avg Completeness", f"{summary['avg_completeness']:.1%}")
            with col3:
                lineages = set(seq_info.pangolin_lineage for seq_info, _ in sequences)
                st.metric("Unique Lineages", len(lineages))
            with col4:
                # Better country counting
                countries = set()
                for seq_info, _ in sequences:
                    if seq_info.country and seq_info.country.strip() and seq_info.country != "Unknown":
                        countries.add(seq_info.country.strip())
                st.metric("Countries", len(countries))
            
            # Lineage distribution
            st.subheader("Lineage Distribution")
            
            lineage_counts = {}
            for seq_info, _ in sequences:
                lineage = seq_info.pangolin_lineage
                lineage_counts[lineage] = lineage_counts.get(lineage, 0) + 1
            
            df_lineages = pd.DataFrame(
                list(lineage_counts.items()), 
                columns=['Lineage', 'Count']
            ).sort_values('Count', ascending=False)
            
            # Charts
            col1, col2 = st.columns([3, 2])
            with col1:
                fig = px.bar(
                    df_lineages.head(10), 
                    x='Count', 
                    y='Lineage',
                    orientation='h',
                    title="Top 10 Lineages"
                )
                fig.update_layout(height=400)
                st.plotly_chart(fig, use_container_width=True)
            
            with col2:
                st.markdown("**Lineage Summary**")
                for i, row in df_lineages.head(10).iterrows():
                    percentage = (row['Count'] / summary['total_sequences']) * 100
                    st.markdown(f"**{row['Lineage']}**: {row['Count']} ({percentage:.1f}%)")
            
            # Quality metrics
            st.subheader("Sequence Quality Overview")
            
            quality_data = []
            for seq_info, _ in sequences:
                if seq_info.quality_stats:
                    quality_data.append({
                        'Lineage': seq_info.pangolin_lineage,
                        'Completeness': seq_info.quality_stats['completeness'],
                        'Gap_Percentage': seq_info.quality_stats['gap_percentage']
                    })
            
            if quality_data:
                df_quality = pd.DataFrame(quality_data)
                
                col1, col2 = st.columns(2)
                with col1:
                    # FIXED: Use nbins instead of bins for Plotly
                    fig = px.histogram(
                        df_quality, 
                        x='Completeness',
                        nbins=20,  # Changed from bins=20 to nbins=20
                        title="Sequence Completeness Distribution"
                    )
                    st.plotly_chart(fig, use_container_width=True)
                
                with col2:
                    fig = px.box(
                        df_quality,
                        y='Gap_Percentage',
                        title="Gap Percentage Distribution"
                    )
                    st.plotly_chart(fig, use_container_width=True)
            
            # Country distribution (if we have country data)
            country_data = []
            for seq_info, _ in sequences:
                if seq_info.country and seq_info.country.strip() and seq_info.country != "Unknown":
                    country_data.append(seq_info.country.strip())
            
            if country_data:
                st.subheader("Geographic Distribution")
                country_counts = pd.Series(country_data).value_counts()
                if len(country_counts) > 1:
                    fig = px.pie(
                        values=country_counts.values,
                        names=country_counts.index,
                        title="Sequences by Country"
                    )
                    st.plotly_chart(fig, use_container_width=True)
                else:
                    st.info("Country information limited - most sequences lack geographic metadata")
        
        else:
            st.warning("No surveillance data found at data/MAIN.fasta")
            st.info("Upload sequences using the 'Sequence Upload & Analysis' module")
            
            # Show example of what the dashboard will display
            st.subheader("Dashboard Preview")
            st.markdown("""
            Once you upload sequences, this dashboard will display:
            - **Total sequences** and quality metrics
            - **Lineage distribution** charts and statistics
            - **Geographic distribution** of variants
            - **Quality assessment** plots and filters
            """)
    
    except Exception as e:
        st.error(f"Dashboard error: {e}")
        st.info("Check system components in sidebar")

def upload_analysis_page():
    """File upload and sequence analysis"""
    st.header("Sequence Upload & Analysis")
    
    # Instructions
    with st.expander("Instructions", expanded=True):
        st.markdown("""
        **Upload SARS-CoV-2 sequences for analysis:**
        
        **Supported formats**: FASTA (.fasta, .fa, .txt)  
        **Flexible headers**: Multiple formats supported including:
        
        **Standard format:**
        ```
        >AccessionID|Title|Length|Date|Organism|Lineage
        >OQ123456|SARS-CoV-2/USA/Sample/2024|29903|2024-01-15|SARS-CoV-2|XEC.2
        ```
        
        **Simple format:**
        ```
        >SequenceName|Lineage
        >Test_Sequence|B.1.1.7
        ```
        
        **Basic format:**
        ```
        >SequenceName
        >Sample_123
        ```
        
        **Analysis includes**: Quality filtering, gap handling, lineage validation, statistics
        """)
    
    # File upload
    uploaded_file = st.file_uploader(
        "Select FASTA file",
        type=['fasta', 'fa', 'txt'],
        help="Upload SARS-CoV-2 sequences in FASTA format"
    )
    
    if uploaded_file:
        st.success(f"File uploaded: {uploaded_file.name} ({uploaded_file.size:,} bytes)")
        
        # Save to temporary file for processing
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp_file:
            try:
                content = str(uploaded_file.read(), 'utf-8')
                tmp_file.write(content)
                tmp_path = tmp_file.name
            except UnicodeDecodeError:
                st.error("File encoding error. Please ensure the file is in UTF-8 format.")
                return
        
        try:
            # Show file preview
            with st.expander("File Preview"):
                lines = content.split('\n')[:20]  # First 20 lines
                for i, line in enumerate(lines):
                    if line.strip():
                        if line.startswith('>'):
                            st.markdown(f"`{i+1:2d}: {line}`")
                        else:
                            st.markdown(f"`{i+1:2d}: {line[:50]}{'...' if len(line) > 50 else ''}`")
                if len(content.split('\n')) > 20:
                    st.markdown("...")
            
            # Analysis options
            st.subheader("Analysis Options")
            col1, col2 = st.columns([3, 1])
            with col1:
                filter_quality = st.checkbox("Filter low-quality sequences", value=True)
                min_completeness = st.slider("Minimum completeness", 0.5, 1.0, 0.8, 0.05)
            
            with col2:
                if st.button("Analyze Sequences", type="primary", use_container_width=True):
                    analyze_uploaded_sequences(tmp_path, filter_quality, min_completeness)
        
        except Exception as e:
            st.error(f"Analysis preparation error: {e}")
        
        finally:
            # Cleanup
            if 'tmp_path' in locals():
                try:
                    os.unlink(tmp_path)
                except:
                    pass

def analyze_uploaded_sequences(file_path: str, filter_quality: bool, min_completeness: float):
    """Analyze uploaded sequences"""
    
    with st.spinner("Analyzing sequences..."):
        try:
            # Load sequences
            sequences = load_sequences(file_path, filter_low_quality=filter_quality, min_completeness=min_completeness)
            
            if not sequences:
                st.error("No valid sequences found in uploaded file")
                st.info("Check the file format and header structure")
                return
            
            # Display results
            st.success(f"Successfully loaded and analyzed {len(sequences)} sequences")
            
            # Summary statistics
            summary = get_sequence_summary(sequences)
            
            st.subheader("Analysis Summary")
            col1, col2, col3, col4 = st.columns(4)
            with col1:
                st.metric("Total Sequences", summary['total_sequences'])
            with col2:
                st.metric("Avg Completeness", f"{summary['avg_completeness']:.1%}")
            with col3:
                st.metric("Avg Gap %", f"{summary['avg_gap_percentage']:.1f}%")
            with col4:
                st.metric("Usable Sequences", summary['usable_sequences'])
            
            # Sequence details
            st.subheader("Sequence Analysis Results")
            
            # Create results DataFrame
            results_data = []
            for seq_info, sequence in sequences:
                results_data.append({
                    'Accession': seq_info.accession,
                    'Lineage': seq_info.pangolin_lineage,
                    'Country': seq_info.country or 'Unknown',
                    'Date': seq_info.collection_date or 'Not specified',
                    'Length': len(sequence),
                    'Completeness': f"{seq_info.quality_stats['completeness']:.1%}" if seq_info.quality_stats else 'N/A',
                    'Gap %': f"{seq_info.quality_stats['gap_percentage']:.1f}%" if seq_info.quality_stats else 'N/A',
                    'Status': 'Usable' if seq_info.is_usable else 'Low Quality'
                })
            
            df_results = pd.DataFrame(results_data)
            st.dataframe(df_results, use_container_width=True)
            
            # Visualizations
            if len(sequences) > 1:
                st.subheader("Analysis Visualizations")
                
                col1, col2 = st.columns(2)
                with col1:
                    # Lineage distribution
                    lineage_counts = df_results['Lineage'].value_counts()
                    fig = px.pie(
                        values=lineage_counts.values,
                        names=lineage_counts.index,
                        title="Lineage Distribution"
                    )
                    st.plotly_chart(fig, use_container_width=True)
                
                with col2:
                    # Quality distribution - FIXED: Use nbins instead of bins
                    df_results['Completeness_Float'] = df_results['Completeness'].str.rstrip('%').astype(float) / 100
                    fig = px.histogram(
                        df_results,
                        x='Completeness_Float',
                        nbins=20,  # Changed from bins=20 to nbins=20
                        title="Quality Distribution"
                    )
                    fig.update_xaxis(title="Completeness")
                    st.plotly_chart(fig, use_container_width=True)
            
            # Download option
            csv_data = df_results.to_csv(index=False)
            st.download_button(
                "Download Results as CSV",
                csv_data,
                "sequence_analysis_results.csv",
                "text/csv",
                use_container_width=True
            )
            
            # Option to run anomaly detection
            if USING_ORIGINAL_MODULES:
                st.subheader("Advanced Threat Analysis")
                if st.button("Run Anomaly Detection", type="secondary", use_container_width=True):
                    run_anomaly_detection_on_upload(sequences)
            else:
                st.info("Advanced anomaly detection available with full module installation")
        
        except Exception as e:
            st.error(f"Analysis failed: {e}")
            st.error("Please check your file format and try again")

def run_anomaly_detection_on_upload(sequences):
    """Run anomaly detection on uploaded sequences"""
    
    st.info("Anomaly detection requires sequence alignment. This feature connects to the full pipeline.")
    st.markdown("**Next steps:**")
    st.markdown("1. Sequences will be aligned against WIV04 reference")
    st.markdown("2. Mutation patterns will be analyzed") 
    st.markdown("3. Pandemic potential scores will be calculated")
    st.markdown("4. Risk assessment report will be generated")
    
    st.warning("Full anomaly detection pipeline integration coming soon!")

def anomaly_detection_page():
    """Anomaly detection results and analysis"""
    st.header("Anomaly Detection")
    
    if not USING_ORIGINAL_MODULES:
        st.warning("Anomaly detection requires the full module suite")
        st.info("Install enhanced_anomaly_detector.py and alignment_module.py to enable this feature")
        return
    
    # Check for existing analysis
    if Path("outputs/alignments/aligned_sequences.fasta").exists():
        st.info("Loading existing anomaly detection results...")
        
        try:
            with st.spinner("Analyzing sequences for pandemic threats..."):
                detector = EnhancedAnomalyDetector("outputs/alignments/aligned_sequences.fasta")
                results = detector.run_analysis()
            
            display_anomaly_results(results)
        
        except Exception as e:
            st.error(f"Anomaly detection failed: {e}")
            st.info("Make sure enhanced_anomaly_detector.py is available")
    
    else:
        st.warning("No aligned sequence data found")
        st.info("Run sequence alignment first using the Upload & Analysis module")
        
        # Show expected analysis capabilities
        show_threat_detection_capabilities()

def show_threat_detection_capabilities():
    """Show threat detection capabilities"""
    st.subheader("Threat Detection Capabilities")
    
    threat_categories = {
        "Mutation Analysis": [
            "PRRA furin cleavage site detection",
            "ACE2 binding mutations (N501Y, Q493E)",
            "Immune escape variants (E484K, R346T)", 
            "Transmissibility markers (P681H, L452R)"
        ],
        "Recombination Detection": [
            "XEC recombinant lineages",
            "Cross-lineage genetic exchange",
            "Unusual insertion patterns",
            "Hybrid variant identification"
        ],
        "Risk Assessment": [
            "Pandemic potential scoring",
            "Geographic spread analysis", 
            "Temporal emergence patterns",
            "Lineage consistency validation"
        ]
    }
    
    for category, features in threat_categories.items():
        with st.expander(category):
            for feature in features:
                st.markdown(f"• {feature}")

def display_anomaly_results(results: List):
    """Display anomaly detection results"""
    st.subheader("Threat Assessment Summary")
    st.success("Anomaly detection completed")

def reports_page():
    """Surveillance reports and exports"""
    st.header("Surveillance Reports")
    
    # Check for existing reports
    reports_dir = Path("outputs/reports")
    if reports_dir.exists():
        reports = list(reports_dir.glob("*.txt"))
        
        if reports:
            st.success(f"Found {len(reports)} surveillance reports")
            
            for report_path in reports:
                with st.expander(f"Report: {report_path.name}"):
                    try:
                        content = report_path.read_text()
                        st.text_area("Report Content", content, height=400)
                        
                        # Download option
                        st.download_button(
                            f"Download {report_path.name}",
                            content,
                            report_path.name,
                            "text/plain"
                        )
                    except Exception as e:
                        st.error(f"Error reading report: {e}")
        else:
            st.info("No reports generated yet")
    else:
        st.warning("Reports directory not found")
        st.info("Reports will be generated after running anomaly detection analysis")

def about_page():
    """About page with system information"""
    st.header("About SARS-CoV-2 Surveillance Agent")
    
    st.markdown("""
    ## System Overview
    
    This surveillance agent provides real-time analysis of SARS-CoV-2 sequences to identify variants
    with pandemic potential. Built for pandemic preparedness research and public health applications.
    
    ## Core Capabilities
    
    **Sequence Analysis**
    - Flexible FASTA file processing with multiple header format support
    - Quality filtering and sequence validation
    - Gap handling and sequence cleaning
    - Metadata extraction and lineage identification
    
    **Phylogenetic Analysis** 
    - MAFFT-based sequence alignment (when available)
    - WIV04 reference genome comparison
    - Mutation pattern detection
    - Lineage assignment and validation
    
    **Threat Assessment** (Full Version)
    - Recombination detection (XEC, JN.1, etc.)
    - High-risk mutation analysis (PRRA, N501Y, Q493E)
    - Pandemic potential scoring
    - Risk level classification
    
    ## Technical Architecture
    
    **Frontend**: Streamlit web application  
    **Backend**: Python with BioPython integration  
    **Alignment**: MAFFT for sequence alignment  
    **Analysis**: Custom anomaly detection algorithms  
    **Visualization**: Plotly for interactive charts  
    **Parser**: Enhanced flexible header format support
    
    ## Development Team
    
    **Developer**: Tobi Lawal (Tlawal5@emory.edu)  
    **Institution**: Emory University
    **Purpose**: Pandemic Preparedness Research
    
    ## Impact & Applications
    
    - **Public Health Surveillance**: Early detection of concerning variants
    - **Research Support**: Automated analysis for virology labs
    - **Pandemic Preparedness**: Risk assessment for public health officials
    - **Educational Tool**: Training platform for bioinformatics analysis
    
    ---
    
    **System Version**: 2.1 Enhanced (Flexible Headers)  
    **Last Updated**: March 2026  
    **License**: MIT License  
    **Repository**: [GitHub](https://github.com/Toby-Dex/sars-cov2-surveillance-agent)
    """)
    
    # System diagnostics
    with st.expander("System Diagnostics"):
        st.subheader("Component Status")
        components = check_system_components()
        
        for component, status in components.items():
            if status:
                st.success(f"{component}: Operational")
            else:
                st.error(f"{component}: Not Available")
        
        st.subheader("File System")
        paths_to_check = [
            "data/MAIN.fasta",
            "temp_ref.fasta", 
            "scripts/sequence_parser.py",
            "scripts/enhanced_anomaly_detector.py",
            "scripts/alignment_module.py",
            "outputs/alignments/",
            "outputs/reports/"
        ]
        
        for path_str in paths_to_check:
            path = Path(path_str)
            if path.exists():
                st.success(f"{path_str}")
            else:
                st.warning(f"{path_str} - Not Found")

if __name__ == "__main__":
    main()
