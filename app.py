#!/usr/bin/env python3
"""
SARS-CoV-2 Surveillance Agent
=======================================================
Integrated surveillance platform with sequence analysis and anomaly detection
Enhanced to work with CSV metadata and improved country extraction
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
        """Extract country code from virus title with enhanced patterns"""
        if not self.title:
            return None
        
        # Multiple patterns for country extraction based on the CSV data
        patterns = [
            r'USA:\s*([^,]+)',           # USA: State format
            r'Japan:\s*([^,]+)',         # Japan: Prefecture format  
            r'Brazil:\s*([^,]+)',        # Brazil: State format
            r'United Kingdom:\s*([^,]+)', # UK: Region format
            r'China:\s*([^,]+)',         # China: Region format
            r'SARS-CoV-2/([^/]+)/',      # Standard SARS-CoV-2/Country/
            r'/([A-Z]{2,3})/',           # /USA/ or /GBR/
            r'hCoV-19/([^/]+)/',         # GISAID format
        ]
        
        for pattern in patterns:
            match = re.search(pattern, self.title)
            if match:
                country = match.group(1).strip()
                # Clean up and standardize country names
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
                    'JPN': 'Japan',
                    'BRA': 'Brazil'
                }
                return country_mapping.get(country, country)
        
        return None

def parse_fasta_header(header: str) -> SequenceInfo:
    """
    Parse FASTA header with enhanced format support based on real data
    
    Enhanced to handle formats seen in the CSV:
    1. Standard 6-field: >AccessionID|Title|Length|Date|Organism|Lineage
    2. Simple 2-field: >SequenceName|Lineage  
    3. CSV-style: Extract from title patterns like "USA: State"
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

def load_csv_metadata(csv_path: str) -> Dict[str, Dict]:
    """Load sequence metadata from CSV file"""
    try:
        df = pd.read_csv(csv_path)
        metadata = {}
        
        for _, row in df.iterrows():
            accession = row['Accession']
            metadata[accession] = {
                'country': row.get('Country', 'Unknown'),
                'lineage': row.get('Pangolin', 'Unknown'),
                'collection_date': row.get('Collection_Date', 'Unknown'),
                'length': row.get('Length', 0),
                'geo_location': row.get('Geo_Location', 'Unknown')
            }
        return metadata
    except Exception as e:
        st.warning(f"Could not load CSV metadata: {e}")
        return {}

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
    
    # Sidebar navigation
    with st.sidebar:
        st.markdown("## Navigation")
        page = st.selectbox("Select Analysis Module", [
            "Dashboard Overview",
            "Data Analysis (CSV Mode)",
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
    elif page == "Data Analysis (CSV Mode)":
        csv_analysis_page()
    elif page == "Sequence Upload & Analysis":
        upload_analysis_page()
    elif page == "Anomaly Detection":
        anomaly_detection_page()
    elif page == "Surveillance Reports":
        reports_page()
    elif page == "About":
        about_page()

def csv_analysis_page():
    """Analyze the uploaded CSV file"""
    st.header("Data Analysis (CSV Mode)")
    
    st.info("Upload the sequences.csv file to analyze the metadata and get accurate statistics")
    
    uploaded_csv = st.file_uploader(
        "Select CSV file",
        type=['csv'],
        help="Upload the sequences CSV file for enhanced analysis"
    )
    
    if uploaded_csv:
        try:
            df = pd.read_csv(uploaded_csv)
            
            st.success(f"Loaded CSV with {len(df)} sequences")
            
            # Display key statistics
            st.subheader("Dataset Overview")
            col1, col2, col3, col4 = st.columns(4)
            
            with col1:
                st.metric("Total Sequences", len(df))
            
            with col2:
                unique_countries = df['Country'].dropna().nunique()
                st.metric("Countries", unique_countries)
            
            with col3:
                unique_lineages = df['Pangolin'].dropna().nunique()
                st.metric("Unique Lineages", unique_lineages)
            
            with col4:
                complete_sequences = len(df[df['Nuc_Completeness'] == 'complete'])
                st.metric("Complete Sequences", complete_sequences)
            
            # Country distribution
            st.subheader("Geographic Distribution")
            country_counts = df['Country'].value_counts()
            
            col1, col2 = st.columns([3, 2])
            with col1:
                fig = px.pie(
                    values=country_counts.values,
                    names=country_counts.index,
                    title="Sequences by Country"
                )
                st.plotly_chart(fig, use_container_width=True)
            
            with col2:
                st.markdown("**Country Breakdown:**")
                for country, count in country_counts.items():
                    percentage = (count / len(df)) * 100
                    st.markdown(f"**{country}**: {count} ({percentage:.1f}%)")
            
            # Lineage distribution  
            st.subheader("Lineage Distribution")
            lineage_counts = df['Pangolin'].value_counts().head(10)
            
            col1, col2 = st.columns([3, 2])
            with col1:
                fig = px.bar(
                    x=lineage_counts.values,
                    y=lineage_counts.index,
                    orientation='h',
                    title="Top 10 Lineages",
                    labels={'x': 'Count', 'y': 'Lineage'}
                )
                fig.update_layout(height=400)
                st.plotly_chart(fig, use_container_width=True)
            
            with col2:
                st.markdown("**Lineage Breakdown:**")
                for lineage, count in lineage_counts.items():
                    percentage = (count / len(df)) * 100
                    st.markdown(f"**{lineage}**: {count} ({percentage:.1f}%)")
            
            # Temporal distribution
            if 'Collection_Date' in df.columns:
                st.subheader("Temporal Distribution")
                df['Collection_Date'] = pd.to_datetime(df['Collection_Date'], errors='coerce')
                df['Year'] = df['Collection_Date'].dt.year
                df['Month'] = df['Collection_Date'].dt.to_period('M')
                
                # Year distribution
                col1, col2 = st.columns(2)
                with col1:
                    year_counts = df['Year'].value_counts().sort_index()
                    fig = px.bar(
                        x=year_counts.index,
                        y=year_counts.values,
                        title="Sequences by Year",
                        labels={'x': 'Year', 'y': 'Count'}
                    )
                    st.plotly_chart(fig, use_container_width=True)
                
                with col2:
                    # Monthly trend for recent data
                    recent_months = df['Month'].value_counts().sort_index().tail(12)
                    fig = px.line(
                        x=recent_months.index.astype(str),
                        y=recent_months.values,
                        title="Monthly Trend (Recent 12 months)",
                        labels={'x': 'Month', 'y': 'Count'}
                    )
                    st.plotly_chart(fig, use_container_width=True)
            
            # Length distribution
            st.subheader("Sequence Length Analysis")
            
            df['Length'] = pd.to_numeric(df['Length'], errors='coerce')
            
            col1, col2 = st.columns(2)
            with col1:
                fig = px.histogram(
                    df,
                    x='Length',
                    nbins=30,
                    title="Sequence Length Distribution"
                )
                st.plotly_chart(fig, use_container_width=True)
            
            with col2:
                st.markdown("**Length Statistics:**")
                st.markdown(f"**Mean**: {df['Length'].mean():.0f} bp")
                st.markdown(f"**Median**: {df['Length'].median():.0f} bp")
                st.markdown(f"**Min**: {df['Length'].min():.0f} bp")
                st.markdown(f"**Max**: {df['Length'].max():.0f} bp")
                st.markdown(f"**Std Dev**: {df['Length'].std():.0f} bp")
            
            # Download processed data
            st.subheader("Download Analysis")
            
            # Create summary report
            summary_data = {
                'Metric': ['Total Sequences', 'Countries', 'Lineages', 'Complete Sequences', 'Mean Length'],
                'Value': [len(df), unique_countries, unique_lineages, complete_sequences, f"{df['Length'].mean():.0f} bp"]
            }
            summary_df = pd.DataFrame(summary_data)
            
            csv_summary = summary_df.to_csv(index=False)
            st.download_button(
                "Download Summary Report",
                csv_summary,
                "surveillance_summary.csv",
                "text/csv"
            )
            
            # Show data table
            with st.expander("Raw Data Preview"):
                st.dataframe(df.head(20), use_container_width=True)
        
        except Exception as e:
            st.error(f"Error analyzing CSV: {e}")

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
                # Enhanced country counting
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
                        nbins=20,
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
            
            # Country distribution from CSV data (if available)
            csv_path = Path("sequences.csv")
            if csv_path.exists():
                try:
                    df_csv = pd.read_csv(csv_path)
                    st.subheader("Geographic Distribution (CSV Data)")
                    
                    country_counts = df_csv['Country'].value_counts()
                    fig = px.pie(
                        values=country_counts.values,
                        names=country_counts.index,
                        title="Sequences by Country (Actual Data)"
                    )
                    st.plotly_chart(fig, use_container_width=True)
                    
                    # Show the real numbers
                    st.markdown("**Actual Country Distribution:**")
                    for country, count in country_counts.items():
                        percentage = (count / len(df_csv)) * 100
                        st.markdown(f"**{country}**: {count} sequences ({percentage:.1f}%)")
                
                except Exception as e:
                    st.warning(f"Could not load CSV data: {e}")
        
        else:
            st.warning("No surveillance data found at data/MAIN.fasta")
            st.info("Upload sequences or CSV data using the analysis modules")
            
            # Show example dashboard
            st.subheader("Dashboard Preview")
            st.markdown("""
            Once you upload data, this dashboard will display:
            - **Total sequences** and quality metrics
            - **Lineage distribution** charts and statistics
            - **Geographic distribution** of variants (5 countries in your dataset)
            - **Quality assessment** plots and filters
            
            **Your Dataset Preview (from CSV)**:
            - 223 sequences from 5 countries
            - Primarily XEC.2 (83.4%) and B.1.1.529 variants
            - USA (83.4%), Japan (7.6%), UK (2.2%), Brazil (1.3%), China (0.4%)
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
        
        **Enhanced country extraction**: Now detects patterns like "USA: State", "Japan: Prefecture"
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
                        nbins=20,
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
            
        except Exception as e:
            st.error(f"Analysis failed: {e}")
            st.error("Please check your file format and try again")

def anomaly_detection_page():
    """Anomaly detection results and analysis"""
    st.header("Anomaly Detection")
    
    if not USING_ORIGINAL_MODULES:
        st.warning("Anomaly detection requires the full module suite")
        st.info("Install enhanced_anomaly_detector.py and alignment_module.py to enable this feature")
        return
    
    st.info("Anomaly detection module would analyze for pandemic threats")

def reports_page():
    """Surveillance reports and exports"""
    st.header("Surveillance Reports")
    st.info("Reports will be generated from surveillance analysis")

def about_page():
    """About page with system information"""
    st.header("About SARS-CoV-2 Surveillance Agent")
    
    st.markdown("""
    ## System Overview
    
    This surveillance agent provides real-time analysis of SARS-CoV-2 sequences to identify variants
    with pandemic potential. Built for pandemic preparedness research and public health applications.
    
    ## Enhanced Features
    
    **CSV Mode**: Direct analysis of sequence metadata from CSV files  
    **Geographic Analysis**: Enhanced country extraction from sequence titles  
    **Flexible Parsing**: Support for multiple FASTA header formats  
    **Quality Assessment**: Comprehensive sequence quality metrics
    
    **Developer**: Tobi Lawal (tlawal5@emory.edu)  
    **Institution**: Emory University  
    **Purpose**: Pandemic Preparedness Research
    
    ---
    
    **System Version**: 2.2 
    **Last Updated**: March 2026  
    **License**: MIT License
    """)

if __name__ == "__main__":
    main()
