#!/usr/bin/env python3
"""
SARS-CoV-2 Surveillance Agent - Streamlit App
=============================================
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
from pathlib import Path
from typing import List, Dict, Optional

# Add scripts directory to path
sys.path.append(str(Path(__file__).parent / "scripts"))

# Import your existing modules
try:
    from sequence_parser import load_sequences, get_sequence_summary, SequenceInfo
    from enhanced_anomaly_detector import EnhancedAnomalyDetector, AnomalyScore
    from alignment_module import AlignmentEngine
except ImportError as e:
    st.error(f"⚠️ Module import error: {e}")
    st.info("Make sure all modules are available in the scripts/ directory")

# Page configuration
st.set_page_config(
    page_title="SARS-CoV-2 Surveillance Agent",
    page_icon="🦠",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS
st.markdown("""
<style>
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
</style>
""", unsafe_allow_html=True)

def main():
    """Main application"""
    
    # Header
    st.markdown("""
    <div class="main-header">
        <h1> SARS-CoV-2 Surveillance Agent</h1>
        <p>Real-time phylogenetic analysis and anomaly detection for pandemic preparedness</p>
        <small>Developed by Tobi Lawal | Emory University</small>
    </div>
    """, unsafe_allow_html=True)
    
    # Sidebar navigation
    with st.sidebar:
        st.markdown("## 🎛️ Navigation")
        page = st.selectbox("Select Analysis Module", [
            "Dashboard Overview",
            "Sequence Upload & Analysis", 
            "⚠️ Anomaly Detection",
            "Surveillance Reports",
            "About"
        ])
        
        st.markdown("---")
        st.markdown("### 🔬 System Status")
        
        # Check system components
        components_status = check_system_components()
        for component, status in components_status.items():
            if status:
                st.success(f"✅ {component}")
            else:
                st.error(f"❌ {component}")
    
    # Route to selected page
    if page == "Dashboard Overview":
        dashboard_page()
    elif page == "Sequence Upload & Analysis":
        upload_analysis_page()
    elif page == "⚠️ Anomaly Detection":
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
    try:
        from sequence_parser import load_sequences
        components["Sequence Parser"] = True
    except:
        components["Sequence Parser"] = False
    
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
                countries = set(seq_info.country for seq_info, _ in sequences if seq_info.country)
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
                    fig = px.histogram(
                        df_quality, 
                        x='Completeness',
                        bins=20,
                        title="Sequence Completeness Distribution"
                    )
                    st.plotly_chart(fig, use_container_width=True)
                
                with col2:
                    fig = px.box(
                        df_quality,
                        y='Gap_Percentage',
                        title="Gap Percentage by Lineage"
                    )
                    st.plotly_chart(fig, use_container_width=True)
        
        else:
            st.warning("⚠️ No surveillance data found at data/MAIN.fasta")
            st.info("Upload sequences using the 'Sequence Upload & Analysis' module")
    
    except Exception as e:
        st.error(f"Dashboard error: {e}")
        st.info("Check system components in sidebar")

def upload_analysis_page():
    """File upload and sequence analysis"""
    st.header("Sequence Upload & Analysis")
    
    # Instructions
    with st.expander("Instructions", expanded=False):
        st.markdown("""
        **Upload SARS-CoV-2 sequences for analysis:**
        
        1. **Supported formats**: FASTA (.fasta, .fa, .txt)
        2. **Expected header format**: `>AccessionID|Title|Length|Date|Organism|Lineage`
        3. **Analysis includes**: Quality filtering, gap handling, lineage validation
        
        **Example header**:
        ```
        >OQ123456|SARS-CoV-2/USA/Example/2024|29903|2024-01-15|SARS-CoV-2|XEC.2
        ```
        """)
    
    # File upload
    uploaded_file = st.file_uploader(
        "Select FASTA file",
        type=['fasta', 'fa', 'txt'],
        help="Upload SARS-CoV-2 sequences in FASTA format"
    )
    
    if uploaded_file:
        st.success(f"✅ File uploaded: {uploaded_file.name} ({uploaded_file.size:,} bytes)")
        
        # Save to temporary file for processing
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp_file:
            content = str(uploaded_file.read(), 'utf-8')
            tmp_file.write(content)
            tmp_path = tmp_file.name
        
        try:
            # Analysis options
            col1, col2 = st.columns([3, 1])
            with col1:
                st.markdown("**Analysis Options:**")
                filter_quality = st.checkbox("Filter low-quality sequences", value=True)
                min_completeness = st.slider("Minimum completeness", 0.5, 1.0, 0.8, 0.05)
            
            with col2:
                if st.button("Analyze Sequences", type="primary"):
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
            sequences = load_sequences(file_path, filter_low_quality=filter_quality)
            
            if not sequences:
                st.error("No valid sequences found in uploaded file")
                return
            
            # Display results
            st.success(f"✅ Loaded {len(sequences)} sequences")
            
            # Summary statistics
            summary = get_sequence_summary(sequences)
            
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
                    'Date': seq_info.collection_date,
                    'Length': len(sequence),
                    'Completeness': f"{seq_info.quality_stats['completeness']:.1%}" if seq_info.quality_stats else 'N/A',
                    'Gap %': f"{seq_info.quality_stats['gap_percentage']:.1f}%" if seq_info.quality_stats else 'N/A',
                    'Status': '✅ Usable' if seq_info.is_usable else '⚠️ Low Quality'
                })
            
            df_results = pd.DataFrame(results_data)
            st.dataframe(df_results, use_container_width=True)
            
            # Download option
            csv_data = df_results.to_csv(index=False)
            st.download_button(
                "📥 Download Results",
                csv_data,
                "sequence_analysis_results.csv",
                "text/csv"
            )
            
            # Option to run anomaly detection
            st.subheader("Advanced Threat Analysis")
            if st.button("Run Anomaly Detection", type="secondary"):
                run_anomaly_detection_on_upload(sequences)
        
        except Exception as e:
            st.error(f"Analysis failed: {e}")

def run_anomaly_detection_on_upload(sequences):
    """Run anomaly detection on uploaded sequences"""
    
    st.info("⚠️ Anomaly detection requires sequence alignment. This feature connects to the full pipeline.")
    st.markdown("**Next steps:**")
    st.markdown("1. Sequences will be aligned against WIV04 reference")
    st.markdown("2. Mutation patterns will be analyzed") 
    st.markdown("3. Pandemic potential scores will be calculated")
    st.markdown("4. Risk assessment report will be generated")
    
    # This would integrate with your alignment_module and enhanced_anomaly_detector
    st.warning("Full anomaly detection pipeline integration coming soon!")

def anomaly_detection_page():
    """Anomaly detection results and analysis"""
    st.header("⚠️ Anomaly Detection")
    
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
        st.warning("⚠️ No aligned sequence data found")
        st.info("Run sequence alignment first using the Upload & Analysis module")
        
        # Show expected analysis capabilities
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

def display_anomaly_results(results: List[AnomalyScore]):
    """Display anomaly detection results"""
    
    # Summary metrics
    critical = len([r for r in results if r.risk_level == "CRITICAL"])
    high = len([r for r in results if r.risk_level == "HIGH"])
    medium = len([r for r in results if r.risk_level == "MEDIUM"])
    low = len(results) - critical - high - medium
    
    st.subheader("🚨 Threat Assessment Summary")
    
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.markdown('<div class="metric-card risk-critical">', unsafe_allow_html=True)
        st.metric("CRITICAL", critical)
        st.markdown('</div>', unsafe_allow_html=True)
    
    with col2:
        st.markdown('<div class="metric-card risk-high">', unsafe_allow_html=True)
        st.metric("HIGH", high) 
        st.markdown('</div>', unsafe_allow_html=True)
    
    with col3:
        st.markdown('<div class="metric-card risk-medium">', unsafe_allow_html=True)
        st.metric("MEDIUM", medium)
        st.markdown('</div>', unsafe_allow_html=True)
    
    with col4:
        st.markdown('<div class="metric-card risk-low">', unsafe_allow_html=True)
        st.metric("LOW", low)
        st.markdown('</div>', unsafe_allow_html=True)
    
    # High-risk alerts
    high_risk = [r for r in results if r.risk_level in ["CRITICAL", "HIGH"]]
    if high_risk:
        st.error(f"🚨 {len(high_risk)} sequences with significant pandemic potential detected")
        
        st.subheader("⚠️ High-Priority Threats")
        for i, result in enumerate(high_risk[:10]):
            with st.expander(f"{i+1}. {result.sequence_id} - {result.risk_level} RISK"):
                col1, col2 = st.columns(2)
                with col1:
                    st.markdown(f"**Lineage**: {result.lineage}")
                    st.markdown(f"**Total Score**: {result.total_score:.3f}")
                    st.markdown(f"**Mutation Score**: {result.mutation_score:.3f}")
                    st.markdown(f"**Recombinant Score**: {result.recombinant_score:.3f}")
                
                with col2:
                    st.markdown("**Detected Threats**:")
                    for flag in result.flags[:5]:
                        st.markdown(f"• {flag}")
    
    else:
        st.success("✅ No critical threats detected in current analysis")
    
    # Risk distribution chart
    st.subheader("📊 Risk Distribution")
    
    risk_data = pd.DataFrame([
        {'Risk Level': 'CRITICAL', 'Count': critical},
        {'Risk Level': 'HIGH', 'Count': high}, 
        {'Risk Level': 'MEDIUM', 'Count': medium},
        {'Risk Level': 'LOW', 'Count': low}
    ])
    
    fig = px.bar(
        risk_data, 
        x='Risk Level', 
        y='Count',
        color='Risk Level',
        color_discrete_map={
            'CRITICAL': '#dc3545',
            'HIGH': '#fd7e14', 
            'MEDIUM': '#ffc107',
            'LOW': '#28a745'
        }
    )
    st.plotly_chart(fig, use_container_width=True)

def reports_page():
    """Surveillance reports and exports"""
    st.header("📈 Surveillance Reports")
    
    # Check for existing reports
    reports_dir = Path("outputs/reports")
    if reports_dir.exists():
        reports = list(reports_dir.glob("*.txt"))
        
        if reports:
            st.success(f"📋 Found {len(reports)} surveillance reports")
            
            for report_path in reports:
                with st.expander(f"📄 {report_path.name}"):
                    try:
                        content = report_path.read_text()
                        st.text_area("Report Content", content, height=400)
                        
                        # Download option
                        st.download_button(
                            f"📥 Download {report_path.name}",
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
    st.header("ℹ️ About SARS-CoV-2 Surveillance Agent")
    
    st.markdown("""
    ## System Overview
    
    This surveillance agent provides real-time analysis of SARS-CoV-2 sequences to identify variants
    with pandemic potential. Built for pandemic preparedness research and public health applications.
    
    ## 🔬 Core Capabilities
    
    **Sequence Analysis**
    - FASTA file processing with quality filtering
    - Gap handling and sequence cleaning
    - Metadata extraction and validation
    
    **Phylogenetic Analysis** 
    - MAFFT-based sequence alignment
    - WIV04 reference genome comparison
    - Mutation pattern detection
    
    **Threat Assessment**
    - Recombination detection (XEC, JN.1, etc.)
    - High-risk mutation analysis (PRRA, N501Y, Q493E)
    - Pandemic potential scoring
    - Risk level classification
    
    ## 📊 Technical Architecture
    
    **Frontend**: Streamlit web application
    **Backend**: Python with BioPython integration
    **Alignment**: MAFFT for sequence alignment  
    **Analysis**: Custom anomaly detection algorithms
    **Visualization**: Plotly for interactive charts
    
    ## 👨‍🔬 Development Team
    
    **Developer**: Tobi Lawal (tobilawal091@gmail.com)  
    **Institution**: Emory University  
    **Purpose**: Pandemic Preparedness Research
    
    ## 📈 Impact & Applications
    
    - **Public Health Surveillance**: Early detection of concerning variants
    - **Research Support**: Automated analysis for virology labs
    - **Pandemic Preparedness**: Risk assessment for public health officials
    - **Educational Tool**: Training platform for bioinformatics analysis
    
    ---
    
    **System Version**: 2.0 Enhanced  
    **Last Updated**: March 2026  
    **License**: MIT License  
    **Repository**: [GitHub](https://github.com/Toby-Dex/sars-cov2-surveillance-agent)
    """)
    
    # System diagnostics
    with st.expander("🔧 System Diagnostics"):
        st.subheader("Component Status")
        components = check_system_components()
        
        for component, status in components.items():
            if status:
                st.success(f"✅ {component}: Operational")
            else:
                st.error(f"❌ {component}: Not Available")
        
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
                st.success(f"✅ {path_str}")
            else:
                st.warning(f"⚠️ {path_str} - Not Found")

if __name__ == "__main__":
    main()
