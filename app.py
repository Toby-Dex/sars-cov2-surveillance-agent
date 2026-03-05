#!/usr/bin/env python3
"""
SARS-CoV-2 Surveillance Agent
============================================
"""

import streamlit as st
import pandas as pd
import sys
from pathlib import Path

# Add scripts directory to path
sys.path.append(str(Path(__file__).parent / "scripts"))

# Page config
st.set_page_config(
    page_title="SARS-CoV-2 Surveillance Agent",
    page_icon="🦠",
    layout="wide"
)

# Main header
st.title("🦠 SARS-CoV-2 Surveillance Agent")
st.markdown("**Real-time phylogenetic analysis and anomaly detection**")
st.markdown("*Developed by Tobi Lawal | Emory University*")

# Sidebar navigation
st.sidebar.title("Navigation")
page = st.sidebar.selectbox("Choose Module", [
    "Dashboard", 
    "Sequence Analysis", 
    "Anomaly Detection"
])

if page == "Dashboard":
    st.header("📊 Dashboard Overview")
    
    try:
        from sequence_parser import load_sequences
        
        if Path("data/MAIN.fasta").exists():
            with st.spinner("Loading surveillance data..."):
                sequences = load_sequences("data/MAIN.fasta")
            
            # Simple metrics
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Total Sequences", len(sequences))
            with col2:
                st.metric("Data Source", "MAIN.fasta")
            with col3:
                st.metric("Status", "✅ Active")
                
            st.success(f"Loaded {len(sequences)} sequences successfully")
        else:
            st.warning("⚠️ No surveillance data found")
            
    except Exception as e:
        st.error(f"Error: {e}")
        
elif page == "Sequence Analysis":
    st.header("Sequence Analysis")
    uploaded_file = st.file_uploader("Upload FASTA file", type=['fasta', 'fa'])
    
    if uploaded_file:
        st.success(f"File uploaded: {uploaded_file.name}")
        
elif page == "Anomaly Detection":
    st.header("⚠️ Anomaly Detection")
    st.info("Identify variants with pandemic potential")

# Footer
st.markdown("---")
st.markdown("*Pandemic Preparedness Research*")
