import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import io
import base64
from src.heatmap_generator import HeatmapGenerator

# Configuration de la page
st.set_page_config(
    page_title="Gene Expression Heatmap Tool",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Titre principal
st.title("🧬 Gene Expression Heatmap Generator")
st.markdown("**Analyze gene expression patterns across developmental stages (L1, L4, D1)**")

# Sidebar pour les paramètres
st.sidebar.header("Parameters")
expr_threshold = st.sidebar.slider(
    "Expression Threshold", 
    min_value=0.0, 
    max_value=10.0, 
    value=0.0, 
    step=0.1,
    help="Values above this threshold are considered 'expressed'"
)

sort_genes = st.sidebar.checkbox(
    "Sort genes by expression pattern", 
    value=True,
    help="Groups genes with similar expression patterns together"
)

# Zone principale pour l'input des gènes
st.header("Enter Gene List")
st.markdown("Enter one gene name per line (e.g., inx-1, inx-7, act-1):")

# Zone de texte pour les gènes
gene_input = st.text_area(
    "Gene Names",
    height=150,
    placeholder="inx-1\ninx-7\ninx-18\nact-1\nunc-54",
    help="Enter gene names separated by line breaks"
)

# Exemple de gènes
if st.button("Load Example Genes"):
    st.session_state.example_genes = "inx-1\ninx-7\ninx-18\nact-1\nunc-54"

if 'example_genes' in st.session_state:
    gene_input = st.text_area(
        "Gene Names",
        value=st.session_state.example_genes,
        height=150,
        key="genes_with_example"
    )

# Bouton pour générer la heatmap
if st.button("Generate Heatmap", type="primary"):
    if not gene_input.strip():
        st.error("Please enter at least one gene name!")
    else:
        # Parse la liste des gènes
        genes = [gene.strip() for gene in gene_input.strip().split('\n') if gene.strip()]
        
        with st.spinner("Generating heatmap..."):
            try:
                # Initialiser le générateur
                generator = HeatmapGenerator(
                    l1_file="data/L1CENGEN.csv",
                    l4_file="data/L4CENGEN.csv", 
                    d1_file="data/D1CENGEN.csv"
                )
                
                # Générer la heatmap
                fig, summary = generator.generate_heatmap(
                    genes, 
                    expr_threshold=expr_threshold,
                    sort_genes_by_pattern=sort_genes
                )
                
                # Afficher les résultats
                col1, col2 = st.columns([2, 1])
                
                with col1:
                    st.subheader("Expression Pattern Heatmap")
                    st.pyplot(fig)
                
                with col2:
                    st.subheader("Summary")
                    st.metric("Total genes requested", len(genes))
                    st.metric("Genes found in data", summary['genes_found'])
                    st.metric("Total cells analyzed", summary['total_cells'])
                    
                    if summary['missing_genes']:
                        st.warning(f"Missing genes: {', '.join(summary['missing_genes'][:5])}")
                        if len(summary['missing_genes']) > 5:
                            st.caption(f"...and {len(summary['missing_genes'])-5} more")
                
                # Boutons de téléchargement
                st.subheader("Download Results")
                col1, col2 = st.columns(2)
                
                with col1:
                    # PNG download
                    img_buffer = io.BytesIO()
                    fig.savefig(img_buffer, format='png', dpi=300, bbox_inches='tight')
                    img_buffer.seek(0)
                    
                    st.download_button(
                        label="Download PNG",
                        data=img_buffer.getvalue(),
                        file_name="gene_heatmap.png",
                        mime="image/png"
                    )
                
                with col2:
                    # PDF download
                    pdf_buffer = io.BytesIO()
                    fig.savefig(pdf_buffer, format='pdf', dpi=300, bbox_inches='tight')
                    pdf_buffer.seek(0)
                    
                    st.download_button(
                        label="Download PDF",
                        data=pdf_buffer.getvalue(),
                        file_name="gene_heatmap.pdf",
                        mime="application/pdf"
                    )
                        
            except Exception as e:
                st.error(f"Error generating heatmap: {str(e)}")
                st.info("Please check that your gene names match the format in the data files.")

# Footer avec informations
st.markdown("---")
st.markdown("""
**Legend:**
- **000**: No expression (none)
- **001**: D1 only  
- **010**: L4 only
- **011**: L4 + D1
- **100**: L1 only
- **101**: L1 + D1  
- **110**: L1 + L4
- **111**: All stages (L1 + L4 + D1)

**Data:** CenGEN single-cell RNA-seq data across developmental stages
""")

# Instructions d'utilisation
with st.expander("How to use"):
    st.markdown("""
    1. **Enter genes**: Type gene names in the text area, one per line
    2. **Adjust parameters**: Use the sidebar to modify expression threshold and sorting
    3. **Generate**: Click the "Generate Heatmap" button
    4. **Download**: Save your results as PNG or PDF files
    
    **Tips:**
    - Gene names should match the format in CenGEN data (e.g., 'inx-1', 'act-1')
    - Higher expression thresholds will be more stringent
    - Sorting by pattern groups similar genes together for easier interpretation
    """)
