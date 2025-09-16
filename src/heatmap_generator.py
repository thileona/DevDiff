import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

class HeatmapGenerator:
    def __init__(self, l1_file, l4_file, d1_file, csv_sep=","):
        """
        Initialize the heatmap generator with data files.
        
        Args:
            l1_file (str): Path to L1 stage CSV file
            l4_file (str): Path to L4 stage CSV file  
            d1_file (str): Path to D1 stage CSV file
            csv_sep (str): CSV separator (default: ",")
        """
        self.l1_file = l1_file
        self.l4_file = l4_file
        self.d1_file = d1_file
        self.csv_sep = csv_sep
        
        # Load data
        self.df_L1 = self._read_gene_by_cell_csv(l1_file)
        self.df_L4 = self._read_gene_by_cell_csv(l4_file)
        self.df_D1 = self._read_gene_by_cell_csv(d1_file)
        
        # Find common genes and cells
        self.genes_all = sorted(set(self.df_L1.index) & set(self.df_L4.index) & set(self.df_D1.index))
        self.cells_union = sorted(set(self.df_L1.columns) | set(self.df_L4.columns) | set(self.df_D1.columns))
        
        # Color palette for patterns
        self.palette = [
            "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3",
            "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3"
        ]
        self.pattern_levels = ["000","001","010","011","100","101","110","111"]
        self.pattern_labels = {
            "000": "none (no stage)",
            "001": "D1 only", 
            "010": "L4 only",
            "011": "L4 + D1",
            "100": "L1 only",
            "101": "L1 + D1",
            "110": "L1 + L4", 
            "111": "L1 + L4 + D1"
        }
        
    def _read_gene_by_cell_csv(self, path):
        """Read and parse gene expression CSV file."""
        df = pd.read_csv(path, sep=self.csv_sep, header=0, dtype=str)
        df.columns = [c.strip() for c in df.columns]
        
        # Detect gene_name column
        lowcols = [c.lower() for c in df.columns]
        if "gene_name" in lowcols:
            idx = lowcols.index("gene_name")
            gene_name_col_name = df.columns[idx]
            gene_id_col = df.columns[0] if df.columns[0] != gene_name_col_name else None
        else:
            if len(df.columns) < 2:
                raise ValueError(f"File {path} must have at least 2 columns")
            gene_id_col = df.columns[0]
            gene_name_col_name = df.columns[1]
        
        # Extract expression columns
        exclude = {gene_name_col_name}
        if gene_id_col is not None:
            exclude.add(gene_id_col)
        expr_cols = [c for c in df.columns if c not in exclude]
        
        # Convert to numeric and set gene names as index
        expr = df[expr_cols].apply(pd.to_numeric, errors='coerce')
        expr.index = df[gene_name_col_name].astype(str).values
        expr.index.name = "gene_name"
        expr = expr[~expr.index.duplicated(keep='first')]
        
        return expr
    
    def _build_binary_matrix(self, expr_df, genes, cells, expr_threshold):
        """Convert expression data to binary matrix (cells x genes)."""
        sub = expr_df.reindex(genes).T  # transpose to cells x genes
        sub = sub.reindex(index=cells)
        sub_filled = sub.fillna(0.0)
        return (sub_filled > expr_threshold).astype(int)
    
    def generate_heatmap(self, gene_list, expr_threshold=0.0, sort_genes_by_pattern=True, max_rows=2000, max_cols=2000):
        """
        Generate heatmap for given gene list.
        
        Args:
            gene_list (list): List of gene names to analyze
            expr_threshold (float): Expression threshold for binary classification
            sort_genes_by_pattern (bool): Whether to sort genes by expression pattern
            max_rows (int): Maximum number of rows to plot
            max_cols (int): Maximum number of columns to plot
            
        Returns:
            tuple: (matplotlib figure, summary dict)
        """
        # Filter genes present in all datasets
        genes_selected = [g for g in gene_list if g in self.genes_all]
        missing_genes = [g for g in gene_list if g not in self.genes_all]
        
        if len(genes_selected) == 0:
            raise ValueError("None of the requested genes are present in all three datasets.")
            
        # Build binary matrices
        bin_L1 = self._build_binary_matrix(self.df_L1, genes_selected, self.cells_union, expr_threshold)
        bin_L4 = self._build_binary_matrix(self.df_L4, genes_selected, self.cells_union, expr_threshold) 
        bin_D1 = self._build_binary_matrix(self.df_D1, genes_selected, self.cells_union, expr_threshold)
        
        # Create pattern arrays
        A = bin_L1.reindex(index=self.cells_union, columns=genes_selected).fillna(0).astype(int).values
        B = bin_L4.reindex(index=self.cells_union, columns=genes_selected).fillna(0).astype(int).values  
        C = bin_D1.reindex(index=self.cells_union, columns=genes_selected).fillna(0).astype(int).values
        
        # Build pattern strings
        strs = np.char.add(
            np.char.add(A.astype(str), B.astype(str)),
            C.astype(str)
        )
        
        # Sort genes by pattern if requested
        if sort_genes_by_pattern:
            gene_L1_any = (A.sum(axis=0) > 0).astype(int)
            gene_L4_any = (B.sum(axis=0) > 0).astype(int) 
            gene_D1_any = (C.sum(axis=0) > 0).astype(int)
            gene_patterns = [f"{l1}{l4}{d1}" for l1, l4, d1 in zip(gene_L1_any, gene_L4_any, gene_D1_any)]
            
            order = {"000":0,"001":1,"010":2,"011":3,"100":4,"101":5,"110":6,"111":7}
            gene_order_idx = sorted(range(len(genes_selected)), key=lambda i: order.get(gene_patterns[i], 0))
            
            genes_selected = [genes_selected[i] for i in gene_order_idx]
            strs = strs[:, gene_order_idx]
        
        # Map patterns to integers
        pattern_to_int = {p:i for i,p in enumerate(self.pattern_levels)}
        vec_map = np.vectorize(lambda s: pattern_to_int.get(s, 0))
        plot_matrix = vec_map(strs)
        
        # Subsample if too large
        n_rows, n_cols = plot_matrix.shape
        row_idx = np.arange(min(n_rows, max_rows))
        col_idx = np.arange(min(n_cols, max_cols))
        
        plot_mat = plot_matrix[np.ix_(row_idx, col_idx)]
        plot_rows = [self.cells_union[i] for i in row_idx]
        plot_cols = [genes_selected[j] for j in col_idx]
        
        # Create figure
        fig = self._create_heatmap_plot(plot_mat, plot_rows, plot_cols)
        
        # Summary statistics
        summary = {
            'genes_found': len(genes_selected),
            'missing_genes': missing_genes,
            'total_cells': len(self.cells_union),
            'genes_plotted': len(plot_cols),
            'cells_plotted': len(plot_rows)
        }
        
        return fig, summary
    
    def _create_heatmap_plot(self, plot_matrix, row_labels, col_labels):
        """Create the matplotlib heatmap figure."""
        # Figure size based on data dimensions
        fig_h = max(8, 0.15 * len(row_labels))
        fig_w = max(12, 0.14 * len(col_labels))
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        
        # Create colormap and plot
        cmap = ListedColormap(self.palette)
        im = ax.imshow(plot_matrix, aspect='auto', cmap=cmap, interpolation='nearest', vmin=0, vmax=7)
        
        # Y-axis (cells)
        ax.set_yticks(np.arange(len(row_labels)))
        ax.set_yticklabels(row_labels, fontsize=6)
        
        # X-axis (genes) - limit labels if too many
        max_x_labels = 50
        if len(col_labels) <= max_x_labels:
            ax.set_xticks(np.arange(len(col_labels)))
            ax.set_xticklabels(col_labels, fontsize=8, rotation=45, ha='right')
        else:
            # Show every nth label
            step = int(np.ceil(len(col_labels) / max_x_labels))
            tick_positions = np.arange(len(col_labels))[::step]
            ax.set_xticks(tick_positions)
            ax.set_xticklabels([col_labels[i] for i in tick_positions], fontsize=8, rotation=45, ha='right')
        
        # Labels and title
        ax.set_xlabel("Genes", fontsize=12)
        ax.set_ylabel("Cells", fontsize=12)
        ax.set_title("Gene Expression Patterns Across Developmental Stages (L1-L4-D1)", fontsize=14, pad=20)
        
        # Colorbar
        cbar = plt.colorbar(im, ax=ax, fraction=0.02, pad=0.02, ticks=np.arange(0,8))
        tick_labels = []
        for i, pat in enumerate(self.pattern_levels):
            tick_labels.append(f"{pat} ({self.pattern_labels[pat]})")
        cbar.ax.set_yticklabels(tick_labels, fontsize=9)
        cbar.set_label('Expression Pattern', fontsize=11)
        
        plt.tight_layout()
        return fig
