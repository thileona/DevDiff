# 🧬 Gene Expression Heatmap Tool

A simple web application to generate gene expression heatmaps across developmental stages (L1, L4, D1) using CenGEN data.

## 🚀 Quick Start

### Option 1: Use the deployed app (Recommended)
Visit: [YOUR_STREAMLIT_URL_HERE] *(Update this after deployment)*

### Option 2: Run locally
```bash
# Clone the repository
git clone https://github.com/YOUR_USERNAME/gene-heatmap-tool.git
cd gene-heatmap-tool

# Install dependencies
pip install -r requirements.txt

# Run the app
streamlit run app.py
```

## 📊 How to Use

1. **Enter genes**: Type gene names in the text area, one per line (e.g., `inx-1`, `act-1`, `unc-54`)
2. **Adjust parameters**: Use the sidebar to modify expression threshold and sorting options
3. **Generate**: Click "Generate Heatmap" 
4. **Download**: Save results as PNG or PDF files

## 🎨 Expression Patterns

The heatmap shows 8 possible expression patterns:

| Pattern | Code | Description |
|---------|------|-------------|
| A | 000 | No expression |
| B | 001 | D1 only |
| C | 010 | L4 only |
| D | 011 | L4 + D1 |
| E | 100 | L1 only |
| F | 101 | L1 + D1 |
| G | 110 | L1 + L4 |
| H | 111 | All stages |

## 📁 Project Structure

```
gene-heatmap-tool/
├── app.py                    # Streamlit application
├── requirements.txt          # Python dependencies
├── data/                     # CSV data files
│   ├── L1CENGEN.csv
│   ├── L4CENGEN.csv
│   └── D1CENGEN.csv
├── src/
│   └── heatmap_generator.py  # Core heatmap generation logic
└── README.md
```

## 🔧 Technical Details

- **Data**: CenGEN single-cell RNA-seq data
- **Framework**: Streamlit for web interface
- **Analysis**: Binary expression classification with customizable thresholds
- **Output**: Interactive heatmaps with downloadable PNG/PDF exports

## 📋 Requirements

- Python 3.8+
- Streamlit
- Pandas
- NumPy  
- Matplotlib

## 🚀 Deployment

### Streamlit Cloud (Recommended)

1. Push your code to GitHub
2. Go to [share.streamlit.io](https://share.streamlit.io)
3. Connect your GitHub repo
4. Deploy with one click!

### Local Development

```bash
# Install dependencies
pip install -r requirements.txt

# Run locally
streamlit run app.py
```

## 📝 Data Format

CSV files should have:
- Genes as rows with a `gene_name` column
- Cells as columns with expression values
- Numeric expression data

## 🤝 Contributing

This is a simple lab tool. Feel free to fork and modify for your needs!

## 📄 License

MIT License - Feel free to use and modify.

---

**Questions?** Contact the lab for support or create an issue on GitHub.
