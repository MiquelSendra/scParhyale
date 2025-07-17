# 🧬 scRNAseq — Parhyale hawaiensis

This project contains the full single-cell RNA sequencing (scRNA-seq) analysis of *Parhyale hawaiensis* embryos. The dataset includes both **intact** and **ablated** embryos sampled at the **beginning**, **middle**, and **end** of the regenerative replacement process. The analysis is implemented in Python using `scanpy`, with preprocessing and alignment in Bash using `STARsolo`.


## 📁 Folder structure

scParhyale/
├── data/                # Raw and processed count matrices
│   ├── raw/             # Filtered STARsolo outputs (barcodes, genes, matrix)
│   └── processed/       # Final AnnData objects (.h5ad)
├── metadata/            # Sample annotations (e.g. library, stage, genotype)
├── notebooks/           # Python notebooks for downstream analysis
├── results/             # QC plots, analysis outputs
├── scripts/             # Bash scripts used for preprocessing (e.g. STARsolo)
├── env/                 # Conda environment setup
├── README.md
└── .gitignore


## ⚙️ Preprocessing summary

### Tools used
- **STAR v2.7.11a** for genome indexing and alignment
- **STARsolo** for single-cell quantification
- **PIPseq T10 v5.0 kit** (3′ scRNA-seq, paired-end 2×50 bp)
- **scanpy** for downstream Python analysis (to be developed)

### Read structure
- **Read 1** contains:
  - 16 bp cell barcode (`CB:1-16`)
  - 12 bp UMI (`UB:17-28`)
- **Read 2** contains the transcript sequence
- Sequencing depth and layout as expected from Illumina NovaSeq 6000


## 🧬 Reference genome indexing

The reference genome was indexed using STAR with the following files:

- Genome FASTA: `genomeV5.fa`
- Gene annotations: `Gene_Models_210210.gtf`

**Command used** (on a high-RAM workstation):

```bash
~/workspace/scRNAseq/tools/STAR-2.7.11a/source/STAR \
  --runThreadN 60 \
  --runMode genomeGenerate \
  --genomeDir ~/workspace/scRNAseq/genome_index \
  --genomeFastaFiles ~/workspace/scRNAseq/genome/genomeV5.fa \
  --sjdbGTFfile ~/workspace/scRNAseq/genome/Gene_Models_210210.gtf \
  --sjdbOverhang 59 \
  --limitGenomeGenerateRAM 480000000000

```

	•	sjdbOverhang = read length - 1 (we used 59 for 60 bp reads)
	•	Output was saved to genome_index/


🎯 Read alignment with STARsolo

Read alignment and quantification were performed using STARsolo in CB_UMI_Complex mode with the following parameters:
	•	--soloType CB_UMI_Complex
	•	--soloCBstart 1 --soloCBlen 16
	•	--soloUMIstart 17 --soloUMIlen 12
	•	--soloFeatures Gene
	•	Filtered matrix output was written to Solo.out/Gene/filtered/ for each library.

To automate the alignment of all libraries (lib01 to lib06), we created a Bash script:

📄 scripts/run_STARsolo_all.sh

This script loops through all FASTQ files and runs STARsolo using the indexed genome and appropriate parameters.



👤 Author

Miquel Sendra
Guignard Lab, IBDM (Marseille)
📧 miquel.sendra-ortola@univ-amu.fr
🔬 Project with Tassos Pavlopoulos and Léo Guignard
