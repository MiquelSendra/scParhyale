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

## Install tools
1. Make project folders (optional but tidy)
```bash
mkdir -p ~/workspace/scRNAseq/tools
cd ~/workspace/scRNAseq/tools

```

2. Download & compile STAR locally
```bash
wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.11a.tar.gz
tar -xvzf 2.7.11a.tar.gz
cd STAR-2.7.11a/source
make STAR
```

This will create the STAR executable in:
```bash
~/workspace/scRNAseq/tools/STAR-2.7.11a/source/STAR

```

You can test it doing 
```bash
./STAR
```

3. Add it to your PATH (optional but useful)

Edit your ~/.bashrc (or ~/.zshrc if using zsh):
```bash
nano ~/.zshrc
```

Add this line at the end:
```bash
export PATH="$HOME/workspace/scRNAseq/tools/STAR-2.7.11a/source:$PATH"
```

then reload:
```bash
source ~/.zshrc
```

4. Test
```bash
STAR --version
```

## 🧬 Reference genome indexing

The reference genome was indexed using STAR with the following files:

- Genome FASTA: `genomeV5.fa`
- Gene annotations:
  - Old `Gene_Models_210210.gtf`
  - New (Matilde Paris) `Annotation_Reg_Embryo_500_cleanOverlap0.5_3_restranded_exp345.gtf`

**Command used** (on a high-RAM workstation):

```bash
~/workspace/scRNAseq/tools/STAR-2.7.11a/source/STAR \
  --runThreadN 60 \
  --runMode genomeGenerate \
  --genomeDir ~/workspace/scRNAseq/genome_index \
  --genomeFastaFiles ~/workspace/scRNAseq/genome/genomeV5.fa \
  --sjdbGTFfile ~/workspace/scRNAseq/genome/Annotation_Reg_Embryo_500_cleanOverlap0.5_3_restranded_exp345.gtf \
  --sjdbOverhang 59 \
  --limitGenomeGenerateRAM 480000000000

```

	•	sjdbOverhang = read length - 1 (we used 59 for 60 bp reads)
	•	Output was saved to genome_index/


##🎯 Read alignment with STARsolo

Read alignment and quantification were performed using STARsolo in CB_UMI_Complex mode with the following parameters:
	•	--soloType CB_UMI_Complex
	•	--soloCBstart 1 --soloCBlen 16
	•	--soloUMIstart 17 --soloUMIlen 12
	•	--soloFeatures Gene
	•	Filtered matrix output was written to Solo.out/Gene/filtered/ for each library.

To automate the alignment of all libraries (lib01 to lib06), we created a Bash script:

📄 scripts/run_STARsolo_all.sh

This script loops through all FASTQ files and runs STARsolo using the indexed genome and appropriate parameters.

Remember to give permissions before running and making sure all paths are correctly specified
```bash
chmod +x ~/workspace/scRNAseq/code/run_STARsolo_all.sh
```

In our workstation it took about 40h to complete with the above specified resources. Run with nohup so that you can disconect from the worktation or close the terminal and it will continue.

```bash
nohup ~/workspace/scRNAseq/code/run_STARsolo_all.sh > ~/workspace/scRNAseq/code/run_log.txt 2>&1 &
```

You can later check the progress...

```bash
tail -f ~/workspace/scRNAseq/code/run_log.txt
```

It will say sth like... processing library 1

##👤 Author

Miquel Sendra
Guignard Lab, IBDM (Marseille)
📧 miquel.sendra-ortola@univ-amu.fr
🔬 Project with Tassos Pavlopoulos and Léo Guignard
