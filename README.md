# Pairwise Sequence Aligner for Homologous Proteins 🧬

A comprehensive Python-based bioinformatics tool built to perform optimal local sequence alignments. Designed to process biological sequence data, evaluate evolutionary relationships, and calculate protein homology with custom substitution matrices and penalties.

## 📌 Overview
In computational biology, simple string matching is insufficient for comparing protein sequences due to evolutionary mutations (substitutions, insertions, and deletions). This tool utilizes the `Biopython` library to find the best local alignment between two protein sequences (e.g., Human and Mouse DDIT3 sequences) using the **PAM250 substitution matrix**. 

## ⚙️ Technical Architecture
* **Language:** Python 3.x
* **Core Dependencies:** `Biopython`
* **Algorithm:** Local pairwise alignment (`pairwise2.align.localds`)
* **Scoring Matrix:** PAM250
* **Penalties:** * Gap Opening: `-9`
  * Gap Extension: `-8.5`

## 🚀 Features
* **Automated FASTA Parsing:** Effortlessly reads standard biological sequence formats (`.fasta`).
* **Biologically Accurate Scoring:** Uses PAM250 to score amino acid matches/mismatches based on evolutionary probability rather than simple text matching.
* **Homology Calculation:** Automatically parses the alignment output to calculate:
  * Sequence Identity (%)
  * Sequence Similarity (%)
  * Number of Gaps per sequence
* **Color-Coded Terminal Feedback:** Provides immediate visual feedback in the CLI based on homology likelihood:
  * 🟢 **Green (>35% ID):** Highly likely homologous.
  * 🟡 **Yellow (20-35% ID):** Potentially homologous (twilight zone).
  * 🔴 **Red (<20% ID):** Unlikely to be homologous.
* **Persistent Reporting:** Exports a detailed, formatted alignment report to `alignment.txt` for further scientific review.

## 📂 Repository Structure
* `sequence_alignment_macOS_eng.py` : Main alignment logic and execution script
* `Homologous-human.fasta` : Input sequence 1 (e.g., DDIT3 Human)
* `Homologous-mouse.fasta` : Input sequence 2 (e.g., DDIT3 Mouse)
* `22a4.fasta` : Additional test sequence (SLC22A4)
* `22a5.fasta` : Additional test sequence (SLC22A5)
* `SEC31_YEAST.fasta` : Additional test sequence (Yeast SEC31)
* `alignment.txt` : Auto-generated output report

## 💻 Installation & Setup

1. **Clone the repository:**
```bash
git clone [https://github.com/yourusername/pairwise-sequence-aligner.git](https://github.com/yourusername/pairwise-sequence-aligner.git)
cd pairwise-sequence-aligner
```

2. **Install dependencies:**
Ensure Python 3 is installed on your machine. Then, install the required `biopython` package:
```bash
pip install biopython
```

3. **Prepare your files:**
Place your input sequences in the root directory as `.fasta` files. If you are using your own sequences, update the file names in the script:
```python
seq1 = SeqIO.read('your_file_1.fasta', 'fasta').seq
seq2 = SeqIO.read('your_file_2.fasta', 'fasta').seq
```

## 🏃‍♂️ Usage
Execute the script via your terminal:
```bash
python sequence_alignment_macOS_eng.py
```

### 📄 Sample Output (`alignment.txt`)
Upon successful execution, the script generates a text file visualizing the sequence matches and a statistical summary:
```text
Score=690
Sequence 1 ID%: 88.17
Sequence 2 ID%: 88.69
Sequence 1 similarity%: 99.41
Sequence 2 similarity%: 100.00
Sequence 1 number of gaps: 0
Sequence 2 number of gaps: 1
```
