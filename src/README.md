# VATAligner

**VATAligner** (Versatile Alignment Tool) is a high-performance, multi-purpose tool designed for DNA and protein sequence alignments. It provides various options for efficient mapping, supporting both short and long nucleotide sequences and protein homology searches.

## Features
- **Fast and efficient**: Multi-threading capabilities for large-scale data processing.
- **Versatile**: Supports DNA and protein alignments.
- **Flexible**: Fine-tune parameters for maximum control over the alignment process.
## Command-line Options

### General Options
| Option               | Description                                   |
|----------------------|-----------------------------------------------|
| `-h`, `--help`        | Show help message.                           |
| `-p`, `--threads`     | Number of CPU threads (default: 1).          |
| `-d`, `--db`          | Specify the database file.                   |
| `-a`, `--vaa`         | VAT alignment archive (vatr) file.           |
| `--dbtype`            | Database type: `nucl` (nucleotide) or `prot` (protein). |

### Makedb Options
| Option              | Description                                      |
|---------------------|--------------------------------------------------|
| `-i`, `--in`        | Input reference file in FASTA format.            |

### Aligner Options
| Option                   | Description                                                                                  |
|--------------------------|----------------------------------------------------------------------------------------------|
| `-q`, `--query`           | Input query file.                                                                            |
| `-k`, `--maxtarget_seqs`  | Maximum number of target sequences to report (default: 25).                                  |
| `--top`                   | Report alignments within this percentage range of top alignment score (default: 98).         |
| `-e`, `--evalue`          | Maximum e-value to report (default: 0.001).                                                  |
| `--min_score`             | Minimum bit score to report alignments (default: 0).                                         |
| `--report_id`             | Minimum identity percentage to report alignments (default: 0).                               |
| `-t`, `--tmpdir`          | Directory for temporary files (default: `/dev/shm`).                                         |
| `--gapopen`               | Gap opening penalty (default: -1, which corresponds to 11 for protein).                      |
| `--gapextend`             | Gap extension penalty (default: -1, which corresponds to 1 for protein).                     |
| `--reward`                | Match reward score for `blastn` only (default: 2).                                           |
| `-S`, `--seed_len`        | Seed length (default: 15 for DNA, 8 for protein).                                            |
| `--penalty`               | Mismatch penalty score for `blastn` only (default: -3).                                      |
| `--match`                 | Match score (default: 5).                                                                    |
| `--mismatch`              | Mismatch score (default: -4).                                                                |
| `--simd_sort`             | Enable SIMD (AVX2) sorting for double-indexing.                                              |
| `--chimera`               | Enable chimera alignment.                                                                    |
| `--circ`                  | Enable circular alignment.                                                                   |
| `--wga`                   | Enable whole-genome alignment.                                                               |
| `--wgs`                   | Enable whole-genome sequencing.                                                              |
| `--splice`                | Enable splice alignments.                                                                    |
| `--dnah`                  | Enable DNA homology search.                                                                  |
| `--avx2`                  | Enable AVX2 hamming distance calculations.                                                   |
| `--spaced`                | Specify spaced seed (default: `null`).                                                       |
| `--matrix`                | Specify scoring matrix for protein alignment (default: `blosum62`).                          |

### Advanced Options
| Option                | Description                                                                                |
|-----------------------|--------------------------------------------------------------------------------------------|
| `-C`, `--max_hits`     | Maximum number of hits to consider for one seed (default: 0).                               |
| `--pre_filter`         | Minimum number of identities for pre-filter hit (default: 0).                               |
| `--xdrop`              | X-drop threshold for ungapped alignment (default: 18).                                      |
| `-X`, `--gapped_xdrop` | X-drop threshold for gapped alignment in bits (default: 18).                                |
| `--ungapped_score`     | Minimum raw alignment score to continue local extension (default: 0).                       |
| `--hit_score`          | Minimum score to keep a tentative alignment (default: 0).                                   |
| `--band`               | Band size for dynamic programming computation (default: 8).                                 |
| `--for_only`           | Enable alignment only on the forward strand.                                                |

## Example Usage

### DNA Alignment

1. **Create a nucleotide database**:
    ```bash
    VAT makevatdb --dbtype nucl --in test_all.fa -d mydb
    ```

2. **Run DNA alignment**:
    ```bash
    VAT dna -d mydb.vatf -q test_reads.fa -a alignment_output
    ```

3. **View the results**:
    ```bash
    VAT view -a alignment_output.vatr -o alignment_output
    vim alignment_output
    ```

### Protein Alignment

1. **Create a protein database**:
    ```bash
    VAT makevatdb --dbtype prot --in protein_ref.fa -d protein_db
    ```

2. **Run protein alignment**:
    ```bash
    VAT protein -d protein_db.vatf -q protein_test.fa -a protein_alignment -p 4
    ```

3. **View the results**:
    ```bash
    VAT view -a protein_alignment.vatr -o protein_alignment
    vim protein_alignment
    ```


