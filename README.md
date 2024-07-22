# VATAligner

Versatile Alignment Tool (VAT) is a fast, multi-purpose aligner for short and long nucleotide sequences mapping and protein homology search. It aims to simplify integrative sequence analysis pipelines and handle sequencing data with previously-unseen characteristics, providing an efficient and high-performance alternative to traditional multi-purpose aligners and special-purpose sequence mapping tools.

## Installation

1. Navigate to the source directory:
    ```sh
    cd src/
    ```
2. Install Boost library:
    ```sh
    ./pre-build_boost
    ```

## Usage

### DNA Alignment
1. Create a nucleotide database:
    ```console
    VAT makevatdb --dbtype nucl --in ../data/test_all.fa -d mt
    ```
2. Perform DNA alignment:
    ```console
    VAT dna -d mt.vatf -q ../data/test_forward_reads.fa -a match
    VAT view -a match.vatr -o match
    vim match
    ```

#### Example with real data:
1. Build the nucleotide database:
    ```console
    make
    ./VAT makevatdb --dbtype nucl --in /home/xuan/test_data/test_data/marine.ref.fn -d marine.ref -b 1
    ```
2. Align DNA sequences using multiple threads:
    ```console
    /usr/bin/time -v ./VAT dna -d marine.ref -q /home/xuan/test_data/test_data/marine.sim.singular.fn/marine.sim.singular.fn -a marine.ref -p 4 -a match
    ```

### Protein Alignment
1. Create a protein database:
    ```console
    VAT makevatdb --dbtype prot --in ../data/pfam_ref.fa -d pfam
    ```
2. Perform protein alignment:
    ```console
    VAT protein -d pfam.vatf -q ../data/pfam_test.fa -a match -p 4
    VAT view -a match.vatr -o match
    ```

For more detailed usage and examples, please refer to our [documentation](https://github.com/xuan13hao/VATAligner.git).
