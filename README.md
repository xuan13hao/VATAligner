# Versatile Alignment Tool 

**VATAligner** (Versatile Alignment Tool) is a fast and efficient multi-purpose sequence aligner. It supports the alignment of both short and long nucleotide sequences, as well as protein homology searches, offering a flexible solution for various sequence analysis needs.

## Features
- Supports both DNA and protein sequence alignments.
- High-performance alignment with multi-threading capabilities.
- Flexible input formats for FASTA/FASTQ sequences.
- Suitable for large-scale sequence analysis in genomics and proteomics.

## Installation

### Prerequisites
- **Boost Library**: Required for VATAligner to function properly.

### Steps to Install

1. Clone the repository and navigate to the `src/` directory:
    ```bash
    git clone https://github.com/xuan13hao/VATAligner.git
    cd VATAligner
    ```

2. Install the Boost library by running the pre-build script (optional):
    ```bash
    cd src
    ./pre-build_boost
    ```

3. Compile the tool:
    ```bash
    mkdir build
    cd build
    cmake ..
    make
    ```

## Performance

**VATAligner** leverages multi-threading to accelerate the alignment process for large datasets. 

## Data Availability

For information on data preparation and access, please refer to the [data_preparation/README.md](data_preparation/README.md) file.

## Help and Options

For a list of available options and command-line flags, refer to the detailed documentation in the [src/README.md](src/README.md) file.
