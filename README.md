# orf-finder
Find open reading frames in all six frames of DNA sequence

## Features
- Scans forward and reverse strands (6 frames)
- Adjustable minimum ORF length
- Outputs coordinates, frame, DNA length, and protein translation
- CSV export

## Installation
    git clone https://github.com/YOUR-USERNAME/orf-finder.git
    cd orf-finder
    pip install -r requirements.txt

## Usage
    python orf_finder.py data/example.fasta --min-aa 30 -o results.csv

## Example output
(paste a few rows of real output here)

## How it works
1. Reads FASTA with Biopython
2. Scans 3 frames on the forward strand and 3 on the reverse complement
3. Each ORF starts at ATG and ends at the first in-frame stop codon
4. Reverse-strand coordinates are converted back to the original sequence

## Testing
    pytest

## Validation
Results were compared against NCBI ORFfinder on [gene/genome name].

## License
MIT