import argparse
import csv
import sys

from Bio import SeqIO
from Bio.Seq import Seq

START = "ATG"
STOPS = {"TAA", "TAG", "TGA"}


def find_orfs(seq, min_aa=30):
    seq = Seq(str(seq).upper())
    length = len(seq)
    orfs = []

    for strand, nuc in ((+1, seq), (-1, seq.reverse_complement())):
        nuc = str(nuc)
        for frame in range(3):
            i = frame
            while i + 3 <= length:
                if nuc[i:i + 3] == START:
                    j = i
                    while j + 3 <= length and nuc[j:j + 3] not in STOPS:
                        j += 3
                    if j + 3 > length:      
                        break
                    end = j + 3           
                    protein = str(Seq(nuc[i:j]).translate())
                    if len(protein) >= min_aa:
                        if strand == 1:
                            start_pos, end_pos = i + 1, end
                        else:               
                            start_pos, end_pos = length - end + 1, length - i
                        orfs.append({
                            "strand": "+" if strand == 1 else "-",
                            "frame": frame + 1,
                            "start": start_pos,
                            "end": end_pos,
                            "length_nt": end - i,
                            "length_aa": len(protein),
                            "protein": protein,
                        })
                    i = end                
                    continue
                i += 3
    return orfs


def main():
    parser = argparse.ArgumentParser(description="Six-frame ORF finder")
    parser.add_argument("fasta", help="input FASTA file")
    parser.add_argument("-m", "--min-aa", type=int, default=30,
                        help="minimum ORF length in amino acids (default: 30)")
    parser.add_argument("-o", "--output", help="write results to this CSV file")
    args = parser.parse_args()

    rows = []
    for record in SeqIO.parse(args.fasta, "fasta"):
        for orf in find_orfs(record.seq, args.min_aa):
            rows.append({"sequence_id": record.id, **orf})

    rows.sort(key=lambda r: r["length_aa"], reverse=True)

    fields = ["sequence_id", "strand", "frame", "start", "end",
              "length_nt", "length_aa", "protein"]
    out = open(args.output, "w", newline="") if args.output else sys.stdout
    writer = csv.DictWriter(out, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)
    if args.output:
        out.close()
        print(f"Found {len(rows)} ORFs. Saved to {args.output}")


if __name__ == "__main__":
    main()