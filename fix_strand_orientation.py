#!/usr/bin/env python3
import argparse
import subprocess
import tempfile
import os
import gzip
import pysam
from Bio.Seq import Seq

def run_minimap2(reference_fasta, input_fastq, sam_output):
    cmd = [
        "minimap2",
        "-ax",  # Output in SAM format
        "lr:hq", # High-quality long-read model
        reference_fasta,
        input_fastq
    ]
    with open(sam_output, "w") as sam_file:
        subprocess.run(cmd, stdout=sam_file, check=True)

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

def reverse_quality(qual):
    return qual[::-1]

def open_output_fastq(path):
    return gzip.open(path, "wt") if path.endswith(".gz") else open(path, "w")

def correct_reads_from_sam(sam_path, output_fastq_path, include_unmapped=False, unmapped_output_path=None):
    samfile = pysam.AlignmentFile(sam_path, "r")
    with open_output_fastq(output_fastq_path) as out_fq:
        unmapped_fq = open_output_fastq(unmapped_output_path) if unmapped_output_path else None
        seen_reads = set() # Only use primary alignment, exclude all others

        for read in samfile:
            if read.query_name in seen_reads:
                continue
            seen_reads.add(read.query_name)

            if read.is_unmapped or read.query_sequence is None:
                if include_unmapped and read.query_sequence:
                    out_fq.write(f"@{read.query_name}_unmapped\n{read.query_sequence}\n+\n{read.qual}\n")
                elif unmapped_fq and read.query_sequence:
                    unmapped_fq.write(f"@{read.query_name}\n{read.query_sequence}\n+\n{read.qual}\n")
                continue

            if read.is_secondary or read.is_supplementary:
                continue

            seq = read.query_sequence
            qual = read.qual

            if read.is_reverse: # Love that this exists
                seq = reverse_complement(seq)
                qual = reverse_quality(qual)

            out_fq.write(f"@{read.query_name}\n{seq}\n+\n{qual}\n")

        if unmapped_fq:
            unmapped_fq.close()

    samfile.close()

def main():
    parser = argparse.ArgumentParser(description="Align FASTQ to FASTA using minimap2 and correct strand orientation.")
    parser.add_argument("input_fastq", help="Input FASTQ (.fastq or .fastq.gz)")
    parser.add_argument("reference_fasta", help="Reference FASTA file")
    parser.add_argument("output_fastq", help="Corrected output FASTQ (.fastq or .fastq.gz)")
    parser.add_argument("--include-unmapped", action="store_true", help="Include unmapped reads in output FASTQ")
    parser.add_argument("--unmapped-output", help="Write unmapped reads to separate FASTQ (.fastq or .fastq.gz)")

    args = parser.parse_args()

    if args.include_unmapped and args.unmapped_output:
        parser.error("Options --include-unmapped and --unmapped-output are mutually exclusive.")

    with tempfile.TemporaryDirectory() as tmpdir:
        sam_path = os.path.join(tmpdir, "aligned.sam")
        print("Running minimap2 alignment...")
        run_minimap2(args.reference_fasta, args.input_fastq, sam_path)
        print("Fixing strand orientation...")
        correct_reads_from_sam(
            sam_path,
            args.output_fastq,
            include_unmapped=args.include_unmapped,
            unmapped_output_path=args.unmapped_output
        )

    print(f"Corrected FASTQ written to: {args.output_fastq}")
    if args.unmapped_output:
        print(f"Unmapped reads written to: {args.unmapped_output}")

if __name__ == "__main__":
    main()
