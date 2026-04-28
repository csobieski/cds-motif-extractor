#!/usr/bin/env python3

import argparse
import os
import re
from collections import defaultdict
from typing import Dict, List, Tuple, Optional

import gffutils
from Bio.Seq import Seq
from Bio.Data import CodonTable


def extract_genome_from_gff3(gff_file: str) -> Dict[str, str]:
    """Extract genome sequences from the ##FASTA section of a GFF3 file."""
    genome_seq: Dict[str, str] = {}
    fasta_flag = False
    seq_id: Optional[str] = None
    sequence_parts: List[str] = []

    with open(gff_file, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")

            if line.startswith("##FASTA"):
                fasta_flag = True
                continue

            if not fasta_flag:
                continue

            if line.startswith(">"):
                if seq_id is not None:
                    genome_seq[seq_id] = "".join(sequence_parts)
                seq_id = line.lstrip(">").split()[0]
                sequence_parts = []
            else:
                sequence_parts.append(line.strip())

    if seq_id is not None:
        genome_seq[seq_id] = "".join(sequence_parts)

    return genome_seq


def find_overlapping_matches(sequence: str, motif: str) -> List[Tuple[int, int]]:
    """Find all overlapping exact matches."""
    spans = []
    i = 0
    motif_len = len(motif)

    while motif_len > 0:
        i = sequence.find(motif, i)
        if i == -1:
            break
        spans.append((i, i + motif_len))
        i += 1

    return spans


def find_regex_matches(sequence: str, pattern: str) -> List[Tuple[int, int]]:
    """Find regex matches."""
    regex = re.compile(pattern)
    return [(m.start(), m.end()) for m in regex.finditer(sequence)]


def build_cds_groups(gff_file: str):
    """Group CDS features by parent gene/transcript."""
    db = gffutils.create_db(
        gff_file,
        dbfn=":memory:",
        force=True,
        keep_order=True,
        merge_strategy="merge",
        sort_attribute_values=True,
    )

    groups = defaultdict(list)

    for cds in db.features_of_type("CDS"):
        parent = cds.attributes.get("Parent", cds.attributes.get("ID", [cds.id]))[0]
        groups[(parent, cds.seqid, cds.strand)].append(cds)

    return groups


def assemble_cds_sequence(parts, genome: Dict[str, str]) -> Tuple[str, int, int]:
    """Assemble joined CDS features and orient sequence according to strand."""
    parts_sorted = sorted(parts, key=lambda f: f.start)
    seqid = parts_sorted[0].seqid
    strand = parts_sorted[0].strand

    pieces = []

    for feature in parts_sorted:
        contig = genome.get(seqid, "")
        if contig:
            pieces.append(contig[feature.start - 1:feature.end])

    cds_seq = "".join(pieces)

    if strand == "-":
        cds_seq = str(Seq(cds_seq).reverse_complement())

    cds_start = min(f.start for f in parts_sorted)
    cds_end = max(f.end for f in parts_sorted)

    return cds_seq, cds_start, cds_end


def extract_matching_cds(
    gff_file: str,
    query_seq: str,
    search_in: str = "dna",
    use_regex: bool = False,
    transl_table: int = 11,
    padding: int = 0,
) -> List[str]:
    """
    Extract CDS regions containing a DNA or protein motif.

    search_in:
        dna     search motif in nucleotide CDS sequence
        protein search motif in translated CDS sequence
    """

    if search_in not in {"dna", "protein"}:
        raise ValueError("search_in must be either 'dna' or 'protein'")

    genome = extract_genome_from_gff3(gff_file)
    cds_groups = build_cds_groups(gff_file)

    basename = os.path.basename(gff_file).replace(".gff3", "")
    fasta_entries = []

    motif = query_seq.upper()

    for (parent, seqid, strand), parts in cds_groups.items():
        if seqid not in genome:
            continue

        cds_seq, cds_start, cds_end = assemble_cds_sequence(parts, genome)

        if search_in == "dna":
            search_sequence = cds_seq.upper()
            unit = "nt"
        else:
            try:
                search_sequence = str(
                    Seq(cds_seq).translate(table=transl_table, to_stop=False)
                ).upper()
            except CodonTable.TranslationError:
                trim = len(cds_seq) % 3
                trimmed_seq = cds_seq[:-trim] if trim else cds_seq
                search_sequence = str(
                    Seq(trimmed_seq).translate(table=transl_table, to_stop=False)
                ).upper()
            unit = "aa"

        if use_regex:
            spans = find_regex_matches(search_sequence, query_seq)
        else:
            spans = find_overlapping_matches(search_sequence, motif)

        if not spans:
            continue

        contig_len = len(genome[seqid])

        if strand == "+":
            region_start = max(1, cds_start - padding)
            region_end = cds_end
        else:
            region_start = cds_start
            region_end = min(contig_len, cds_end + padding)

        region_seq = genome[seqid][region_start - 1:region_end]

        if strand == "-":
            region_seq = str(Seq(region_seq).reverse_complement())

        attrs = parts[0].attributes
        name_bits = []

        for key in ("gene", "Name", "product", "locus_tag"):
            if key in attrs:
                name_bits.extend(attrs[key])

        gene_name = "|".join(name_bits) if name_bits else parent
        match_str = ",".join(f"{start}-{end}" for start, end in spans)

        header = (
            f">{basename}|parent={parent}|name={gene_name}"
            f"|contig={seqid}|coords={cds_start}-{cds_end}|strand={strand}"
            f"|search={search_in}{'(regex)' if use_regex else ''}"
            f"|motif={query_seq}|matches_{unit}={match_str}|padding={padding}"
        )

        fasta_entries.append(f"{header}\n{region_seq}\n")

    return fasta_entries


def extract_from_directory(
    gff_dir: str,
    query_seq: str,
    output_file: str,
    search_in: str = "dna",
    use_regex: bool = False,
    padding: int = 0,
    transl_table: int = 11,
) -> Tuple[List[str], List[str]]:

    gff_files = [
        os.path.join(gff_dir, filename)
        for filename in os.listdir(gff_dir)
        if filename.endswith(".gff3")
    ]

    present = []
    missing = []

    with open(output_file, "w") as out_f:
        for gff_file in sorted(gff_files):
            basename = os.path.basename(gff_file)

            try:
                entries = extract_matching_cds(
                    gff_file=gff_file,
                    query_seq=query_seq,
                    search_in=search_in,
                    use_regex=use_regex,
                    transl_table=transl_table,
                    padding=padding,
                )

                if entries:
                    present.append(basename)
                    for entry in entries:
                        out_f.write(entry + "\n")
                else:
                    missing.append(basename)

            except Exception as exc:
                missing.append(f"{basename} (error: {exc})")

    print(f"Search motif: {query_seq}")
    print(f"Search mode: {search_in}")
    print(f"Output FASTA: {output_file}\n")

    if present:
        print("Files with at least one match:")
        for filename in present:
            print(f"  - {filename}")

    if missing:
        print("\nFiles with no match or error:")
        for filename in missing:
            print(f"  - {filename}")

    print("\nFinished.")

    return present, missing


def main():
    parser = argparse.ArgumentParser(
        description="Extract CDS sequences from GFF3 files based on DNA or protein motif matches."
    )

    parser.add_argument(
        "-i", "--input-dir",
        required=True,
        help="Directory containing GFF3 files with embedded FASTA sections."
    )

    parser.add_argument(
        "-q", "--query",
        required=True,
        help="DNA or protein motif to search."
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output FASTA file."
    )

    parser.add_argument(
        "--search-in",
        choices=["dna", "protein"],
        default="dna",
        help="Search motif in nucleotide CDS sequence or translated protein sequence. Default: dna."
    )

    parser.add_argument(
        "--regex",
        action="store_true",
        help="Interpret query as a regular expression."
    )

    parser.add_argument(
        "--padding",
        type=int,
        default=0,
        help="Number of upstream nucleotides to include before the CDS. Default: 0."
    )

    parser.add_argument(
        "--translation-table",
        type=int,
        default=11,
        help="NCBI translation table. Default: 11, bacterial."
    )

    args = parser.parse_args()

    extract_from_directory(
        gff_dir=args.input_dir,
        query_seq=args.query,
        output_file=args.output,
        search_in=args.search_in,
        use_regex=args.regex,
        padding=args.padding,
        transl_table=args.translation_table,
    )


if __name__ == "__main__":
    main()
