#!/usr/bin/env python
"""Run the real XYalign commands inside the CI image (Python 2.7 or 3)."""
from __future__ import print_function

import argparse
import csv
import math
import os
import random
import subprocess
import sys
import tarfile

import pysam


CHROMOSOMES = ("chr1", "chrX", "chrY")
CHROMOSOME_LENGTH = 20000
WINDOW_SIZE = 5000
READ_LENGTH = 100
FRAGMENT_LENGTH = 300


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def run(command):
    print("[smoke] Running: {}".format(" ".join(command)))
    sys.stdout.flush()
    subprocess.check_call(command)


def write_fixture(directory):
    """Create a FASTA and sorted, indexed paired-end BAM with read groups."""
    rng = random.Random(7319)
    sequences = dict((chrom, "".join(rng.choice("ACGT")
                                    for _ in range(CHROMOSOME_LENGTH)))
                     for chrom in CHROMOSOMES)
    reference = os.path.join(directory, "reference.fa")
    with open(reference, "w") as handle:
        for chrom in CHROMOSOMES:
            handle.write(">{}\n".format(chrom))
            sequence = sequences[chrom]
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")
    pysam.faidx(reference)

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": chrom, "LN": CHROMOSOME_LENGTH}
               for chrom in CHROMOSOMES],
        "RG": [{"ID": "smoke", "SM": "smoke", "PL": "ILLUMINA",
                "LB": "smoke", "PU": "smoke"}],
    }
    unsorted_bam = os.path.join(directory, "unsorted.bam")
    read_names = set()
    with pysam.AlignmentFile(unsorted_bam, "wb", header=header) as bam:
        for reference_id, chrom in enumerate(CHROMOSOMES):
            for window in range(CHROMOSOME_LENGTH // WINDOW_SIZE):
                # Every analysis window has reads. Counts vary between windows.
                for number in range(5 + reference_id + window):
                    start = window * WINDOW_SIZE + 300 + number * 37
                    mate_start = start + FRAGMENT_LENGTH - READ_LENGTH
                    name = "{}_{}_{}".format(chrom, window, number)
                    read_names.add(name)
                    for is_second in (False, True):
                        position = mate_start if is_second else start
                        read = pysam.AlignedSegment()
                        read.query_name = name
                        # SAM stores the reverse-strand read in reference order.
                        read.query_sequence = sequences[chrom][
                            position:position + READ_LENGTH]
                        read.flag = 147 if is_second else 99
                        read.reference_id = reference_id
                        read.reference_start = position
                        read.mapping_quality = 60
                        read.cigar = [(0, READ_LENGTH)]
                        read.next_reference_id = reference_id
                        read.next_reference_start = start if is_second else mate_start
                        read.template_length = -FRAGMENT_LENGTH if is_second else FRAGMENT_LENGTH
                        read.query_qualities = pysam.qualitystring_to_array("I" * READ_LENGTH)
                        read.set_tag("RG", "smoke")
                        read.set_tag("NM", 0)
                        bam.write(read)
    bam_path = os.path.join(directory, "input.bam")
    pysam.sort("-o", bam_path, unsorted_bam)
    pysam.index(bam_path)
    os.remove(unsorted_bam)
    chromosomes = os.path.join(directory, "chromosomes.txt")
    with open(chromosomes, "w") as handle:
        handle.write("\n".join(CHROMOSOMES) + "\n")
    mask = os.path.join(directory, "mask.bed")
    with open(mask, "w") as handle:
        handle.write("chrX\t100\t200\n")
    return reference, bam_path, chromosomes, mask, read_names


def check_file(path):
    require(os.path.isfile(path) and os.path.getsize(path) > 0,
            "Missing or empty output: {}".format(path))


def report_members(directory):
    archive = os.path.join(directory, "reports.tar.gz")
    check_file(archive)
    check_file(os.path.join(directory, "runner.log"))
    with tarfile.open(archive, "r:gz") as handle:
        members = [member.name for member in handle.getmembers() if member.isfile()]
    require(bool(members), "The report archive has no files")
    return members


def check_stats(directory):
    for filename in ("chromosome_depth.tsv", "chromosome_mapq.tsv"):
        path = os.path.join(directory, filename)
        check_file(path)
        with open(path) as handle:
            rows = list(csv.reader(handle, delimiter="\t"))
        require(len(rows) == 4, "Expected one header and three chromosome rows")
        require(set(row[0] for row in rows[1:]) == set(CHROMOSOMES),
                "Statistics do not contain all fixture chromosomes")
        for row in rows[1:]:
            values = [float(value) for value in row[1:]]
            require(values and all(value > 0 and not math.isinf(value)
                                   and not math.isnan(value) for value in values),
                    "Statistics contain a missing, zero, or non-finite value")


def check_prepared_references(directory):
    archive = os.path.join(directory, "references.tar")
    check_file(archive)
    check_file(os.path.join(directory, "runner.log"))
    with tarfile.open(archive) as handle:
        members = [member for member in handle.getmembers() if member.isfile()]
        fasta_members = [member for member in members
                         if member.name.endswith((".fa", ".fasta"))]
        require(len(fasta_members) == 2, "Expected both prepared reference FASTAs")
        y_sequences = []
        for member in fasta_members:
            content = handle.extractfile(member).read().decode("ascii")
            contigs = {}
            name = None
            for line in content.splitlines():
                if line.startswith(">"):
                    name = line[1:].split()[0]
                    contigs[name] = ""
                else:
                    contigs[name] += line
            require(set(contigs) == set(CHROMOSOMES), "Prepared reference lost contigs")
            require(contigs["chrX"][100:200].upper() == "N" * 100,
                    "The supplied reference mask was not applied")
            y_sequences.append(contigs["chrY"].upper())
        require(sum(sequence == "N" * CHROMOSOME_LENGTH
                    for sequence in y_sequences) == 1,
                "Expected exactly one reference with chromosome Y masked")
    return archive


def check_remapped(directory, y_status, read_names):
    bam_path = os.path.join(directory, "remapped.bam")
    check_file(bam_path)
    check_file(bam_path + ".bai")
    report_members(directory)
    run(["samtools", "quickcheck", "-v", bam_path])
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        require(tuple(bam.references) == CHROMOSOMES,
                "Remapped BAM changed the reference dictionary")
        require(bool(bam.header.to_dict().get("RG")), "Remapped BAM lost read groups")
        reads = list(bam.fetch(until_eof=True))
    primary = [read for read in reads if not read.is_secondary and not read.is_supplementary]
    require(len(primary) == 2 * len(read_names), "Remapping lost or duplicated primary reads")
    require(set(read.query_name for read in primary) == read_names,
            "Remapping changed the read names")
    require(all(read.is_paired for read in primary), "Remapping lost paired-end flags")
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        require(sum(1 for _ in bam.fetch("chr1")) > 0,
                "Remapped BAM lost the autosomal alignments")
        require(sum(1 for _ in bam.fetch("chrX")) > 0,
                "Chromosome X reads did not remap")
        y_count = sum(1 for _ in bam.fetch("chrY"))
    require((y_count > 0) if y_status == "present" else (y_count == 0),
            "Chromosome Y mappings do not match the selected reference")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", default="/opt/sexchromosome/scripts/xyalign_runner.py")
    parser.add_argument("--work-dir", required=True)
    args = parser.parse_args()
    directory = os.path.abspath(args.work_dir)
    if not os.path.isdir(directory):
        os.makedirs(directory)
    reference, bam_path, chromosomes, mask, read_names = write_fixture(directory)
    runner = [sys.executable, args.runner]
    reference_args = ["--reference", reference, "--reference-index", reference + ".fai"]
    bam_args = ["--bam", bam_path, "--bam-index", bam_path + ".bai"]

    # Also check that analysis works when WDL does not supply an optional FAI.
    basic_output = os.path.join(directory, "analyze")
    run(runner + ["analyze"] + bam_args + ["--reference", reference,
        "--sample-id", "smoke", "--chromosomes-file", chromosomes,
        "--x-chromosome", "chrX", "--y-chromosome", "chrY", "--cpus", "2",
        "--output-dir", basic_output])
    check_stats(basic_output)
    report_members(basic_output)

    characterized_output = os.path.join(directory, "characterize")
    run(runner + ["analyze"] + bam_args + reference_args + [
        "--sample-id", "smoke", "--chromosomes-file", chromosomes,
        "--x-chromosome", "chrX", "--y-chromosome", "chrY", "--cpus", "2",
        "--output-dir", characterized_output, "--characterize",
        "--window-size", str(WINDOW_SIZE), "--num-permutations", "10",
        "--num-bootstraps", "10"])
    check_stats(characterized_output)
    members = report_members(characterized_output)
    for suffix in ("permutation_results.txt", "ks_results.txt", "bootstrap_results.txt"):
        require(any(name.endswith(suffix) for name in members),
                "Characterization report missing: {}".format(suffix))

    prepared_output = os.path.join(directory, "prepare")
    run(runner + ["prepare-reference"] + reference_args + [
        "--reference-mask", mask, "--y-chromosome", "chrY", "--cpus", "2",
        "--output-dir", prepared_output])
    prepared = check_prepared_references(prepared_output)
    for status in ("present", "absent"):
        output = os.path.join(directory, "remap_" + status)
        run(runner + ["remap"] + bam_args + reference_args + [
            "--prepared-references", prepared, "--sample-id", "smoke",
            "--x-chromosome", "chrX", "--y-chromosome", "chrY",
            "--y-status", status, "--cpus", "2", "--output-dir", output,
            "--java-memory", "1g"])
        check_remapped(output, status, read_names)
    print("[smoke] Analysis, characterization, reference preparation, and both remap branches passed")


if __name__ == "__main__":
    main()
