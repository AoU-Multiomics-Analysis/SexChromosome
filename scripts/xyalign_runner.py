#!/usr/bin/env python
"""Validate task-local inputs and run unchanged XYalign commands."""
from __future__ import print_function

import argparse
import contextlib
import logging
import os
import re
import shutil
import subprocess
import sys
import tarfile

LOG = logging.getLogger('xyalign_runner')
TOKEN = re.compile(r'^[A-Za-z0-9_][A-Za-z0-9_.:-]*$')
INDEX_SUFFIXES = ['.fai', '.dict', '.amb', '.ann', '.bwt', '.pac', '.sa']


def require_file(path):
    if not path or re.match(r'^[A-Za-z][A-Za-z0-9+.-]*://', path):
        raise ValueError('Localization error: supply a readable local file: {0}'.format(path))
    path = os.path.abspath(path)
    if not os.path.isfile(path) or not os.access(path, os.R_OK):
        raise ValueError('Localization error: file is not readable: {0}'.format(path))
    if os.path.getsize(path) == 0:
        raise ValueError('Input file is empty: {0}'.format(path))
    return path


def safe_token(value, label):
    if not TOKEN.match(value) or value in ('.', '..'):
        raise ValueError('Invalid {0}: use letters, numbers, underscores, dots, colons, or hyphens'.format(label))
    return value


def stage_file(source, target, copy=False):
    source = require_file(source)
    if os.path.lexists(target):
        raise ValueError('Task-local input already exists: {0}'.format(target))
    if copy:
        shutil.copyfile(source, target)
    else:
        os.symlink(source, target)
    return target


def run(command):
    # Each argument is separate. Upstream receives only controlled relative paths.
    LOG.info('Run: %s', repr(command))
    subprocess.check_call(command)


def read_chromosomes(path, x_chromosome, y_chromosome):
    with open(require_file(path)) as handle:
        chromosomes = [line.strip() for line in handle if line.strip()]
    if len(chromosomes) != len(set(chromosomes)):
        raise ValueError('Chromosome names must be unique')
    for chrom in chromosomes:
        safe_token(chrom, 'chromosome name')
    if x_chromosome == y_chromosome:
        raise ValueError('X and Y chromosome names must differ')
    if x_chromosome not in chromosomes or y_chromosome not in chromosomes:
        raise ValueError('The chromosome list must contain X and Y')
    if not set(chromosomes) - {x_chromosome, y_chromosome}:
        raise ValueError('The chromosome list must contain at least one autosome')
    return chromosomes


def stage_reference(args):
    stage_file(args.reference, 'reference.fa')
    if args.reference_index:
        stage_file(args.reference_index, 'reference.fa.fai', copy=True)
    else:
        run(['samtools', 'faidx', 'reference.fa'])
    os.utime('reference.fa.fai', None)
    import pysam
    with pysam.FastaFile('reference.fa') as reference:
        contigs = list(zip(reference.references, reference.lengths))
    for name, length in contigs:
        safe_token(name, 'reference contig')
        if length < 1:
            raise ValueError('Reference contigs must have positive lengths')
    if args.y_chromosome not in dict(contigs):
        raise ValueError('Y chromosome is absent from the reference')
    return contigs


def stage_bam(args, reference_contigs, chromosomes):
    stage_file(args.bam, 'input.bam')
    stage_file(args.bam_index, 'input.bam.bai', copy=True)
    os.utime('input.bam.bai', None)
    run(['samtools', 'quickcheck', '-v', 'input.bam'])
    import pysam
    with pysam.AlignmentFile('input.bam', 'rb') as bam:
        bam.check_index()
        if bam.header.to_dict().get('HD', {}).get('SO') != 'coordinate':
            raise ValueError('The BAM must declare coordinate sort order')
        contigs = list(zip(bam.references, bam.lengths))
        if contigs != reference_contigs:
            raise ValueError('BAM and reference contig names, lengths, and order must match')
        if not set(chromosomes).issubset(set(bam.references)):
            raise ValueError('An analysis chromosome is absent from the BAM')
        # Upstream inserts read-group values into shell strings during remapping.
        if args.mode == 'remap':
            for group in bam.header.to_dict().get('RG', []):
                for value in group.values():
                    if re.search(r"[^A-Za-z0-9_.:/ +,=-]", str(value)):
                        raise ValueError('A read-group value contains unsupported shell characters')


def validate_mask(path, contigs):
    lengths = dict(contigs)
    with open(require_file(path)) as handle:
        for line in handle:
            if not line.strip() or line.startswith(('#', 'track ', 'browser ')):
                continue
            fields = line.split()
            if len(fields) < 3 or fields[0] not in lengths:
                raise ValueError('The reference mask must contain valid BED intervals')
            start, end = int(fields[1]), int(fields[2])
            if not 0 <= start < end <= lengths[fields[0]]:
                raise ValueError('A reference mask interval is outside its contig')


def archive_reports(output, paths):
    with tarfile.open(output, 'w:gz') as archive:
        for path in paths:
            if os.path.exists(path):
                archive.add(path, arcname=path)


def extract_references(path, destination):
    allowed = {'references/' + prefix + '.fa' + suffix
               for prefix in ('xx', 'xy') for suffix in [''] + INDEX_SUFFIXES}
    with tarfile.open(require_file(path), 'r:') as archive:
        members = archive.getmembers()
        if not members or any(not member.isfile() or member.name not in allowed for member in members):
            raise ValueError('Invalid prepared-reference archive')
        if {member.name for member in members} != allowed:
            raise ValueError('Prepared-reference archive is incomplete')
        # Members are regular files with fixed names. No directories or links are accepted.
        for member in members:
            target = os.path.join(destination, member.name)
            if not os.path.isdir(os.path.dirname(target)):
                os.makedirs(os.path.dirname(target))
            with contextlib.closing(archive.extractfile(member)) as source, open(target, 'wb') as output:
                shutil.copyfileobj(source, output)
    for prefix in ('xx', 'xy'):
        for suffix in INDEX_SUFFIXES:
            os.utime(os.path.join(destination, 'references', prefix + '.fa' + suffix), None)


def analyze(args):
    chromosomes = read_chromosomes(args.chromosomes_file, args.x_chromosome, args.y_chromosome)
    contigs = stage_reference(args)
    stage_bam(args, contigs, chromosomes)
    common = ['--bam', 'input.bam', '--ref', 'reference.fa', '--sample_id', args.sample_id,
              '--cpus', str(args.cpus), '--chromosomes'] + chromosomes
    run(['xyalign', '--CHROM_STATS', '--output_dir', 'chrom-stats'] + common)
    for source, target in [('depth', 'chromosome_depth.tsv'), ('mapq', 'chromosome_mapq.tsv')]:
        path = 'chrom-stats/results/{0}_chrom_stats_{1}.txt'.format(args.sample_id, source)
        shutil.copyfile(require_file(path), os.path.join(args.output_dir, target))
    reports = ['chrom-stats']
    if args.characterize:
        command = ['xyalign', '--CHARACTERIZE_SEX_CHROMS', '--output_dir', 'characterization',
                   '--x_chromosome', args.x_chromosome, '--y_chromosome', args.y_chromosome,
                   '--window_size', str(args.window_size), '--num_permutations', str(args.num_permutations),
                   '--num_bootstraps', str(args.num_bootstraps)] + common
        if not args.call_variants:
            command += ['--platypus_calling', 'none', '--no_variant_plots']
        run(command)
        require_file('characterization/bed/{0}_full_dataframe_depth_mapq_preprocessing.csv'.format(args.sample_id))
        reports.append('characterization')
    archive_reports(os.path.join(args.output_dir, 'reports.tar.gz'), reports)


def prepare_reference(args):
    contigs = stage_reference(args)
    os.mkdir('references')
    command = ['xyalign', '--PREPARE_REFERENCE', '--ref', 'reference.fa',
               '--output_dir', 'preparation', '--cpus', str(args.cpus),
               '--y_chromosome', args.y_chromosome, '--xx_ref_out', 'references/xx.fa',
               '--xy_ref_out', 'references/xy.fa', '--bwa_index', 'True']
    if args.reference_mask:
        validate_mask(args.reference_mask, contigs)
        stage_file(args.reference_mask, 'mask.bed')
        command += ['--reference_mask', 'mask.bed']
    run(command)
    with tarfile.open(os.path.join(args.output_dir, 'references.tar'), 'w') as archive:
        for prefix in ('xx', 'xy'):
            for suffix in [''] + INDEX_SUFFIXES:
                path = 'references/' + prefix + '.fa' + suffix
                archive.add(require_file(path), arcname=path)


def remap(args):
    contigs = stage_reference(args)
    stage_bam(args, contigs, [args.x_chromosome, args.y_chromosome])
    extract_references(args.prepared_references, '.')
    import pysam
    for prefix in ('xx', 'xy'):
        with pysam.FastaFile('references/' + prefix + '.fa') as reference:
            if list(zip(reference.references, reference.lengths)) != contigs:
                raise ValueError('Prepared reference does not match the original reference dictionary')
    command = ['xyalign', '--REMAPPING', '--bam', 'input.bam', '--ref', 'reference.fa',
               '--xx_ref_in', 'references/xx.fa', '--xy_ref_in', 'references/xy.fa',
               '--sample_id', args.sample_id, '--output_dir', 'remapping',
               '--x_chromosome', args.x_chromosome, '--y_chromosome', args.y_chromosome,
               '--cpus', str(args.cpus), '--xmx', args.java_memory,
               '--shufflesh_path', 'shuffle.sh',
               '--y_present' if args.y_status == 'present' else '--y_absent']
    if args.single_end:
        command.append('--single_end')
    run(command)
    source = 'remapping/bam/{0}.merged.bam'.format(args.sample_id)
    require_file(source)
    run(['samtools', 'quickcheck', '-v', source])
    target = os.path.join(args.output_dir, 'remapped.bam')
    shutil.move(source, target)
    # Check every alignment. Upstream does not propagate all subprocess failures.
    with pysam.AlignmentFile(target, 'rb') as bam:
        for record in bam.fetch(until_eof=True):
            pass
    run(['samtools', 'index', '-@', str(args.cpus), target])
    with pysam.AlignmentFile(target, 'rb') as bam:
        bam.check_index()
    archive_reports(os.path.join(args.output_dir, 'reports.tar.gz'), ['remapping/logfiles'])


def make_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest='mode')
    for name in ('analyze', 'prepare-reference', 'remap'):
        sub = subparsers.add_parser(name)
        sub.add_argument('--reference', required=True)
        sub.add_argument('--reference-index')
        sub.add_argument('--y-chromosome', default='chrY')
        sub.add_argument('--cpus', type=int, default=2)
        sub.add_argument('--output-dir', required=True)
        if name != 'prepare-reference':
            sub.add_argument('--bam', required=True)
            sub.add_argument('--bam-index', required=True)
            sub.add_argument('--sample-id', required=True)
            sub.add_argument('--x-chromosome', default='chrX')
        if name == 'analyze':
            sub.add_argument('--chromosomes-file', required=True)
            sub.add_argument('--characterize', action='store_true')
            sub.add_argument('--remapping-requested', action='store_true')
            sub.add_argument('--y-status', default='unknown')
            sub.add_argument('--call-variants', action='store_true')
            sub.add_argument('--window-size', type=int, default=5000)
            sub.add_argument('--num-permutations', type=int, default=10000)
            sub.add_argument('--num-bootstraps', type=int, default=10000)
        if name == 'prepare-reference':
            sub.add_argument('--reference-mask')
        if name == 'remap':
            sub.add_argument('--prepared-references', required=True)
            sub.add_argument('--y-status', required=True, choices=['present', 'absent'])
            sub.add_argument('--single-end', action='store_true')
            sub.add_argument('--java-memory', default='8g')
    return parser


def main(argv=None):
    parser = make_parser()
    args = parser.parse_args(argv)
    if args.mode is None:
        parser.error('A subcommand is required')
    # Resolve localized inputs before changing directory.
    for name in ('reference', 'reference_index', 'bam', 'bam_index', 'chromosomes_file',
                 'reference_mask', 'prepared_references'):
        value = getattr(args, name, None)
        if value is not None:
            setattr(args, name, require_file(value))
    safe_token(args.y_chromosome, 'Y chromosome')
    if hasattr(args, 'sample_id'):
        safe_token(args.sample_id, 'sample ID')
        safe_token(args.x_chromosome, 'X chromosome')
        if args.x_chromosome == args.y_chromosome:
            raise ValueError('X and Y chromosome names must differ')
    if args.cpus < 1:
        raise ValueError('CPU count must be positive')
    if args.mode == 'analyze':
        if args.remapping_requested and args.y_status not in ('present', 'absent'):
            raise ValueError('Remapping requires an explicit present or absent Y status')
        if min(args.window_size, args.num_permutations, args.num_bootstraps) < 1:
            raise ValueError('Window size and replicate counts must be positive')
        if args.call_variants and not args.characterize:
            raise ValueError('Variant calling requires --characterize')
    if args.mode == 'remap' and not re.match(r'^[1-9][0-9]*[mMgG]$', args.java_memory):
        raise ValueError('Java memory must use a positive integer with m or g')
    args.output_dir = os.path.abspath(args.output_dir)
    os.makedirs(args.output_dir)
    work = os.path.join(args.output_dir, 'work')
    os.mkdir(work)
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    for handler in [logging.StreamHandler(), logging.FileHandler(os.path.join(args.output_dir, 'runner.log'))]:
        handler.setFormatter(formatter)
        root_logger.addHandler(handler)
    LOG.info('Start %s; upstream XYalign is unchanged', args.mode)
    run(['xyalign', '--version'])
    previous = os.getcwd()
    try:
        os.chdir(work)
        {'analyze': analyze, 'prepare-reference': prepare_reference, 'remap': remap}[args.mode](args)
        LOG.info('Completed %s', args.mode)
    finally:
        os.chdir(previous)


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        LOG.error('%s', error)
        sys.exit(1)
