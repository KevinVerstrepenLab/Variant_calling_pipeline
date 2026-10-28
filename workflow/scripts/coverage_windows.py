#!/usr/bin/env python3
"""Prepare fixed-size BED windows and reference GC/ambiguity metadata."""
import argparse
import csv
import gzip
import json
from pathlib import Path


def fasta_sequences(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    name, sequence = None, []
    with opener(path, 'rt') as handle:
        for line in handle:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(sequence).upper()
                name, sequence = line[1:].split()[0], []
            else:
                sequence.append(line.strip())
        if name is not None:
            yield name, ''.join(sequence).upper()


def write_windows(reference, bed, metadata, settings):
    window = int(settings['window_bp'])
    count = 0
    Path(bed).parent.mkdir(parents=True, exist_ok=True)
    with open(bed, 'w') as regions, open(metadata, 'w') as table:
        writer = csv.writer(table, delimiter='\t')
        writer.writerow(['chromosome', 'start', 'end', 'gc', 'n_fraction', 'eligible'])
        for name, sequence in fasta_sequences(reference):
            if len(sequence) < settings['min_contig_bp'] or name in settings['exclude_contigs']:
                continue
            for start in range(0, len(sequence), window):
                chunk = sequence[start:start+window]
                end = start + len(chunk)
                canonical = sum(chunk.count(base) for base in 'ACGT')
                n_fraction = 1 - canonical / len(chunk)
                gc = (chunk.count('G') + chunk.count('C')) / canonical if canonical else 0
                regions.write(f'{name}\t{start}\t{end}\n')
                writer.writerow([name, start, end, gc, n_fraction, int(n_fraction <= settings['max_n_fraction'])])
                count += 1
    if not count:
        raise ValueError('No reference windows remain; check min_contig_bp and exclude_contigs')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', required=True)
    parser.add_argument('--bed', required=True)
    parser.add_argument('--metadata', required=True)
    parser.add_argument('--settings', required=True)
    args = parser.parse_args()
    write_windows(args.reference, args.bed, args.metadata, json.loads(args.settings))


if __name__ == '__main__':
    main()
