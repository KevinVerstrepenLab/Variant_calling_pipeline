#!/usr/bin/env python3
"""Screen fixed-window coverage for local relative gains/losses, without ploidy calls."""
import argparse
import csv
import gzip
import json
from pathlib import Path

import numpy as np


def weighted_median(values, weights):
    values, weights = np.asarray(values), np.asarray(weights)
    order = np.argsort(values)
    return float(values[order][np.searchsorted(np.cumsum(weights[order]), weights.sum()/2)])


def gc_correct(depth, gc, weights, training, settings):
    if not settings['gc_correct']:
        return depth.copy(), 'disabled'
    train = np.flatnonzero(training & (depth > 0))
    if len(train) < settings['min_gc_windows'] or np.ptp(gc[train]) < 0.02:
        return depth.copy(), 'skipped: insufficient windows or GC variation'
    groups = np.array_split(train[np.argsort(gc[train])], max(1, min(settings['gc_bins'], len(train)//10)))
    x = np.array([np.median(gc[group]) for group in groups])
    y = np.array([weighted_median(depth[group], weights[group]) for group in groups])
    # Collapse tied GC positions before interpolation.
    unique = np.unique(x)
    y = np.array([np.median(y[x == value]) for value in unique])
    trend = np.interp(gc, unique, y)
    neutral = weighted_median(trend[training], weights[training])
    return depth / (trend / neutral), 'applied: GC-bin median trend'


def smooth_contiguous(values, chromosomes, starts, ends, eligible, width):
    result = np.full(len(values), np.nan)
    first = 0
    while first < len(values):
        if not eligible[first]:
            first += 1
            continue
        last = first + 1
        while (last < len(values) and eligible[last] and chromosomes[last] == chromosomes[first]
               and starts[last] == ends[last-1]):
            last += 1
        half = width // 2
        for i in range(first, last):
            result[i] = np.median(values[max(first, i-half):min(last, i+half+1)])
        first = last
    return result


def candidate_regions(chromosomes, starts, ends, smoothed, settings):
    states = np.zeros(len(smoothed), dtype=int)
    states[smoothed >= settings['log2_threshold']] = 1
    states[smoothed <= -settings['log2_threshold']] = -1
    regions = []
    first = 0
    while first < len(states):
        last = first + 1
        while (last < len(states) and states[last] == states[first]
               and chromosomes[last] == chromosomes[first] and starts[last] == ends[last-1]):
            last += 1
        if states[first] and last-first >= settings['min_region_windows']:
            regions.append([chromosomes[first], int(starts[first]), int(ends[last-1]),
                            'gain' if states[first] == 1 else 'loss', last-first,
                            float(np.median(smoothed[first:last]))])
        first = last
    return regions


def analyze(regions_path, metadata_path, settings):
    with open(metadata_path) as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    chrom = np.array([row['chromosome'] for row in rows])
    starts = np.array([int(row['start']) for row in rows])
    ends = np.array([int(row['end']) for row in rows])
    gc = np.array([float(row['gc']) for row in rows])
    eligible = np.array([row['eligible'] == '1' for row in rows])
    depth = []
    with gzip.open(regions_path, 'rt') as handle:
        for i, line in enumerate(handle):
            parts = line.split()
            if i >= len(rows) or (parts[0], int(parts[1]), int(parts[2])) != (chrom[i], starts[i], ends[i]):
                raise ValueError('Coverage windows do not match reference metadata')
            depth.append(float(parts[3]))
    depth = np.array(depth)
    if len(depth) != len(rows) or not np.isfinite(depth).all() or (depth < 0).any():
        raise ValueError('Missing, invalid or negative window depths')
    baseline = eligible.copy()
    if settings['baseline_contigs']:
        missing = set(settings['baseline_contigs']) - set(chrom)
        if missing:
            raise ValueError(f'Baseline contigs absent from selected windows: {sorted(missing)}')
        baseline &= np.isin(chrom, settings['baseline_contigs'])
    if not baseline.any():
        raise ValueError('No eligible baseline windows remain')
    lengths = ends-starts
    raw_baseline = weighted_median(depth[baseline], lengths[baseline])
    if raw_baseline <= 0:
        raise ValueError('Median baseline depth is zero; no stable coverage baseline')
    corrected, correction = gc_correct(depth, gc, lengths, baseline, settings)
    corrected_baseline = weighted_median(corrected[baseline], lengths[baseline])
    pseudo = settings['pseudocount']
    raw_log2 = np.log2((depth+pseudo)/(raw_baseline+pseudo))
    log2 = np.log2((corrected+pseudo)/(corrected_baseline+pseudo))
    log2[~eligible] = np.nan
    smoothed = smooth_contiguous(log2, chrom, starts, ends, eligible, settings['smooth_windows'])
    regions = candidate_regions(chrom, starts, ends, smoothed, settings)
    qc = {'baseline_depth': raw_baseline, 'corrected_baseline_depth': corrected_baseline,
          'low_depth': bool(raw_baseline < settings['min_baseline_depth']),
          'gc_correction': correction, 'windows': len(rows), 'eligible_windows': int(eligible.sum()),
          'candidate_regions': len(regions), 'settings': settings,
          'interpretation': 'Relative coverage screen; no absolute ploidy or statistical CNV calls'}
    return rows, depth, corrected, raw_log2, log2, smoothed, regions, qc


def write_results(prefix, sample, result, settings):
    rows, depth, corrected, raw_log2, log2, smoothed, regions, qc = result
    prefix = str(prefix)
    Path(prefix).parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(prefix+'.windows.tsv.gz', 'wt') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(list(rows[0]) + ['mean_depth', 'corrected_depth', 'raw_log2_ratio', 'log2_ratio', 'smoothed_log2_ratio'])
        for i, row in enumerate(rows):
            writer.writerow(list(row.values()) + [depth[i], corrected[i], raw_log2[i], log2[i], smoothed[i]])
    with open(prefix+'.candidate_regions.tsv', 'w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['chromosome', 'start_0based', 'end_exclusive', 'direction', 'windows', 'median_smoothed_log2_ratio'])
        writer.writerows(regions)
    plot_qc = draw_plot(prefix, sample, rows, depth, raw_log2, log2, smoothed, qc, settings)
    Path(prefix+'.qc.json').write_text(json.dumps(dict(qc, sample=sample, plot=plot_qc), indent=2)+'\n')


def display_layout(rows, settings):
    """Filter and wrap reference-order contigs without changing analysis windows."""
    lengths = {}
    for row in rows:
        chromosome = row['chromosome']
        lengths[chromosome] = max(lengths.get(chromosome, 0), int(row['end']))
    minimum = settings.get('plot_min_contig_bp')
    minimum = 50000 if minimum is None else minimum
    labels = [name for name, length in lengths.items() if length >= minimum]
    if not labels:
        raise ValueError(f'No sequences meet plot_min_contig_bp={minimum}; lower this display threshold')
    per_row = settings.get('chromosomes_per_row') or 16
    bands = [labels[start:start+per_row] for start in range(0, len(labels), per_row)]
    return bands, lengths


def display_limits(values, quantiles, anchors=(), minimum_span=1.0, nonnegative=False):
    """Trim only the display range; inward rank rounding excludes isolated extremes."""
    finite = np.sort(np.asarray(values)[np.isfinite(values)])
    if not len(finite):
        finite = np.array([0.])
    if quantiles is None:
        low, high = float(finite[0]), float(finite[-1])
    else:
        lower = int(np.ceil((len(finite)-1)*quantiles[0]))
        upper = int(np.floor((len(finite)-1)*quantiles[1]))
        if lower > upper:
            lower, upper = 0, len(finite)-1
        low, high = float(finite[lower]), float(finite[upper])
    for anchor in anchors:
        low, high = min(low, anchor), max(high, anchor)
    span = max(high-low, minimum_span)
    middle = (low+high)/2
    low, high = middle-0.55*span, middle+0.55*span
    if nonnegative:
        low = max(0, low)
    return float(low), float(high)


def scatter_with_overflow(ax, x, y, limits, color, size, label=None):
    """Show out-of-range observations as edge triangles rather than squeezing axes."""
    low, high = limits
    valid = np.isfinite(y)
    inside = valid & (y >= low) & (y <= high)
    colors = np.asarray(color) if not isinstance(color, str) else color
    ax.scatter(x[inside], y[inside], s=size,
               c=colors[inside] if not isinstance(colors, str) else colors,
               label=label, rasterized=True)
    counts = {}
    for name, mask, marker, level in (
            ('below', valid & (y < low), 'v', low+0.015*(high-low)),
            ('above', valid & (y > high), '^', high-0.015*(high-low))):
        counts[name] = int(mask.sum())
        if counts[name]:
            ax.scatter(x[mask], np.full(counts[name], level), marker=marker, s=16,
                       c=colors[mask] if not isinstance(colors, str) else colors, rasterized=True)
    return counts


def draw_plot(prefix, sample, rows, depth, raw_log2, log2, smoothed, qc, settings):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    bands, lengths = display_layout(rows, settings)
    chromosomes = np.array([row['chromosome'] for row in rows])
    shown = np.isin(chromosomes, [name for band in bands for name in band])
    quantiles = settings.get('outlier_quantiles', [0.01, 0.99])
    depth_limits = display_limits(depth[shown], quantiles, anchors=[qc['baseline_depth']],
                                  minimum_span=max(1, qc['baseline_depth']*0.2), nonnegative=True)
    log_limits = display_limits(np.concatenate([raw_log2[shown], log2[shown], smoothed[shown]]),
                                quantiles, anchors=[-settings['log2_threshold'], 0, settings['log2_threshold']],
                                minimum_span=1.2)
    fig = plt.figure(figsize=(settings['width'], settings['height']*len(bands)), constrained_layout=True)
    grid = fig.add_gridspec(len(bands), 1, hspace=0.08)
    plot_qc = {'displayed_contigs': [name for band in bands for name in band],
               'hidden_contigs': list(dict.fromkeys(chromosomes[~shown].tolist())),
               'sequence_rows': len(bands),
               'plot_min_contig_bp': 50000 if settings.get('plot_min_contig_bp') is None else settings['plot_min_contig_bp'],
               'depth_ylim': depth_limits, 'log2_ylim': log_limits, 'rows': []}
    for band_index, labels in enumerate(bands):
        tracks = grid[band_index].subgridspec(2, 1, hspace=0.04)
        top = fig.add_subplot(tracks[0])
        bottom = fig.add_subplot(tracks[1], sharex=top)
        offsets, centers, total = {}, [], 0
        for chromosome in labels:
            offsets[chromosome] = total
            centers.append((total+lengths[chromosome]/2)/1e6)
            total += lengths[chromosome]
        indices = np.flatnonzero(np.isin(chromosomes, labels))
        x = np.array([(offsets[rows[i]['chromosome']] + (int(rows[i]['start'])+int(rows[i]['end']))/2)/1e6
                      for i in indices])
        depth_overflow = scatter_with_overflow(top, x, depth[indices], depth_limits, '0.45', 3)
        raw_overflow = scatter_with_overflow(bottom, x, raw_log2[indices], log_limits, '0.8', 3, 'Raw')
        values = log2[indices]
        colors = np.where(values >= settings['log2_threshold'], '#d73027',
                          np.where(values <= -settings['log2_threshold'], '#4575b4', '0.5'))
        norm_overflow = scatter_with_overflow(bottom, x, values, log_limits, colors, 5, 'Normalized')
        for chromosome in labels:
            mask = chromosomes[indices] == chromosome
            trace = smoothed[indices][mask].copy()
            trace[(trace < log_limits[0]) | (trace > log_limits[1])] = np.nan
            bottom.plot(x[mask], trace, color='black', linewidth=0.9)
        top.axhline(qc['baseline_depth'], color='black', linestyle='--', linewidth=0.8)
        top.set_ylabel('Mean depth (reads/base)')
        genomic_axis = top.secondary_xaxis('top')
        genomic_axis.set_xlabel(f'Sequence row {band_index+1}: cumulative position (Mb)')
        top.tick_params(axis='x', labelbottom=False, bottom=False)
        for value in (-settings['log2_threshold'], 0, settings['log2_threshold']):
            bottom.axhline(value, color='0.4', linestyle='--', linewidth=0.6)
        for ax, limits in ((top, depth_limits), (bottom, log_limits)):
            for offset in list(offsets.values())[1:]:
                ax.axvline(offset/1e6, color='0.85', linewidth=0.6)
            ax.grid(axis='y', alpha=0.15)
            ax.set_xlim(0, total/1e6)
            ax.set_ylim(*limits)
        for ax, count in ((top, sum(depth_overflow.values())),
                          (bottom, sum(raw_overflow.values())+sum(norm_overflow.values()))):
            if count:
                ax.text(0.01, 0.97, f'{count} points outside range (triangles)', transform=ax.transAxes,
                        va='top', fontsize=8, bbox={'facecolor': 'white', 'alpha': 0.8, 'edgecolor': 'none'})
        bottom.set_ylabel('log2 coverage / sample baseline')
        bottom.set_xticks(centers, labels, rotation=60, ha='right')
        bottom.set_xlabel('Reference contigs in genome order')
        bottom.legend(loc='upper right', fontsize=8)
        plot_qc['rows'].append({'contigs': labels, 'depth_outside_range': depth_overflow,
                                'raw_log2_outside_range': raw_overflow, 'normalized_outside_range': norm_overflow})
    status = 'LOW DEPTH' if qc['low_depth'] else 'relative coverage screen'
    fig.suptitle(f"{sample} | {status} | {settings['window_bp']:,} bp windows | "
                 f"baseline {qc['baseline_depth']:.1f}x | GC: {qc['gc_correction']}", fontsize=11)
    for suffix in ('.png', '.pdf'):
        fig.savefig(str(prefix)+suffix, dpi=settings['dpi'])
    plt.close(fig)
    return plot_qc


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('regions', 'metadata', 'settings', 'prefix', 'sample'):
        parser.add_argument('--'+key, required=True)
    args = parser.parse_args()
    settings = json.loads(args.settings)
    result = analyze(args.regions, args.metadata, settings)
    write_results(args.prefix, args.sample, result, settings)


if __name__ == '__main__':
    main()
