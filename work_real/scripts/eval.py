#!/usr/bin/env python3
"""Evaluate read-to-locus mappings against a fixed mapping-truth population.

The mapping-truth table contains: read_name, hap, contig, start, end. A mapping
agrees when it is on the truth contig and its interval IoU meets --ov.
Opposite MAT/PAT contigs remain errors; reciprocal haplotype PAFs distinguish
homologous-locus errors from other-locus errors.

All mapping-truth reads are scored. Select one or more of the three precomputed
stratifications with --strata: alignment, hapmer, or variant.

Usage:
  eval.py --mapping-truth chr6.reads.confident.mapping_truth.tsv \
      --strata hapmer --strata alignment name=mapping.tsv ...
"""
import argparse
import collections
import os
import re


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
WORK_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STRATIFICATION_FILES = {
    "alignment": os.path.join(WORK_DIR, "eval", "alignment_strata.tsv"),
    "hapmer": os.path.join(WORK_DIR, "eval", "hapmer_strata.tsv"),
    "variant": os.path.join(WORK_DIR, "eval", "variant_strata.tsv"),
}
STRATIFICATION_ORDER = {
    "alignment": ("no_map", "dup", "hard", "med", "easy"),
    "hapmer": (
        "no_hapmer", "ambiguous", "disagree", "vhard", "hard", "med", "easy"
    ),
    "variant": ("no_variant", "variant"),
}


def load_mapping_truth(path):
    mapping_truth = {}
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            mapping_truth[p[0]] = (p[1], p[2], int(p[3]), int(p[4]))
    return mapping_truth


def load_map(path, ncol):
    data = {}
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= ncol:
                data[p[0]] = p
    return data


def load_mapping(path):
    """Load one representative mapping per read.

    PAF produced with --secondary=no can still contain multiple primary-tagged
    records for a split alignment.  The mapping TSV lacks query spans and
    alignment scores, so select the record with the longest reference span and
    retain the first record on an exact tie.
    """
    mappings = {}
    best_spans = {}
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 5:
                start, end = int(p[3]), int(p[4])
                span = end - start
                if p[0] not in best_spans or span > best_spans[p[0]]:
                    mappings[p[0]] = (p[1], p[2], start, end)
                    best_spans[p[0]] = span
    return mappings


def interval_iou(s1, e1, s2, e2):
    inter = max(0, min(e1, e2) - max(s1, s2))
    union = max(e1, e2) - min(s1, s2)
    return inter / union if union > 0 else 0.0


def load_liftover(paths):
    """Load primary PAF alignments as query-to-target coordinate blocks."""
    lifts = collections.defaultdict(list)
    for path in paths:
        with open(path) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) < 12:
                    continue
                tp = next((t[5:] for t in p[12:] if t.startswith("tp:A:")), None)
                if tp is not None and tp != "P":
                    continue
                cigar = next((t[5:] for t in p[12:] if t.startswith("cg:Z:")), None)
                if cigar is None:
                    continue
                qname, qstart, qend = p[0], int(p[2]), int(p[3])
                strand = p[4]
                tname, tstart, tend = p[5], int(p[7]), int(p[8])
                qpos = qstart
                tpos = tstart if strand == "+" else tend
                blocks = []
                for length_s, op in CIGAR_RE.findall(cigar):
                    length = int(length_s)
                    if op in "M=X":
                        if strand == "+":
                            blocks.append((qpos, qpos + length, tpos,
                                           tpos + length, strand))
                            tpos += length
                        else:
                            blocks.append((qpos, qpos + length, tpos - length,
                                           tpos, strand))
                            tpos -= length
                        qpos += length
                    elif op in "IS":
                        qpos += length
                    elif op in "DN":
                        tpos += length if strand == "+" else -length
                lifts[(qname, tname)].append((qstart, qend, blocks))
    return lifts


def project_intervals(records, start, end, min_coverage):
    """Project [start,end) through each sufficiently covering PAF record."""
    projections = []
    span = end - start
    if span <= 0:
        return projections
    for qstart, qend, blocks in records:
        if end <= qstart or start >= qend:
            continue
        pieces = []
        mapped = 0
        for bqs, bqe, bts, bte, strand in blocks:
            if bqe <= start:
                continue
            if bqs >= end:
                break
            os = max(start, bqs)
            oe = min(end, bqe)
            if os >= oe:
                continue
            mapped += oe - os
            if strand == "+":
                ts = bts + os - bqs
                te = bts + oe - bqs
            else:
                ts = bte - (oe - bqs)
                te = bte - (os - bqs)
            pieces.append((ts, te))
        if mapped / span >= min_coverage:
            projections.append((min(s for s, _e in pieces),
                                max(e for _s, e in pieces)))
    return projections


def classify(truth, mapping, lifts, ov):
    """Classify agreement and the type of disagreement."""
    if mapping is None:
        return "unmapped"
    _m_hap, m_ctg, m_s, m_e = mapping
    _t_hap, t_ctg, t_s, t_e = truth
    if m_ctg == t_ctg:
        if interval_iou(t_s, t_e, m_s, m_e) >= ov:
            return "agreement"
        return "other_locus"
    projections = project_intervals(lifts.get((m_ctg, t_ctg), []),
                                    m_s, m_e, ov)
    if not projections:
        return "other_locus"
    if any(interval_iou(t_s, t_e, s, e) >= ov for s, e in projections):
        return "homologous_locus"
    return "other_locus"


def parse_elapsed(value):
    """Convert [[hours:]minutes:]seconds from /usr/bin/time to seconds."""
    try:
        fields = [float(field) for field in value.split(":")]
    except ValueError:
        return None
    if not fields or len(fields) > 3:
        return None
    return sum(field * (60 ** power)
               for power, field in enumerate(reversed(fields)))


def parse_time(path):
    """Return (CPU seconds, wall seconds, maximum RSS GB) for a successful run."""
    if not path or not os.path.exists(path):
        return None
    user = system = wall = rss = status = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("User time"):
                user = float(line.split(":")[-1])
            elif line.startswith("System time"):
                system = float(line.split(":")[-1])
            elif line.startswith("Elapsed (wall clock) time"):
                wall = parse_elapsed(line.rsplit("):", 1)[-1].strip())
            elif line.startswith("Maximum resident set size"):
                rss = float(line.split(":")[-1]) / 1e6
            elif line.startswith("Exit status"):
                status = int(line.split(":")[-1])
    if (user is None or system is None or wall is None or rss is None
            or status != 0):
        return None
    cpu = user + system if user is not None and system is not None else None
    return cpu, wall, rss


def aggregate_times(paths):
    """Aggregate sequential timed steps: sum CPU/wall and take peak RSS."""
    timings = [parse_time(path) for path in paths]
    if any(timing is None for timing in timings):
        return None
    return (sum(timing[0] for timing in timings),
            sum(timing[1] for timing in timings),
            max(timing[2] for timing in timings))


def preprocessing_time_groups(time_dir):
    """Return preprocessing groups matching the current stage 05/06 outputs."""
    himapper_dir = os.path.join(time_dir, "himapper")
    def himapper_paths(*names):
        return [os.path.join(himapper_dir, name) for name in names]

    common = ("ref-profile.time", "reads-profile.time")
    return {
        "himapper/ref-all_read-seed": himapper_paths(
            *common, "reads-vs-ref-profile.time", "read-seeds.time"),
        "himapper/ref-seed_read-all": himapper_paths(
            *common, "ref-vs-reads-profile.time", "ref-seeds.time"),
        "himapper/ref-seed_read-seed": himapper_paths(
            *common,
            "reads-vs-ref-profile.time", "ref-vs-reads-profile.time",
            "ref-seeds.time", "read-seeds.time"),
        "mm2": [os.path.join(time_dir, "mm2", "index.time")],
        "wm2": [
            os.path.join(time_dir, "wm2", "ref-meryl.time"),
            os.path.join(time_dir, "wm2", "repetitive-kmers.time"),
        ],
    }


def mapping_family(mapping_path, time_dir):
    """Infer the preprocessing family from a mapping path under time_dir."""
    mapping_dir = os.path.abspath(os.path.dirname(mapping_path))
    time_dir = os.path.abspath(time_dir)
    try:
        relative = os.path.relpath(mapping_dir, time_dir)
    except ValueError:
        return None
    if relative == os.pardir or relative.startswith(os.pardir + os.sep):
        return None
    groups = preprocessing_time_groups(time_dir)
    parts = relative.split(os.sep)
    family = parts[0]
    if family == "himapper" and len(parts) > 1:
        family = os.path.join(family, parts[1])
    return family if family in groups else None


def score(categories, population, strata, stratum):
    counts = collections.Counter()
    for name in population:
        if strata is not None and strata.get(name) != stratum:
            continue
        counts[categories[name]] += 1
    return sum(counts.values()), counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mapping-truth", required=True)
    ap.add_argument(
        "--strata", choices=tuple(STRATIFICATION_FILES), required=True,
        action="append",
        help="stratification axis; repeat to include multiple axes")
    ap.add_argument(
        "--strata-only", action="store_true",
        help="report named strata without the aggregate all stratum")
    ap.add_argument("--liftover-paf", action="append", default=[],
                    help="query-to-target haplotype PAF for disagreement classification")
    ap.add_argument(
        "--ov", type=float, default=0.5,
        help="minimum interval IoU and liftover coverage [0.5]")
    ap.add_argument(
        "--time-dir", default=WORK_DIR,
        help=f"work directory containing time files [{WORK_DIR}]")
    ap.add_argument("tools", nargs="+", help="name=mapping.tsv ...")
    args = ap.parse_args()

    mapping_truth = load_mapping_truth(args.mapping_truth)
    population = list(mapping_truth)
    axes = list(dict.fromkeys(args.strata))
    stratifications = []
    for axis in axes:
        strata_path = STRATIFICATION_FILES[axis]
        strata = load_map(strata_path, 2)
        strata = {name: fields[1] for name, fields in strata.items()}
        missing_strata = [name for name in population if name not in strata]
        if missing_strata:
            ap.error(
                f"{strata_path} lacks {len(missing_strata)} mapping-truth reads")
        stratum_values = {strata[name] for name in population}
        expected_order = STRATIFICATION_ORDER[axis]
        unexpected = sorted(stratum_values - set(expected_order))
        ordered_strata = list(expected_order) + unexpected
        stratifications.append((axis, strata_path, strata, ordered_strata))
    lifts = load_liftover(args.liftover_paf)

    tools = []
    for spec in args.tools:
        name, path = spec.split("=", 1)
        mappings = load_mapping(path)
        categories = {
            read: classify(
                mapping_truth[read], mappings.get(read), lifts, args.ov)
            for read in population
        }
        tools.append((name, path, categories))

    print(f"# mapping truth: {args.mapping_truth}\n")

    def rate(a, b):
        return f"{100 * a / b:.3f}" if b else "N/A"

    def report_matrix(axis, strata, ordered_strata):
        columns = []
        for stratum in ordered_strata:
            row_ns = []
            agreements = []
            metrics = []
            for _name, _path, categories in tools:
                n, counts = score(categories, population, strata, stratum)
                row_ns.append(n)
                agreements.append(counts["agreement"])
                metrics.append((
                    rate(counts["agreement"], n),
                    f"{counts['homologous_locus']:,}",
                    f"{counts['other_locus']:,}",
                    f"{counts['unmapped']:,}",
                ))
            if len(set(row_ns)) != 1:
                raise RuntimeError(
                    f"inconsistent population size among tools for {axis}/{stratum}")
            eligible = [
                i for i, (name, _path, _categories) in enumerate(tools)
                if name != "wm2"
            ]
            if row_ns[0] and eligible:
                best = max(agreements[i] for i in eligible)
                metrics = [
                    (
                        agreement + (
                            "*" if i in eligible and agreements[i] == best
                            else ""),
                        homologous,
                        other,
                        unmapped,
                    )
                    for i, (agreement, homologous, other, unmapped)
                    in enumerate(metrics)
                ]
            columns.append((stratum, f"{row_ns[0]:,}", metrics))

        def unmark(value):
            return value[:-1] if value.endswith("*") else value

        metric_widths = [
            max(
                len(unmark(metrics[i][j]))
                for _stratum, _n, metrics in columns
                for i in range(len(tools)))
            for j in range(4)
        ]
        labels = [
            f"{stratum} ({n})"
            for stratum, n, _metrics in columns
        ]
        formatted_columns = []
        for _stratum, _n, metrics in columns:
            formatted_columns.append([
                f"{unmark(agreement):>{metric_widths[0]}}"
                f"{'*' if agreement.endswith('*') else ' '}"
                f"[{homologous:>{metric_widths[1]}}/"
                f"{other:>{metric_widths[2]}}/"
                f"{unmapped:>{metric_widths[3]}}]"
                for agreement, homologous, other, unmapped in metrics
            ])

        tool_width = max(
            len("tool"), *(len(name) for name, _path, _categories in tools))
        column_widths = [
            max(
                len(label),
                *(len(cell) for cell in formatted_columns[i]))
            for i, label in enumerate(labels)
        ]

        print(f"=== accuracy: {axis} ===")
        print(
            f"{'tool':<{tool_width}} "
            + " ".join(
                f"{label:>{width}}"
                for label, width in zip(labels, column_widths)))
        for tool_index, (name, _path, _categories) in enumerate(tools):
            cells = [
                formatted_columns[column_index][tool_index]
                for column_index in range(len(columns))
            ]
            print(
                f"{name:<{tool_width}} "
                + " ".join(
                    f"{cell:>{width}}"
                    for cell, width in zip(cells, column_widths)))
        print()

    print("# %accuracy [homologous_locus/other_locus/unmapped]")
    print("# * = highest %accuracy in each stratum (wm2 excluded)\n")
    if not args.strata_only:
        report_matrix("overall", None, ("all",))
    for axis, _strata_path, strata, ordered_strata in stratifications:
        report_matrix(axis, strata, ordered_strata)

    if args.time_dir:
        prep_groups = preprocessing_time_groups(args.time_dir)

        def seconds(number):
            return f"{int(number + 0.5):,}" if number is not None else "N/A"

        def rss(number):
            return f"{number:,.2f}" if number is not None else "N/A"

        timing_rows = []
        for name, path, _categories in tools:
            map_timing = parse_time(
                os.path.join(os.path.dirname(path), "map.time"))
            family = mapping_family(path, args.time_dir)
            prep_timing = (
                aggregate_times(prep_groups[family]) if family else None)
            if prep_timing is not None and map_timing is not None:
                total_timing = (
                    prep_timing[0] + map_timing[0],
                    prep_timing[1] + map_timing[1],
                    max(prep_timing[2], map_timing[2]),
                )
            else:
                total_timing = None
            prep_cpu, prep_wall, prep_rss = (
                prep_timing if prep_timing is not None
                else (None, None, None))
            map_cpu, map_wall, map_rss = (
                map_timing if map_timing is not None
                else (None, None, None))
            total_cpu, total_wall, peak_rss = (
                total_timing if total_timing is not None
                else (None, None, None))

            timing_rows.append((
                name,
                (seconds(total_cpu), seconds(prep_cpu), seconds(map_cpu)),
                (seconds(total_wall), seconds(prep_wall), seconds(map_wall)),
                (rss(peak_rss), rss(prep_rss), rss(map_rss)),
                (total_cpu, total_wall, peak_rss),
            ))

        minimum_totals = []
        for metric in range(3):
            values = [
                row[4][metric] for row in timing_rows
                if row[4][metric] is not None
            ]
            minimum_totals.append(min(values) if values else None)

        value_widths = [
            [
                max(len(row[metric][part]) for row in timing_rows)
                for part in range(3)
            ]
            for metric in range(1, 4)
        ]
        formatted_timing_rows = []
        for row in timing_rows:
            cells = [
                f"{row[metric][0]:>{value_widths[metric - 1][0]}}"
                f"{'*' if (minimum_totals[metric - 1] is not None and row[4][metric - 1] == minimum_totals[metric - 1]) else ' '}"
                f"[{row[metric][1]:>{value_widths[metric - 1][1]}}/"
                f"{row[metric][2]:>{value_widths[metric - 1][2]}}]"
                for metric in range(1, 4)
            ]
            formatted_timing_rows.append((row[0], *cells))

        headers = ("tool", "CPU time (s)", "wall time (s)", "peak RSS (GB)")
        widths = [
            max(
                len(headers[i]),
                *(len(row[i]) for row in formatted_timing_rows))
            for i in range(len(headers))
        ]
        print("=== resources ===")
        print("# time / memory [prep/map]")
        print("# * = lowest total in each resource column")
        print("  ".join(
            f"{header:<{width}}" if i == 0 else f"{header:>{width}}"
            for i, (header, width) in enumerate(zip(headers, widths))))
        for row in formatted_timing_rows:
            print("  ".join(
                f"{cell:<{width}}" if i == 0 else f"{cell:>{width}}"
                for i, (cell, width) in enumerate(zip(row, widths))))


if __name__ == "__main__":
    main()
