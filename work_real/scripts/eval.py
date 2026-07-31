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
    mapping_truth, _excluded, _total = load_mapping_truth_filtered(path)
    return mapping_truth


def load_mapping_truth_filtered(path, exclusion_index=None):
    """Load truth records, optionally dropping any interval touching a mask."""
    mapping_truth = {}
    excluded = 0
    total = 0
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            record = (p[1], p[2], int(p[3]), int(p[4]))
            total += 1
            if (
                exclusion_index is not None
                and truth_overlaps_exclusion(record, exclusion_index)
            ):
                excluded += 1
                continue
            mapping_truth[p[0]] = record
    return mapping_truth, excluded, total


def merge_bed_records(records):
    """Return sorted, same-contig overlap/book-end merged BED3 records."""
    merged = []
    for contig, start, end in sorted(records):
        if start < 0 or end <= start:
            raise ValueError(f"invalid BED interval: {contig}:{start}-{end}")
        if (
            merged
            and merged[-1][0] == contig
            and start <= merged[-1][2]
        ):
            old_contig, old_start, old_end = merged[-1]
            merged[-1] = (old_contig, old_start, max(old_end, end))
        else:
            merged.append((contig, start, end))
    return tuple(merged)


def load_exclusion_beds(paths):
    """Load one or more BEDs as a single 0-based, half-open interval union."""
    records = []
    for path in paths:
        with open(path) as fh:
            for line_number, line in enumerate(fh, 1):
                if (
                    not line.strip()
                    or line.startswith(("#", "track", "browser"))
                ):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    raise ValueError(
                        f"{path}:{line_number}: fewer than 3 BED columns")
                try:
                    start, end = int(fields[1]), int(fields[2])
                except ValueError as error:
                    raise ValueError(
                        f"{path}:{line_number}: non-integer BED coordinate"
                    ) from error
                records.append((fields[0], start, end))
    if paths and not records:
        raise RuntimeError("truth-exclusion BEDs contain no intervals")
    merged = merge_bed_records(records)
    index = collections.defaultdict(list)
    for contig, start, end in merged:
        index[contig].append((start, end))
    return dict(index), merged


def truth_overlaps_exclusion(truth_record, exclusion_index):
    """Whether a truth interval has at least one base in the exclusion union."""
    _hap, contig, start, end = truth_record
    return any(
        start < mask_end and end > mask_start
        for mask_start, mask_end in exclusion_index.get(contig, ())
    )


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
    """Return named preprocessing stages for current stage 05/06 outputs."""
    himapper_dir = os.path.join(time_dir, "himapper")
    def himapper_paths(*names):
        return [os.path.join(himapper_dir, name) for name in names]

    dual_profile = os.path.join(himapper_dir, "reads-dual-profile.time")
    if os.path.exists(dual_profile):
        # FastKdual produces reads.prof, reads.ktab, and reads.ref.prof from
        # one partition/sort.  Keep the legacy layout as a fallback for runs
        # prepared before FastKdual was introduced.
        ref_and_reads = ("ref-profile.time", "reads-dual-profile.time")
        read_seed_profiles = ref_and_reads
        ref_seed_profiles = (*ref_and_reads, "ref-vs-reads-profile.time")
        both_seed_profiles = ref_seed_profiles
    else:
        ref_and_reads = ("ref-profile.time", "reads-profile.time")
        read_seed_profiles = (*ref_and_reads, "reads-vs-ref-profile.time")
        ref_seed_profiles = (*ref_and_reads, "ref-vs-reads-profile.time")
        both_seed_profiles = (
            *ref_and_reads,
            "reads-vs-ref-profile.time",
            "ref-vs-reads-profile.time",
        )

    return {
        "himapper/ref-all_read-seed": {
            "FastK": himapper_paths(*read_seed_profiles),
            "seed": himapper_paths("read-seeds.time"),
        },
        "himapper/ref-seed_read-all": {
            "FastK": himapper_paths(*ref_seed_profiles),
            "seed": himapper_paths("ref-seeds.time"),
        },
        "himapper/ref-seed_read-seed": {
            "FastK": himapper_paths(*both_seed_profiles),
            "seed": himapper_paths("ref-seeds.time", "read-seeds.time"),
        },
        "mm2": {
            "prep": [os.path.join(time_dir, "mm2", "index.time")],
        },
        "wm": {
            "prep": [
                os.path.join(time_dir, "wm", "ref-meryl.time"),
                os.path.join(time_dir, "wm", "repetitive-kmers.time"),
            ],
        },
        "wm2": {
            "prep": [
                os.path.join(time_dir, "wm2", "ref-meryl.time"),
                os.path.join(time_dir, "wm2", "repetitive-kmers.time"),
            ],
        },
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
        "--strata-file", action="append", default=[], metavar="AXIS=TSV",
        help="override a stratification TSV path; repeat by axis")
    ap.add_argument(
        "--exclude-truth-bed", action="append", default=[], metavar="BED",
        help=(
            "exclude reads whose 0-based half-open truth interval overlaps "
            "BED by at least one base; repeated BEDs are unioned"))
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

    exclusion_index, exclusion_records = load_exclusion_beds(
        args.exclude_truth_bed)
    mapping_truth, excluded_truth, source_truth = load_mapping_truth_filtered(
        args.mapping_truth,
        exclusion_index if exclusion_records else None,
    )
    population = list(mapping_truth)
    strata_files = dict(STRATIFICATION_FILES)
    for specification in args.strata_file:
        if "=" not in specification:
            ap.error(
                f"--strata-file must be AXIS=TSV: {specification}")
        axis, path = specification.split("=", 1)
        if axis not in STRATIFICATION_FILES:
            ap.error(
                f"unknown --strata-file axis {axis!r}; choose from "
                + ", ".join(STRATIFICATION_FILES))
        if not path:
            ap.error(f"empty --strata-file path for axis {axis!r}")
        strata_files[axis] = path
    axes = list(dict.fromkeys(args.strata))
    stratifications = []
    for axis in axes:
        strata_path = strata_files[axis]
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

    print(f"# mapping truth: {args.mapping_truth}")
    if exclusion_records:
        for path in args.exclude_truth_bed:
            print(f"# truth exclusion mask: {path}")
        print(
            "# truth exclusion rule: any >=1-bp truth-interval overlap "
            "(0-based, half-open)")
        print(
            f"# truth exclusion union: {len(exclusion_records):,} intervals, "
            f"{sum(end - start for _contig, start, end in exclusion_records):,} bp")
        print(
            f"# truth population: original={source_truth:,} "
            f"excluded={excluded_truth:,} evaluated={len(mapping_truth):,}")
    print()

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
            stage_paths = prep_groups[family] if family else {"prep": []}
            stage_labels = tuple(stage_paths)
            stage_timings = tuple(
                aggregate_times(paths) if paths else None
                for paths in stage_paths.values()
            )
            prep_paths = [
                stage_path
                for paths in stage_paths.values()
                for stage_path in paths
            ]
            prep_timing = aggregate_times(prep_paths) if prep_paths else None
            if prep_timing is not None and map_timing is not None:
                total_timing = (
                    prep_timing[0] + map_timing[0],
                    prep_timing[1] + map_timing[1],
                    max(prep_timing[2], map_timing[2]),
                )
            else:
                total_timing = None
            timing_rows.append((
                name,
                stage_labels + ("map",),
                stage_timings + (map_timing,),
                total_timing,
            ))

        minimum_totals = []
        for metric in range(3):
            values = [
                row[3][metric] for row in timing_rows
                if row[3] is not None and row[3][metric] is not None
            ]
            minimum_totals.append(min(values) if values else None)

        formatted_timing_values = []
        for row in timing_rows:
            metrics = []
            for metric, formatter in enumerate((seconds, seconds, rss)):
                total_value = row[3][metric] if row[3] is not None else None
                total = formatter(total_value)
                parts = tuple(
                    formatter(timing[metric]) if timing is not None else "N/A"
                    for timing in row[2]
                )
                marker = (
                    "*" if (
                        minimum_totals[metric] is not None
                        and total_value == minimum_totals[metric]
                    ) else " "
                )
                metrics.append((total, marker, row[1], parts))
            formatted_timing_values.append((row[0], *metrics))

        value_widths = []
        for metric in range(1, 4):
            total_width = max(
                len(row[metric][0]) for row in formatted_timing_values)
            part_widths = {}
            for row in formatted_timing_values:
                _total, _marker, labels, parts = row[metric]
                widths = part_widths.setdefault(labels, [0] * len(parts))
                for part_index, part in enumerate(parts):
                    widths[part_index] = max(
                        widths[part_index], len(part))
            detail_width = max(
                2 + sum(widths) + max(0, len(widths) - 1)
                for widths in part_widths.values()
            )
            value_widths.append((total_width, part_widths, detail_width))

        formatted_timing_rows = []
        for row in formatted_timing_values:
            cells = []
            for metric in range(1, 4):
                total, marker, labels, parts = row[metric]
                total_width, part_widths, detail_width = (
                    value_widths[metric - 1])
                widths = part_widths[labels]
                detail = "[" + "/".join(
                    f"{part:>{width}}"
                    for part, width in zip(parts, widths)
                ) + "]"
                cells.append(
                    f"{total:>{total_width}}{marker}"
                    f"{detail:<{detail_width}}")
            formatted_timing_rows.append((row[0], *cells))

        headers = ("tool", "CPU time (s)", "wall time (s)", "peak RSS (GB)")
        widths = [
            max(
                len(headers[i]),
                *(len(row[i]) for row in formatted_timing_rows))
            for i in range(len(headers))
        ]
        print("=== resources ===")
        print("# himapper time / memory [FastK/seed/map]")
        print("# Other tools time / memory [prep/map]")
        print("# * = lowest total in each resource column")
        print("  ".join(
            f"{header:<{width}}" if i == 0 else f"{header:>{width}}"
            for i, (header, width) in enumerate(zip(headers, widths))))
        for row in formatted_timing_rows:
            print("  ".join(
                f"{cell:<{width}}" if i == 0 else f"{cell:>{width}}"
                for i, (cell, width) in enumerate(zip(row, widths))).rstrip())


if __name__ == "__main__":
    main()
