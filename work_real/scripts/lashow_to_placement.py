#!/usr/bin/env python3
"""Convert himapper LAshow output to original FASTA scaffold coordinates.

LAshow identifies DAM sequences as 1-based ``scaffold.contig`` pairs.  A FASTA
scaffold containing N runs is represented by multiple sequences in a DAM.
Forward LAshow coordinates already include the sequence's ``fpulse`` offset.
For reverse hits, however, LAshow assumes that the whole scaffold was reverse
complemented, whereas himapper reverse complements each N-delimited DAM contig.
``--ref-layout`` supplies the metadata needed to correct that difference.

By default, reconstruct DALIGNER START/NEXT chains, join compatible nearby
chains (100,000 bp distance, 1% overlap), and emit the path with the largest
union of aligned query intervals for each read:
  columns: read_name  hap  ctg  start  end  strand
This shared default applies to both ONT and HiFi. --no-merge-chains restores
the historical largest-query-bounding-span chain selection.
With --all, emit every ORIGINAL START/NEXT chain in LAshow input order,
without between-chain joining:
  columns: read_name  hap  ctg  start  end  strand  ref_bounding_span

Coordinates are zero-based, half-open FASTA-scaffold coordinates.  The final
--all column is ``end - start`` of the bounded reference interval, including
any internal gaps between alignment pieces.  It is not aligned bases or query
span.  Unmapped reads are omitted, and --all does not mark which chain the
default mode would select.

hap = mat/pat from the contig-name suffix (_MATERNAL/_PATERNAL).

Usage:
  DBshow -n ref.dam > ref.dam.layout
  LAshow reads.dam ref.dam reads.ref.las | \\
    lashow_to_placement.py [--all | --no-merge-chains] \\
      [--merge-distance 100000] [--merge-overlap-percent 1] \\
      --ref-layout ref.dam.layout \\
      <reads.names> <ctg.names> > mapping.tsv
  (reads.names / ctg.names: one name per line, in DAM/fasta order)
"""
import argparse
import re
import sys

DEFAULT_MERGE_DISTANCE = 100000
DEFAULT_MERGE_OVERLAP_PERCENT = 1

NUM = r"[\d,]+"
LINE = re.compile(
    r"^[>+\- ]+"
    r"(" + NUM + r")\.(\d+)\s+"           # A-id . A-ctg   (read index)
    r"(" + NUM + r")\.(\d+)\s+"           # B-id . B-ctg   (contig index)
    r"([nc])\s+"                            # n=fwd, c=rev
    r"[<\[]\s*(" + NUM + r")\s*\.\.\s*"
    r"(" + NUM + r")\s*[>\]]\s+x\s+"        # A coords (ignored)
    r"[<\[]\s*(" + NUM + r")\s*\.\.\s*"
    r"(" + NUM + r")\s*[>\]]"               # B coords (ref start..end)
)
LAYOUT = re.compile(
    r"^>(.*?)\s+::\s+Contig\s+(\d+)"
    r"\[\s*(" + NUM + r")\s*,\s*(" + NUM + r")\s*\]\s*$"
)


def num(s):
    return int(s.replace(",", ""))


def hap_of(ctg):
    if ctg.endswith("_MATERNAL"):
        return "mat"
    if ctg.endswith("_PATERNAL"):
        return "pat"
    return "NA"


def load_names(path):
    with open(path) as fh:
        return [ln.strip().split()[0] for ln in fh if ln.strip()]


def load_ref_layout(path, scaffolds):
    """Load ``DBshow -n`` output keyed by zero-based scaffold/contig IDs."""
    scaffold_ids = {}
    for index, name in enumerate(scaffolds):
        if name in scaffold_ids:
            raise ValueError(f"duplicate scaffold name in name list: {name}")
        scaffold_ids[name] = index

    contigs = {}
    scaffold_lengths = {}
    with open(path) as fh:
        for line_number, line in enumerate(fh, 1):
            match = LAYOUT.match(line.rstrip("\n"))
            if match is None:
                if line.strip():
                    raise ValueError(
                        f"{path}:{line_number}: invalid DBshow layout line"
                    )
                continue
            scaffold_name = match.group(1).split()[0]
            if scaffold_name not in scaffold_ids:
                raise ValueError(
                    f"{path}:{line_number}: unknown scaffold {scaffold_name}"
                )
            scaffold_index = scaffold_ids[scaffold_name]
            contig_index = int(match.group(2))
            fpulse = num(match.group(3))
            end = num(match.group(4))
            if fpulse < 0 or end < fpulse:
                raise ValueError(
                    f"{path}:{line_number}: invalid contig interval "
                    f"[{fpulse},{end}]"
                )
            key = (scaffold_index, contig_index)
            if key in contigs:
                raise ValueError(
                    f"{path}:{line_number}: duplicate scaffold/contig {key}"
                )
            contigs[key] = (fpulse, end - fpulse)
            scaffold_lengths[scaffold_index] = max(
                scaffold_lengths.get(scaffold_index, 0), end
            )

    missing = [
        scaffolds[index]
        for index in range(len(scaffolds))
        if index not in scaffold_lengths
    ]
    if missing:
        raise ValueError(
            f"{path}: no DAM layout records for scaffold(s): "
            + ", ".join(missing)
        )
    for scaffold_index, scaffold in enumerate(scaffolds):
        records = sorted(
            (contig_index, fpulse, contig_length)
            for (index, contig_index), (fpulse, contig_length) in contigs.items()
            if index == scaffold_index
        )
        indices = [record[0] for record in records]
        if indices != list(range(len(records))):
            raise ValueError(
                f"{path}: non-contiguous DAM contig indices for {scaffold}: "
                f"{indices}"
            )
        if records[0][1] != 0:
            raise ValueError(
                f"{path}: unsupported leading-N/all-N scaffold {scaffold}: "
                f"contig 0 starts at {records[0][1]}"
            )
        previous_end = 0
        for contig_index, fpulse, contig_length in records:
            if fpulse < previous_end:
                raise ValueError(
                    f"{path}: overlapping/out-of-order DAM contig "
                    f"{scaffold}.{contig_index}"
                )
            if contig_length == 0 and contig_index != records[-1][0]:
                raise ValueError(
                    f"{path}: zero-length non-final DAM contig "
                    f"{scaffold}.{contig_index}"
                )
            previous_end = fpulse + contig_length
    return contigs, scaffold_lengths


def parse_alignment(
    match,
    line_number,
    reads,
    scaffolds,
    contigs,
    scaffold_lengths,
):
    """Parse one LAshow row and restore scaffold reference coordinates."""
    read_index = num(match.group(1)) - 1
    reference_index = num(match.group(3)) - 1
    contig_index = int(match.group(4)) - 1
    if not (
        0 <= read_index < len(reads)
        and 0 <= reference_index < len(scaffolds)
    ):
        raise ValueError(
            f"LAshow line {line_number}: read/reference index is out of range"
        )

    strand = "+" if match.group(5) == "n" else "-"
    query_start, query_end = sorted(
        (num(match.group(6)), num(match.group(7)))
    )
    reference_start, reference_end = sorted(
        (num(match.group(8)), num(match.group(9)))
    )
    key = (reference_index, contig_index)
    if key not in contigs:
        raise ValueError(
            f"LAshow line {line_number}: missing DAM layout for "
            f"{reference_index + 1}.{contig_index + 1}"
        )
    fpulse, contig_length = contigs[key]
    if strand == "-":
        delta = (
            2 * fpulse
            + contig_length
            - scaffold_lengths[reference_index]
        )
        reference_start += delta
        reference_end += delta
    contig_end = fpulse + contig_length
    if (
        reference_start < fpulse
        or reference_end > contig_end
        or reference_end < reference_start
    ):
        raise ValueError(
            f"LAshow line {line_number}: corrected interval "
            f"[{reference_start},{reference_end}] is outside DAM contig "
            f"{reference_index + 1}.{contig_index + 1} interval "
            f"[{fpulse},{contig_end}]"
        )

    return {
        "read_index": read_index,
        "reference_index": reference_index,
        "contig_index": contig_index,
        "strand": strand,
        "query_start": query_start,
        "query_end": query_end,
        "reference_start": reference_start,
        "reference_end": reference_end,
    }


def bound_chain(chain):
    """Validate and return the bounding placement of one START/NEXT chain."""
    if not chain:
        raise ValueError("cannot bound an empty alignment chain")
    first = chain[0]
    for previous, segment in zip(chain, chain[1:]):
        if (
            segment["read_index"] != first["read_index"]
            or segment["reference_index"] != first["reference_index"]
            or segment["contig_index"] != first["contig_index"]
            or segment["strand"] != first["strand"]
        ):
            raise ValueError(
                "one START/NEXT chain mixes a read, reference, DAM contig, "
                "or strand"
            )
        if segment["query_start"] < previous["query_start"]:
            raise ValueError("query coordinates are not monotone within chain")
        if (
            first["strand"] == "+"
            and segment["reference_start"] < previous["reference_start"]
        ):
            raise ValueError(
                "forward reference coordinates are not monotone within chain"
            )
        if (
            first["strand"] == "-"
            and segment["reference_end"] > previous["reference_end"]
        ):
            raise ValueError(
                "reverse reference coordinates are not monotone within chain"
            )

    bounded = dict(first)
    bounded["query_start"] = min(
        segment["query_start"] for segment in chain
    )
    bounded["query_end"] = max(segment["query_end"] for segment in chain)
    bounded["reference_start"] = min(
        segment["reference_start"] for segment in chain
    )
    bounded["reference_end"] = max(
        segment["reference_end"] for segment in chain
    )
    return bounded


def choose_chain(chains):
    """Historical selector: maximum query span, first exact tie."""
    if not chains:
        raise ValueError("cannot choose from an empty chain collection")
    best = bound_chain(chains[0])
    best_span = best["query_end"] - best["query_start"]
    for chain in chains[1:]:
        bounded = bound_chain(chain)
        span = bounded["query_end"] - bounded["query_start"]
        if span > best_span:
            best = bounded
            best_span = span
    return best


def query_interval_union(intervals):
    """Count aligned query bases once, without rewarding internal gaps."""
    merged = []
    for start, end in sorted(intervals):
        if end <= start:
            raise ValueError("empty/reversed aligned query interval")
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return tuple(merged)


def choose_merged_chain(chains, distance=DEFAULT_MERGE_DISTANCE,
                        overlap_percent=DEFAULT_MERGE_OVERLAP_PERCENT):
    """Select the exact maximum aligned-query-union path of compatible chains.

    Both query and oriented-reference endpoints must advance on the same read,
    reference, DAM contig and strand. Signed gaps gq/gr satisfy
    max(abs(gq), abs(gr), abs(gq-gr)) <= distance. Overlap on each axis is at
    most overlap_percent of the shorter chain's bounding span on that axis.
    Ties prefer fewer chains, then the earliest original chain-index sequence.
    Truth coordinates, gap penalties and search-state cutoffs are not used.

    Keep this implementation self-contained: experiment snapshots copy this
    script as a single file. The historical START/NEXT groups remain intact.
    """
    if distance < 0 or not 0 <= overlap_percent < 100:
        raise ValueError("invalid chain merge distance or overlap percentage")
    if not chains:
        raise ValueError("cannot choose from an empty chain collection")
    bounded = [bound_chain(chain) for chain in chains]
    if len(chains) == 1:
        return bounded[0]
    queries = [query_interval_union(
        (s["query_start"], s["query_end"]) for s in chain
    ) for chain in chains]
    oriented = [
        (b["reference_start"], b["reference_end"]) if b["strand"] == "+"
        else (-b["reference_end"], -b["reference_start"])
        for b in bounded
    ]

    def compatible(i, j):
        a, b = bounded[i], bounded[j]
        if any(a[key] != b[key] for key in (
            "read_index", "reference_index", "contig_index", "strand"
        )):
            return False
        aus, aue = oriented[i]
        bus, bue = oriented[j]
        if not (a["query_start"] < b["query_start"]
                and a["query_end"] < b["query_end"]
                and aus < bus and aue < bue):
            return False
        gq = b["query_start"] - a["query_end"]
        gr = bus - aue
        return (
            max(abs(gq), abs(gr), abs(gq - gr)) <= distance
            and 100 * max(0, -gq) <= overlap_percent * min(
                a["query_end"] - a["query_start"],
                b["query_end"] - b["query_start"])
            and 100 * max(0, -gr) <= overlap_percent * min(aue - aus, bue - bus)
        )

    order = sorted(range(len(chains)), key=lambda i: (
        bounded[i]["query_start"], bounded[i]["query_end"], i))
    states = [{} for _ in chains]
    best_path = best_key = None
    for j in order:
        states[j][queries[j]] = (j,)
        for i in order:
            if not compatible(i, j):
                continue
            for query, path in states[i].items():
                merged = query_interval_union((*query, *queries[j]))
                candidate = (*path, j)
                previous = states[j].get(merged)
                if previous is None or (len(candidate), candidate) < (len(previous), previous):
                    states[j][merged] = candidate
        # Equal unions ending at the same chain have identical future gains;
        # retaining the preferred path is exact, not a beam approximation.
        for query, path in states[j].items():
            key = (-sum(end - start for start, end in query), len(path), path)
            if best_key is None or key < best_key:
                best_key, best_path = key, path
    selected = [bounded[i] for i in best_path]
    result = dict(selected[0])
    for key in ("query_start", "reference_start"):
        result[key] = min(b[key] for b in selected)
    for key in ("query_end", "reference_end"):
        result[key] = max(b[key] for b in selected)
    return result


def emit_placement(placement, reads, scaffolds, *, include_span=False):
    """Write one bounded chain placement."""
    scaffold = scaffolds[placement["reference_index"]]
    fields = [
        reads[placement["read_index"]],
        hap_of(scaffold),
        scaffold,
        str(placement["reference_start"]),
        str(placement["reference_end"]),
        placement["strand"],
    ]
    if include_span:
        fields.append(
            str(placement["reference_end"] - placement["reference_start"])
        )
    print("\t".join(fields))


def emit_read(read_index, chains, reads, scaffolds, *, emit_all,
              merge_distance=DEFAULT_MERGE_DISTANCE,
              merge_overlap_percent=DEFAULT_MERGE_OVERLAP_PERCENT):
    """Emit either every chain or the representative chain for one read."""
    if emit_all:
        placements = [bound_chain(chain) for chain in chains]
    elif merge_distance is None:
        placements = [choose_chain(chains)]
    else:
        placements = [choose_merged_chain(chains, merge_distance, merge_overlap_percent)]
    for placement in placements:
        if placement["read_index"] != read_index:
            raise ValueError("one read group contains a foreign chain segment")
        emit_placement(
            placement,
            reads,
            scaffolds,
            include_span=emit_all,
        )


def convert_chains(
    lines,
    reads,
    scaffolds,
    contigs,
    scaffold_lengths,
    *,
    emit_all=False,
    merge_distance=DEFAULT_MERGE_DISTANCE,
    merge_overlap_percent=DEFAULT_MERGE_OVERLAP_PERCENT,
):
    """Stream sorted LAshow rows and emit bounded candidate chains."""
    current_read = None
    chains = []
    current_chain = None

    for line_number, line in enumerate(lines, 1):
        match = LINE.match(line)
        if match is None:
            continue
        segment = parse_alignment(
            match,
            line_number,
            reads,
            scaffolds,
            contigs,
            scaffold_lengths,
        )
        read_index = segment["read_index"]
        if current_read is not None and read_index < current_read:
            raise ValueError(
                f"LAshow line {line_number}: reads are not in sorted order"
            )
        if current_read is None or read_index != current_read:
            if current_read is not None:
                emit_read(
                    current_read,
                    chains,
                    reads,
                    scaffolds,
                    emit_all=emit_all,
                    merge_distance=merge_distance,
                    merge_overlap_percent=merge_overlap_percent,
                )
            current_read = read_index
            chains = []
            current_chain = None

        flag = line.lstrip()[0]
        if flag in ">+":
            current_chain = [segment]
            chains.append(current_chain)
        elif flag == "-":
            if current_chain is None:
                raise ValueError(
                    f"LAshow line {line_number}: NEXT record without START"
                )
            current_chain.append(segment)
        else:
            raise ValueError(
                f"LAshow line {line_number}: unknown chain flag {flag!r}"
            )

    if current_read is not None:
        emit_read(
            current_read,
            chains,
            reads,
            scaffolds,
            emit_all=emit_all,
            merge_distance=merge_distance,
            merge_overlap_percent=merge_overlap_percent,
        )


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--all",
        action="store_true",
        help=(
            "emit every START/NEXT candidate chain in input order, with a "
            "reference-bounding-span column"
        ),
    )
    ap.add_argument(
        "--no-merge-chains", action="store_true",
        help="restore historical maximum-query-bounding-span single-chain selection",
    )
    ap.add_argument(
        "--merge-distance", type=int, default=DEFAULT_MERGE_DISTANCE,
        help="maximum absolute query/reference gap and gap difference, in bp [100000]",
    )
    ap.add_argument(
        "--merge-overlap-percent", type=int, default=DEFAULT_MERGE_OVERLAP_PERCENT,
        help="maximum overlap as percent of the shorter chain span on each axis [1]",
    )
    ap.add_argument(
        "--ref-layout",
        required=True,
        help="output of 'DBshow -n ref.dam'",
    )
    ap.add_argument("reads_names", help="read names, one per line, in DAM/fasta order")
    ap.add_argument("ctg_names", help="contig names, one per line, in DAM/fasta order")
    args = ap.parse_args()
    if args.merge_distance < 0:
        ap.error("--merge-distance must be nonnegative")
    if not 0 <= args.merge_overlap_percent < 100:
        ap.error("--merge-overlap-percent must be between 0 and 99")
    reads = load_names(args.reads_names)
    ctgs = load_names(args.ctg_names)
    contigs, scaffold_lengths = load_ref_layout(args.ref_layout, ctgs)

    convert_chains(
        sys.stdin,
        reads,
        ctgs,
        contigs,
        scaffold_lengths,
        emit_all=args.all,
        merge_distance=None if args.no_merge_chains else args.merge_distance,
        merge_overlap_percent=args.merge_overlap_percent,
    )


if __name__ == "__main__":
    main()
