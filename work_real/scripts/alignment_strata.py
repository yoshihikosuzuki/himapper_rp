#!/usr/bin/env python3
"""Estimate mapping difficulty by self-mapping mapping-truth reference spans.

For each read, extract its assigned reference span and map it back to the diploid
panel while retaining secondary alignments.

Output:
  <read_name> <level> <n_nm0> <nearest_other_NM>

mapping_truth.tsv columns: read_name, hap, contig, start, end

Usage:
  alignment_strata.py <panel.fasta> <mapping_truth.tsv> <out.tsv>
      [--threads 8] [--near 2] [--med 10]
"""
import argparse
import collections
import os
import shlex
import subprocess
import sys
import tempfile


NO_OTHER_HIT = 999999


def parse_ref(fa):
    """{contig_name: sequence}.  Name = first whitespace token of the header."""
    seqs, name, cur = {}, None, []
    with open(fa) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(cur)
                name = line[1:].split()[0]
                cur = []
            else:
                cur.append(line.strip())
        if name is not None:
            seqs[name] = "".join(cur)
    return seqs


def parse_mapping_truth(mapping_truth_tsv):
    """{read_name: (ctg, start, end)}."""
    mapping_truth = {}
    with open(mapping_truth_tsv) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 5:
                mapping_truth[p[0]] = (p[2], int(p[3]), int(p[4]))
    return mapping_truth


def update_hit_stats(stats, nm):
    """Add one full-length-hit NM without retaining every hit."""
    n_nm0, nearest_nonzero = stats
    if nm == 0:
        n_nm0 += 1
    elif nearest_nonzero is None or nm < nearest_nonzero:
        nearest_nonzero = nm
    return n_nm0, nearest_nonzero


def classify_hit_stats(stats, near, med):
    """Reproduce the original sorted-NM classification from compact stats."""
    n_nm0, nearest_nonzero = stats
    if n_nm0 == 0:
        nearest = (
            nearest_nonzero
            if nearest_nonzero is not None
            else NO_OTHER_HIT
        )
        level = "no_map"
    elif n_nm0 >= 2:
        # After dropping one exact self hit, another NM=0 hit remains.
        nearest = 0
        level = "dup"
    else:
        nearest = (
            nearest_nonzero
            if nearest_nonzero is not None
            else NO_OTHER_HIT
        )
        if nearest <= near:
            level = "hard"
        elif nearest <= med:
            level = "med"
        else:
            level = "easy"
    return level, n_nm0, nearest


def selfmap_full_length_nm(
    panel,
    seqs,
    mapping_truth,
    out,
    threads,
    tmp_dir=None,
    index_batch_size=None,
    minimap2_log=None,
):
    """Cut assigned spans, self-map to the panel, and return full-length-hit NMs."""
    temporary_root = tmp_dir or os.path.dirname(os.path.abspath(out))
    os.makedirs(temporary_root, exist_ok=True)
    minimap2_log = minimap2_log or out + ".minimap2.log"
    os.makedirs(os.path.dirname(os.path.abspath(minimap2_log)), exist_ok=True)

    with tempfile.TemporaryDirectory(
        prefix="alignment_strata.", dir=temporary_root
    ) as temporary_directory:
        spans = os.path.join(temporary_directory, "spans.fa")
        id2name = {}
        with open(spans, "w") as fh:
            for i, (name, (c, s, e)) in enumerate(mapping_truth.items()):
                if c in seqs and 0 <= s < e <= len(seqs[c]):
                    sid = f"span_{i}"
                    id2name[sid] = name
                    fh.write(f">{sid}\n{seqs[c][s:e]}\n")
        print(
            f"alignment_strata: self-mapping {len(id2name)} assigned spans ...",
            file=sys.stderr)

        command = [
            "minimap2",
            "-t", str(threads),
            "-c",
            "--secondary=yes",
            "-p0.1",
            "-N50",
        ]
        if index_batch_size:
            command.extend(("-I", index_batch_size))
        command.extend((panel, spans))
        print(
            "alignment_strata: " + shlex.join(command),
            file=sys.stderr,
        )

        hit_stats = collections.defaultdict(lambda: (0, None))
        with open(minimap2_log, "w") as log_handle:
            process = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=log_handle,
                text=True,
            )
            assert process.stdout is not None
            try:
                for line in process.stdout:
                    col = line.rstrip("\n").split("\t")
                    if len(col) < 12:
                        raise RuntimeError(
                            f"malformed minimap2 PAF line: {line[:200]!r}")
                    sid = col[0]
                    if sid not in id2name:
                        continue
                    qlen, qs, qe = int(col[1]), int(col[2]), int(col[3])
                    nm = next(
                        (
                            int(tag[5:])
                            for tag in col[12:]
                            if tag.startswith("NM:i:")
                        ),
                        None,
                    )
                    if nm is not None and (qe - qs) >= 0.99 * qlen:
                        name = id2name[sid]
                        hit_stats[name] = update_hit_stats(
                            hit_stats[name], nm)
            except BaseException:
                process.terminate()
                process.wait()
                raise
            finally:
                process.stdout.close()
            return_code = process.wait()
        if return_code != 0:
            raise subprocess.CalledProcessError(return_code, command)

        return dict(hit_stats), id2name


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("panel")
    ap.add_argument("mapping_truth")
    ap.add_argument("out")
    ap.add_argument("--threads", type=int, default=8,
                    help="number of minimap2 threads [8]")
    ap.add_argument("--near", type=int, default=2, help="hard NM threshold")
    ap.add_argument("--med", type=int, default=10, help="med NM threshold")
    ap.add_argument(
        "--tmp-dir",
        help="directory for the potentially large extracted-span FASTA")
    ap.add_argument(
        "--index-batch-size",
        help="minimap2 -I value; use 8G for the 6-Gbp diploid panel")
    ap.add_argument(
        "--minimap2-log",
        help="persistent minimap2 stderr path [OUT.minimap2.log]")
    args = ap.parse_args()

    seqs = parse_ref(args.panel)
    mapping_truth = parse_mapping_truth(args.mapping_truth)
    hit_stats, id2name = selfmap_full_length_nm(
        args.panel,
        seqs,
        mapping_truth,
        args.out,
        args.threads,
        tmp_dir=args.tmp_dir,
        index_batch_size=args.index_batch_size,
        minimap2_log=args.minimap2_log,
    )

    counts = collections.Counter()
    output_directory = os.path.dirname(os.path.abspath(args.out))
    os.makedirs(output_directory, exist_ok=True)
    descriptor, temporary_output = tempfile.mkstemp(
        prefix=f".{os.path.basename(args.out)}.",
        dir=output_directory,
        text=True,
    )
    try:
        with os.fdopen(descriptor, "w") as fh:
            for _sid, name in id2name.items():
                lvl, n0, nearest = classify_hit_stats(
                    hit_stats.get(name, (0, None)),
                    args.near,
                    args.med,
                )
                counts[lvl] += 1
                fh.write(
                    f"{name}\t{lvl}\t{n0}\t"
                    f"{nearest if nearest < NO_OTHER_HIT else -1}\n"
                )
        os.replace(temporary_output, args.out)
    except BaseException:
        try:
            os.unlink(temporary_output)
        except FileNotFoundError:
            pass
        raise
    print("alignment_strata: " + "  ".join(
        f"{k}={counts[k]}"
        for k in ["easy", "med", "hard", "dup", "no_map"]),
        file=sys.stderr)


if __name__ == "__main__":
    main()
