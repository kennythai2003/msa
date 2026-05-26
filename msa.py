from _align3_core import ffi, lib
import time
import os
import sys
import tracemalloc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def align3(u, v, w):
    max_len = len(u) + len(v) + len(w) + 1
    aln_len = ffi.new('int *')
    perfect = ffi.new('int *')
    out_u = ffi.new('char[]', max_len)
    out_v = ffi.new('char[]', max_len)
    out_w = ffi.new('char[]', max_len)
    score = lib.msa_align3(u.encode(), len(u),
                           v.encode(), len(v),
                           w.encode(), len(w),
                           aln_len, perfect, out_u, out_v, out_w)
    aln_u = ffi.string(out_u).decode()
    aln_v = ffi.string(out_v).decode()
    aln_w = ffi.string(out_w).decode()
    return score, aln_len[0], perfect[0], aln_u, aln_v, aln_w


DATASETS = {
    "msa1": ["NM_000558.5.fna", "NM_008218.2.fna", "NM_013096.2.fna"],
    "msa2": ["NM_178850.3.fna", "NM_001030004.3.fna", "NM_001287184.2.fna"],
    "msa3": ["NM_010019.4.fna", "NM_014326.5.fna", "NM_001243563.1.fna"],
    "msa4": ["NM_000457.6.fna", "NM_008261.3.fna", "XM_016937951.3.fna"],
    "msa5": ["NM_000492.4.fna", "NM_021050.2.fna", "NM_031506.1.fna"],
}


def read_fasta(filepath):
    seq = ""
    with open(filepath) as f:
        for line in f:
            if not line.startswith(">"):
                seq += line.strip()
    return seq.upper()


def find_conserved_regions(aln_u, aln_v, aln_w, min_run=5):
    regions = []
    run_start = None
    run_len = 0
    for col in range(len(aln_u)):
        a, b, c = aln_u[col], aln_v[col], aln_w[col]
        if a == b == c and a != '-':
            if run_start is None:
                run_start = col
            run_len += 1
        else:
            if run_len >= min_run:
                regions.append((run_start, run_start + run_len - 1, run_len))
            run_start = None
            run_len = 0
    if run_len >= min_run:
        regions.append((run_start, run_start + run_len - 1, run_len))
    return regions


if __name__ == "__main__":

    for dataset, files in DATASETS.items():
        print(f"\nRunning dataset: {dataset}")
        print(f"Files: {files}\n")

        folder = dataset
        seqs = []
        for fname in files:
            path = os.path.join(folder, fname)
            seq = read_fasta(path)
            seqs.append(seq)
            print(f"  {fname}: {len(seq)} bp")

        u, v, w = seqs
        print(f"\nStarting alignment ({len(u)} x {len(v)} x {len(w)})...")

        tracemalloc.start()
        start = time.time()
        score, aln_len, perfect, aln_u, aln_v, aln_w = align3(u, v, w)
        elapsed = time.time() - start
        _, peak_bytes = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        print(f"Optimal SP score: {score}")
        print(f"Alignment length: {aln_len} columns")
        print(f"Perfect match columns: {perfect}")
        print(f"Running time: {elapsed:.2f} seconds")
        print(f"Peak memory (Python): {peak_bytes / 1024 / 1024:.2f} MB")

        regions = find_conserved_regions(aln_u, aln_v, aln_w, min_run=5)
        if regions:
            print(
                f"Conserved regions (>= 5 consecutive matches): {len(regions)}")
            for start_col, end_col, length in regions[:5]:
                print(
                    f"  col {start_col}–{end_col}  ({length} bp)  {aln_u[start_col:end_col+1]}")
            if len(regions) > 5:
                print(f"  ... and {len(regions)-5} more")
