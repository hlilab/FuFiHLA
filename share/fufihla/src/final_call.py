#!/usr/bin/env python3
import sys
import re
import argparse
import gzip

###########################################################
# Call Closest Allele
###########################################################

def parse_gene_annot(file_path):
    """
    Parses the gene annotation file and returns a dict mapping allele
    names (e.g. "HLA-C*03:04:18") to a list of CDS intervals (1-based inclusive).
    """
    gene_data = {}
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split()
            info = {k: v for k, v in (fld.split('=',1) for fld in fields if '=' in fld)}
            if 'allele' in info and 'CDS' in info:
                allele = info['allele']
                cds_intervals = []
                for token in info['CDS'].split(','):
                    if '..' in token:
                        try:
                            s,e = token.split('..')
                            cds_intervals.append((int(s), int(e)))
                        except ValueError:
                            pass
                gene_data[allele] = cds_intervals
    return gene_data

def check_exon_coverage(cds_intervals, ref_start, ref_end, threshold=1.0):
    for s,e in cds_intervals:
        length = e - s + 1
        ov = max(0, min(e, ref_end) - max(s, ref_start) + 1)
        if ov / length < threshold:
            return False
    return True

def parse_cs(cs, ref_offset=1):
    token_pattern = re.compile(r'(:\d+)|(\*[A-Za-z]{2})|([\+\-][A-Za-z]+)')
    pos = ref_offset
    muts = []
    for m in token_pattern.finditer(cs):
        tok = m.group(0)
        if tok.startswith(':'):
            pos += int(tok[1:])
        elif tok.startswith('*'):
            muts.append(('sub', pos, tok))
            pos += 1
        elif tok.startswith('+'):
            muts.append(('ins', pos, tok))
        elif tok.startswith('-'):
            length = len(tok) - 1
            muts.append(('del', pos, tok))
            pos += length
    return muts

PENALTY_CDS     = 9999
PENALTY_NONCDS  = 1

def mutation_penalty(mutations, cds_intervals):
    pen = 0
    for _, pos, _ in mutations:
        in_cds = any(s <= pos <= e for s,e in cds_intervals)
        pen += (PENALTY_CDS if in_cds else PENALTY_NONCDS)
    return pen

def has_cds_mutation(mutations, cds_intervals):
    for _, pos, _ in mutations:
        if any(s <= pos <= e for s,e in cds_intervals):
            return True
    return False

def calculate_pid(cs):
    matches = total = 0
    token_pattern = re.compile(r'(:\d+)|(\*[A-Za-z]{2})|([\+\-][A-Za-z]+)')
    for m in token_pattern.finditer(cs):
        tok = m.group(0)
        if tok.startswith(':'):
            l = int(tok[1:])
            matches += l
            total   += l
        elif tok.startswith('*'):
            total   += 1
        elif tok.startswith('-'):
            total   += (len(tok) - 1)
    pid = (matches / total * 100) if total else 0.0
    return matches, total, pid

def find_csstr(line):
    m = re.search(r'cs:Z:(\S+)', line)
    return m.group(1) if m else ''

def find_nm(line):
    m = re.search(r'\bNM:i:(\d+)\b', line)
    return int(m.group(1)) if m else None

def opengz(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)

###########################################################
# CDS perfect-match search
###########################################################

def load_cds_perfect_hits(paf_cds_gz):
    """
    CDS PAF from:
      minimap2 -x splice:hq --ds --end-bonus 10 consensus.fa ref.cds.fa.gz

    We want alignments where the *consensus query* is fully covered AND NM:i:0:
      qstart == 0
      qend   == qlen
      NM:i:0

    Return: dict[rec_seq] = cds_allele (target name)
    """
    hits = {}
    with opengz(paf_cds_gz, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            ll = line.split("\t")
            if len(ll) < 12:
                continue
            cds_allele = ll[0]
            rec_seq = ll[5]
            try:
                qlen = int(ll[1])
                qstart = int(ll[2])
                qend = int(ll[3])
            except ValueError:
                continue

            nm = find_nm(line)
            if nm is None:
                continue

            if qstart == 0 and qend == qlen and nm == 0:
                # keep first perfect hit if multiple
                if rec_seq not in hits:
                    hits[rec_seq] = cds_allele
    return hits

###########################################################
# Main
###########################################################

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene-annot", required=True, help="gene_annot.info")
    ap.add_argument("--paf-gene", required=True, help="PAF (complete allele mapping), gz ok")
    ap.add_argument("--paf-cds", required=True, help="PAF (CDS splice mapping), gz ok")
    args = ap.parse_args()

    gene_data = parse_gene_annot(args.gene_annot)
    cds_perfect = load_cds_perfect_hits(args.paf_cds)

    da = {}

    # Read complete-allele PAF
    with opengz(args.paf_gene, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            ll = line.split('\t')
            if len(ll) < 10:
                continue

            allele  = ll[5]
            rec_seq = ll[0]

            try:
                ref_len  = float(ll[6])
                ref_s    = int(ll[7])
                ref_e    = int(ll[8])
                aln_len  = ref_e - ref_s
            except (ValueError, IndexError):
                continue

            # Skip short alignments
            if ref_len > 0 and aln_len / ref_len < 0.5:
                continue

            # If we have CDS intervals, enforce exon coverage
            if allele in gene_data:
                cds = gene_data[allele]
                if not check_exon_coverage(cds, ref_s + 1, ref_e, threshold=0.95):
                    continue

            ref_start_1based = ref_s + 1
            csstr = find_csstr(line)
            mutations = parse_cs(csstr, ref_offset=ref_start_1based)

            matches, total_bases, pid = calculate_pid(csstr)

            nm = find_nm(line)

            # Compute penalty + detect CDS mutation 
            cds_mut = False
            if allele in gene_data:
                pen = mutation_penalty(mutations, gene_data[allele])
                cds_mut = has_cds_mutation(mutations, gene_data[allele])
            else:
                pen = csstr.count('*') + csstr.count('+') + csstr.count('-')

            da.setdefault(rec_seq, []).append([allele, aln_len, pen, pid, nm, cds_mut, line])

    # Pick best per reconstructed seq
    for rec in da:
        best = sorted(da[rec], key=lambda x: (x[2], -x[1]))[0]
        allele, aln_len, pen, pid, nm, cds_mut, raw = best

        # Rule:
        # - if NM==0 OR no mutation in CDS: keep complete allele
        # - else if CDS has a full-length perfect match (NM==0) -> output CDS allele
        #   but if cds allele is the same as complete allele string, keep complete
        # - else keep complete
        chosen = allele

        if not ((nm == 0) or (cds_mut is False)):
            cds_hit = cds_perfect.get(rec)
            if cds_hit is not None:
                if cds_hit == allele:
                    chosen = allele
                else:
                    chosen = cds_hit
            else:
                chosen = allele

        print(f"{chosen}\t{aln_len}\t{pen}\t{pid:.2f}\t@{raw}")

if __name__ == '__main__':
    main()

