#!/usr/bin/env python3
import sys, os, re, gzip
import glob, subprocess
from collections import defaultdict

if len(sys.argv) != 6:
    sys.stderr.write(
        "Usage: call_variants_bcf.py "
        "<template_list> <filtered_paf.gz> "
        "<reads_fa.gz> <all_allele_seq.fa.gz> <outdir>\n"
    )
    sys.exit(1)

template_file    = sys.argv[1]
filtered_paf_gz  = sys.argv[2]
reads_fa         = sys.argv[3]
allele_refs_fa   = sys.argv[4]
OUTDIR           = sys.argv[5]

SAMPLE = os.path.basename(OUTDIR)
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
KEEP_MAJOR = os.path.join(SCRIPT_DIR, "keep_major.py")




ds_regex = re.compile(r'ds:Z:(.*)')
err_regex = re.compile(r'([+\-])([ACGTNacgtn\[\]]+)')
###########################################################
# Load Template Alleles from sys.argv[1]
###########################################################
def load_templates(template_file):
    template_alleles = {}
    with open(template_file, 'rt') as fp:
        for line in fp:
            allele = line.strip()
            if not allele:
                continue
            gene = allele.split('*')[0]
            template_alleles.setdefault(gene, []).append(allele)
    return template_alleles

###########################################################
# Candidate Alignment Processing from s4_filtered.paf.gz
###########################################################
# Read candidate alignments from sys.argv[2] (a gzip file).
# Extract alignment information from template allele 

dat = {}
template_alleles = load_templates(sys.argv[1])
with gzip.open(sys.argv[2], 'rt') as fp:
    for line in fp:
        line = line.strip()
        if not line:
            continue
        llst = line.split('\t')
        # In s4_filtered.paf.gz:
        # Field[0]: candidate allele (template allele)
        # Field[1]: allele length.
        # Field[5]: read identifier.
        # Field[7]: ref_start.
        # Field[8]: ref_end.
        read = llst[5]
        allele = llst[0]
        gene = allele.split('*')[0]
        # Only keep candidate alignments if the allele is in the template set.
        if allele not in template_alleles.get(gene, []):
            continue
        allele_len = int(llst[1])
        ref_start = llst[7]
        ref_end = llst[8]
        match_len = int(ref_end) - int(ref_start)
        if match_len <= 0:
            continue
        nm_field = llst[12].split(':')
        if len(nm_field) < 3:
            continue
        nm = int(nm_field[2])
        m = ds_regex.search(line)
        if not m:
            continue
        ds_str = m.group(1)
        # Count raw mismatch symbols.
        mismatch = ds_str.count('*') + ds_str.count('+') + ds_str.count('-')
        # Ignore mismatch for bracketed events less than or equal to 2 bp -- possible sequencing error
        for sign, chunk in err_regex.findall(ds_str):
            stripped = re.sub(r'[\[\]]', '', chunk)
            if len(stripped) <= 2:
                mismatch -= 1
        mismatch_rate = round(mismatch * 1.0 / match_len, 7)
        coverage = match_len * 1.0 / allele_len
        # Only keep candidate alignments with mismatch_rate < 0.03.
        if mismatch_rate < 0.03:
            dat.setdefault(read, []).append([
                gene,          # [0]
                allele,        # [1]
                match_len,     # [2]
                mismatch,      # [3]
                nm,            # [4]
                coverage,      # [5]
                mismatch_rate, # [6]
                ref_start,     # [7] (string)
                ref_end,       # [8] (string)
                ds_str         # [9]
            ])

###########################################################
# Read Assignment by Mismatch in the Overlapped Mapping Interval
###########################################################
def subtract_mismatches_from_start(ds_str, length):
    events = parse_ds_string_full(ds_str, ref_start=0)
    mismatch_sub = 0
    nm_sub = 0
    region_start = 0
    region_end = length
    for e in events:
        e_start = e['ref_start']
        e_end = e['ref_end']
        if e_end <= region_start:
            continue
        if e_start >= region_end and e['type'] != 'insertion':
            break
        overlap_len = min(e_end, region_end) - max(e_start, region_start)
        if e['type'] == 'mismatch':
            mismatch_sub += 1
            nm_sub += 1
        elif e['type'] == 'insertion':
            if e_start < region_end:
                stripped = re.sub(r'[\[\]]', '', e['all_bases'])
                if len(stripped) > 0:
                    mismatch_sub += 1
                nm_sub += len(stripped)
        elif e['type'] == 'deletion':
            if overlap_len > 0:
                stripped = re.sub(r'[\[\]]', '', e['all_bases'])
                if len(stripped) > 0:
                    mismatch_sub += 1
                nm_sub += len(stripped)
    return mismatch_sub, nm_sub

def subtract_mismatches_from_end(ds_str, length):
    events = parse_ds_string_full(ds_str, ref_start=0)
    if not events:
        return 0, 0
    mismatch_sub = 0
    nm_sub = 0
    ref_alignment_length = max(e['ref_end'] for e in events)
    region_end = ref_alignment_length
    region_start = max(0, ref_alignment_length - length)
    for e in events:
        e_start = e['ref_start']
        e_end = e['ref_end']
        if e_end <= region_start:
            continue
        if e_start >= region_end and e['type'] != 'insertion':
            continue
        overlap_len = min(e_end, region_end) - max(e_start, region_start)
        if e['type'] == 'mismatch':
            mismatch_sub += 1
            nm_sub += 1
        elif e['type'] == 'insertion':
            if e_start < region_end:
                stripped = re.sub(r'[\[\]]', '', e['all_bases'])
                if len(stripped) > 0:
                    mismatch_sub += 1
                nm_sub += len(stripped)
        elif e['type'] == 'deletion':
            if overlap_len > 0:
                stripped = re.sub(r'[\[\]]', '', e['all_bases'])
                if len(stripped) > 0:
                    mismatch_sub += 1
                nm_sub += len(stripped)
    return mismatch_sub, nm_sub

ds_pattern = re.compile(r"""
    :(\d+)
  | \*([ACGTNacgtn])([ACGTNacgtn])
  | \+([ACGTNacgtn\[\]]+)
  | -([ACGTNacgtn\[\]]+)
""", re.VERBOSE)

def parse_ds_string_full(ds_str, ref_start=0):
    events = []
    ref_pos = ref_start
    tokens = ds_pattern.findall(ds_str)
    for (match_len, snp_from, snp_to, ins_chunk, del_chunk) in tokens:
        if match_len:
            length = int(match_len)
            events.append({
                'type': 'match',
                'ref_start': ref_pos,
                'ref_end': ref_pos + length
            })
            ref_pos += length
            continue
        if snp_from and snp_to:
            events.append({
                'type': 'mismatch',
                'ref_start': ref_pos,
                'ref_end': ref_pos + 1
            })
            ref_pos += 1
            continue
        if ins_chunk:
            events.append({
                'type': 'insertion',
                'ref_start': ref_pos,
                'ref_end': ref_pos,
                'all_bases': ins_chunk
            })
            continue
        if del_chunk:
            length = len(del_chunk)
            events.append({
                'type': 'deletion',
                'ref_start': ref_pos,
                'ref_end': ref_pos + length,
                'all_bases': del_chunk
            })
            ref_pos += length
            continue
    return events

def csstr_to_var(csstr, fro, to) :
    cslst = re.split(r'(:|\*|\+|\-|\~)', csstr)[1:]
    cslst_L = int(len(cslst)/2)
    cur_pos = fro
    out = [fro, to]
    for i in range(cslst_L) :
        tag = cslst[i*2]
        signal = cslst[i*2+1]
        if tag == ':' :
            cur_pos += int(signal)
        elif tag == '*' :
            out.append(str(cur_pos)+':'+'*'+ signal)
            cur_pos += 1
        elif tag == '+' :
            out.append(str(cur_pos)+':'+'+'+ signal)
        elif tag == '-' :
            out.append(str(cur_pos)+':'+'-'+ signal)
            cur_pos += len(signal)
    return(out)

best_for_read = {}

for read, cand_list in dat.items():
    if not cand_list:
        continue

    # group candidates by gene
    by_gene = {}
    for c in cand_list:
        by_gene.setdefault(c[0], []).append(c)

    per_gene_best_sets = {}  # gene -> set(tuple(candidate))

    for gene, cand in by_gene.items():
        if not cand:
            continue

        allele_best = cand[0]
        best_set = {tuple(allele_best)}

        for allele in cand[1:]:
            if allele[1] == allele_best[1]:
                continue

            mismatch0  = allele_best[3]
            nm0        = allele_best[4]
            start_best = int(allele_best[7])
            end_best   = int(allele_best[8])
            ds_str0    = allele_best[9]

            mismatch1  = allele[3]
            nm1        = allele[4]
            start_curr = int(allele[7])
            end_curr   = int(allele[8])
            ds_str1    = allele[9]

            # Adjust start
            if start_curr > start_best:
                change_start = start_curr - start_best
                sub_m0, sub_nm0 = subtract_mismatches_from_start(ds_str0, change_start)
                mismatch0 -= sub_m0
                nm0       -= sub_nm0
            elif start_curr < start_best:
                change_start = start_best - start_curr
                sub_m1, sub_nm1 = subtract_mismatches_from_start(ds_str1, change_start)
                mismatch1 -= sub_m1
                nm1       -= sub_nm1

            # Adjust end
            if end_curr > end_best:
                change_end = end_curr - end_best
                sub_m1, sub_nm1 = subtract_mismatches_from_end(ds_str1, change_end)
                mismatch1 -= sub_m1
                nm1       -= sub_nm1
            elif end_curr < end_best:
                change_end = end_best - end_curr
                sub_m0, sub_nm0 = subtract_mismatches_from_end(ds_str0, change_end)
                mismatch0 -= sub_m0
                nm0       -= sub_nm0

            if mismatch1 < mismatch0 or (mismatch1 == mismatch0 and nm1 < nm0):
                allele_best = allele
                best_set = {tuple(allele_best)}
            elif mismatch1 == mismatch0 and nm1 == nm0:
                best_set.add(tuple(allele))

        per_gene_best_sets[gene] = best_set

    best_for_read[read] = {
        gene: [(al[1], al[5], al[3], al[4]) for al in best_set]
        for gene, best_set in per_gene_best_sets.items()
    }

########################################################################
# Write read->allele mapping
########################################################################


with open(os.path.join(OUTDIR, "read_assignment.txt"), "w") as f:
    for read, per_gene in best_for_read.items():
        allele_names = []
        for gene, best_alleles in per_gene.items():
            allele_names.extend([allele for allele, _, _, _ in best_alleles])
        f.write(f"{read}\t{','.join(allele_names)}\n")

allele_to_reads = defaultdict(set)
allele_support  = defaultdict(float)

for read_id, per_gene in best_for_read.items():
    if not per_gene:
        continue

    for gene, gene_winners in per_gene.items():
        if not gene_winners:
            continue
        k = len(gene_winners)
        frac = 1.0 / k
        for allele_name, _, _, _ in gene_winners:
            allele_to_reads[allele_name].add(read_id)
            allele_support[allele_name] += frac

with open(os.path.join(OUTDIR, "allele_to_reads.txt"), "w") as out:
    for allele, reads in allele_to_reads.items():
        out.write(f"{allele}\t{','.join(sorted(reads))}\n")


###########################################################
# Consensus Sequence Reconstruction based on Templates
###########################################################
INPUT_READS_FA = reads_fa
REF_ALLELES_FA = allele_refs_fa
THREADS          = 4
MM2OPTS          = ["-ax","asm20","-B2","--ds",f"-t{THREADS}","--end-bonus=20"]
GENE_LIST_PATH = os.environ.get("GENE_LIST") 
with open(GENE_LIST_PATH) as f: 
    GENES = [line.strip() for line in f if line.strip()]
HOMOZYGOUS_THRESHOLD = float(os.environ.get("FUFIHLA_HOM_THRESHOLD", "4"))
print(f"Homozygous ratio threshold used: {HOMOZYGOUS_THRESHOLD}", file=sys.stderr)

subdirs = ["lists","fas","bam","vcf","consensus","paf","vcf_log"]
for sub in subdirs:
    os.makedirs(os.path.join(OUTDIR, sub), exist_ok=True)


def first3_fields(allele: str) -> str:
    return ":".join(allele.split(":")[:3])

# Detect and force homozygosity when one allele has ≥4× more reads
gene_alleles = {}
for gene in GENES:
    # Find all called alleles for this gene
    alleles = sorted(a for a in allele_to_reads if a.startswith(gene + "*"))
    if len(alleles) == 2:
        r0, r1 = alleles[0], alleles[1]
        n0 = allele_support.get(r0, 0.0)
        n1 = allele_support.get(r1, 0.0)
        if n1 == 0.0 and n0 > 0.0:
            gene_alleles[gene] = [r0, r0]
        elif n0 == 0.0 and n1 > 0.0:
            gene_alleles[gene] = [r1, r1]
        elif n0 >= HOMOZYGOUS_THRESHOLD * n1:
            gene_alleles[gene] = [r0, r0]
        elif n1 >= HOMOZYGOUS_THRESHOLD * n0:
            gene_alleles[gene] = [r1, r1]
        else:
            gene_alleles[gene] = [r0, r1]
    else:
        gene_alleles[gene] = alleles

new_allele_list = []
known_allele_list = []

for gene in GENES:
    alleles = gene_alleles.get(gene, [])
    if len(alleles) != 2:
        print(f"Warning: expected 2 alleles for {gene}, found {len(alleles)} → {alleles}", file=sys.stderr)


    # If first three fields are the same but fourth field different
    same_first3 = False
    template_allele = None
    merged_reads = None

    if len(alleles) == 2:
        a0, a1 = alleles[0], alleles[1]
        same_first3 = (a0 != a1) and (first3_fields(a0) == first3_fields(a1))
        #same_first3 = (first3_fields(a0) == first3_fields(a1))
        if same_first3:
            # Deterministic choice of template
            template_allele = min(a0, a1)
            merged_reads = set(allele_to_reads.get(a0, set())) | set(allele_to_reads.get(a1, set()))
    
    if same_first3:
        allele = template_allele
        allele_safe = allele.replace(":", "_")
        reads_set = merged_reads

        list_file    = os.path.join(OUTDIR, "lists", f"{gene}_{allele_safe}.list")
        reads_fa_out = os.path.join(OUTDIR, "fas",  f"{gene}_{allele_safe}_reads.fa")
        ref_normal   = os.path.join(OUTDIR, "fas",  f"{allele}.fa")
        ref_fa       = os.path.join(OUTDIR, "fas",  f"{allele_safe}.fa")
        bam          = os.path.join(OUTDIR, "bam",  f"{gene}_shared.bam")
        vcf          = os.path.join(OUTDIR, "vcf",  f"{gene}_shared.vcf.gz")
        vcf_log      = os.path.join(OUTDIR, "vcf_log", f"{gene}_shared.vcf")

        cons_fa_h1   = os.path.join(OUTDIR, "consensus", f"{gene}_asm1.fa")
        cons_fa_h2   = os.path.join(OUTDIR, "consensus", f"{gene}_asm2.fa")

        with open(list_file, "w") as lf:
            lf.write("\n".join(sorted(reads_set)) + "\n")

        with open(reads_fa_out, "w") as rf:
            subprocess.run(["seqtk", "subseq", INPUT_READS_FA, list_file],
                           stdout=rf, check=True)

        ref_list = os.path.join(OUTDIR, "lists", f"{gene}_{allele_safe}_ref.list")
        with open(ref_list, "w") as rl:
            rl.write(allele_safe + "\n")

        with open(ref_fa, "w") as rf2:
            cmd = (
                f"gzip -dc {allele_refs_fa} | "
                f"sed 's/:/_/g' | "
                f"seqtk subseq - {ref_list}"
            )
            subprocess.run(cmd, shell=True, stdout=rf2, check=True)

        with open(ref_normal, "w") as rfn:
            subprocess.run(f"sed '1s/_/:/g' {ref_fa}",
                           shell=True, check=True, stdout=rfn)

        subprocess.run(["samtools", "faidx", ref_fa], check=True)

        p1 = subprocess.Popen(["minimap2"] + MM2OPTS + [ref_fa, reads_fa_out],
                              stdout=subprocess.PIPE)
        p2 = subprocess.Popen(["samtools", "sort", "-o", bam, f"-@{THREADS}", "-m2g"],
                              stdin=p1.stdout)
        p1.stdout.close()
        p2.communicate()
        subprocess.run(["samtools", "index", bam], check=True)

        cl = subprocess.Popen(["longcallD", "call", ref_fa, bam],
                              stdout=subprocess.PIPE)

        tee = subprocess.Popen(["tee", vcf_log],
                               stdin=cl.stdout, stdout=subprocess.PIPE)
        cl.stdout.close()

        with open(vcf, "wb") as vf:
            subprocess.run(["bgzip", "-c"],
                           stdin=tee.stdout, stdout=vf, check=True)
        tee.stdout.close()

        subprocess.run(["bcftools", "index", "-t", vcf], check=True)

        with gzip.open(vcf, "rt") as vaf:
            has_variant = any(line.strip() and not line.startswith('#') for line in vaf)

        for hap, cons_fa in [(1, cons_fa_h1), (2, cons_fa_h2)]:
            with open(cons_fa, "w") as cf:
                p1 = subprocess.Popen(["bcftools", "consensus", "-H", str(hap), "-f", ref_fa, vcf],
                                      stdout=subprocess.PIPE, text=True)
                p2 = subprocess.Popen(["sed", fr"s/^>/>cons_h{hap}_/"],
                                      stdin=p1.stdout, stdout=cf, text=True)
                p1.stdout.close()
                p2.communicate()

        if has_variant:
            new_allele_list.extend([cons_fa_h1, cons_fa_h2])
        else:
            known_allele_list.extend([cons_fa_h1, cons_fa_h2])

        continue

    else:
        for idx, allele in enumerate(alleles, start=1):
            allele_safe = allele.replace(":", "_")
            list_file = os.path.join(OUTDIR,    "lists",      f"{allele_safe}.list")
            reads_fa_out = os.path.join(OUTDIR, "fas", f"{allele_safe}_reads.fa")
            ref_normal   = os.path.join(OUTDIR, "fas", f"{allele}.fa")
            ref_fa    = os.path.join(OUTDIR, "fas", f"{allele_safe}.fa")
            bam       = os.path.join(OUTDIR,    "bam",        f"{gene}_asm{idx}.bam")
            vcf       = os.path.join(OUTDIR,    "vcf",        f"{gene}_asm{idx}.vcf.gz")
            cons_fa   = os.path.join(OUTDIR,    "consensus",  f"{gene}_asm{idx}.fa")

            # 1) Write the list of reads
            with open(list_file, "w") as lf:
                lf.write("\n".join(sorted(allele_to_reads.get(allele, []))) + "\n")

            # 2) Extract reads FASTA
            reads_fa_out = os.path.join(OUTDIR, "fas", f"{gene}_{allele_safe}_reads.fa")
            with open(reads_fa_out, "w") as rf:
                subprocess.run(
                    ["seqtk", "subseq", INPUT_READS_FA, list_file],
                    stdout=rf, check=True
                )

            # 3) Extract reference allele sequence via seqtk:
            # 3a) Write a one‑line “ref_list” of the sanitized (":" to "_") allele name
            ref_list = os.path.join(OUTDIR, "lists", f"{gene}_{allele_safe}_ref.list")
            with open(ref_list, "w") as rl:
                rl.write(allele_safe + "\n")


            # 3b) Decompress the big allele FASTA, sanitize headers, and subseq
            with open(ref_fa, "w") as rf2:
                cmd = (
                    f"gzip -dc {allele_refs_fa} | "
                    f"sed 's/:/_/g' | "
                    f"seqtk subseq - {ref_list}"
                )
                subprocess.run(cmd, shell=True, stdout=rf2, check=True)


            with open(ref_normal, "w") as rfn:
                subprocess.run(
                    f"sed '1s/_/:/g' {ref_fa}",
                    shell=True,
                    check=True,
                    stdout=rfn
                )


            # 3c) Index the per‑allele FASTA for downstream tools
            subprocess.run(["samtools", "faidx", ref_fa], check=True)       

            # 4) Map reads → allele, sort & index
            p1 = subprocess.Popen(
                ["minimap2"] + MM2OPTS + [ref_fa, reads_fa_out],
                stdout=subprocess.PIPE
            )
            p2 = subprocess.Popen(
                ["samtools","sort", "-o", bam, f"-@{THREADS}", "-m2g"],
                stdin=p1.stdout
            )
            p1.stdout.close()
            p2.communicate()
            subprocess.run(["samtools","index", bam], check=True)

            # 5) Call variants & index using longcallD by Yan Gao 
            # "https://github.com/yangao07/longcallD"
            
            cl = subprocess.Popen(
                ["longcallD", "call", ref_fa, bam],
                stdout=subprocess.PIPE
            )

            # write raw VCF for log purposes
            vcf_log = os.path.join(OUTDIR, "vcf_log", f"{gene}_asm{idx}.vcf")
            tee = subprocess.Popen(
                ["tee", vcf_log],
                stdin=cl.stdout,
                stdout=subprocess.PIPE
            )
            cl.stdout.close()

            km = subprocess.Popen(
                ["python3", KEEP_MAJOR],
                stdin=tee.stdout,
                stdout=subprocess.PIPE
            )
            tee.stdout.close()

            # compress & write final VCF
            with open(vcf, "wb") as vf:
                subprocess.run(
                    ["bgzip","-c"],
                    stdin=km.stdout,
                    stdout=vf,
                    check=True
                )
            km.stdout.close()

            # index it
            subprocess.run(["bcftools","index","-t", vcf], check=True)

            # 6a) Check whether the vcf is empty; if so, skip consensus/mapping
            with gzip.open(vcf, "rt") as vaf:
                has_variant = any(line.strip() and not line.startswith('#') for line in vaf)
            if has_variant:
                new_allele_list.append(cons_fa)
            else:
                known_allele_list.append(cons_fa) 

            # 6b) Build consensus
            with open(cons_fa, "w") as cf:
                p1 = subprocess.Popen(
                    ["bcftools", "consensus", "-f", ref_fa, vcf],
                    stdout=subprocess.PIPE,
                    text=True
                )
                p2 = subprocess.Popen(
                    ["sed", fr"s/^>/>cons_h{idx}_/"],
                    stdin=p1.stdout,
                    stdout=cf,
                    text=True
                )
                p1.stdout.close()
                p2.communicate()         
 
with open(os.path.join(OUTDIR, "new_allele.fa"), "w") as out:
    for fname in new_allele_list:
        with open(fname, "r") as infile:
            out.write(infile.read())

def normalize_header(header_line: str) -> str:
    token = header_line.strip().split()[0]
    allele = token.lstrip('>')

    for prefix in ("cons_h1_", "cons_h2_", "cons_"):
        if allele.startswith(prefix):
            allele = allele[len(prefix):]
            break

    return allele.replace('_', ':')

with open(os.path.join(OUTDIR, "known_allele.paf"), "w") as out_paf:
    for cons_fa in known_allele_list:
        # Read the header from the consensus file to get the allele name
        with open(cons_fa) as fh:
            header = next(l for l in fh if l.startswith(">"))
        allele_name = normalize_header(header)
        allele_fa = os.path.join(OUTDIR, "fas", f"{allele_name}.fa")

        # Map assembled consensus (query) to the known allele (reference)
        subprocess.run(
            ["minimap2", "-c", "--cs", allele_fa, cons_fa],
            stdout=out_paf,
            check=True
        )
