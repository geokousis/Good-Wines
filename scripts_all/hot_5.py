#!/usr/bin/env python3
import sys
import gzip
import bisect
from tqdm import tqdm

# Optionally use cyvcf2 for VCF parsing.
try:
    from cyvcf2 import VCF
    USE_CYVCF = True
except ImportError:
    USE_CYVCF = False

########################################
# Helper Functions
########################################

def convert_gt(gt_str):
    """
    Convert a genotype string to a numeric code:
      - "0/0" or "0|0" -> "0"
      - "0/1" or "1/0" -> "1"
      - "1/1" -> "2"
    Returns None if missing.
    """
    gt_str = gt_str.replace("|", "/")
    alleles = gt_str.split("/")
    if len(alleles) < 2 or "." in alleles:
        return None
    if alleles[0] == "0" and alleles[1] == "0":
        return "0"
    elif alleles[0] == "1" and alleles[1] == "1":
        return "2"
    else:
        return "1"

########################################
# Parsing Functions
########################################

def parse_tab(file):
    """
    Parse a tab-delimited file representing a simplified VCF.
    Assumes:
      - Column 1: chromosome
      - Column 2: position (integer)
      - Columns 3+: genotype information.
    (For tab-delimited input we assume all events are SNVs.)
    
    Returns:
      events_by_chr: {chrom: [(pos, pattern)]}
      variant_info_by_chr: {chrom: [(pos, sample_id, is_indel)]} with is_indel always False.
      sample_count: number of samples.
    """
    events_by_chr = {}
    variant_info_by_chr = {}
    sample_count = None
    for line in file:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        chrom = parts[0]
        try:
            pos = int(parts[1])
        except ValueError:
            continue
        genotypes = parts[2:]
        if sample_count is None:
            sample_count = len(genotypes)
        elif sample_count != len(genotypes):
            sys.stderr.write("Warning: Inconsistent number of sample columns.\n")
        try:
            numeric = [int(x) for x in genotypes]
        except Exception as e:
            sys.stderr.write(f"Error converting genotype values: {e}\n")
            continue
        pattern = "-".join(str(x) for x in numeric)
        events_by_chr.setdefault(chrom, []).append((pos, pattern))
        nonzero_samples = [i for i, val in enumerate(numeric) if val != 0]
        sample_id = nonzero_samples[0] if len(nonzero_samples)==1 else -1
        variant_info_by_chr.setdefault(chrom, []).append((pos, sample_id, False))
    for chrom in events_by_chr:
        events_by_chr[chrom].sort(key=lambda x: x[0])
    for chrom in variant_info_by_chr:
        variant_info_by_chr[chrom].sort(key=lambda x: x[0])
    if sample_count is None:
        sample_count = 0
    return events_by_chr, variant_info_by_chr, sample_count

def parse_vcf_manual(file):
    """
    Manually parse a VCF/gVCF file.
    
    For each record:
      - Compute a genotype pattern (using convert_gt) across samples.
      - For hotspot detection, include only SNVs (records where the variant is not an indel and ALT is not "*")
        in events_by_chr.
      - Always record the record in variant_info_by_chr with an is_indel flag.
      - Also record every record's position in entries_by_chr (for continuous coverage checking).
      - Print a warning if a genotype is missing.
      
    Returns:
      events_by_chr: {chrom: [(pos, pattern)]} for SNVs only.
      variant_info_by_chr: {chrom: [(pos, sample_id, is_indel)]} for all records.
      entries_by_chr: {chrom: sorted list of positions} for every record.
      sample_count: number of samples.
    """
    events_by_chr = {}
    variant_info_by_chr = {}
    entries_by_chr = {}
    sample_count = 0
    for raw_line in file:
        line = raw_line.decode('utf-8') if isinstance(raw_line, bytes) else raw_line
        if not line or line.startswith("##"):
            continue
        if line.startswith("#"):
            cols = line.strip().split()
            if len(cols) > 8:
                sample_count = len(cols[9:])
            continue
        parts = line.strip().split('\t')
        if len(parts) < 5:
            continue
        chrom, pos_str, _id, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
        try:
            pos = int(pos_str)
        except ValueError:
            continue
        # Record every entry for continuous coverage.
        entries_by_chr.setdefault(chrom, []).append(pos)
        
        # If ALT indicates a reference call, record as reference.
        if alt in ("<NON_REF>", "."):
            variant_info_by_chr.setdefault(chrom, []).append((pos, 0, False))
            continue
        
        # Determine if record is an indel or "*" event.
        is_indel = (len(ref) != len(alt)) or (alt == "*")
        
        # For hotspot detection, include only SNVs.
        if not is_indel:
            pattern_parts = []
            unique_alt_samples = []
            if sample_count <= 1:
                sample_field = parts[9] if len(parts) > 9 else "."
                if sample_field in (".", "./."):
                    sys.stderr.write(f"Warning: Missing genotype at {chrom}:{pos}\n")
                    continue
                pat = convert_gt(sample_field)
                if pat is None:
                    continue
                pattern_parts.append(pat)
                unique_alt_samples.append(0)
                alt_sample_id = 0
            else:
                if len(parts) < 9:
                    continue
                fmt = parts[8].split(':')
                try:
                    gt_index = fmt.index("GT")
                except ValueError:
                    gt_index = 0
                for i in range(sample_count):
                    sample_field = parts[9+i] if 9+i < len(parts) else "."
                    if sample_field in (".", "./."):
                        sys.stderr.write(f"Warning: Missing genotype for sample {i} at {chrom}:{pos}\n")
                        pattern_parts.append("NA")
                        continue
                    gt_val = convert_gt(sample_field)
                    pattern_parts.append(gt_val if gt_val is not None else "NA")
                    if gt_val is not None and gt_val != "0":
                        unique_alt_samples.append(i)
                if len(unique_alt_samples) == 0:
                    alt_sample_id = 0
                elif len(unique_alt_samples) == 1:
                    alt_sample_id = unique_alt_samples[0]
                else:
                    alt_sample_id = -1
            pattern = "-".join(pattern_parts)
            events_by_chr.setdefault(chrom, []).append((pos, pattern))
        else:
            # For indels, do not add to events_by_chr, but ensure alt_sample_id is defined.
            alt_sample_id = 0
        # Always record the record.
        variant_info_by_chr.setdefault(chrom, []).append((pos, alt_sample_id, is_indel))
    for chrom in events_by_chr:
        events_by_chr[chrom].sort(key=lambda x: x[0])
    for chrom in variant_info_by_chr:
        variant_info_by_chr[chrom].sort(key=lambda x: x[0])
    for chrom in entries_by_chr:
        entries_by_chr[chrom].sort()
    return events_by_chr, variant_info_by_chr, entries_by_chr, sample_count

def parse_vcf(file):
    """
    Wrapper for VCF parsing.
    Uses cyvcf2 if available; otherwise, uses manual parsing.
    """
    if USE_CYVCF:
        return parse_vcf_cyvcf(file)
    else:
        return parse_vcf_manual(file)

def parse_vcf_cyvcf(file):
    """
    Parse VCF using cyvcf2.
    
    For each record:
      - Compute a genotype pattern using convert_gt.
      - For hotspot detection, include only SNVs.
      - Record every record in variant_info_by_chr (with is_indel flag) and in entries_by_chr.
    """
    events_by_chr = {}
    variant_info_by_chr = {}
    entries_by_chr = {}
    vcf = VCF(file)
    sample_count = len(vcf.samples)
    for rec in tqdm(vcf, desc="Parsing VCF with cyvcf2", unit="rec"):
        chrom = rec.CHROM
        pos = rec.POS  # 1-based
        entries_by_chr.setdefault(chrom, []).append(pos)
        if rec.ALT is None or rec.ALT[0] in (".", "<NON_REF>"):
            variant_info_by_chr.setdefault(chrom, []).append((pos, 0, False))
            continue
        alt = rec.ALT[0]
        is_indel = (len(rec.REF) != len(alt)) or (alt == "*")
        pattern_parts = []
        alt_samples = []
        for i, gt in enumerate(rec.gt_bases):
            if gt in (".", "./."):
                sys.stderr.write(f"Warning: Missing genotype for sample {i} at {chrom}:{pos}\n")
                pattern_parts.append("NA")
                continue
            val = convert_gt(gt)
            pattern_parts.append(val if val is not None else "NA")
            if val is not None and val != "0":
                alt_samples.append(i)
        pattern = "-".join(pattern_parts)
        alt_sample_id = alt_samples[0] if len(alt_samples)==1 else -1
        if not is_indel:
            events_by_chr.setdefault(chrom, []).append((pos, pattern))
        variant_info_by_chr.setdefault(chrom, []).append((pos, alt_sample_id, is_indel))
    for chrom in events_by_chr:
        events_by_chr[chrom].sort(key=lambda x: x[0])
    for chrom in variant_info_by_chr:
        variant_info_by_chr[chrom].sort(key=lambda x: x[0])
    for chrom in entries_by_chr:
        entries_by_chr[chrom].sort()
    return events_by_chr, variant_info_by_chr, entries_by_chr, sample_count

########################################
# Hotspot & Flank Functions
########################################

def detect_diverse_hotspots(events_by_chr, window=250, min_unique=3):
    """
    Detect hotspots characterized by high diversity of genotype patterns.
    
    For each chromosome, iterate over events (each is (pos, pattern)). For each starting event,
    look ahead within a given window (e.g., 250 bp) and collect the unique genotype patterns.
    If the number of unique patterns is at least min_unique, record the hotspot.
    
    Returns a list of hotspots with keys:
      'chrom', 'start', 'end', 'unique_count', and 'patterns' (list of unique patterns).
    """
    hotspots = []
    for chrom, events in tqdm(events_by_chr.items(), desc="Detecting diverse hotspots", unit="chrom"):
        n = len(events)
        for i in range(n):
            pos_i, pattern_i = events[i]
            cluster_patterns = {pattern_i}
            cluster_start = pos_i
            cluster_end = pos_i
            for j in range(i+1, n):
                pos_j, pattern_j = events[j]
                if pos_j > pos_i + window:
                    break
                cluster_patterns.add(pattern_j)
                cluster_end = pos_j
            if len(cluster_patterns) >= min_unique:
                hotspots.append({
                    "chrom": chrom,
                    "start": cluster_start,
                    "end": cluster_end,
                    "unique_count": len(cluster_patterns),
                    "patterns": list(cluster_patterns)
                })
    return hotspots

def adjust_hotspot_window(hotspot, window_size, flank_size=20):
    """
    Center the hotspot (using its core boundaries) and adjust it to a fixed window plus flanks.
    The new region spans (window_size + 2*flank_size) bp centered on the core hotspot.
    """
    center = (hotspot["start"] + hotspot["end"]) // 2
    final_length = window_size + 2 * flank_size
    new_start = center - final_length // 2
    new_end = new_start + final_length
    hotspot["start"] = new_start
    hotspot["end"] = new_end
    return hotspot

def region_is_continuously_covered(chrom, start, end, entries_by_chr):
    """
    Check if every base in the region [start, end] is covered by an entry in the gVCF.
    No missing bases are allowed.
    """
    if chrom not in entries_by_chr:
        return False
    all_positions = entries_by_chr[chrom]
    for pos in range(start, end + 1):
        idx = bisect.bisect_left(all_positions, pos)
        if idx >= len(all_positions) or all_positions[idx] != pos:
            return False
    return True

def flank_is_clean(chrom, start, end, variant_info_by_chr, allowed_variants=1):
    """
    Check if the flank region [start, end] is continuously covered and "clean."
    "Clean" means:
      - Every base in the region must have a corresponding record in variant_info_by_chr (i.e. no missing bases).
      - In the flank, at most allowed_variants records are variant calls (i.e. where sample_id != 0),
        and no record is an indel.
    """
    if chrom not in variant_info_by_chr:
        return False
    records = variant_info_by_chr[chrom]
    idx = bisect.bisect_left(records, (start, -float('inf'), False))
    variant_count = 0
    for pos in range(start, end + 1):
        if idx < len(records) and records[idx][0] == pos:
            _, sample_id, is_indel = records[idx]
            if is_indel:
                return False
            if sample_id != 0:
                variant_count += 1
                if variant_count > allowed_variants:
                    return False
            idx += 1
        else:
            # No record found for this base; region is not continuously covered.
            return False
    return True

def region_has_no_indels(chrom, start, end, variant_info_by_chr):
    """
    Check if the entire region [start, end] has no indel (or '*' event) records.
    """
    if chrom not in variant_info_by_chr:
        return True
    records = variant_info_by_chr[chrom]
    idx = bisect.bisect_left(records, (start, -float('inf'), False))
    while idx < len(records) and records[idx][0] <= end:
        _, _, is_indel = records[idx]
        if is_indel:
            return False
        idx += 1
    return True

########################################
# Main Routine
########################################

def main():
    if len(sys.argv) < 2:
        print("Usage: python hot_3.py <input_file> [window_size]")
        sys.exit(1)
    input_path = sys.argv[1]
    # window_size here refers to the core window size for output adjustment.
    window_size = 250  
    if len(sys.argv) >= 3:
        try:
            window_size = int(sys.argv[2])
        except Exception as e:
            sys.stderr.write(f"Invalid window size provided, using default 250 bp: {e}\n")
    
    # Determine input type: VCF/gVCF (including gzipped) or tab-delimited.
    is_vcf = input_path.lower().endswith(('.vcf', '.g.vcf', '.gvcf', '.vcf.gz'))
    try:
        if is_vcf:
            f = gzip.open(input_path, 'rt') if input_path.lower().endswith('.gz') else open(input_path, 'r')
            events_by_chr, variant_info_by_chr, entries_by_chr, sample_count = parse_vcf(f)
            f.close()
        else:
            f = open(input_path, 'r')
            events_by_chr, variant_info_by_chr, sample_count = parse_tab(f)
            # For tab-delimited, assume every event is covered.
            entries_by_chr = {chrom: [pos for pos, _ in events] for chrom, events in events_by_chr.items()}
            f.close()
    except Exception as e:
        sys.stderr.write(f"Error reading input file: {e}\n")
        sys.exit(1)
    
    # Detect diverse (pattern) hotspots.
    hotspots = detect_diverse_hotspots(events_by_chr, window=250, min_unique=3)
    
    # Adjust each hotspot to a fixed window (core window_size plus flanks).
    final_hotspots = [adjust_hotspot_window(hs, window_size) for hs in hotspots]
    
    # Filter out hotspots where the entire region (core+flanks) is not continuously covered.
    if not all(region_is_continuously_covered(hs["chrom"], hs["start"], hs["end"], entries_by_chr) for hs in final_hotspots):
        final_hotspots = [hs for hs in final_hotspots if region_is_continuously_covered(hs["chrom"], hs["start"], hs["end"], entries_by_chr)]
    
    # Filter out hotspots where the flanking regions (outer 20 bp on each side) are not "clean."
    final_hotspots = [hs for hs in final_hotspots if 
                      flank_is_clean(hs["chrom"], hs["start"], hs["start"]+20-1, variant_info_by_chr, allowed_variants=1)
                      and flank_is_clean(hs["chrom"], hs["end"]-20+1, hs["end"], variant_info_by_chr, allowed_variants=1)]
    
    # Additionally, filter out hotspots if any record in the entire region is an indel.
    final_hotspots = [hs for hs in final_hotspots if region_has_no_indels(hs["chrom"], hs["start"], hs["end"], variant_info_by_chr)]
    
    # Output in BED format (0-based start, 1-based end).
    for entry in final_hotspots:
        bed_start = max(entry['start'] - 1, 0)
        bed_end = entry['end']
        print(f"{entry['chrom']}\t{bed_start}\t{bed_end}\t{entry['unique_count']}\tTrue\t{'|'.join(entry['patterns'])}")

if __name__ == "__main__":
    main()
