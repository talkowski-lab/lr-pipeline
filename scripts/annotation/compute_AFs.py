#!/usr/bin/env python3

# Copyright (c) 2009-2026, Broad Institute, Inc. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the conditions in GATK-SV's
# BSD 3-Clause license are met. See:
# https://github.com/broadinstitute/gatk-sv/blob/main/LICENSE.TXT
#
# Derived from src/sv-pipeline/05_annotation/scripts/compute_AFs.py of broadinstitute/gatk-sv (branch kj_project_gnomad_lr)
# https://github.com/broadinstitute/gatk-sv/blob/main/src/sv-pipeline/05_annotation/scripts/compute_AFs.py

"""
Compute allele frequencies for sex & population combinations given an SV VCF
"""

import argparse
import sys
from collections import Counter

import pysam

ALLOWED_POPS = set(["afr", "ami", "amr", "asj", "eas", "fin", "mid", "nfe", "sas", "oth"])


def check_motifs_canonically_unique(motifs):
    """
    Check if motifs are canonically unique.
    Two motifs are NOT canonically unique if they have the same character counts and one can be created from the
    other by rotation.
    """
    for i in range(len(motifs)):
        for j in range(i + 1, len(motifs)):
            motif1 = motifs[i]
            motif2 = motifs[j]
            if len(motif1) != len(motif2):
                continue
            if Counter(motif1) != Counter(motif2):
                continue
            doubled = motif1 + motif1
            if motif2 in doubled:
                return False
    return True


def create_pop_dict(popfile):
    """
    Makes dictionary of sample-population pairs.
    """
    pop_dict = {}
    for sample in popfile:
        pop_dict[sample.split("\t")[0]] = sample.split("\t")[1]
    return pop_dict


def get_genotype_allele_value(gt, allele_idx, field_values, cast=None):
    """
    Map a per-genotype-allele value back to a specific allele index.
    """
    if field_values is None:
        return None

    values = field_values.split(",") if isinstance(field_values, str) else list(field_values)
    for gt_allele, raw_value in zip(gt, values):
        if gt_allele != allele_idx or raw_value is None or raw_value == ".":
            continue
        return cast(raw_value) if cast is not None else raw_value

    return None


def resolve_lps_alleles(carriers, num_alleles):
    """
    Assign an LPS value to each allele index from carriers' per-haplotype LPS values, without relying on their order.
    The LPS table follows the TRGT VCF's haplotype order, which needn't match the callset GT, so a value is only
    assigned when a carrier pins it down. Unresolved alleles are None.
    """
    lps = {}
    # A carrier whose haplotypes share one LPS value fixes every allele it carries
    for gt, values in carriers:
        if len(set(values)) == 1:
            for allele_idx in gt:
                lps.setdefault(allele_idx, values[0])

    # A heterozygous carrier with one allele resolved fixes the other, repeated until nothing changes
    changed = True
    while changed:
        changed = False
        for gt, values in carriers:
            if len(gt) != 2:
                continue
            for allele_idx, other_idx in (gt, gt[::-1]):
                if allele_idx not in lps and other_idx in lps and lps[other_idx] in values:
                    remaining = list(values)
                    remaining.remove(lps[other_idx])
                    lps[allele_idx] = remaining[0]
                    changed = True

    return [lps.get(allele_idx) for allele_idx in range(num_alleles)]


def is_biallelic(record):
    """
    Check if a record has at most one alternate allele.
    """
    return len(record.alleles) <= 2


def load_par_regions(par_bed):
    """
    Read pseudoautosomal regions from a BED file as (chrom, start, end) tuples.
    """
    regions = []
    if par_bed is None:
        return regions
    with open(par_bed) as bed:
        for line in bed:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3 or line.startswith(("#", "track", "browser")):
                continue
            regions.append((fields[0], int(fields[1]), int(fields[2])))
    return regions


def in_par(record, par_regions):
    """
    Check if variant overlaps pseudoautosomal region, matching bedtools intersect overlap semantics.
    """
    sstart, send = sorted([record.start, record.stop])
    for chrom, pstart, pend in par_regions:
        if chrom != record.chrom:
            continue
        if sstart == send:
            if pstart <= sstart <= pend:
                return True
        elif sstart < pend and send > pstart:
            return True
    return False


def update_sex_freqs(record, pop=None):
    """
    Recompute allele frequencies for variants on sex chromosomes outside of PARs.
    """
    if pop is not None:
        m_prefix = "_".join([pop, "XY"])
        f_prefix = "_".join([pop, "XX"])
    else:
        m_prefix = "XY"
        f_prefix = "XX"

    m_an = record.info.get("AN_" + m_prefix, 0)
    m_ac = sum(record.info.get("AC_" + m_prefix, (0,)))

    f_an = record.info.get("AN_" + f_prefix, 0)
    f_ac = sum(record.info.get("AC_" + f_prefix, (0,)))

    adj_an = m_an + f_an
    adj_ac = m_ac + f_ac
    if adj_an > 0:
        adj_af = adj_ac / adj_an
    else:
        adj_af = 0

    if pop is None:
        record.info["AN"] = adj_an
        record.info["AC"] = (adj_ac,)
        record.info["AF"] = (adj_af,)
    else:
        record.info["AN_" + pop] = adj_an
        record.info["AC_" + pop] = (adj_ac,)
        record.info["AF_" + pop] = (adj_af,)

    return record


def gather_allele_freqs(record, samples, males_set, females_set, par_regions, pop_dict, pops, sex_chroms, lps_dict):
    """
    Wrapper to compute allele frequencies for all sex & population pairings.
    """
    # Add PAR annotation to record if applicable
    if record.chrom in sex_chroms and len(par_regions) > 0:
        if in_par(record, par_regions):
            rec_in_par = True
            record.info["par"] = True
        else:
            rec_in_par = False
    else:
        rec_in_par = False

    # Get allele frequencies for all populations
    calc_allele_freq(record, samples, lps_dict=lps_dict)
    if len(males_set) > 0:
        if record.chrom in sex_chroms and not rec_in_par:
            calc_allele_freq(record, males_set, prefix="XY", hemi=True, lps_dict=lps_dict)
        else:
            calc_allele_freq(record, males_set, prefix="XY", lps_dict=lps_dict)
    if len(females_set) > 0:
        calc_allele_freq(record, females_set, prefix="XX", lps_dict=lps_dict)

    # Adjust global allele frequencies on sex chromosomes, if famfile provided
    if record.chrom in sex_chroms and not rec_in_par and is_biallelic(record) and len(males_set) + len(females_set) > 0:
        update_sex_freqs(record)

    # Get allele frequencies per population
    if len(pops) > 0:
        for pop in pops:
            pop_samps = [s for s in samples if pop_dict.get(s, None) == pop]
            calc_allele_freq(record, pop_samps, prefix=pop, lps_dict=lps_dict)
            if len(males_set) > 0:
                if record.chrom in sex_chroms and not rec_in_par:
                    calc_allele_freq(
                        record, list([s for s in pop_samps if s in males_set]), prefix=pop + "_XY", hemi=True, lps_dict=lps_dict
                    )
                else:
                    calc_allele_freq(
                        record, list([s for s in pop_samps if s in males_set]), prefix=pop + "_XY", lps_dict=lps_dict
                    )
            if len(females_set) > 0:
                calc_allele_freq(record, list([s for s in pop_samps if s in females_set]), prefix=pop + "_XX", lps_dict=lps_dict)

            # Adjust per-pop allele frequencies on sex chromosomes, if famfile provided
            if record.chrom in sex_chroms and not rec_in_par and is_biallelic(record) and len(males_set) + len(females_set) > 0:
                update_sex_freqs(record, pop=pop)

        # Get grpmax AF biallelic sites only
        if is_biallelic(record) and (record.alts and len(record.alts) == 1):
            AFs = [record.info["AF_{0}".format(pop)][0] for pop in pops]
            grpmax = max(AFs)
            record.info["AF_grpmax"] = grpmax

            grp_label = pops[AFs.index(grpmax)]
            record.info["grpmax"] = grp_label

            record.info["AC_grpmax"] = int(record.info["AC_{0}".format(grp_label)][0])
            record.info["AN_grpmax"] = int(record.info["AN_{0}".format(grp_label)])
            record.info["nhomalt_grpmax"] = int(record.info["nhomalt_{0}".format(grp_label)])

    return record


def calc_allele_freq(record, samples, prefix=None, hemi=False, lps_dict=None):
    """
    Computes allele frequencies for a single record based on a list of samples.
    """
    if is_biallelic(record):
        # Count alleles & genotypes
        GTs = [record.samples[s]["GT"] for s in samples]
        AC = 0
        AN = 0
        n_alt_count_0 = 0
        n_alt_count_1 = 0
        n_alt_count_2 = 0
        n_gts_with_gt_0_alts = 0  # Used specifically for hemizygous sites
        for GT in GTs:
            AN += len([allele for allele in GT if allele != "." and allele is not None])
            AC += len([allele for allele in GT if allele != "." and allele != 0 and allele is not None])
            if GT == (0, 0):
                n_alt_count_0 += 1
            if (
                len([allele for allele in GT if allele == 0 and allele != "." and allele is not None]) == 1
                and len([allele for allele in GT if allele != 0 and allele != "." and allele is not None]) == 1
            ):
                n_alt_count_1 += 1
                n_gts_with_gt_0_alts += 1
            if len([allele for allele in GT if allele != 0 and allele != "." and allele is not None]) == 2:
                n_alt_count_2 += 1
                n_gts_with_gt_0_alts += 1

        # Adjust hemizygous allele number and allele count, if optioned
        if hemi:
            AN = AN // 2
            # For hemizygous sites, AC must be the sum of all non-reference *genotypes*, not alleles
            AC = n_gts_with_gt_0_alts

        # Calculate allele frequency
        if AN > 0:
            AF = AC / AN
            AF = round(AF, 6)
        else:
            AF = 0

        # Add AN, AC, and AF to INFO field
        record.info[("AN" if prefix is None else "AN_" + prefix)] = AN
        record.info[("AC" if prefix is None else "AC_" + prefix)] = (AC,)
        record.info[("AF" if prefix is None else "AF_" + prefix)] = (AF,)

        # Add NCR to INFO field
        n_samples = len(samples)
        if n_samples > 0:
            if hemi:
                ncr = 1 - (AN / n_samples)
            else:
                ncr = 1 - (AN / (2 * n_samples))
            ncr = round(ncr, 6)
        else:
            ncr = 1.0
        record.info[("NCR" if prefix is None else "NCR_" + prefix)] = ncr

        # Calculate genotype frequencies
        n_bi_genos = n_alt_count_0 + n_alt_count_1 + n_alt_count_2
        if n_bi_genos > 0:
            freq_homref = n_alt_count_0 / n_bi_genos
            freq_het = n_alt_count_1 / n_bi_genos
            freq_homalt = n_alt_count_2 / n_bi_genos
        else:
            freq_homref = 0
            freq_het = 0
            freq_homalt = 0
        if hemi:
            freq_hemialt = freq_het + freq_homalt

        # Add n_bi_genos, n_homref, nhet, nhomalt, freq_homref, freq_het, and freq_homalt to INFO field
        record.info["n_bi_genos" + ("_" + prefix if prefix else "")] = n_bi_genos
        if hemi:
            record.info["nhemiref" + ("_" + prefix if prefix else "")] = n_alt_count_0
            record.info["nhemialt" + ("_" + prefix if prefix else "")] = n_gts_with_gt_0_alts
            record.info["freq_hemiref" + ("_" + prefix if prefix else "")] = freq_homref
            record.info["freq_hemialt" + ("_" + prefix if prefix else "")] = freq_hemialt
        record.info["nhomref" + ("_" + prefix if prefix else "")] = n_alt_count_0
        record.info["nhet" + ("_" + prefix if prefix else "")] = n_alt_count_1
        record.info["nhomalt" + ("_" + prefix if prefix else "")] = n_alt_count_2
        record.info["freq_homref" + ("_" + prefix if prefix else "")] = freq_homref
        record.info["freq_het" + ("_" + prefix if prefix else "")] = freq_het
        record.info["freq_homalt" + ("_" + prefix if prefix else "")] = freq_homalt

    else:
        # Get all sample GTs and filter out missing genotypes
        GTs_raw = [(s, record.samples[s]["GT"]) for s in samples]
        valid_GTs = [(s, gt) for s, gt in GTs_raw if gt is not None and any(a is not None for a in gt)]

        # Flatten all alleles into a single list of indices
        all_alleles = [int(a) for s, gt in valid_GTs for a in gt if a is not None and a != "."]
        allele_counts = Counter(all_alleles)

        # Calculate AN, AC and AF for each allele
        AN = len(all_alleles)
        num_alleles = len(record.alts) + 1
        allele_indices = list(range(num_alleles))
        ac_allele = [int(allele_counts.get(k, 0)) for k in allele_indices]
        af_allele = [round(v / AN, 6) if AN > 0 else 0 for v in ac_allele]

        # Add AN, AC, and AF to INFO field
        record.info[("AN" if prefix is None else "AN_" + prefix)] = AN
        record.info[("AC" if prefix is None else "AC_" + prefix)] = tuple(ac_allele[1:])
        record.info[("AF" if prefix is None else "AF_" + prefix)] = tuple(af_allele[1:])

        # Add NCR to INFO field
        n_samples = len(samples)
        if n_samples > 0:
            ncr = 1 - (AN / (2 * n_samples))
            ncr = round(ncr, 6)
        else:
            ncr = 1.0
        record.info[("NCR" if prefix is None else "NCR_" + prefix)] = ncr

        if prefix is None:
            # Calculate AP_allele for each allele
            ap_allele = []
            for allele_idx in allele_indices:
                ap_val = None
                for sample_name, gt in valid_GTs:
                    ap_val = get_genotype_allele_value(
                        gt,
                        allele_idx,
                        record.samples[sample_name].get("AP"),
                        cast=float,
                    )
                    if ap_val is not None:
                        break
                ap_allele.append(ap_val if ap_val is not None else 0)
            record.info["AP_allele"] = tuple(ap_allele)

            # Calculate MC_allele and LPS_allele for single motif sites
            motifs_field = record.info.get("MOTIFS")
            if isinstance(motifs_field, (tuple, list)):
                motifs = list(motifs_field)
            else:
                motifs = motifs_field.split(",")
            if len(motifs) == 1 and check_motifs_canonically_unique(motifs):
                mc_allele = []
                for allele_idx in allele_indices:
                    mc_val = None
                    for sample_name, gt in valid_GTs:
                        mc_val = get_genotype_allele_value(
                            gt,
                            allele_idx,
                            record.samples[sample_name].get("MC"),
                            cast=int,
                        )
                        if mc_val is not None:
                            break
                    mc_allele.append(mc_val if mc_val is not None else 0)
                record.info["MC_allele"] = tuple(mc_allele)

                lps_key = (record.info.get("TRID"), motifs[0])
                if lps_dict is not None and lps_key in lps_dict:
                    carriers = []
                    for sample_name, gt in valid_GTs:
                        values = lps_dict[lps_key].get(sample_name)
                        if values is not None and None not in gt and len(values) == len(gt):
                            carriers.append((gt, values))
                    record.info["LPS_allele"] = tuple(resolve_lps_alleles(carriers, num_alleles))

    return record


def info_line(field_id, number, field_type, description):
    """
    Format a VCF INFO header line.
    """
    return f'##INFO=<ID={field_id},Number={number},Type={field_type},Description="{description}">'


def group_info_lines(suffix, label, ac_description, hemi):
    """
    Build the INFO header lines for one sex, population or population-sex group.
    """
    lines = [
        info_line(f"AN_{suffix}", "1", "Integer", f"Total number of {label} alleles genotyped."),
        info_line(f"AC_{suffix}", "A", "Integer", f"{ac_description} {label} alleles observed."),
        info_line(f"AF_{suffix}", "A", "Float", f"{label} allele frequency."),
        info_line(f"NCR_{suffix}", "1", "Float", f"{label} no-call rate."),
        info_line(
            f"n_bi_genos_{suffix}",
            "1",
            "Integer",
            f"Total number of {label} samples with complete genotypes (biallelic sites only).",
        ),
        info_line(
            f"nhomref_{suffix}",
            "1",
            "Integer",
            f"Number of {label} samples with homozygous reference genotypes (biallelic sites only).",
        ),
        info_line(
            f"nhet_{suffix}", "1", "Integer", f"Number of {label} samples with heterozygous genotypes (biallelic sites only)."
        ),
        info_line(
            f"nhomalt_{suffix}",
            "1",
            "Integer",
            f"Number of {label} samples with homozygous alternate genotypes (biallelic sites only).",
        ),
        info_line(
            f"freq_homref_{suffix}", "1", "Float", f"{label} homozygous reference genotype frequency (biallelic sites only)."
        ),
        info_line(f"freq_het_{suffix}", "1", "Float", f"{label} heterozygous genotype frequency (biallelic sites only)."),
        info_line(
            f"freq_homalt_{suffix}", "1", "Float", f"{label} homozygous alternate genotype frequency (biallelic sites only)."
        ),
    ]
    if hemi:
        lines.extend(
            [
                info_line(
                    f"nhemiref_{suffix}",
                    "1",
                    "Integer",
                    f"Number of {label} samples with hemizygous reference genotypes (biallelic sites only).",
                ),
                info_line(
                    f"nhemialt_{suffix}",
                    "1",
                    "Integer",
                    f"Number of {label} samples with hemizygous alternate genotypes (biallelic sites only).",
                ),
                info_line(
                    f"freq_hemiref_{suffix}",
                    "1",
                    "Float",
                    f"{label} hemizygous reference genotype frequency (biallelic sites only).",
                ),
                info_line(
                    f"freq_hemialt_{suffix}",
                    "1",
                    "Float",
                    f"{label} hemizygous alternate genotype frequency (biallelic sites only).",
                ),
            ]
        )
    return lines


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("vcf", help='Input vcf. Also accepts "stdin" and "-".')
    parser.add_argument("fout", help='Output vcf. Also accepts "stdout" and "-".')
    parser.add_argument(
        "-p",
        "--popfile",
        help="Two-column file of samples & population assignments. If omitted, population-specific AFs are not computed.",
    )
    parser.add_argument(
        "-f", "--famfile", help="Input .fam file (used for sex-specific AFs). If omitted, sex-specific AFs are not computed."
    )
    parser.add_argument("-l", "--lpsfile", help="TSV file of LPS values per sample (for multiallelic sites).")
    parser.add_argument("--par", help="BED file of pseudoautosomal regions (used for sex-specific AFs).")
    args = parser.parse_args()

    # Define constants
    sex_chroms = ["chrX", "chrY"]

    # Open input VCF
    if args.vcf in "- stdin".split():
        vcf = pysam.VariantFile(sys.stdin)
    else:
        vcf = pysam.VariantFile(args.vcf)

    # Get list of all samples in vcf
    samples_list = list(vcf.header.samples)

    # Get lists of males and females, if famfile provided
    if args.famfile is not None:
        sexes = ["XY", "XX"]
        famfile = [line.rstrip("\n") for line in open(args.famfile)]
        males_set = set([line.split("\t")[1] for line in famfile if line.split("\t")[4] == "1"])
        males_set = set(s for s in samples_list if s in males_set)
        females_set = set([line.split("\t")[1] for line in famfile if line.split("\t")[4] == "2"])
        females_set = set(s for s in samples_list if s in females_set)
    else:
        sexes = []
        males_set = set()
        females_set = set()

    # Get dictionary of populations, if popfile provided
    if args.popfile is not None:
        popfile = [line.rstrip("\n") for line in open(args.popfile)]
        pop_dict = create_pop_dict(popfile)
        pops = list(set(pop_dict.values()))
        pops = sorted([p for p in pops if p != "."])
        for label in pops:
            if label not in ALLOWED_POPS:
                raise ValueError(f"Invalid label: '{label}'.")
    else:
        pop_dict = {}
        pops = []

    # Get LPS values per (TRID, motif), keeping only the loci in the VCF
    lps_dict = None
    if args.lpsfile is not None:
        vcf_trids = set(rec.info.get("TRID") for rec in vcf.fetch())
        vcf.reset()
        lps_dict = {}
        # Older concatenated tables were sorted with their header, leaving it last rather than first
        lines = [line.rstrip("\n").rstrip("\r") for line in open(args.lpsfile) if line.strip()]
        header = next(line for line in lines if line.startswith("trid\t")).split("\t")
        sample_names = header[2:]
        for line in lines:
            fields = line.split("\t")
            if len(fields) < 2 or fields[0] == "trid":
                continue
            trid, motif = fields[0], fields[1]
            if trid not in vcf_trids:
                continue
            sample_values = fields[2:]
            lps_dict[(trid, motif)] = {}
            for i, sample_name in enumerate(sample_names):
                if i < len(sample_values):
                    val_str = sample_values[i].strip()
                    # A missing haplotype can't be placed without knowing the order, so drop the whole call
                    values = [v.strip() for v in val_str.split(",")]
                    if val_str == "" or "." in values:
                        lps_dict[(trid, motif)][sample_name] = None
                    else:
                        lps_dict[(trid, motif)][sample_name] = [int(v) for v in values]
                else:
                    lps_dict[(trid, motif)][sample_name] = None

    # Get PAR regions
    par_regions = load_par_regions(args.par)

    # Define base fields
    INFO_ADD = [
        info_line("AN", "1", "Integer", "Total number of alleles genotyped."),
        info_line("AC", "A", "Integer", "Number of alleles observed."),
        info_line("AF", "A", "Float", "Allele frequency."),
        info_line("NCR", "1", "Float", "No-call rate."),
        info_line("n_bi_genos", "1", "Integer", "Total number of samples with complete genotypes (biallelic sites only)."),
        info_line("nhomref", "1", "Integer", "Number of samples with homozygous reference genotypes (biallelic sites only)."),
        info_line("nhet", "1", "Integer", "Number of samples with heterozygous genotypes (biallelic sites only)."),
        info_line("nhomalt", "1", "Integer", "Number of samples with homozygous alternate genotypes (biallelic sites only)."),
        info_line("freq_homref", "1", "Float", "Homozygous reference genotype frequency (biallelic sites only)."),
        info_line("freq_het", "1", "Float", "Heterozygous genotype frequency (biallelic sites only)."),
        info_line("freq_homalt", "1", "Float", "Homozygous alternate genotype frequency (biallelic sites only)."),
        info_line("AP_allele", ".", "Float", "Allele purity for each allele index (multiallelic sites only)."),
        info_line("MC_allele", ".", "Integer", "Motif count for each allele index (multiallelic sites only)."),
    ]
    if args.lpsfile is not None:
        INFO_ADD.append(
            info_line(
                "LPS_allele",
                ".",
                "Integer",
                "Longest polymer sequence for each allele index, missing where no carrier resolves it (multiallelic sites only).",
            )
        )
    if len(par_regions) > 0:
        INFO_ADD.append(info_line("par", "0", "Flag", "Variant overlaps pseudoautosomal region."))

    # Define sex fields
    for sex in sexes:
        INFO_ADD.extend(group_info_lines(sex, sex, "Number of", sex == "XY"))

    # Define pop fields
    if len(pops) > 0:
        INFO_ADD.extend(
            [
                info_line("AN_grpmax", "1", "Integer", "Allele number for the grpmax population (biallelic sites only."),
                info_line("AC_grpmax", "1", "Integer", "Allele count for the grpmax population (biallelic sites only."),
                info_line("AF_grpmax", "1", "Float", "Maximum allele frequency across any population (biallelic sites only."),
                info_line("grpmax", "1", "String", "Population label with maximum allele frequency (biallelic sites only)."),
                info_line(
                    "nhomalt_grpmax",
                    "1",
                    "Integer",
                    "Number of homozygous-alternate genotypes in the grpmax population (biallelic sites only).",
                ),
            ]
        )
    for pop in pops:
        INFO_ADD.extend(group_info_lines(pop, pop, "Number of non-reference", False))
        for sex in sexes:
            INFO_ADD.extend(group_info_lines("_".join((pop, sex)), " ".join((pop, sex)), "Number of non-reference", sex == "XY"))

    # Add new INFO fields to header
    for line in INFO_ADD:
        vcf.header.add_line(line)

    # Prep output VCF
    if args.fout in "- stdout".split():
        fout = pysam.VariantFile(sys.stdout, "w", header=vcf.header)
    else:
        fout = pysam.VariantFile(args.fout, "w", header=vcf.header)

    # Write allele frequencies for each record to new VCF
    for r in vcf.fetch():
        newrec = gather_allele_freqs(r, samples_list, males_set, females_set, par_regions, pop_dict, pops, sex_chroms, lps_dict)
        fout.write(newrec)

    # Close VCF
    fout.close()


if __name__ == "__main__":
    main()
