#!/usr/bin/env python3
# Derived from scripts/makeigvpesr.py of yuliamostovoy/variant-interpretation
# https://github.com/yuliamostovoy/variant-interpretation/blob/main/scripts/makeigvpesr.py

import argparse
import os

import numpy as np
import pandas as pd

parser = argparse.ArgumentParser("makeigvsplit_trio.py")
parser.add_argument("-v", "--varfile", type=str, help="variant file including CHR, POS, END and SVID")
# parser.add_argument('-n', '--nestedrepeats', type=str, help='nested repeats sequences')
# parser.add_argument('-s', '--simplerepeats', type=str, help='simple repeats sequences')
# parser.add_argument('-e', '--emptytrack', type=str, help='empty track')
# parser.add_argument('-f', '--fasta', type=str, help='reference sequences')
# parser.add_argument('sample', type=str, help='name of sample to make igv on')
parser.add_argument("-fam_id", "--fam_id", type=str, help="family to plot")
parser.add_argument("-p", "--ped", type=str, help="ped file")
# parser.add_argument('cram_list', type=str, help='comma separated list of all cram files to run igv on')
parser.add_argument("-samples", "--samples", type=str, help="List of all samples to run igv on")
parser.add_argument("-crams", "--crams", type=str, help="File of all cram files to run igv on")
parser.add_argument("-o", "--outdir", type=str, help="output folder")
parser.add_argument("-b", "--buff", type=str, help="length of buffer to add around variants", default=500)
parser.add_argument("-c", "--chromosome", type=str, help="name of chromosome to make igv on", default="all")
parser.add_argument("-i", "--igvfile", type=str, help="name of chromosome to make igv on", default="all")
parser.add_argument("-bam", "--bamfiscript", type=str, help="name of chromosome to make igv on", default="all")
parser.add_argument("-m", "--igvmaxwindow", type=str, help="max length of SV to appear in IGV", default=10e10)
parser.add_argument(
    "--long_read",
    dest="long_read",
    action="store_true",
    help="render for long reads: omit viewaspairs (paired-end display is meaningless for long reads)",
)
parser.add_argument(
    "--genome", type=str, help="reference genome fasta to load in IGV (needed for BAM/arbitrary references)", default=None
)
parser.add_argument(
    "--status_labels",
    dest="status_labels",
    action="store_true",
    help="label/order each BAM track by affected+carrier status (long-read); "
    "requires the varfile 6th column to hold the per-variant carrier sample list",
)
parser.add_argument(
    "--annotation_beds",
    nargs="*",
    default=[],
    help="optional reference BED files to load as IGV feature tracks (e.g. segdups, N-gaps)",
)
parser.add_argument("--annotation_names", nargs="*", default=[], help="track labels for --annotation_beds, in the same order")
parser.add_argument(
    "--genes", type=str, default=None, help="optional gene annotation file (gtf/gff3/bed/refGene) to load as an IGV gene track"
)

args = parser.parse_args()


buff = int(args.buff)
# fasta = args.fasta
varfile = args.varfile
pedigree = args.ped
fam_id = args.fam_id
igv_max_window = args.igvmaxwindow
long_read = args.long_read
genome = args.genome
status_labels = args.status_labels
annotation_beds = args.annotation_beds
annotation_names = args.annotation_names
genes = args.genes


def status_label(sample, affected, carriers):
    aff, car = sample in affected, sample in carriers
    if aff and car:
        return "AFFECTED_CARRIER"
    if car:
        return "CARRIER"
    if aff:
        return "AFFECTED"
    return "unaffected"


def status_rank(sample, affected, carriers):
    aff, car = sample in affected, sample in carriers
    return 0 if (aff and car) else (1 if car else (2 if aff else 3))


def sample_from_bam(path):
    b = os.path.basename(path)
    return b[:-4] if b.endswith(".bam") else os.path.splitext(b)[0]


outstring = os.path.basename(varfile)[0:-4]
bamdir = "pe_bam"
outdir = args.outdir
igvfile = args.igvfile
bamfiscript = args.bamfiscript
###################################

# crams = args.crams
chromosome = args.chromosome
# nested_repeats = args.nestedrepeats
# simple_repeats = args.simplerepeats
# empty_track = args.emptytrack


def ped_info_readin(ped_file):
    out = {}
    fin = open(ped_file)
    for line in fin:
        pin = line.strip().split()
        if pin[1] not in out.keys():
            out[pin[1]] = [pin[1]]
        if not (pin[2]) == 0:
            out[pin[1]].append(pin[2])
        if not (pin[3]) == 0:
            out[pin[1]].append(pin[3])
    fin.close()
    return out


def cram_info_readin(cram_file):
    out = {}
    fin = open(cram_file)
    for line in fin:
        pin = line.strip().split()
        if pin[0] not in out.keys():
            out[pin[0]] = pin[1:]
    fin.close()
    return out


# ped_info = ped_info_readin(args.ped)
# cram_info = cram_info_readin(args.cram_list)

# If file inputs
cram_colnames = colnames = ["cram"]
cram = pd.read_csv(args.crams, sep="\t", names=cram_colnames, header=None).replace(np.nan, "", regex=True)
cram_list = cram["cram"].tolist()
# cram_list = [c.replace('gs://', '/cromwell_root/') for c in cram_list_]

# sample_colnames = colnames=[ 'samples']
# sample = pd.read_csv(args.samples, sep='\t', names= sample_colnames, header=None).replace(np.nan, '', regex=True)
# samples_list = sample['samples'].tolist()

samples_list = args.samples.split(",")
# cram_list=args.crams.split(',')
mydict = {key: value for key, value in zip(samples_list, cram_list)}
ped = pd.read_csv(pedigree, sep="\t", header=None).replace(np.nan, "", regex=True)
ped = ped.iloc[:, :6]
ped.columns = ["FamilyID", "IndividualID", "FatherID", "MotherID", "Sex", "Affected"]
ped["FatherID"] = ped["FatherID"].astype(str)
ped["MotherID"] = ped["MotherID"].astype(str)
ped.Affected = pd.to_numeric(ped.Affected)
affected_samples = set(ped.loc[ped["Affected"] == 2, "IndividualID"].astype(str))
# status mode reorders per-variant below; skip this whole-family reorder
for sample_id in samples_list if not status_labels else []:
    if ped.loc[(ped["IndividualID"] == sample_id)]["Affected"].iloc[0] == 2:
        if (ped.loc[(ped["IndividualID"] == sample_id)]["MotherID"].iloc[0] != "0") | (
            ped.loc[(ped["IndividualID"] == sample_id)]["FatherID"].iloc[0] != "0"
        ):
            print(sample_id)
            proband_cram_file = mydict[sample_id]
            cram_list.remove(proband_cram_file)
            cram_list.insert(0, proband_cram_file)
        else:
            affected_cram_file = mydict[sample_id]
            cram_list.remove(affected_cram_file)
            cram_list.insert(1, affected_cram_file)
print(cram_list)


def build_tracks(carriers):
    """Ordered [(load_arg, squish_name)] for a carrier set, creating symlinks as a side
    effect in status mode. Track label + order depend only on affected status (family-level,
    constant) and this carrier set, so the returned list is identical for any two variants
    that share a carrier set -- which is what lets us load a group's tracks exactly once."""
    tracks = []
    if status_labels:
        ordered = sorted(
            cram_list, key=lambda c: (status_rank(sample_from_bam(c), affected_samples, carriers), cram_list.index(c))
        )
        for c in ordered:
            s = sample_from_bam(c)
            link = status_label(s, affected_samples, carriers) + "." + s + ".bam"
            if not os.path.lexists(link):
                os.symlink(os.path.abspath(c), link)
            for idx in (c + ".bai", (c[:-4] if c.endswith(".bam") else c) + ".bai"):
                if os.path.exists(idx):
                    if not os.path.lexists(link + ".bai"):
                        os.symlink(os.path.abspath(idx), link + ".bai")
                    break
            tracks.append((link, link))
    else:
        for c in cram_list:
            tracks.append((c, os.path.basename(c)))
    return tracks


def carrier_key(dat):
    """Group key: the carrier set (status mode) determines the track config. Non-status mode
    uses one fixed config for the whole family, so everything collapses to a single group."""
    if not status_labels:
        return frozenset()
    return frozenset(dat[5].split(",")) if len(dat) > 5 and dat[5] else frozenset()


# Read all variants, then bucket them by carrier set. Each variant still produces its own
# PNG(s) named by ID, so regrouping the render order does not change any output file.
groups = {}
group_order = []
with open(varfile, "r") as f:
    for line in f:
        dat = line.rstrip().split("\t")
        Chr = dat[0]
        if not chromosome == "all" and not Chr == chromosome:
            continue
        k = carrier_key(dat)
        if k not in groups:
            groups[k] = []
            group_order.append(k)
        groups[k].append(dat)

with open(bamfiscript, "w") as h:
    h.write("#!/bin/bash\n")
    h.write("set -e\n")
    h.write("mkdir -p {}\n".format(bamdir))
    h.write("mkdir -p {}\n".format(outdir))
    with open(igvfile, "w") as g:
        for k in group_order:
            carriers = set(k)
            # new session per carrier group; load its BAM tracks + the (shared) annotation and
            # gene tracks exactly ONCE, then render every variant in the group without reloading.
            g.write("new\n")
            if genome:
                g.write("genome " + genome + "\n")
            # keep annotation/gene feature tracks at a standard (collapsed) height
            g.write("preference EXPAND_FEATURE_TRACKS false\n")
            # cap read depth (~100x) so deep long-read pileups don't exhaust the JVM heap
            g.write("preference SAM.DOWNSAMPLE_READS true\n")
            g.write("preference SAM.SAMPLING_READ_LIMIT 100\n")

            # track the loaded read files so we can re-squish just their alignment tracks
            # after a blanket collapse (see snapshot block)
            loaded_bams = []
            for load_arg, squish_name in build_tracks(carriers):
                g.write("load " + load_arg + "\n")
                loaded_bams.append(squish_name)

            # reference annotation tracks (segdups, N-gaps, ...) so IGV-only variants
            # still show them; symlink to the given name for a clean track label
            for i, abed in enumerate(annotation_beds):
                ext = ".bed.gz" if abed.endswith(".bed.gz") else (os.path.splitext(abed)[1] or ".bed")
                if annotation_names and i < len(annotation_names):
                    alink = annotation_names[i] + ext
                    if not os.path.lexists(alink):
                        os.symlink(os.path.abspath(abed), alink)
                    g.write("load " + alink + "\n")
                else:
                    g.write("load " + os.path.abspath(abed) + "\n")
            # gene annotation track
            if genes:
                g.write("load " + os.path.abspath(genes) + "\n")

            for dat in groups[k]:
                Chr = dat[0]
                Start = int(dat[1])
                End = int(dat[2])
                ID = dat[3]
                Length = End - Start

                Length_total = int(Length + (Length) * 1.5)

                # long-read reads are noisy at the per-read indel level: hide sub-5bp
                # indels for larger variants, but show them when the variant itself is small
                if long_read:
                    if Length > 50:
                        g.write("preference SAM.HIDE_SMALL_INDEL true\n")
                        g.write("preference SAM.SMALL_INDEL_BP_THRESHOLD 5\n")
                        # also stop labeling sub-5bp insertions (SAM.LARGE_INSERTIONS_THRESOLD
                        # is IGV's actual, misspelled key) so hidden indels don't leave markers
                        g.write("preference SAM.LARGE_INSERTIONS_THRESOLD 5\n")
                    else:
                        g.write("preference SAM.HIDE_SMALL_INDEL false\n")
                        g.write("preference SAM.LARGE_INSERTIONS_THRESOLD 1\n")

                if Length_total < int(igv_max_window):
                    if Length_total < 1000:
                        Start_Buff = int(Start - 500)
                        End_Buff = int(End + 500)
                    else:
                        Start_Buff = int(Start - (Length * 0.25))
                        End_Buff = int(End + (Length * 0.25))
                    g.write("preference SAM.SHOW_CENTER_LINE false\n")
                    g.write("goto " + Chr + ":" + str(Start_Buff) + "-" + str(End_Buff) + "\n")
                    g.write("region " + Chr + ":" + str(Start) + "-" + str(End) + "\n")
                    if Length <= 50:
                        g.write("sort base\n")
                    if not long_read:
                        g.write("viewaspairs\n")
                    g.write("collapse\n")
                    for _bam in loaded_bams:
                        g.write("squish " + _bam + " Alignments\n")
                    g.write("snapshotDirectory " + outdir + "\n")
                    g.write("snapshot " + fam_id + "_" + ID + ".png\n")
                else:
                    # split view: one snapshot per breakpoint, each centered on the junction
                    # with a single center line (not a double-edged region-of-interest marker)
                    g.write("preference SAM.SHOW_CENTER_LINE true\n")
                    g.write("goto " + Chr + ":" + str(Start - buff) + "-" + str(Start + buff) + "\n")
                    if not long_read:
                        g.write("viewaspairs\n")
                    g.write("collapse\n")
                    for _bam in loaded_bams:
                        g.write("squish " + _bam + " Alignments\n")
                    g.write("snapshotDirectory " + outdir + "\n")
                    g.write("snapshot " + fam_id + "_" + ID + ".left.png\n")
                    g.write("goto " + Chr + ":" + str(End - buff) + "-" + str(End + buff) + "\n")
                    if not long_read:
                        g.write("viewaspairs\n")
                    g.write("collapse\n")
                    for _bam in loaded_bams:
                        g.write("squish " + _bam + " Alignments\n")
                    g.write("snapshotDirectory " + outdir + "\n")
                    g.write("snapshot " + fam_id + "_" + ID + ".right.png\n")
        g.write("exit\n")
