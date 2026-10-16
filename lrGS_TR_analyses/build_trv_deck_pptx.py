#!/usr/bin/env python3
"""Build the TRV / STR summary deck as a PowerPoint file from the figures made by make_trv_deck_figures.py.

Slide text mirrors TRV_STR_deck.html; numbers come from that run's summary.json (quoted, not recomputed).
Figures are embedded as images; full slide notes go into the speaker notes.

Usage:
  build_trv_deck_pptx.py --figures-dir deck/figures --out deck/TRV_STR_deck.pptx
"""
import argparse
import os

from PIL import Image
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.oxml.ns import qn
from pptx.util import Inches, Pt

NAVY = RGBColor(0x14, 0x2A, 0x4A)
INK = RGBColor(0x1D, 0x21, 0x28)
INK2 = RGBColor(0x54, 0x5B, 0x66)
ACCENT = RGBColor(0x1B, 0x5A, 0xA6)
WARN = RGBColor(0xB5, 0x54, 0x1F)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
ICE = RGBColor(0xC9, 0xDB, 0xF2)
RULE = RGBColor(0xD8, 0xDC, 0xE1)
FONT = "Arial"
SLIDE_W, SLIDE_H = 13.333, 7.5
MARGIN = 0.5


def set_bg(slide, color):
    fill = slide.background.fill
    fill.solid()
    fill.fore_color.rgb = color


def textbox(slide, x, y, w, h, paras, size=14, color=INK, bullets=False, bold_first=False, name=None):
    """paras: list of str or (str, dict) with optional keys bold, color, size, italic."""
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    if name:
        tb.name = name
    tf = tb.text_frame
    tf.word_wrap = True
    for side in ("left", "right", "top", "bottom"):
        setattr(tf, f"margin_{side}", Inches(0))
    for i, p in enumerate(paras):
        text, opt = (p, {}) if isinstance(p, str) else p
        para = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        if bullets:
            pPr = para._p.get_or_add_pPr()
            pPr.set("marL", str(Inches(0.22)))
            pPr.set("indent", str(-Inches(0.22)))
            bu_clr = pPr.makeelement(qn("a:buClr"), {})
            srgb = bu_clr.makeelement(qn("a:srgbClr"), {"val": "1B5AA6"})
            bu_clr.append(srgb)
            pPr.append(bu_clr)
            bu = pPr.makeelement(qn("a:buChar"), {"char": "•"})
            pPr.append(bu)
        para.space_after = Pt(opt.get("space_after", 6 if bullets else 2))
        run = para.add_run()
        run.text = text
        f = run.font
        f.name = FONT
        f.size = Pt(opt.get("size", size))
        f.bold = opt.get("bold", bold_first and i == 0)
        f.italic = opt.get("italic", False)
        f.color.rgb = opt.get("color", color)
    return tb


def add_picture(slide, path, x, y, w=None, h=None):
    iw, ih = Image.open(path).size
    if w is not None and h is not None:
        if w / h > iw / ih:
            w = h * iw / ih
        else:
            h = w * ih / iw
    elif w is not None:
        h = w * ih / iw
    else:
        w = h * iw / ih
    pic = slide.shapes.add_picture(path, Inches(x), Inches(y), Inches(w), Inches(h))
    pic.line.color.rgb = RULE
    pic.line.width = Pt(0.75)
    return w, h


def header(slide, eyebrow, title, num, total):
    textbox(slide, MARGIN, 0.32, SLIDE_W - 2 * MARGIN, 0.3, [(eyebrow.upper(), {"color": ACCENT, "bold": True})], size=11,
            name="Eyebrow")
    textbox(slide, MARGIN, 0.62, SLIDE_W - 2 * MARGIN, 0.8, [(title, {"bold": True})], size=24, color=INK, name="Title")
    textbox(slide, SLIDE_W - MARGIN - 1.2, SLIDE_H - 0.38, 1.2, 0.25, [(f"{num} / {total}", {"color": INK2})], size=10,
            name="Slide number").text_frame.paragraphs[0].alignment = 3


def notes(slide, text):
    slide.notes_slide.notes_text_frame.text = text


def table(slide, x, y, w, rows, col_w, size=11, header_color=INK2):
    shape = slide.shapes.add_table(len(rows), len(rows[0]), Inches(x), Inches(y), Inches(w), Inches(0.3 * len(rows)))
    tbl = shape.table
    for j, cw in enumerate(col_w):
        tbl.columns[j].width = Inches(cw)
    for i, row in enumerate(rows):
        for j, val in enumerate(row):
            cell = tbl.cell(i, j)
            cell.fill.background()
            cell.margin_left = cell.margin_right = Inches(0.06)
            cell.margin_top = cell.margin_bottom = Inches(0.03)
            tf = cell.text_frame
            tf.text = str(val)
            p = tf.paragraphs[0]
            p.alignment = 1 if j == 0 else 3
            f = p.runs[0].font
            f.name = FONT
            f.size = Pt(size)
            f.bold = i == 0
            f.color.rgb = header_color if i == 0 else INK
    tbl.first_row = False
    tbl.horz_banding = False
    return shape


MOTIF_SLIDES = [
    (1, "STR motif 1 bp · homopolymers · TRExplorer, hprc_hgsvc", "Coding homopolymers are 5–6× less "
                                                                  "variable than non-coding ones",
     ["1.09M loci, only 260 coding (0.02%): homopolymers are strongly depleted from coding sequence.",
      "Non-ref AF: coding 0.047 vs UTR 0.23, intronic 0.25, intergenic 0.27.",
      "Any length change: 84% of coding loci vs 97% elsewhere; every homopolymer change is a frameshift.",
      "Obs/exp 0.66–0.91 in every decile except 9 (9–46 loci per decile). Decile 9 is one locus, ZNF814 (104 obs vs 2 exp)."]),
    (2, "STR motif 2 bp · TRExplorer, hprc_hgsvc", "Coding dinucleotide repeats rarely change length",
     ["292k loci; 274 coding, 3.4k UTR.",
      "Non-ref AF: coding 0.11 vs 0.23–0.26 non-coding (2.4× lower).",
      "Change ≥ 2 bp: 24% of coding loci vs 73% elsewhere; coding median max |Δ| is 0 bp.",
      "Obs/exp swings 0.6–1.6 with 11–35 loci per decile; no trend resolvable. Outliers: MUC5AC, PRDM9."]),
    (3, "STR motif 3 bp · TRExplorer, hprc_hgsvc", "Coding trinucleotide repeats tolerate in-frame "
                                                   "change but are constrained at 4–8 copies",
     ["293k loci; 5,039 coding (1.7%), the largest coding set of any motif.",
      "Non-ref AF: coding 0.091 vs UTR 0.126, non-coding 0.16–0.17; widest gap at 4–8 copies (coding 25–60% of intergenic).",
      "Change ≥ 3 bp: 40% coding vs 37–38% non-coding, so in-frame codon changes are tolerated; "
      "frameshift alleles are 20–100× depleted.",
      "Obs/exp: deciles 1–3 0.81–0.89, deciles 4–8 0.82–1.18. Decile 9 spike (2.0) is ABCA7 and FTCD.",
      "Long repeats (≥ 9 copies) are enriched in constrained genes (ZFHX3, ARID1B, HOXA13) and remain polymorphic."]),
    (4, "STR motif 4 bp · TRExplorer, hprc_hgsvc", "Coding tetranucleotide repeats are rare and rarely change length",
     ["189k loci; 110 coding.",
      "Non-ref AF: coding 0.081 vs 0.15 UTR and 0.20–0.21 non-coding.",
      "Change ≥ 4 bp (one motif, shifts the frame): 14% of coding loci vs 62% elsewhere.",
      "Obs/exp: 2–17 loci per decile, too few to read a trend."]),
    (5, "STR motif 5 bp · TRExplorer, hprc_hgsvc", "Only 24 coding pentanucleotide loci",
     ["68k loci; 24 coding, 808 UTR.",
      "Non-ref AF: coding 0.10 vs 0.16 UTR and 0.19–0.20 non-coding.",
      "Change ≥ 5 bp: 25% of coding loci vs 54–57% elsewhere.",
      "Obs/exp not interpretable with 1–6 loci per decile."]),
    (6, "STR motif 6 bp · TRExplorer, hprc_hgsvc", "Coding hexanucleotide repeats change by whole motifs, which keep the frame",
     ["27k loci; 249 coding.",
      "Non-ref AF: coding 0.089 vs 0.14 UTR and 0.17–0.18 non-coding.",
      "Change ≥ 6 bp (two codons, in frame): 77% of coding loci vs 51–52% non-coding; coding hexamers vary by whole units.",
      "Obs/exp flat at 0.73–1.11 across deciles (14–43 loci each)."]),
]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--figures-dir", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()
    fig = lambda name: os.path.join(args.figures_dir, f"{name}.png")  # noqa: E731

    prs = Presentation()
    prs.slide_width, prs.slide_height = Inches(SLIDE_W), Inches(SLIDE_H)
    blank = prs.slide_layouts[6]
    total = 13
    n = 0

    # 1 title
    s = prs.slides.add_slide(blank)
    n += 1
    set_bg(s, NAVY)
    textbox(s, 0.8, 1.2, 11.5, 0.4, [("GNOMAD LONG-READ · TANDEM REPEATS", {"bold": True, "color": ICE})], size=13)
    textbox(s, 0.8, 1.7, 11.5, 1.5, [("Tandem repeat variation in two long-read cohorts", {"bold": True, "color": WHITE})],
            size=40)
    textbox(s, 0.8, 3.25, 11.0, 0.8, [("Catalog, data quality, mutability by repeat size and motif, and signs of selection on "
                                       "coding STRs", {"color": ICE})], size=18)
    textbox(s, 0.8, 4.6, 11.5, 1.6, [
        ("Cohorts: hprc_hgsvc (292 samples, AN 584) · AoU_I (1,027 samples)", {"color": WHITE}),
        ("Input: gnomAD_LR.{cohort}.vep_parsed.annotated.bed.gz", {"color": WHITE}),
        ("GRCh38 · genes GENCODE v39 (gnomAD-SV r3) · constraint gnomAD v4.1.1 (LOEUF)", {"color": WHITE}),
        ("Sections: 1 overview · 2 TRV catalog · 3 AoU ±1 bp alleles · 4 mutation rate · 5 STRs by motif and constraint",
         {"color": ICE})], size=14)
    notes(s, "Two cohorts: hprc_hgsvc (292 samples) and AoU_I (1,027 samples). Sections 4-5 use hprc_hgsvc only because "
             "AoU_I repeat alleles carry ±1 bp errors (section 3).")

    # 2 overview
    s = prs.slides.add_slide(blank)
    n += 1
    header(s, "1 · Variant overview", "AoU_I has 1.4× more SNVs and 1.9× more TRV sites than hprc_hgsvc", n, total)
    w, h = add_picture(s, fig("s1_variant_overview"), MARGIN, 1.6, w=7.4)
    caption = "PASS sites per variant class; indel vs SV split at |allele_length| 50 bp."
    textbox(s, MARGIN, 1.6 + h + 0.1, 7.4, 0.3, [(caption, {"color": INK2})], size=10)
    rows = [["PASS sites", "hprc_hgsvc", "AoU_I"],
            ["SNV", "42,165,343", "59,535,077"], ["Deletion <50 bp", "4,634,840", "6,343,851"],
            ["Insertion <50 bp", "4,923,508", "6,145,024"], ["Tandem repeat (TRV)", "2,446,683", "4,598,792"],
            ["Insertion SV ≥50 bp", "197,236", "540,382"], ["Deletion SV ≥50 bp", "105,447", "178,784"],
            ["Duplication", "6,132", "9,512"]]
    table(s, 8.3, 1.6, 4.5, rows, [2.1, 1.2, 1.2], size=11)
    textbox(s, 8.3, 4.35, 4.5, 2.6, [
        "SNVs and short indels grow with cohort size (3.5× more samples in AoU_I).",
        "TRVs nearly double and insertion SVs almost triple in AoU_I (ins:del SV 3.0 vs 1.9): repeat-rich classes that need the "
        "QC on slide 4.",
        "Almost all TRV sites are PASS (361 and 382 non-PASS)."], size=12, bullets=True)
    notes(s, "Counts from count_variant_overview.sh over the full annotated bed files. TRV non-PASS sites are all "
             "LOW_COVERAGE_REGION.")

    # 3 TRV catalog
    s = prs.slides.add_slide(blank)
    n += 1
    header(s, "2 · TRV catalog", "Most TRVs are short, homopolymer or dinucleotide, and non-coding", n, total)
    w, h = add_picture(s, fig("s2_trv_counts"), MARGIN, 1.55, w=SLIDE_W - 2 * MARGIN)
    textbox(s, MARGIN, 1.55 + h + 0.08, 12.3, 0.3, [("Left: genic context (coding = TR span inside one "
                                                     "CDS block). Middle: repeat "
                                                     "span from TRID. Right: shortest TRID motif.", {"color": INK2})], size=10)
    y = 1.55 + h + 0.5
    textbox(s, MARGIN, y, 6.0, SLIDE_H - y - 0.5, [
        "2.45M TRV sites in hprc_hgsvc, 4.60M in AoU_I; 85–88% from the TRExplorer catalog, the rest Vamos.",
        "Genic context is the same in both: intergenic 59%, intronic 39%, UTR 1.3–1.4%."], size=13, bullets=True)
    textbox(s, 6.85, y, 6.0, SLIDE_H - y - 0.5, [
        "Coding is 0.4% in hprc_hgsvc but 0.9% in AoU_I (8.7k vs 40k sites inside CDS).",
        "Median span 13 bp (hprc) / 11 bp (AoU). AoU_I has 4.3× more 3-bp-motif sites but only 1.2× more homopolymers: "
        "the ±1 bp issue."], size=13, bullets=True)
    notes(s, "Homopolymers dominate (1.23M in hprc_hgsvc). 90% of TRVs are under 43 bp.")

    # 4 pm1
    s = prs.slides.add_slide(blank)
    n += 1
    header(s, "3 · Data quality", "AoU_I repeat alleles carry scattered ±1 bp errors", n, total)
    w, h = add_picture(s, fig("s3_pm1_issue"), MARGIN, 1.55, w=SLIDE_W - 2 * MARGIN)
    y = 1.55 + h + 0.3
    textbox(s, MARGIN, y, 7.2, SLIDE_H - y - 0.5, [
        "hprc_hgsvc: multiple-of-3 share of coding length changes rises 50% → 88% with AC; ±1 bp falls to 8%.",
        "AoU_I: ±1 bp share rises 46% → 88% and multiple-of-3 falls below chance (4–13%). 99.8% of coding loci carry a ±1 bp "
        "allele (≈1,200 per sample vs ≈26).",
        "Edits are scattered inside the repeat as separate low-AC alleles, inflating ALT alleles per locus 3–5×.",
        ("Consequence: sections 4–5 use hprc_hgsvc only.", {"bold": True, "color": WARN})], size=12, bullets=True)
    textbox(s, 8.1, y, 4.7, 0.3, [("GBX1 (chr7-151167152), REF CGGCAGTGGCGGCGGCGGCGGC", {"bold": True})], size=11)
    rows = [["Cohort", "AC", "Δ", "Allele"], ["AoU_I", "729", "0", "CCGCAGTGG…"], ["AoU_I", "7", "−1", "CG[-]CAGTGG…"],
            ["AoU_I", "4", "+1", "…GGCGG[G]C"], ["AoU_I", "3", "−1", "…GGCG[-]CGGC"], ["hprc", "185", "0", "CCGCAGTGG…"],
            ["hprc", "1", "+9", "+3 codons"]]
    table(s, 8.1, y + 0.35, 4.7, rows, [0.9, 0.7, 0.6, 2.5], size=10)
    notes(s, "In coding repeats, length changes that are a multiple of 3 should rise with allele "
             "frequency. They do in hprc_hgsvc; "
             "in AoU_I the ±1 bp share rises instead. 385,720 distinct ±1 bp alleles across 39,750 "
             "coding loci inside CDS. Only 4% "
             "of edits sit in homopolymer runs >= 4. Common alleles agree between cohorts; the noise is a cloud of 1-bp variants "
             "around each real allele. GBX1: AoU_I has 42 ALT alleles (27 are ±1 bp); hprc_hgsvc has 2. Next: per-sample ±1 bp "
             "counts vs platform, coverage and read quality; collapse alleles to motif units.")

    # 5 mutation rate
    s = prs.slides.add_slide(blank)
    n += 1
    header(s, "4 · Mutation rate vs size and motif · hprc_hgsvc",
           "Variability rises with repeat length; at fixed copy number, longer motifs vary more", n, total)
    w, h = add_picture(s, fig("s4_mutation_rate"), MARGIN, 1.5, h=5.55)
    x = MARGIN + w + 0.3
    textbox(s, x, 1.5, SLIDE_W - MARGIN - x, 5.5, [
        "STR vs size: ALT alleles 2–3 at 10 bp → 30–40 at 50–75 bp; non-ref AF 0.15 → 0.7–0.95.",
        "Fixed bp: shorter motifs vary more (30–40 bp: 32 / 18 / 7–9 alleles for 1 / 2 / 3–6 bp).",
        "Fixed copies: longer motifs vary more (10–14 copies: AF 0.52 → 0.76 from 2 to 5 bp).",
        "VNTR: AF 0.1–0.2 below 100 bp, 0.5–0.8 above 1 kb; motif class barely matters at fixed size.",
        ("Catalogs differ several-fold at equal size, so STR panels use TRExplorer only.", {"color": INK2})],
        size=12, bullets=True)
    notes(s, "Rows: non-ref AF (sum AC / AN), distinct ALT alleles (AC > 0), max |ALT − REF| (median, "
             "IQR). PASS, single-component "
             "loci, >= 2 copies. STR panels: TRExplorer (1.96M loci). VNTR panels: all catalogs (154k "
             "loci, 87% Vamos). Homopolymers "
             "sit below dinucleotides at equal copies; low-copy bins mostly reflect each catalog's minimum locus size.")

    # 6 model
    s = prs.slides.add_slide(blank)
    n += 1
    header(s, "5 · Constraint model · hprc_hgsvc coding TRVs",
           "An expected-variability model finds mild depletion in the most LoF-intolerant genes", n, total)
    w, h = add_picture(s, fig("s5_model"), MARGIN, 1.55, w=SLIDE_W - 2 * MARGIN)
    y = 1.55 + h + 0.3
    textbox(s, MARGIN, y, 6.0, SLIDE_H - y - 0.5, [
        ("Model", {"bold": True, "space_after": 4}),
        "Poisson GLM on 8,502 coding TRVs (inside CDS, PASS, single-component, ≥ 2 copies); response = distinct ALT alleles.",
        "Predictors: log size + log size², motif class × catalog, motif class × log size, GC, purity.",
        "Deviance explained 0.41; purity strongest (β 0.97, p 5×10⁻⁵⁸); GC n.s. (p 0.37)."], size=12, bullets=True)
    textbox(s, 6.85, y, 6.0, SLIDE_H - y - 0.5, [
        ("Obs/exp by LOEUF decile (10 equal bins of 17,418 MANE genes)", {"bold": True, "space_after": 4}),
        "Deciles 1–3: 0.90, 0.80, 0.85 (upper CI ≤ 1.03); deciles 4–8: 0.92–1.14.",
        "Decile 9 peak (1.34) is driven by VNTR-like alleles in a few genes (ABCA7: 540 obs vs 7.4 exp).",
        ("Trained on the same coding loci: obs/exp is relative to the average coding TRV.", {"color": WARN})],
        size=12, bullets=True)
    notes(s, "Calibration in 20 quantile bins. Per-gene boxes show genes with >= 3 expected ALT alleles. "
             "Pooled obs/exp = sum obs / "
             "sum exp per decile with 500 gene bootstraps. Multi-gene loci count toward each gene. Most genes have 1–4 loci.")

    # 7-12 motifs
    for motif, eyebrow, title, bullets in MOTIF_SLIDES:
        s = prs.slides.add_slide(blank)
        n += 1
        header(s, f"5 · {eyebrow}", title, n, total)
        w, h = add_picture(s, fig(f"s5_motif{motif}"), MARGIN, 1.6, h=4.7)
        x = MARGIN + w + 0.3
        textbox(s, x, 1.6, SLIDE_W - MARGIN - x, 5.4, bullets, size=12, bullets=True)
        notes(s, "Panels: a loci by genic context; b non-ref AF and c distinct ALT alleles vs copy "
                 "number by context; d max |ALT − REF| "
                 "ECDF; e per-gene log2 obs/exp (genes with >= 3 expected) and f pooled obs/exp by LOEUF "
                 "decile (numbers = coding "
                 "loci per decile). Coding = locus inside one CDS block. " + " ".join(bullets))

    # 14 summary
    s = prs.slides.add_slide(blank)
    n += 1
    set_bg(s, NAVY)
    textbox(s, MARGIN, 0.45, 12.3, 0.3, [("SUMMARY", {"bold": True, "color": ICE})], size=11)
    textbox(s, MARGIN, 0.75, 12.3, 0.7, [("What we learned and what to check next", {"bold": True, "color": WHITE})], size=28)
    textbox(s, MARGIN, 1.7, 6.0, 0.4, [("Findings", {"bold": True, "color": ICE})], size=16)
    textbox(s, MARGIN, 2.15, 6.0, 4.9, [
        ("TRVs: 2.45M sites (hprc_hgsvc) and 4.60M (AoU_I), mostly short non-coding homopolymers and dinucleotides.",
         {"color": WHITE}),
        ("AoU_I repeat alleles carry ±1 bp errors that inflate ALT alleles, 3-bp sites and coding frameshifts.",
         {"color": WHITE}),
        ("Variability scales with repeat length and copy number; at equal copies longer motifs vary more.", {"color": WHITE}),
        ("Coding STRs have 1.8–5.8× lower non-ref AF than intergenic for every motif 1–6 bp; frame-shifting changes are the "
         "most depleted.", {"color": WHITE}),
        ("LoF-intolerant genes (LOEUF deciles 1–3) show 10–20% fewer ALT alleles than expected; weak and noisy per gene.",
         {"color": WHITE})], size=13, bullets=True)
    textbox(s, 6.85, 1.7, 6.0, 0.4, [("Next steps", {"bold": True, "color": ICE})], size=16)
    textbox(s, 6.85, 2.15, 6.0, 3.6, [
        ("Per-sample ±1 bp counts in AoU_I vs platform, coverage and read quality; collapse alleles to motif units.",
         {"color": WHITE}),
        ("Model frameshift and in-frame alleles separately; AF-weighted or AC ≥ 2 outcomes.", {"color": WHITE}),
        ("Train the expectation on non-coding loci matched for size, motif and purity.", {"color": WHITE}),
        ("Repeat on AoU_I after cleanup to gain power for 4–6 bp motifs.", {"color": WHITE})], size=13, bullets=True)
    textbox(s, 6.85, 6.2, 6.0, 0.6, [("github.com/talkowski-lab/lr-pipeline · xz_analyses · lrGS_TR_analyses/",
                                      {"color": ICE})], size=10)
    notes(s, "Scripts: count_variant_overview.sh, build_trv_locus_table.py, make_trv_deck_figures.py, build_trv_deck_pptx.py.")

    assert n == total, n
    prs.save(args.out)
    print(f"wrote {args.out} ({n} slides)")


if __name__ == "__main__":
    main()
