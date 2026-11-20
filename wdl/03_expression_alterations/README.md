# 03_expression_alterations

Workflows for analyses of expression alterations caused by long-read variants.

- `PhasedHetPerGene.wdl`: phased heterozygous coding (LoF, missense, synonymous) and UTR (5', 3') variants per gene and
  sample, plus a per-gene summary (samples with >= 1 / 3 / 5 of each category). Used to estimate how easily Iso-Seq /
  Kinnex reads of each gene can be phased (a read covering a phased het in exons / UTRs can be assigned to a haplotype).
