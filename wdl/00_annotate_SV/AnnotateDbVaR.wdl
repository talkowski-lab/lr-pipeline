version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateDbVaR {
    input {
        File vcf
        File vcf_idx
        File dbvar_vcf
        File dbvar_vcf_idx
        Array[String] contigs
        String prefix

        Int min_length

        Int? records_per_shard

        Int del_breakpoint_window = 500
        Float del_reciprocal_overlap = 0.7
        Float del_size_similarity = 0.7

        Int dup_breakpoint_window = 500
        Float dup_reciprocal_overlap = 0.7
        Float dup_size_similarity = 0.7

        Int ins_breakpoint_window = 200
        Float ins_reciprocal_overlap = 0.0
        Float ins_size_similarity = 0.5

        String utils_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_subset_dbvar_vcf
        RuntimeAttr? runtime_attr_convert_symbolic
        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_annotate
        RuntimeAttr? runtime_attr_concat_shards
        RuntimeAttr? runtime_attr_concat
    }

    Boolean single_contig = length(contigs) == 1

    scatter (contig in contigs) {
        call Helpers.SubsetVcfByArgs {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                include_args = "abs(INFO/allele_length) >= ~{min_length} && (INFO/allele_type != \"trv\")",
                extra_args = if single_contig then "" else "--regions " + contig,
                prefix = "~{prefix}.~{contig}.filtered",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }

        call Helpers.ConvertToSymbolic {
            input:
                vcf = SubsetVcfByArgs.subset_vcf,
                vcf_idx = SubsetVcfByArgs.subset_vcf_idx,
                move_all_dups = false,
                prefix = "~{prefix}.~{contig}.symbolic",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_convert_symbolic
        }

        if (!single_contig) {
            call Helpers.SubsetVcfToContig as SubsetDbVaRVcf {
                input:
                    vcf = dbvar_vcf,
                    vcf_idx = dbvar_vcf_idx,
                    contig = contig,
                    prefix = "~{prefix}.~{contig}.dbvar",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_dbvar_vcf
            }
        }

        File contig_dbvar_vcf = select_first([SubsetDbVaRVcf.subset_vcf, dbvar_vcf])
        File contig_dbvar_vcf_idx = select_first([SubsetDbVaRVcf.subset_vcf_idx, dbvar_vcf_idx])

        if (defined(records_per_shard)) {
            call Helpers.ShardVcfByRecords {
                input:
                    vcf = ConvertToSymbolic.processed_vcf,
                    vcf_idx = ConvertToSymbolic.processed_vcf_idx,
                    records_per_shard = select_first([records_per_shard]),
                    prefix = "~{prefix}.~{contig}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_shard
            }
        }

        Array[File] shard_array = select_first([ShardVcfByRecords.shards, []])
        Array[File] shard_idx_array = select_first([ShardVcfByRecords.shard_idxs, []])
        Array[File] vcfs_to_process = if length(shard_array) > 0 then shard_array else [ConvertToSymbolic.processed_vcf]
        Array[File] vcf_idxs_to_process = if length(shard_idx_array) > 0 then shard_idx_array else [ConvertToSymbolic.processed_vcf_idx]

        scatter (i in range(length(vcfs_to_process))) {
            call AnnotateDbVaRIds {
                input:
                    vcf = vcfs_to_process[i],
                    vcf_idx = vcf_idxs_to_process[i],
                    original_vcf = SubsetVcfByArgs.subset_vcf,
                    original_vcf_idx = SubsetVcfByArgs.subset_vcf_idx,
                    dbvar_vcf = contig_dbvar_vcf,
                    dbvar_vcf_idx = contig_dbvar_vcf_idx,
                    del_breakpoint_window = del_breakpoint_window,
                    del_reciprocal_overlap = del_reciprocal_overlap,
                    del_size_similarity = del_size_similarity,
                    dup_breakpoint_window = dup_breakpoint_window,
                    dup_reciprocal_overlap = dup_reciprocal_overlap,
                    dup_size_similarity = dup_size_similarity,
                    ins_breakpoint_window = ins_breakpoint_window,
                    ins_reciprocal_overlap = ins_reciprocal_overlap,
                    ins_size_similarity = ins_size_similarity,
                    prefix = "~{prefix}.~{contig}.shard_~{i}.dbvar_annotated",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_annotate
            }
        }

        if (defined(records_per_shard)) {
            call Helpers.ConcatTsvs as ConcatShards {
                input:
                    tsvs = AnnotateDbVaRIds.annotations_tsv,
                    sort_output = true,
                    prefix = "~{prefix}.~{contig}.dbvar_annotated",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_concat_shards
            }
        }

        File final_annotations_tsv = select_first([ConcatShards.concatenated_tsv, AnnotateDbVaRIds.annotations_tsv[0]])
    }

    if (!single_contig) {
        call Helpers.ConcatTsvs as MergeAnnotations {
            input:
                tsvs = final_annotations_tsv,
                sort_output = false,
                prefix = "~{prefix}.dbvar_annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat
        }
    }

    output {
        File annotations_tsv_dbvar = select_first([MergeAnnotations.concatenated_tsv, final_annotations_tsv[0]])
    }
}

task AnnotateDbVaRIds {
    input {
        File vcf
        File vcf_idx
        File original_vcf
        File original_vcf_idx
        File dbvar_vcf
        File dbvar_vcf_idx
        Int del_breakpoint_window
        Float del_reciprocal_overlap
        Float del_size_similarity
        Int dup_breakpoint_window
        Float dup_reciprocal_overlap
        Float dup_size_similarity
        Int ins_breakpoint_window
        Float ins_reciprocal_overlap
        Float ins_size_similarity
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        if [[ "~{vcf_idx}" != "~{vcf}.tbi" ]]; then
            ln -sf "~{vcf_idx}" "~{vcf}.tbi"
        fi
        if [[ "~{dbvar_vcf_idx}" != "~{dbvar_vcf}.tbi" ]]; then
            ln -sf "~{dbvar_vcf_idx}" "~{dbvar_vcf}.tbi"
        fi

        bcftools query -f '%ID\t%REF\t%ALT\n' ~{original_vcf} > original_ref_alt.tsv

        python3 <<'EOF'
import pysam

DEL_SIZE_SIM = float(~{del_size_similarity})
DEL_REC_OVL = float(~{del_reciprocal_overlap})
DEL_BP_WIN = int(~{del_breakpoint_window})
DUP_SIZE_SIM = float(~{dup_size_similarity})
DUP_REC_OVL = float(~{dup_reciprocal_overlap})
DUP_BP_WIN = int(~{dup_breakpoint_window})
INS_SIZE_SIM = float(~{ins_size_similarity})
INS_REC_OVL = float(~{ins_reciprocal_overlap})
INS_BP_WIN = int(~{ins_breakpoint_window})

PARAMS = {
    'DEL': {'allowed': {'DEL', 'CNV'}, 'size_sim': DEL_SIZE_SIM, 'rec_ovl': DEL_REC_OVL, 'bp_win': DEL_BP_WIN},
    'DUP': {'allowed': {'DUP', 'CNV'}, 'size_sim': DUP_SIZE_SIM, 'rec_ovl': DUP_REC_OVL, 'bp_win': DUP_BP_WIN},
    'INS': {'allowed': {'INS'}, 'size_sim': INS_SIZE_SIM, 'rec_ovl': INS_REC_OVL, 'bp_win': INS_BP_WIN},
}

def parse_query_len(rec):
    raw = rec.info.get('allele_length')
    if isinstance(raw, (list, tuple)):
        raw = raw[0]
    try:
        return abs(int(raw)) if raw is not None else 0
    except (ValueError, TypeError):
        return 0

def parse_dbvar_len(rec):
    raw = rec.info.get('SVLEN')
    if raw is not None:
        if isinstance(raw, (list, tuple)):
            raw = raw[0]
        try:
            return abs(int(raw))
        except (ValueError, TypeError):
            pass
    end = rec.info.get('END')
    if end is not None:
        try:
            return abs(int(end) - (rec.start + 1))
        except (ValueError, TypeError):
            pass
    return 0

def passes_del_dup(qs, qe, ql, cs, ce, cl, size_sim, rec_ovl, bp_win):
    if ql == 0 or cl == 0:
        return False
    if min(ql, cl) / max(ql, cl) < size_sim:
        return False
    overlap = max(0, min(qe, ce) - max(qs, cs))
    if overlap / min(ql, cl) < rec_ovl:
        return False
    if abs(qs - cs) + abs(qe - ce) > bp_win:
        return False
    return True

def passes_ins(qs, qe, ql, cs, ce, cl, size_sim, rec_ovl, bp_win):
    if ql == 0 or cl == 0:
        return False
    if min(ql, cl) / max(ql, cl) < size_sim:
        return False
    overlap = max(0, min(qe, ce) - max(qs, cs))
    if overlap / min(ql, cl) < rec_ovl:
        return False
    if abs(qs - cs) + abs(qe - ce) > bp_win:
        return False
    return True

original_ref_alt = {}
with open("original_ref_alt.tsv") as f:
    for line in f:
        parts = line.rstrip('\n').split('\t')
        if len(parts) >= 3:
            original_ref_alt[parts[0]] = (parts[1], parts[2])

dbvar = pysam.VariantFile("~{dbvar_vcf}")
query = pysam.VariantFile("~{vcf}")

with open("~{prefix}.annotations.unsorted.tsv", 'w') as out:
    for rec in query:
        # Structure search criteria
        svtype = rec.info.get('SVTYPE')
        if svtype not in PARAMS:
            continue
        params = PARAMS[svtype]
        chrom = rec.chrom
        if chrom not in dbvar.header.contigs:
            continue
        
        qs = rec.start
        qe = rec.stop
        ql = parse_query_len(rec)
        var_id = rec.id if rec.id else '.'
        orig_ref, orig_alt = original_ref_alt.get(var_id, (rec.ref if rec.ref else '.', rec.alts[0] if rec.alts else '.'))
        bp_win = params['bp_win']
        size_sim = params['size_sim']
        region_start = max(0, qs - bp_win)
        region_end = qe + bp_win

        # Look for the best (minimum breakpoint distance) canonical match
        best_match = None
        best_bp_dist = float('inf')
        cands = dbvar.fetch(chrom, region_start, region_end)
        for cand in cands:
            ctype = cand.info.get('SVTYPE')
            if ctype is None or not any(a in ctype for a in params['allowed']):
                continue
            cs = cand.start
            ce = cand.stop
            cl = parse_dbvar_len(cand)
            cand_id = cand.id if cand.id else cand.info.get('DBVARID', '.')
            if isinstance(cand_id, (list, tuple)):
                cand_id = cand_id[0]

            if svtype == 'INS':
                ok = passes_ins(qs, qe, ql, cs, ce, cl, size_sim, params['rec_ovl'], bp_win)
            else:
                ok = passes_del_dup(qs, qe, ql, cs, ce, cl, size_sim, params['rec_ovl'], bp_win)
            if ok:
                bp_dist = abs(qs - cs) + abs(qe - ce)
                if bp_dist < best_bp_dist:
                    best_bp_dist = bp_dist
                    best_match = str(cand_id)

        # Fallback for DUPs by looking for a match with an INS
        if svtype == 'DUP' and best_match is None:
            orig_pos_info = rec.info.get('ORIGINAL_POS', None)
            if orig_pos_info is not None:
                if isinstance(orig_pos_info, (list, tuple)):
                    orig_pos_info = orig_pos_info[0]
                ins_qs = int(orig_pos_info) - 1
                ins_params = PARAMS['INS']
                ins_bp_win = ins_params['bp_win']
                ins_size_sim = ins_params['size_sim']
                ins_region_start = max(0, ins_qs - ins_bp_win)
                ins_region_end = ins_qs + 1 + ins_bp_win
                for cand in dbvar.fetch(chrom, ins_region_start, ins_region_end):
                    if cand.info.get('SVTYPE') != 'INS':
                        continue
                    cs = cand.start
                    ce = cand.stop
                    cl = parse_dbvar_len(cand)
                    cand_id = cand.id if cand.id else cand.info.get('DBVARID', '.')
                    if isinstance(cand_id, (list, tuple)):
                        cand_id = cand_id[0]
                    if passes_ins(ins_qs, ins_qs + 1, ql, cs, ce, cl, ins_size_sim, ins_params['rec_ovl'], ins_bp_win):
                        bp_dist = abs(ins_qs - cs) + abs(ins_qs + 1 - ce)
                        if bp_dist < best_bp_dist:
                            best_bp_dist = bp_dist
                            best_match = str(cand_id)

        # Add best match to output
        if best_match:
            out_pos = rec.info.get('ORIGINAL_POS', None)
            if out_pos is not None:
                if isinstance(out_pos, (list, tuple)):
                    out_pos = out_pos[0]
                out_pos = int(out_pos)
            else:
                out_pos = qs + 1
            out.write(f"{chrom}\t{out_pos}\t{orig_ref}\t{orig_alt}\t{var_id}\t{best_match}\n")

query.close()
dbvar.close()
EOF

        sort -k1,1 -k2,2n ~{prefix}.annotations.unsorted.tsv > ~{prefix}.annotations.tsv
    >>>

    output {
        File annotations_tsv = "~{prefix}.annotations.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size([vcf, original_vcf, dbvar_vcf], "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
