#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// importing modules
include {GERMLINE_SV; GERMLINE_SV_PHASED} from './modules/germline_call'
include {SOMATIC_SV; SOMATIC_SV_PHASED} from './modules/somatic_call'
include {POPULATION_SV} from './modules/population_call'
include {COMPRESS_VCF; INDEX_VCF; FILTER_VCF; VCF_TO_FASTA} from './modules/vcf_processing'
include {COMPRESS_VCF as COMPRESS_VCF_FILTERED; INDEX_VCF as INDEX_VCF_FILTERED} from './modules/vcf_processing'
include {REPEATMASK; REFINE_SVA_ANNOTATIONS} from './modules/repeat_mask'

workflow read_metadata {
       take: 
       metadata
       
       main:
           bam_files = Channel.fromPath(metadata, checkIfExists: true)
                .splitCsv(header: true)
                .map { row -> 
                    // Validating required fields
                    if (!row.sample_name) error "Missing sample_name in row: ${row}"
                    if (!row.sample_id) error "Missing sample_id in row: ${row}"
                    if (!row.donor_id) error "Missing donor_id in row: ${row}"
                    if (!row.project) error "Missing project in row: ${row}"
                    if (!row.condition) error "Missing condition in row: ${row}"
                    if (!row.brain_region) error "Missing brain_region in row: ${row}"
                    if (!row.age) error "Missing age in row: ${row}"
                    if (!row.path) error "Missing path in row: ${row}"
                    if (!row.bam_file) error "Missing bam_file in row: ${row}"
                    if (!row.index_file) error "Missing index_file in row: ${row}"

                    def meta = [sample_name: row.sample_name, sample_id: row.sample_id, donor_id: row.donor_id, dataset: row.project, condition: row.condition, brain_region: row.brain_region, age: row.age, sex: row.sex, pmi: row.pmi]
                    def bam_file = file(row.path + '/' + row.bam_file, checkIfExists: true)
                    def index_file = file(row.path + '/' + row.index_file, checkIfExists: true)
                    tuple(meta, bam_file, index_file)
                }
       
       emit: 
       bam_files
   }

// Modified structural_variant_call workflow to include the common_ref step. This allows for flexibility in using different reference genomes and tags for structural variant calling.
workflow structural_variant_call {
    take:
    ref_tag
    reference
    bam_files

    main:
    bam_files_ref = bam_files
                        .map { it -> tuple(it[0], it[1], it[2], ref_tag, reference) }
    
    GERMLINE_SV(bam_files_ref)
    SOMATIC_SV(bam_files_ref)
    
    emit:
    germline_vcf = GERMLINE_SV.out.vcf_files
    germline_snf = GERMLINE_SV.out.snf_files
    somatic_vcf = SOMATIC_SV.out.vcf_files
    somatic_snf = SOMATIC_SV.out.snf_files
}

workflow structural_variant_call_phased {
    take:
    ref_tag
    reference
    bam_files

    main:
    bam_files_ref = bam_files
                        .map { it -> tuple(it[0], it[1], it[2], ref_tag, reference) }
    
    GERMLINE_SV_PHASED(bam_files_ref)
    SOMATIC_SV_PHASED(bam_files_ref)
    
    emit:
    germline_vcf = GERMLINE_SV_PHASED.out.vcf_files
    germline_snf = GERMLINE_SV_PHASED.out.snf_files
    somatic_vcf = SOMATIC_SV_PHASED.out.vcf_files
    somatic_snf = SOMATIC_SV_PHASED.out.snf_files
}

workflow population_variant_call {
    take:
    population_tag
    ref_tag
    reference
    snf_files
    
    main:
    def expected_ref_tag = ref_tag
    population_in = snf_files
                        .map {  ch_ref_tag, sv_type, snf_file -> 
                            if (ch_ref_tag != ref_tag) {
                                error "Reference tag mismatch: The snf files were generated with reference tag '${ch_ref_tag}', but the population call is being run with reference tag '${ref_tag}'. Please ensure that the correct reference is being used."
                            }
                            tuple(population_tag, ref_tag, reference, sv_type, snf_file)
                        }
    POPULATION_SV(population_in)

    population_sv = POPULATION_SV.out.vcf_files
    
    emit:
    population_sv
}

// Requirement 1: the last element of the tuple in the input channel should be a vcf file. Accepts only 1 vcf file per tuple.
// Requirement 2: the first element of the tuple should be unique for each tuple in the input channel. For example, the metadata element from germline or somatic channels, or the population tag from the population channel.
// The following workflow takes any vcf file containing output channel, compresses and indexes the vcf file and emits the indexed vcf file. 
// The remaining contents of the channel are preserved and emitted as well. 
// Input: tuple(non-vcf elements..., vcf_file) Output: tuple(non-vcf elements..., compressed_vcf_file, indexed_vcf_file)
// This workflow is useful for indexing vcf files generated from the structural variant calling workflows.
workflow index_vcf {
    take:
    vcf_ch

    main:
    vcf_ch_split = vcf_ch
                        .multiMap {vcf_tuple ->
                            def key = vcf_tuple[0..2]
                            non_vcf_elements: tuple(key, *vcf_tuple[0..-2])
                            vcf: tuple(key, vcf_tuple[-1])
                        }
    COMPRESS_VCF(vcf_ch_split.vcf)
    COMPRESS_VCF.out.compressed_vcf_files | INDEX_VCF

    indexed_vcf_ch = vcf_ch_split.non_vcf_elements
                        .join(INDEX_VCF.out.indexed_vcf_files, by: 0)
                        .map {it -> it[1..-1]}

    emit:
    indexed_vcf_ch
}

workflow index_vcf_filtered {
    take:
    vcf_ch

    main:
    vcf_ch_split = vcf_ch
                        .multiMap {vcf_tuple ->
                            def key = vcf_tuple[0..2]
                            non_vcf_elements: tuple(key, *vcf_tuple[0..-2])
                            vcf: tuple(key, vcf_tuple[-1])
                        }
    COMPRESS_VCF_FILTERED(vcf_ch_split.vcf)
    COMPRESS_VCF_FILTERED.out.compressed_vcf_files | INDEX_VCF_FILTERED

    indexed_vcf_ch = vcf_ch_split.non_vcf_elements
                        .join(INDEX_VCF_FILTERED.out.indexed_vcf_files, by: 0)
                        .map {it -> it[1..-1]}

    emit:
    indexed_vcf_ch
}

workflow filter_vcf {
    take:
    vcf_ch

    main:
    vcf_ch_split = vcf_ch
                        .multiMap {vcf_tuple ->
                            def key = vcf_tuple[0..2]
                            non_vcf_elements: tuple(key, *vcf_tuple[0..-2])
                            vcf: tuple(key, vcf_tuple[-1])
                        }
    FILTER_VCF(vcf_ch_split.vcf)
    

    filtered_vcf_ch = vcf_ch_split.non_vcf_elements
                        .join(FILTER_VCF.out.filtered_vcf_files, by: 0)
                        .map {it -> it[1..-1]}

    emit:
    filtered_vcf_ch
}

workflow te_annot_sv {
    take:
    te_lib
    te_lib_tag
    vcf_ch

    main:
    vcf_ch_split = vcf_ch
                        .multiMap {vcf_tuple ->
                            def key = vcf_tuple[0..2]
                            non_vcf_elements: tuple(key, *vcf_tuple[0..-3])
                            vcf: tuple(key, vcf_tuple[-2], vcf_tuple[-1])
                        }
    VCF_TO_FASTA(vcf_ch_split.vcf)

    REPEATMASK(te_lib, te_lib_tag, VCF_TO_FASTA.out.sv_fasta_files)

    te_annot_sv_ch = vcf_ch_split.non_vcf_elements
                        .join(vcf_ch_split.vcf, by: 0)
                        .join(REPEATMASK.out.repeatmask_out, by: 0)
                        .map {it -> it[1..-1]}

    emit:
    te_annot_sv_ch
}

workflow refine_sva {
    take:
    ref_tag
    reference
    mast2_fa
    hexamer_fa
    te_annot_sv_ch

    main:

    xtea_in = te_annot_sv_ch
                        .multiMap {  te_ann_tuple -> 
                            if (te_ann_tuple[1] != ref_tag) {
                                error "Reference tag mismatch: The files were generated with reference tag '${te_ann_tuple[1]}', but the population call is being run with reference tag '${ref_tag}'. Please ensure that the correct reference is being used."
                            }

                            def key = te_ann_tuple[0..2]
                            other_elements: tuple(key, *te_ann_tuple[0..-3])
                            in_elements: tuple(key, reference, hexamer_fa, mast2_fa, te_ann_tuple[-2], te_ann_tuple[-1])
                        }

    REFINE_SVA_ANNOTATIONS(xtea_in.in_elements)

    refined_sva_ch = xtea_in.other_elements
                        .join(REFINE_SVA_ANNOTATIONS.out.refined_sva_out, by: 0)
                        .map {it -> it[1..-1]}

    emit:
    refined_sva_ch
}

workflow sv_sva_pipeline {
    
    take:
    metadata
    ref_tag
    reference
    mast2_fa
    hexamer_fa
    population_tag
    te_lib
    te_lib_tag

    main:
    read_metadata(metadata)
    
    bam_files = read_metadata.out.bam_files
    
    structural_variant_call(ref_tag, reference, bam_files)

    snf_files_ch = structural_variant_call.out.germline_snf
                     .mix(structural_variant_call.out.somatic_snf)
                     .map { meta_data, refs_tag, sv_type, snf_file -> tuple(refs_tag, sv_type, snf_file) }
                     .groupTuple(by: [0, 1])

    population_sv_ch = population_variant_call(population_tag, ref_tag, reference, snf_files_ch)

    unfiltered_indexed_vcf_ch = index_vcf(population_sv_ch)
    
    filtered_indexed_vcf_ch = filter_vcf(population_sv_ch) | index_vcf_filtered

    te_annot_sv_ch = te_annot_sv(te_lib, te_lib_tag, filtered_indexed_vcf_ch)

    refine_sva(ref_tag, reference, mast2_fa, hexamer_fa, te_annot_sv_ch)

}

workflow sv_sva_phased_pipeline {
    
    take:
    metadata
    ref_tag
    reference
    mast2_fa
    hexamer_fa
    population_tag
    te_lib
    te_lib_tag

    main:
    read_metadata(metadata)
    
    bam_files = read_metadata.out.bam_files
    
    structural_variant_call_phased(ref_tag, reference, bam_files)

    snf_files_ch = structural_variant_call_phased.out.germline_snf
                     .mix(structural_variant_call_phased.out.somatic_snf)
                     .map { meta_data, refs_tag, sv_type, snf_file -> tuple(refs_tag, sv_type, snf_file) }
                     .groupTuple(by: [0, 1])

    population_sv_ch = population_variant_call(population_tag, ref_tag, reference, snf_files_ch)

    unfiltered_indexed_vcf_ch = index_vcf(population_sv_ch)
    
    filtered_indexed_vcf_ch = filter_vcf(population_sv_ch) | index_vcf_filtered

    te_annot_sv_ch = te_annot_sv(te_lib, te_lib_tag, filtered_indexed_vcf_ch)

    refine_sva(ref_tag, reference, mast2_fa, hexamer_fa, te_annot_sv_ch)

}

workflow {
    def ref = file(params.reference, checkIfExists: true)
    
    if (!(params.phased instanceof Boolean)) {
        error "Invalid value for params.phased: '${params.phased}'. It must be either true or false."
    }
    
    if (params.phased)
        sv_sva_phased_pipeline(params.metadata, params.reference_tag, ref, params.mast2_fa, params.hexamer_fa, params.population_tag, params.te_library, params.te_library_tag)
    else
        sv_sva_pipeline(params.metadata, params.reference_tag, ref, params.mast2_fa, params.hexamer_fa, params.population_tag, params.te_library, params.te_library_tag)
    
}