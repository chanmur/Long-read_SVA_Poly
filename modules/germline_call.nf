process GERMLINE_SV {
    tag "GermlineSV ${metadata.sample_id}"
    label 'STRUCTURAL_VARIANT_CALLING'

    publishDir( 
        path: "${params.out_dir}/01_structural_variants/germline/${metadata.condition}",  
        mode: 'copy',
        saveAs: { filename ->
                if (filename.endsWith('.vcf')) {
                    "vcf/${filename}"
                }
                else if (filename.endsWith('.snf')) {
                    "snf/${filename}"
                }
                else {
                    null
                }
            }
    )

    input:
    tuple val(metadata), path(bam_file), path(bai_file), val(ref_tag), path(reference)

    output:
    val metadata.sample_id, emit: sample_ids
    tuple val(metadata), val(ref_tag), val("germline"), path("${metadata.sample_name}_${ref_tag}_germline.vcf"), emit: vcf_files
    tuple val(metadata), val(ref_tag), val("germline"), path("${metadata.sample_name}_${ref_tag}_germline.snf"), emit: snf_files

    script:
    """
    sniffles --input ${bam_file} --vcf ${metadata.sample_name}_${ref_tag}_germline.vcf --snf ${metadata.sample_name}_${ref_tag}_germline.snf --reference ${reference} --minsvlen 10 --threads 28
    """
}

process GERMLINE_SV_PHASED {
    tag "GermlineSV ${metadata.sample_id}"
    label 'STRUCTURAL_VARIANT_CALLING'

    publishDir( 
        path: "${params.out_dir}/01_structural_variants/germline/${metadata.condition}",  
        mode: 'copy',
        saveAs: { filename ->
                if (filename.endsWith('.vcf')) {
                    "vcf/${filename}"
                }
                else if (filename.endsWith('.snf')) {
                    "snf/${filename}"
                }
                else {
                    null
                }
            }
    )

    input:
    tuple val(metadata), path(bam_file), path(bai_file), val(ref_tag), path(reference)

    output:
    val metadata.sample_id, emit: sample_ids
    tuple val(metadata), val(ref_tag), val("germline_phased"), path("${metadata.sample_name}_${ref_tag}_germline_phased.vcf"), emit: vcf_files
    tuple val(metadata), val(ref_tag), val("germline_phased"), path("${metadata.sample_name}_${ref_tag}_germline_phased.snf"), emit: snf_files

    script:
    """
    sniffles --input ${bam_file} --vcf ${metadata.sample_name}_${ref_tag}_germline_phased.vcf --snf ${metadata.sample_name}_${ref_tag}_germline_phased.snf --reference ${reference} --minsvlen 10 --phase --threads 28
    """
}