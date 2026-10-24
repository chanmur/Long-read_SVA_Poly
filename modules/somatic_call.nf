process SOMATIC_SV {
    tag "SomaticSV ${metadata.sample_id}"
    label 'STRUCTURAL_VARIANT_CALLING'

    publishDir (
        path: "${params.out_dir}/01_structural_variants/somatic/${metadata.condition}",  
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
    tuple val(metadata), val(ref_tag), val("somatic"), path("${metadata.sample_name}_${ref_tag}_somatic.vcf"), emit: vcf_files
    tuple val(metadata), val(ref_tag), val("somatic"), path("${metadata.sample_name}_${ref_tag}_somatic.snf"), emit: snf_files

    script:
    """
    sniffles --input ${bam_file} --vcf ${metadata.sample_name}_${ref_tag}_somatic.vcf --snf ${metadata.sample_name}_${ref_tag}_somatic.snf --reference ${reference} --minsvlen 10 --threads 28 --mosaic
    """
}

process SOMATIC_SV_PHASED {
    tag "SomaticSV ${metadata.sample_id}"
    label 'STRUCTURAL_VARIANT_CALLING'

    publishDir (
        path: "${params.out_dir}/01_structural_variants/somatic/${metadata.condition}",  
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
    tuple val(metadata), val(ref_tag), val("somatic_phased"), path("${metadata.sample_name}_${ref_tag}_somatic_phased.vcf"), emit: vcf_files
    tuple val(metadata), val(ref_tag), val("somatic_phased"), path("${metadata.sample_name}_${ref_tag}_somatic_phased.snf"), emit: snf_files

    script:
    """
    sniffles --input ${bam_file} --vcf ${metadata.sample_name}_${ref_tag}_somatic_phased.vcf --snf ${metadata.sample_name}_${ref_tag}_somatic_phased.snf --reference ${reference} --mosaic --phase --minsvlen 10 --threads 28 
    """
}