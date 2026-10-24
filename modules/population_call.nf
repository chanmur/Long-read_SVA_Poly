process POPULATION_SV {
    tag "PopulationSV $population_tag"
    label 'STRUCTURAL_VARIANT_CALLING'

    publishDir "${params.out_dir}/01_structural_variants/population",  
        mode: 'copy',
        pattern: "*.vcf"

    input:
    tuple val(population_tag), val(ref_tag), path(reference), val(sv_type), path(snf_files)


    output:
    val population_tag, emit: population_tags
    val ref_tag, emit: ref_tags
    val sv_type, emit: sv_types
    tuple val(population_tag), val(ref_tag), val(sv_type), path("${population_tag}_${ref_tag}_${sv_type}.vcf"), emit: vcf_files

    script:
    """
    sniffles --input ${snf_files} --vcf ${population_tag}_${ref_tag}_${sv_type}.vcf --reference ${reference} --allow-overwrite --threads 28 --output-rnames --combine-max-inmemory-results 1000
    """
}