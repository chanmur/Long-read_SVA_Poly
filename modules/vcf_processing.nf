process COMPRESS_VCF {
    tag "COMPRESS_VCF $vcf_file"
    label 'VCF_PROCESSING'

    publishDir "${params.out_dir}/01_structural_variants/population",  
        mode: 'copy',
        pattern: "${vcf_file}.gz"

    input:
    tuple val(key), path(vcf_file)


    output:
    val key, emit: keys
    tuple val(key), path("${vcf_file}.gz"), emit: compressed_vcf_files

    script:
    """
    bgzip -c -@ 5 ${vcf_file} > ${vcf_file}.gz
    """
}

process INDEX_VCF {
    tag "INDEX_VCF $compressed_vcf_file"
    label 'VCF_PROCESSING'

    publishDir "${params.out_dir}/01_structural_variants/population",  
        mode: 'copy',
        pattern: "${compressed_vcf_file}.tbi"

    input:
    tuple val(key), path(compressed_vcf_file)


    output:
    val key, emit: keys
    tuple val(key), path(compressed_vcf_file), path("${compressed_vcf_file}.tbi"), emit: indexed_vcf_files

    script:
    """
    tabix -p vcf ${compressed_vcf_file}
    """
}

process FILTER_VCF {
    tag "FILTER_VCF $vcf_file"
    label 'VCF_PROCESSING'
    label 'BCFTOOLS'

    publishDir "${params.out_dir}/01_structural_variants/population",  
        mode: 'copy',
        pattern: "${vcf_file.simpleName}_filtered.vcf"

    input:
    tuple val(key), path(vcf_file)


    output:
    val key, emit: keys
    tuple val(key), path("${vcf_file.simpleName}_filtered.vcf"), emit: filtered_vcf_files

    script:
    def output_vcf = "${vcf_file.simpleName}_filtered.vcf"
    """
    bcftools filter -i '(INFO/SVTYPE == "INS" || INFO/SVTYPE == "DEL") && INFO/SVLEN != "." && ABS(INFO/SVLEN) >= 50 && ALT[0] != "<INS>" && ALT[0] != "<DEL>"' -Ov -o ${output_vcf} ${vcf_file}
    """
}

process VCF_TO_FASTA {
    tag "VCF_TO_FASTA $vcf_file"
    label 'VCF_PROCESSING'
    label 'BCFTOOLS'

    publishDir "${params.out_dir}/01_structural_variants/population",  
        mode: 'copy',
        pattern: "${vcf_file.simpleName}.fa"

    input:
    tuple val(key), path(vcf_file), path(vcf_index)


    output:
    val key, emit: keys
    tuple val(key), path("${vcf_file.simpleName}.fa"), emit: sv_fasta_files

    script:
    def output_fa = "${vcf_file.simpleName}.fa"
    """
    bcftools query -i 'INFO/ILEN>0' -f '>%ID\\n%ALT\\n' ${vcf_file} > ${output_fa}
    bcftools query -i 'INFO/ILEN<0' -f '>%ID\\n%REF\\n' ${vcf_file} >> ${output_fa}
    """
}