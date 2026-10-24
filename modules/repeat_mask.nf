process REPEATMASK {
    tag "RepeatMask $sv_fasta"
    label 'TE_ANNOTATIONS'

    publishDir "${params.out_dir}/02_te_annotation",  
        mode: 'copy',
        pattern: "${outdir}"

    input:
    val(te_lib)
    val(te_lib_tag)
    tuple val(key), path(sv_fasta)

    output:
    val(key), emit: keys
    tuple val(key), path(sv_fasta), path("${outdir}"), emit: repeatmask_outdir
    tuple val(key), path(sv_fasta), path("${outdir}/${outfile_prefix}.fa.out"), emit: repeatmask_out
    tuple val(key), path(sv_fasta), path("${outdir}/${outfile_prefix}.fa.tbl"), emit: repeatmask_tbl, optional: true
    tuple val(key), path(sv_fasta), path("${outdir}/${outfile_prefix}.fa.masked"), emit: repeatmask_masked, optional: true
    

    script:
    outdir = "${sv_fasta.simpleName}_${te_lib_tag}"
    outfile_prefix = "${sv_fasta.simpleName}_${te_lib_tag}"
    def repeatmask_input = "${outfile_prefix}.fa"
    
    """
    mkdir ${outdir}
    ln -s -- "${sv_fasta}" "${repeatmask_input}"
    RepeatMasker -lib /repeat_library/${te_lib} -pa 22 -dir ${outdir} ${repeatmask_input}
    """
}

process REFINE_SVA_ANNOTATIONS {
    tag "REFINE_SVA_ANNOTATIONS $rmsk_sv_out"
    label 'TE_ANNOTATIONS'

    publishDir "${params.out_dir}/02_te_annotation",  
        mode: 'copy',
        pattern: "${outfile}"

    input:
    tuple val(key), path(reference), path(hexamer_fa), path(MAST2_fa), path(sv_fasta), path(rmsk_sv_out)


    output:
    val(key), emit: keys
    tuple val(key), path(sv_fasta), path(rmsk_sv_out), path("${outfile}"), emit: refined_sva_out

    script:
    script_dir = "/sw/pkg/xTea/v0.1.0-long/xtea_long"
    work_dir = "xtea_workdir"
    outfile = "${sv_fasta.simpleName}_xtea_sva.txt"
    """
    python ${script_dir}/classify_TE_insertion_from_rmsk.py --seq ${sv_fasta} --rmsk ${rmsk_sv_out} --ref ${reference} --hexam ${hexamer_fa} --mast2 ${MAST2_fa} -n 22 -p ${work_dir} -o ${outfile}
    """
}