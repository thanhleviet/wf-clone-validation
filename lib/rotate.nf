process ROTATE {
    // publishDir  
   
    tag {contigs.simpleName}
    
    errorStrategy = 'ignore'

    cpus 4
    memory '16.GB'

    input:
        path(contigs)
        env STATUS
    output:
        path("${contigs.baseName}.fasta"), emit: polished
        env STATUS, emit: status
    script:
    """
    plannotate batch -i ${contigs} -f ${contigs.baseName} -b ${params.db_directory} -h
    convert.py rotate ${contigs.baseName}.gbk
    """
}