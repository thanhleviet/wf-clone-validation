#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.fofn = "/share/nick_plasmids/canu_input.full.fofn"
params.out_dir = "/share/nick_plasmids/results.nov/compare"
params.db_directory = "/share/nick_plasmids/wf-clone-validation-db"

ch_input0 = Channel.fromPath(params.fofn)
                  .splitCsv(header: true)
                  .map {row -> tuple(row.sample_id, row.ref, row.query)}


process ROTATE {
    tag {sample_id}
    
    conda '/home/ubuntu/miniconda3/envs/pLannotate'
    
    cpus 4
    
    memory '16.GB'
    
    input:
    tuple val(sample_id), path(ref), path(query)
    output:
    tuple val(sample_id), path(ref), path(query)
    script:
    """
    plannotate batch -i ${ref} -f ${sample_id}_ref -b ${params.db_directory} -h
    convert.py rotate ${sample_id}_ref.gbk
    mv ${sample_id}_ref.final.fasta ${ref}
    """
}


process VCF {
    publishDir "${params.out_dir}/map", mode: "copy"
    
    tag {sample_id}
    
    cpus 8
    
    memory '12.GB'

    input:
    tuple val(sample_id), path("ref.fa"), path("query.fa")
    output:
    path("${sample_id}.*")

    script:
    """
    minimap2 -c --cs ref.fa query.fa | sort -k6,6 -k8,8n >  ${sample_id}.asm
    k8 ${projectDir}/bin/paftools.js stat ${sample_id}.asm > ${sample_id}.stat
    k8 ${projectDir}/bin/paftools.js call -g 10 -f ref.fa -L100 ${sample_id}.asm > ${sample_id}.vcf
    k8 ${projectDir}/bin/paftools.js call  -g 10 -L100 ${sample_id}.asm > ${sample_id}.txt
    """
}

process LASTAL {
    publishDir "${params.out_dir}/lastal/${sample_id}/", mode: "copy" 
    
    tag {sample_id}
    errorStrategy 'ignore'
    cpus 8
    
    memory '12.GB' 

    input:
    tuple val(sample_id), path("ref.fa"), path("query.fa")
    output:
    path("${sample_id}.{maf,png,html,tab}")
    
    script:
    """
    lastdb ${sample_id} ref.fa
    lastal -j4 ${sample_id} query.fa > ${sample_id}.maf
    last-dotplot ${sample_id}.maf ${sample_id}.png
    maf-convert tab  ${sample_id}.maf > ${sample_id}.tab
    maf-convert html  ${sample_id}.maf > ${sample_id}.html
    """
}

process MUSCLE {
    publishDir "${params.out_dir}/muscle", mode: "copy"  
    
    tag {sample_id}
    
    cpus 4
    memory '16.GB'
    errorStrategy 'ignore'

    input:
        tuple val(sample_id), path("ref.fa"), path("query.fa")
    output:
        path("${sample_id}.aln")
    script:
    """
    cat ref.fa query.fa > in.fa
    muscle -in in.fa -out ${sample_id}.aln
    """
}

process PLANNOTATE {
    publishDir "${params.out_dir}/plannotate", mode: "copy" 

    tag {sample_id}
    
    conda '/home/ubuntu/miniconda3/envs/pLannotate'

    cpus 8

    input:
    tuple val(sample_id), path(ref), path(contigs)
    output:
    tuple val(sample_id), path("${ref.baseName}.gbk"), path("${contigs.baseName}.gbk"), emit: gbk
    path("*.html"), emit: html
    script:
    """
    plannotate batch -i ${ref} -f ${ref.baseName} -b ${params.db_directory} -h
    plannotate batch -i ${contigs} -f ${contigs.baseName} -b ${params.db_directory} -h
    """
}

process CLINKER {
    publishDir "${params.out_dir}/clinker", mode: "copy" 
    
    tag {sample_id}
    
    conda '/home/ubuntu/miniconda3/envs/clinker'

    cpus 8

    input:
    tuple val(sample_id), path(ref), path(contigs)
    
    output:
    path("*.html")
    
    script:
    """
    convert.py fix-feature-type $ref
    mv new_${ref} $ref
    
    convert.py fix-feature-type $contigs
    mv new_${contigs} $contigs
    clinker $ref $contigs -p ${sample_id}.html
    """
}

workflow {
    ch_input = ROTATE(ch_input0)
    VCF(ch_input)
    LASTAL(ch_input)
    MUSCLE(ch_input)

    PLANNOTATE(ch_input)
    CLINKER(PLANNOTATE.out.gbk)
}