nextflow.enable.dsl=2

fasta_ch = file("./data_test/REF/*.fasta")
fastq_ch = Channel.fromFilePairs("./data_test/FASTQ/*_{1,2}.fastq.gz").view()

process trim {
    input:
        tuple val(label), path(fastq)
		
    output:
        path "${label}"

    publishDir "TRIM", mode: 'copy'

    shell:
    """
	trim_galore --paired !{fastq[0]} !{fastq[1]} -o !{label}
    """
}

process buildIndex {
    input:
        path fasta

    output:
        path "transcripts_index"

    publishDir "INDEX", mode: 'copy'

    shell:
    """
    salmon index -t !{fasta} -i transcripts_index
    """
}

process quant {
    input:
        path index
        tuple val(label), path(fastq)

    output:
        path "${label}"

    publishDir "EXPRESSION", mode: 'copy'

    shell:
    """
    salmon quant -i !{index} -l A -1 !{fastq[0]} -2 !{fastq[1]} -o !{label}
    """
}


process splice {
    input:
		        

    output:
       

    publishDir 

    shell:
    """
    multipleFieldSelection.py 
    """
}

workflow {

    trim(fastq_ch)
    buildIndex(fasta_ch)
    quant(buildIndex.out, trim.out)
    splice()

}
