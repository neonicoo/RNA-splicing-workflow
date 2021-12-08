nextflow.enable.dsl=2

fasta_ch = file("./data_test/REF/*.fasta")
fastq_ch = Channel.fromFilePairs("./data_test/FASTQ/*_{1,2}.fastq.gz").view()

process trim {
    input:
        tuple val(label), path(fastq)
		
    output:
        path "${label}"

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


process generateEvents {
    input:
		path quant
		path annotation

    output:
       path "iso_tpm.txt"

    publishDir "SUPPA", mode: 'copy'

    shell:
    """
    suppa.py generateEvents -i !{annotation} -o ensemble.events -e SE SS MX RI FL -f ioe
    awk '
		FNR==1 && NR!=1 { while (/^<header>/) getline; }
		1 {print}
	'./SUPPA/*.ioe > ./SUPPA/ensemble.events.ioe
    """
}



workflow {

    trim(fastq_ch)
    buildIndex(fasta_ch)
    quant(buildIndex.out, trim.out)
    splice()

}
