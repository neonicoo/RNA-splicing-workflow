nextflow.enable.dsl=2

params.reads = null
params.ref = null
params.formatscript = null 
params.annot = null

gtf_file = file(params.annot)
R_file = file(params.formatscript)
fasta_file = file(params.ref)
fastq_ch = Channel.fromFilePairs(params.reads + "*_{1,2}.{fastq,fq}.gz").view()

process trim {
    input:
        tuple val(label), path(fastq)
		
    output:
        tuple val("${label}"), file("${label}/${label}_1_val_1.fq.gz"), file("${label}/${label}_2_val_2.fq.gz")
	
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

    shell:
    """
    salmon index -t !{fasta} -i transcripts_index
    """
}

process quant {
    input:
        path index
        tuple val(label), path(fastq1), path(fastq2)

    output:
        path "${label}"

    shell:
    """
    salmon quant -i !{index} -l A -1 !{fastq1} -2 !{fastq2} -o !{label}
    """
}

process formatTPM {
	input:
		path quant

	output:
		file "iso_tpm_formatted.txt"

	shell:
	"""
	multipleFieldSelection.py -i !{quant}/quant.sf -k 1 -f 4 -o iso_tpm.txt
	Rscript !{R_file} iso_tpm.txt
	"""
}


process generateEvents {
    output:
       file "transcripts.events.ioe"

    shell:
    """
    suppa.py generateEvents -i !{gtf_file} -o transcripts.events -e SE SS MX RI FL -f ioe
    awk '
        FNR==1 && NR!=1 { while (/^<header>/) getline; }
        1 {print}
    ' *.ioe > transcripts.events.ioe
    """
}

process psiPerEvent{
	input:
		file transcripts
		file isoTPM
		
	output:
		file "TRA2_events.psi"
	
	publishDir "SPLICING", mode: 'copy'
	
	shell:
	"""
	suppa.py psiPerEvent -i !{transcripts} -e !{isoTPM} -o TRA2_events
	"""
}



workflow {
    trim(fastq_ch)
    buildIndex(fasta_file)
    quant(buildIndex.out, trim.out)
    formatTPM(quant.out)
    generateEvents()
    psiPerEvent(generateEvents.out, formatTPM.out)
}
