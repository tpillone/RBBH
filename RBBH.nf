

process gbk2faa {
	publishDir 'data/faa', mode: 'copy', overwrite: true

	//container "$params.chlamdb_container"

	cpus 1

    echo true 

	input:
	file genome_gbk from Channel.fromPath(params.gbk_dir + "/*.gbff")

	output:
	file "*.faa" into faa_files

	script:
	"""
        #!/usr/bin/env python
        import RBBH 
        print("OK depart!")
        RBBH.gbk2faa("$genome_gbk", 
                     "${genome_gbk.baseName}")
	"""
}

faa_files.collect().into{
    faa_list1
    faa_list2
}


process ssearch_all_vs_all {
    publishDir 'data/ssearch', mode: 'copy', overwrite: true

    container "$params.fasta_container"

    cpus 2

    input:
      each faa1 from faa_list1
      each faa2 from faa_list2

    output:
        file "*.tab" into ssearch_results
    script:
    """
    ssearch36 -T ${task.cpus} -m 8 -s BL50 -b 10 -d 0 -E 1e-5 $faa1 $faa2 -z 1 > ${faa1.baseName}_vs_${faa2.baseName}.tab
    """
}


process RBBH_table {
	publishDir 'data/RBBH_tables', mode: 'copy', overwrite: true

	//container "$params.chlamdb_container"

	cpus 1

    echo true 

	input:
	file genome_gbk_files from Channel.fromPath(params.gbk_dir + "/*.gbff").collect()
    file ssearch_result_files from ssearch_results.collect()

	output:
	file "*RBBH.tab" into rbbh_tables

	script:
	"""
        #!/usr/bin/env python
        import RBBH 
        print("OK depart!")
        RBBH.RBBH_table("$genome_gbk_files".split(" "), 
                       "$ssearch_result_files".split(" "))
	"""
}