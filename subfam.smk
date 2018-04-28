SAMPLE_NAMES = ['3.6k.fam01151.fa']

rule all:
	input:
		expand("{sample}_results/output/subfamily/orf2subfamily.tsv", sample = SAMPLE_NAMES)
		
rule createDB_mmseqs2:
	input:
		"data/{sample}"
	output:
		"{sample}_results/database/mmseqs/{sample}.mmseqsDB"
	run:
		shell("mmseqs createdb {input} {output}")

rule cluster_mmseqs2:
	input:
		db = "{sample}_results/database/mmseqs/{sample}.mmseqsDB",
	output:
		cluster = "{sample}_results/database/mmseqs/{sample}.mmseqsDB_clu",
		tmp = "{sample}_results/database/mmseqs/tmp"
	threads: 6
	run:
		shell("mkdir {output.tmp}")
		shell("mmseqs cluster {input.db} {output.cluster} {output.tmp} --threads {threads} -s 6 -c 0.7 --cov-mode 0 --max-seqs 5000 -e 0.001 --cluster-mode 0")

rule createtsv_mmseqs2:
	input:
		db = "{sample}_results/database/mmseqs/{sample}.mmseqsDB",
		cluster = "{sample}_results/database/mmseqs/{sample}.mmseqsDB_clu"
	output:
		"{sample}_results/database/mmseqs/{sample}.mmseqsDB_clu.tsv"
	shell:
		"mmseqs createtsv {input.db} {input.db} {input.cluster} {output}"


rule reportSubfamily:
	input:
		"{sample}_results/database/mmseqs/{sample}.mmseqsDB_clu.tsv",
		"data/{sample}"
	output:
		"{sample}_results/output/subfamily/orf2subfamily.tsv",
		"{sample}_results/output/subfamily/config.json"	
	shell:
		"/home/meheurap/snakemake/proteinCluster/creatingSubfamilyReport.py {input} {output}"

