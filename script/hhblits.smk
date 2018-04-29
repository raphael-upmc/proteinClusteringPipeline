configfile : "3.6k.fam01151.fa_results/output/subfamily/config.json"

SAMPLE_NAMES = ['3.6k.fam01151.fa']

rule all:
	input:
		expand("{sample}_results/database/hhblits/hhm/{clusterId}.hhm",clusterId=config["clusters"], sample = SAMPLE_NAMES),
		expand("{sample}_results/database/hhblits/hhr/{clusterId}.hhr",clusterId=config["clusters"], sample = SAMPLE_NAMES),		
		expand("{sample}_results/database/hhblits/db", sample = SAMPLE_NAMES),
		expand("{sample}_results/database/hhblits/allvsall.hhr", sample = SAMPLE_NAMES)		

rule mafft:
	input:
		"{sample}_results/output/subfamily/fasta/{clusterId}.fa"
	output:
		"{sample}_results/database/hhblits/mafft/{clusterId}.mafft"
	run:
		shell("mafft --auto {input} > {output}")

rule a3m:
	input:
		"{sample}_results/database/hhblits/mafft/{clusterId}.mafft"
	output:
		"{sample}_results/database/hhblits/mafft/{clusterId}.a3m"
	run:
     		shell("reformat.pl fas a3m {input} {output} -M 50")

rule reformat_a3m:
	input:
		"{sample}_results/database/hhblits/mafft/{clusterId}.a3m"
	output:
		"{sample}_results/database/hhblits/mafft/{clusterId}_renamed.a3m"     
	run:
     		shell("/home/meheurap/snakemake/proteinCluster/addingNameToA3m.py {input} {output}")


rule add_ss:
	input:
		"{sample}_results/database/hhblits/mafft/{clusterId}_renamed.a3m"     
	output:
		"{sample}_results/database/hhblits/a3m/{clusterId}.a3m"
	run:
     		shell("addss.pl {input} {output} -a3m")

rule hhm:
	input:
		"{sample}_results/database/hhblits/a3m/{clusterId}.a3m"
	output:
		"{sample}_results/database/hhblits/hhm/{clusterId}.hhm"
	run:		
		shell("hhmake -i {input} -o {output}")

rule makehhblitsdb :
	input:
		expand("{sample}_results/database/hhblits/a3m/{clusterId}.a3m",clusterId=config["clusters"], sample = SAMPLE_NAMES),
		expand("{sample}_results/database/hhblits/hhm/{clusterId}.hhm",clusterId=config["clusters"], sample = SAMPLE_NAMES)
	output:
		"{sample}_results/database/hhblits/db_a3m.ffdata",
		"{sample}_results/database/hhblits/db_a3m.ffindex",
		"{sample}_results/database/hhblits/db_hhm.ffdata",
		"{sample}_results/database/hhblits/db_hhm.ffindex",
		"{sample}_results/database/hhblits/db_cs219.ffdata",
		"{sample}_results/database/hhblits/db_cs219.ffindex",
		touch("{sample}_results/database/hhblits/db")
	run:
		shell("/home/meheurap/snakemake/proteinCluster/makehhblitsdb.py {output}")

rule runhhblits:
	input:
		hhm = "{sample}_results/database/hhblits/hhm/{clusterId}.hhm",
		db = "{sample}_results/database/hhblits/db"
	output:
		"{sample}_results/database/hhblits/hhr/{clusterId}.hhr"
	run:
		shell("hhblits -i {input.hhm} -o {output} -d {input.db}  -v 0 -p 50 -E 0.001 -z 1 -Z 32000 -B 0 -b 0 -n 2 -cpu 1")

rule parsehhblitsresults:
	input:
		expand("{sample}_results/database/hhblits/hhr/{clusterId}.hhr",clusterId=config["clusters"], sample = SAMPLE_NAMES)
	output:
		"{sample}_results/database/hhblits/allvsall.hhr"
	run:
		shell("/home/meheurap/snakemake/proteinCluster/parsingHhblitsResults.py {output}")