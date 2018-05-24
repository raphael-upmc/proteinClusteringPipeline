import sys,os,re

directory = config["directory"] 

rule all:
	input:
#		expand("{directory}/hhblits/hhm/{clusterId}.hhm",clusterId=config["clusters"],directory=config["directory"]),
		expand("{directory}/hhblits/hhr/{clusterId}.hhr",clusterId=config["clusters"],directory=config["directory"]),
#		expand("{directory}/hhblits/db",directory=config["directory"]),

rule mafft:
	input:
		"{directory}/subfamiliesFasta/{clusterId}.fa"
	output:
		"{directory}/hhblits/mafft/{clusterId}.mafft"
	run:
		shell("mafft --auto {input} > {output} 2>/dev/null")


rule hhfilter:
	input:
		"{directory}/hhblits/mafft/{clusterId}.mafft"
	output:
		"{directory}/hhblits/mafft/{clusterId}.a3m"
	run:
		print(config["clusters"])
		print(wildcards.clusterId)
		print(config["clusters"][wildcards.clusterId])		
		if int(config["clusters"][wildcards.clusterId]) > 100 :
		   shell("hhfilter -M 50 -diff 100 -i {input} -o {output}")
		   shell("reformat.pl a3m a3m {output} {output} -M 50 -r;")
		else:
		   shell("reformat.pl fas a3m {input} {output} -M 50 -r;")
		   		 		  

rule reformat_a3m:
	input:
		"{directory}/hhblits/mafft/{clusterId}.a3m"
	output:
		"{directory}/hhblits/mafft/{clusterId}_renamed.a3m"     
	run:
     		shell("/home/meheurap/proteinClusteringPipeline/scripts/addingNameToA3m.py {input} {output}")


rule add_ss:
	input:
		"{directory}/hhblits/mafft/{clusterId}_renamed.a3m"     
	output:
		"{directory}/hhblits/a3m/{clusterId}.a3m"
	run:
     		shell("addss.pl {input} {output} -a3m")

rule hhm:
	input:
		"{directory}/hhblits/a3m/{clusterId}.a3m"
	output:
		"{directory}/hhblits/hhm/{clusterId}.hhm"
	run:		
		shell("hhmake -add_cons -M 50 -i {input} -o {output}")

rule makehhblitsdb :
	input:
		expand("{directory}/hhblits/a3m/{clusterId}.a3m",clusterId=config["clusters"],directory=config["directory"]),
		expand("{directory}/hhblits/hhm/{clusterId}.hhm",clusterId=config["clusters"],directory=config["directory"])
	output:
		"{directory}/hhblits/db_a3m.ffdata",
		"{directory}/hhblits/db_a3m.ffindex",
		"{directory}/hhblits/db_hhm.ffdata",
		"{directory}/hhblits/db_hhm.ffindex",
		"{directory}/hhblits/db_cs219.ffdata",
		"{directory}/hhblits/db_cs219.ffindex",
		touch("{directory}/hhblits/db")
	run:
		shell("/home/meheurap/proteinClusteringPipeline/scripts/makehhblitsdb.py {output}")

rule runhhblits:
	input:
		hhm = "{directory}/hhblits/hhm/{clusterId}.hhm",
		db = "{directory}/hhblits/db"
	output:
		"{directory}/hhblits/hhr/{clusterId}.hhr"
	run:
		shell("hhblits -i {input.hhm} -o {output} -d {input.db}  -v 0 -p 50 -E 0.001 -z 1 -Z 32000 -B 0 -b 0 -n 2 -cpu 1")

