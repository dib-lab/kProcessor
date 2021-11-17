

urls=dict([ l.strip().split("\t") for l in  open("urls")])


rule all:
    input:	urls.keys()


rule download:
    output: "{sample}"

    params:  url=lambda wildcards: urls[f"{wildcards.sample}"]
    conda: "env.yaml"
    log: "{sample}.download.log"
    shell:
       	"""
		folderName=$( basename {wildcards.sample} )
#		echo $folderName
#		mkdir -p  $folderName	
	        axel -n 10 -o {output} {params.url} &> {log}
        """


	 