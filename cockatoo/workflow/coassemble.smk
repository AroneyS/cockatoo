#############
### Setup ###
#############
ruleorder: single_assembly > query_processing > singlem_appraise

import os

output_dir = os.path.abspath("coassemble")
logs_dir = output_dir + "/logs"

mapped_reads_1 = {read: output_dir + f"/mapping/{read}_unmapped.1.fq.gz" for read in config["reads_1"]}
mapped_reads_2 = {read: output_dir + f"/mapping/{read}_unmapped.2.fq.gz" for read in config["reads_2"]}

def get_reads(wildcards):
    version = wildcards.version
    if version == "":
        return config["reads_1"].values()
    elif version == "unmapped_":
        return mapped_reads_1.values()
    else:
        raise ValueError("Version should be empty or unmapped")

def get_cat(wildcards):
    return "zcat" if [r for r in get_reads(wildcards)][0].endswith(".gz") else "cat"

rule all:
    input:
        output_dir + "/target/elusive_clusters.tsv",
        output_dir + "/commands/coassemble_commands.sh",
        output_dir + "/commands/recover_commands.sh",
        output_dir + "/summary.tsv",

rule summary:
    input:
        elusive_clusters=output_dir + "/target/elusive_clusters.tsv",
        read_size=output_dir + "/unmapped_read_size.csv" if config["assemble_unmapped"] else [],
    output:
        summary=output_dir + "/summary.tsv",
    script:
        "scripts/summarise_coassemblies.py"

#####################
### SingleM reads ###
#####################
rule singlem_pipe_reads:
    input:
        reads_1=lambda wildcards: config["reads_1"][wildcards.read],
        reads_2=lambda wildcards: config["reads_2"][wildcards.read] if config["reads_2"] else [],
    output:
        output_dir + "/pipe/{read}_read.otu_table.tsv"
    log:
        logs_dir + "/pipe/{read}_read.log"
    params:
        singlem_metapackage=config["singlem_metapackage"],
        reverse_arg="--reverse " if config["reads_2"] else "",
    conda:
        "env/singlem.yml"
    shell:
        "singlem pipe "
        "--forward {input.reads_1} "
        "{params.reverse_arg} {input.reads_2} "
        "--otu-table {output} "
        "--metapackage {params.singlem_metapackage} "
        "&> {log}"

#######################
### SingleM genomes ###
#######################
rule genome_transcripts:
    input:
        lambda wildcards: config["genomes"][wildcards.genome],
    output:
        output_dir + "/transcripts/{genome}_protein.fna"
    log:
        logs_dir + "/transcripts/{genome}_protein.log"
    params:
        prodigal_meta = "-p meta" if config["prodigal_meta"] else ""
    conda:
        "env/prodigal.yml"
    shell:
        "prodigal "
        "-i {input} "
        "-d {output} "
        "{params.prodigal_meta} "
        "&> {log} "

rule singlem_pipe_genomes:
    input:
        output_dir + "/transcripts/{genome}_protein.fna"
    output:
        output_dir + "/pipe/{genome}_bin.otu_table.tsv"
    log:
        logs_dir + "/pipe/{genome}_bin.log"
    params:
        singlem_metapackage=config["singlem_metapackage"]
    conda:
        "env/singlem.yml"
    shell:
        "singlem pipe "
        "--forward {input} "
        "--otu-table {output} "
        "--metapackage {params.singlem_metapackage} "
        "&> {log}"

rule singlem_summarise_genomes:
    input:
        expand(output_dir + "/pipe/{genome}_bin.otu_table.tsv", genome=config["genomes"])
    output:
        output_dir + "/summarise/bins_summarised.otu_table.tsv"
    log:
        logs_dir + "/summarise/genomes.log"
    params:
        singlem_metapackage=config["singlem_metapackage"]
    conda:
        "env/singlem.yml"
    shell:
        "singlem summarise "
        "--input-otu-tables {input} "
        "--output-otu-table {output} "
        "--exclude-off-target-hits "
        "--metapackage {params.singlem_metapackage} "
        "&> {log}"

########################
### SingleM appraise ###
########################
rule singlem_appraise:
    input:
        reads=expand(output_dir + "/pipe/{read}_read.otu_table.tsv", read=config["reads_1"]),
        bins=output_dir + "/summarise/bins_summarised.otu_table.tsv",
    output:
        unbinned=output_dir + "/appraise/unbinned.otu_table.tsv",
        binned=output_dir + "/appraise/binned.otu_table.tsv",
    log:
        logs_dir + "/appraise/appraise.log"
    params:
        sequence_identity=config["appraise_sequence_identity"],
        singlem_metapackage=config["singlem_metapackage"],
    conda:
        "env/singlem.yml"
    shell:
        "singlem appraise "
        "--metagenome-otu-tables {input.reads} "
        "--genome-otu-tables {input.bins} "
        "--metapackage {params.singlem_metapackage} "
        "--output-unaccounted-for-otu-table {output.unbinned} "
        "--output-binned-otu-table {output.binned} "
        "--imperfect "
        "--sequence-identity {params.sequence_identity} "
        "--output-found-in "
        "&> {log}"

###################################
### SingleM query (alternative) ###
###################################
rule query_processing:
    input:
        pipe_reads=expand(output_dir + "/pipe/{read}_read.otu_table.tsv", read=config["reads_1"]),
        query_reads=expand(output_dir + "/query/{read}_query.otu_table.tsv", read=config["reads_1"]),
    output:
        unbinned=output_dir + "/appraise/unbinned.otu_table.tsv",
        binned=output_dir + "/appraise/binned.otu_table.tsv",
    log:
        logs_dir + "/query/processing.log"
    params:
        sequence_identity=config["appraise_sequence_identity"],
        window_size=60,
    script:
        "scripts/query_processing.py"

#####################################
### Single assembly (alternative) ###
#####################################
rule single_assembly:
    input:
        reads=expand(output_dir + "/pipe/{read}_read.otu_table.tsv", read=config["reads_1"]),
    output:
        unbinned=output_dir + "/appraise/unbinned.otu_table.tsv" if config["single_assembly"] else [],
        binned=output_dir + "/appraise/binned.otu_table.tsv" if config["single_assembly"] else [],
    log:
        logs_dir + "/appraise/appraise.log"
    script:
        "scripts/single_assembly.py"

######################
### Target elusive ###
######################
rule count_bp_reads:
    input:
        get_reads
    output:
        output_dir + "/{version,.*}read_size.csv"
    params:
        names=list(config["reads_1"].keys()),
        cat=get_cat,
    threads:
        8
    shell:
        "parallel -k -j {threads} "
        "echo -n {{1}}, '&&' "
        "{params.cat} {{2}} '|' sed -n 2~4p '|' tr -d '\"\n\"' '|' wc -m "
        "::: {params.names} :::+ {input} "
        "> {output}"

rule target_elusive:
    input:
        unbinned=output_dir + "/appraise/unbinned.otu_table.tsv"
    output:
        output_edges=output_dir + "/target/elusive_edges.tsv",
        output_targets=output_dir + "/target/targets.tsv",
    log:
        logs_dir + "/elusive_targets.log"
    params:
        min_coassembly_coverage=config["min_coassembly_coverage"],
        taxa_of_interest=config["taxa_of_interest"],
    script:
        "scripts/target_elusive.py"

rule cluster_graph:
    input:
        elusive_edges=output_dir + "/target/elusive_edges.tsv",
        read_size=output_dir + "/read_size.csv",
    output:
        elusive_clusters=output_dir + "/target/elusive_clusters.tsv"
    log:
        logs_dir + "/cluster_graph.log"
    params:
        max_coassembly_size=config["max_coassembly_size"],
        num_coassembly_samples=config["num_coassembly_samples"],
        max_coassembly_samples=config["max_coassembly_samples"],
        max_recovery_samples=config["max_recovery_samples"],
    script:
        "scripts/cluster_graph.py"

#####################################
### Map reads to matching genomes ###
#####################################
rule collect_genomes:
    input:
        appraise_binned=output_dir + "/appraise/binned.otu_table.tsv",
    output:
        temp(output_dir + "/mapping/{read}_reference.fna"),
    params:
        genomes=config["genomes"],
        sample="{read}",
    script:
        "scripts/collect_reference_bins.py"

rule map_reads:
    input:
        reads_1=lambda wildcards: config["reads_1"][wildcards.read],
        reads_2=lambda wildcards: config["reads_2"][wildcards.read] if config["reads_2"] else [],
        genomes=output_dir + "/mapping/{read}_reference.fna",
    output:
        dir=temp(directory(output_dir + "/mapping/{read}_coverm")),
        unmapped_bam=temp(output_dir + "/mapping/{read}_unmapped.bam"),
        reads_1=output_dir + "/mapping/{read}_unmapped.1.fq.gz",
        reads_2=output_dir + "/mapping/{read}_unmapped.2.fq.gz" if config["reads_2"] else "",
    log:
        logs_dir + "/mapping/{read}.log",
    params:
        genomes="{read}_reference.fna",
        reads_1=lambda wildcards: os.path.basename(config["reads_1"][wildcards.read]),
        coverm_log=logs_dir + "/mapping/{read}_coverm.log",
        view_log=logs_dir + "/mapping/{read}_view.log",
        fastq_log=logs_dir + "/mapping/{read}_fastq.log",
        reverse_arg="-2" if config["reads_2"] else "",
    threads:
        8
    conda:
        "env/coverm.yml"
    shell:
         "coverm make "
         "-r {input.genomes} "
         "-1 {input.reads_1} "
         "{params.reverse_arg} {input.reads_2} "
         "-o {output.dir} "
         "-t {threads} "
         "&> {params.coverm_log} "
         "&& "
         "samtools view "
         "-@ $(({threads} - 1)) "
         "-b -f12 {output.dir}/{params.genomes}.{params.reads_1}.bam "
         "2> {params.view_log} "
         "> {output.unmapped_bam} "
         "&& "
         "samtools fastq "
         "-@ $(({threads} - 1)) "
         "{output.unmapped_bam} "
         "-1 {output.reads_1} "
         "{params.reverse_arg} {output.reads_2} "
         "-0 /dev/null "
         "-s /dev/null "
         "-n "
         "&> {params.fastq_log} "

rule finish_mapping:
    input:
        mapped_reads_1.values(),
        mapped_reads_2.values(),
    output:
        output_dir + "/mapping/done"
    shell:
        "touch {output}"

##############################
### Create Aviary commands ###
##############################
rule aviary_commands:
    input:
        output_dir + "/mapping/done" if config["assemble_unmapped"] else [],
        elusive_clusters=output_dir + "/target/elusive_clusters.tsv",
    output:
        coassemble_commands=output_dir + "/commands/coassemble_commands.sh",
        recover_commands=output_dir + "/commands/recover_commands.sh"
    threads:
        64
    params:
        reads_1 = mapped_reads_1 if config["assemble_unmapped"] else config["reads_1"],
        reads_2 = mapped_reads_2 if config["assemble_unmapped"] else config["reads_2"],
        dir=output_dir,
        memory=config["aviary_memory"],
        threads=config["aviary_threads"],
    log:
        logs_dir + "/aviary_commands.log"
    script:
        "scripts/aviary_commands.py"
