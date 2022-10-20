##############################
### aviary_commands.py ###
##############################
# Author: Samuel Aroney

import pandas as pd

#################
### Functions ###
#################
def produce_command(row):
    coassembly_sample_names = row["samples"]
    coassembly_samples_1 = [snakemake.params.reads_1[sample] for sample in coassembly_sample_names]
    coassembly_samples_2= [snakemake.params.reads_2[sample] for sample in coassembly_sample_names] if snakemake.params.reads_2 else []

    all_samples_names = row["recover_samples"]
    all_samples_1 = [snakemake.params.reads_1[sample] for sample in all_samples_names]
    all_samples_2 = [snakemake.params.reads_2[sample] for sample in all_samples_names] if snakemake.params.reads_2 else []

    reverse_arg = "-2" if snakemake.params.reads_2 else ""

    coassembly = row["coassembly"]
    output_dir = snakemake.params.dir + "/coassemble"
    threads = snakemake.params.threads
    memory = snakemake.params.memory

    aviary_assemble = ("aviary assemble "
    f"-1 {' '.join(coassembly_samples_1)} "
    f"{reverse_arg} {' '.join(coassembly_samples_2)} "
    f"--output {output_dir}/{coassembly}/assemble "
    f"-n {threads} "
    f"-m {memory} "
    f"&> {output_dir}/logs/{coassembly}_assemble.log ")

    aviary_recover = ("aviary recover "
    f"--assembly {output_dir}/{coassembly}/assemble/assembly/final_contigs.fasta "
    f"-1 {' '.join(all_samples_1)} "
    f"{reverse_arg} {' '.join(all_samples_2)} "
    f"--output {output_dir}/{coassembly}/recover "
    f"-n {threads} "
    f"-m {memory} "
    f"&> {output_dir}/logs/{coassembly}_recover.log ")

    return pd.Series([aviary_assemble, aviary_recover])

################
### Pipeline ###
################
coassemblies = pd.read_csv(snakemake.input.elusive_clusters, sep="\t")

if len(coassemblies.index) == 0:
    with open(snakemake.output.coassemble_commands, "w") as f:
        pass
    with open(snakemake.output.recover_commands, "w") as f:
        pass
    print("No coassemblies to perform")
    exit(0)

# Produce coassembly_commands
coassemblies["samples"] = coassemblies["samples"].apply(lambda item: item.split(","))
coassemblies["recover_samples"] = coassemblies["recover_samples"].apply(lambda item: item.split(","))
coassemblies["coassembly"] = coassemblies.reset_index()["index"].apply(lambda x: "coassembly_" + str(x))
coassemblies[["assemble", "recover"]] = coassemblies.apply(produce_command, axis=1)

coassemblies["assemble"].to_csv(snakemake.output.coassemble_commands, sep="\t", index=False, header=False)
coassemblies["recover"].to_csv(snakemake.output.recover_commands, sep="\t", index=False, header=False)
