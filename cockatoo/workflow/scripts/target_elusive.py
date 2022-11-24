#########################
### target_elusive.py ###
#########################
# Author: Samuel Aroney

import os
import pandas as pd
from itertools import combinations

def pipeline(
    unbinned,
    MIN_COASSEMBLY_COVERAGE = 10,
    TAXA_OF_INTEREST = ""):

    if len(unbinned) == 0:
        return unbinned.rename(columns={"found_in": "target"}), pd.DataFrame(columns=["taxa_group", "weight", "target_ids", "sample1", "sample2"])

    TAXA_LEVEL_SEP = "; "
    taxa_levels = {
        "d": 1,
        "p": 2,
        "c": 3,
        "o": 4,
        "f": 5,
        "g": 6,
        "s": 7,
    }
    try:
        TAXA_LEVEL_OF_INTEREST = taxa_levels[TAXA_OF_INTEREST[0]]
    except IndexError:
        TAXA_LEVEL_OF_INTEREST = 0

    # Group hits by sequence within genes and number to form targets
    unbinned.drop(["found_in"], axis=1, errors="ignore", inplace=True)
    unbinned["target"] = unbinned.groupby(["gene", "sequence"]).ngroup()
    unbinned["target"] = unbinned["target"].astype(str)

    taxonomy = unbinned.groupby("target")["taxonomy"].first()

    # Within each target cluster, find pairs of samples with combined coverage > MIN_COVERAGE
    def find_pairs(df):
        if len(df) < 2: return []

        def sum_coverages(sample1, sample2, df):
            print(df.loc[sample1]["coverage"])
            print(df.loc[sample2]["coverage"])
            return df.loc[sample1]["coverage"] + df.loc[sample2]["coverage"]

        edges = [frozenset([df.loc[sample1]["sample"], df.loc[sample2]["sample"]])
                for (sample1, sample2) in combinations(df.index, 2)
                if sum_coverages(sample1, sample2, df) > MIN_COASSEMBLY_COVERAGE]

        return edges

    sample_pairs = (unbinned
        .groupby("target")
        .apply(find_pairs)
        .to_frame("sample_pairs")
        .explode("sample_pairs")
        .dropna(subset = "sample_pairs")
        .join(taxonomy)
        .reset_index()
        )

    # Select group of interest (corresponding taxa level in brackets)
    # e.g. Root (0) or d__Bacteria (1) or p__Acidobacteriota (2) or c__Acidobacteriae (3) or o__Bryobacterales (4)
    # or f__Bryobacteraceae (5) or g__Solibacter (6) or s__Solibacter sp003136215 (7)

    # Or select level of interest (corresponding taxa level in brackets) - produce a graph for each group in level
    # e.g. Root (0) or Domain (1) or Phylum (2) or Class (3) or Order (4) or Family (5) or Genus (6) or Species (7)
    def get_taxa_group(taxonomy, taxa_level = TAXA_LEVEL_OF_INTEREST, taxa_level_sep = TAXA_LEVEL_SEP):
        try:
            return taxonomy.split(taxa_level_sep)[taxa_level]
        except IndexError:
            return None

    sample_pairs["taxa_group"] = sample_pairs["taxonomy"].apply(get_taxa_group)
    sample_pairs = sample_pairs.dropna(subset = "taxa_group")
    if len(TAXA_OF_INTEREST) > 1:
        sample_pairs = sample_pairs[sample_pairs["taxa_group"] == TAXA_OF_INTEREST]

    # Create weighted graph with nodes as samples and edges weighted by the number of clusters supporting that co-assembly (networkx)
    # Output sparse matrix with taxa_group/sample1/sample2/edge weight (number of supporting clusters with coverage > threshold)
    sparse_edges = (sample_pairs
        .groupby(["taxa_group", "sample_pairs"])["target"]
        .agg(["count", lambda x: ",".join(sorted(x))])
        .reset_index()
        )
    sparse_edges.rename(columns = {"count": "weight", "<lambda_0>": "target_ids"}, inplace=True)

    if len(sparse_edges) == 0:
        sparse_edges[["sample1", "sample2"]] = None
    else:
        sparse_edges[["sample1", "sample2"]] = sparse_edges.sample_pairs.apply(lambda x: pd.Series([i for i in x]))

    sparse_edges.drop("sample_pairs", axis = 1, inplace = True)
    return unbinned, sparse_edges

if __name__ == "__main__":
    os.environ["MODIN_CPUS"] = str(snakemake.threads)
    import modin.pandas as pd

    MIN_COASSEMBLY_COVERAGE = snakemake.params.min_coassembly_coverage
    TAXA_OF_INTEREST = snakemake.params.taxa_of_interest
    unbinned_path = snakemake.input.unbinned
    targets_path = snakemake.output.output_targets
    edges_path = snakemake.output.output_edges

    unbinned = pd.read_csv(unbinned_path, sep="\t")

    targets, edges = pipeline(
        unbinned,
        MIN_COASSEMBLY_COVERAGE=MIN_COASSEMBLY_COVERAGE,
        TAXA_OF_INTEREST=TAXA_OF_INTEREST
        )
    targets.to_csv(targets_path, sep="\t", index=False)
    edges.to_csv(edges_path, sep="\t", index=False)
