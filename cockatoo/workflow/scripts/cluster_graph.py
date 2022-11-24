########################
### cluster_graph.py ###
########################
# Author: Samuel Aroney

import pandas as pd
import re
import networkx as nx
from networkx.algorithms import community
from networkx.algorithms import components
from operator import itemgetter
import itertools

def process_community(c, graph):
    umbrella = set()
    if len(c) <= MAX_RECOVERY_SAMPLES and len(c) >= MIN_COASSEMBLY_SAMPLES:
        umbrella = c
    if len(c) > MAX_COASSEMBLY_SAMPLES: return pd.DataFrame(), umbrella
    if len(c) < MIN_COASSEMBLY_SAMPLES: return pd.DataFrame(), umbrella

    subgraphview = nx.subgraph_view(
        graph,
        filter_node = lambda x: x in c
        )

    total_size = sum([data["read_size"] for _,data in subgraphview.nodes(data=True)])
    if MAX_COASSEMBLY_SIZE and total_size > MAX_COASSEMBLY_SIZE: return pd.DataFrame(), umbrella

    edgeview = subgraphview.edges(data=True)
    total_weight = sum([data["weight"] for _,_,data in edgeview])
    total_targets = len(set([i for l in [data["target_ids"].split(",") for _,_,data in edgeview] for i in l]))

    df = pd.DataFrame({
        "samples": ",".join(c),
        "length": int(len(c)),
        "total_weight": total_weight,
        "total_targets": total_targets,
        "total_size": total_size,
        },
        index = [0])
    return df, umbrella

def pipeline(
        elusive_edges,
        read_size,
        MAX_COASSEMBLY_SIZE=None,
        MAX_COASSEMBLY_SAMPLES=2,
        MIN_COASSEMBLY_SAMPLES=2,
        MAX_RECOVERY_SAMPLES=20):

    if len(elusive_edges) == 0:
        return pd.DataFrame(columns=["samples", "length", "total_weight", "total_targets", "total_size", "recover_samples", "coassembly"])

    elusive_edges["sample1"] = elusive_edges["sample1"].apply(lambda x: re.sub("\.1$", "", x))
    elusive_edges["sample2"] = elusive_edges["sample2"].apply(lambda x: re.sub("\.1$", "", x))

    clusters = pd.DataFrame(columns=["samples", "length", "total_weight", "total_targets", "total_size", "recover_samples"])

    # Create weighted graph and cluster with Girvan-Newman algorithm, removing edges from lightest to heaviest
    graph = nx.from_pandas_edgelist(elusive_edges, source="sample1", target="sample2", edge_attr=True)
    node_attributes = {k: {"read_size": v} for k,v in zip(read_size["sample"], read_size["read_size"])}
    nx.set_node_attributes(graph, node_attributes)

    def lightest(graph):
        u, v, _ = min(graph.edges(data="weight"), key=itemgetter(2))
        return (u, v)

    comp = community.girvan_newman(graph, most_valuable_edge=lightest)
    first_community = [[c for c in components.connected_components(graph)]]
    umbrellas = []
    for iteration in itertools.chain(first_community, comp):
        communities = tuple(sorted(c) for c in iteration)
        for c in communities:
            if len(c) <= MAX_RECOVERY_SAMPLES and len(c) >= MIN_COASSEMBLY_SAMPLES and not c in umbrellas:
                umbrellas.append(c)
            if len(c) > MAX_COASSEMBLY_SAMPLES: continue
            if len(c) < MIN_COASSEMBLY_SAMPLES: continue

            subgraphview = nx.subgraph_view(
                graph,
                filter_node = lambda x: x in c
                )

            total_size = sum([data["read_size"] for _,data in subgraphview.nodes(data=True)])
            if MAX_COASSEMBLY_SIZE and total_size > MAX_COASSEMBLY_SIZE: continue

            edgeview = subgraphview.edges(data=True)
            total_weight = sum([data["weight"] for _,_,data in edgeview])
            total_targets = len(set([i for l in [data["target_ids"].split(",") for _,_,data in edgeview] for i in l]))

            suitable_indices = [i for i,u in enumerate(umbrellas) if all([s in u for s in c])]
            suitable_umbrellas = [umbrellas[i] for i in suitable_indices]
            best_umbrella = suitable_umbrellas[0]

            df = pd.DataFrame({
                "samples": ",".join(c),
                "length": int(len(c)),
                "total_weight": total_weight,
                "total_targets": total_targets,
                "total_size": total_size,
                "recover_samples": ",".join(best_umbrella)
                },
                index = [0])
            clusters = pd.concat([clusters, df], ignore_index=True)

    clusters.drop_duplicates(inplace=True)

    def find_top_samples(row):
        samples = row["samples"].split(",")
        recover_samples = row["recover_samples"].split(",")
        if len(recover_samples) >= MAX_RECOVERY_SAMPLES:
            return ",".join(recover_samples)

        relevant_edges = elusive_edges[(elusive_edges["sample1"].isin(samples)) != (elusive_edges["sample2"].isin(samples))].copy()
        if len(relevant_edges) == 0:
            return ",".join(recover_samples)

        relevant_edges["other_sample"] = relevant_edges.apply(lambda x: x["sample1"] if x["sample1"] not in samples else x["sample2"], axis=1)
        relevant_edges = relevant_edges[~relevant_edges["other_sample"].isin(recover_samples)]

        if MAX_COASSEMBLY_SAMPLES > 1:
            relevant_targets = elusive_edges[(elusive_edges["sample1"].isin(samples)) & (elusive_edges["sample2"].isin(samples))].copy()
            relevant_targets = set([i for l in [l.split(",") for l in relevant_targets["target_ids"].to_list()] for i in l])
            relevant_edges["target_ids"] = relevant_edges["target_ids"].apply(lambda x: [i for i in x.split(",") if i in relevant_targets])
        else:
            relevant_edges["target_ids"] = relevant_edges["target_ids"].apply(lambda x: x.split(","))

        relevant_samples = relevant_edges.groupby("other_sample").agg({"target_ids": lambda x: set(itertools.chain(*x))}).reset_index()
        relevant_samples["length"] = relevant_samples["target_ids"].apply(lambda x: len(x))
        relevant_samples.sort_values(by="length", ascending=False, inplace=True)

        recover_samples += relevant_samples.head(MAX_RECOVERY_SAMPLES - len(recover_samples))["other_sample"].to_list()
        return ",".join(sorted(recover_samples))

    clusters["recover_samples"] = clusters[["samples", "recover_samples"]].apply(find_top_samples, axis=1)

    clusters.sort_values(by=["total_targets", "samples"], ascending=False, inplace=True)
    clusters.reset_index(drop=True, inplace=True)
    clusters["coassembly"] = (clusters.reset_index()["index"].apply(lambda x: "coassembly_" + str(x)))

    return clusters

if __name__ == "__main__":
    MAX_COASSEMBLY_SIZE = snakemake.params.max_coassembly_size * 10**9 if snakemake.params.max_coassembly_size else None
    MAX_COASSEMBLY_SAMPLES = snakemake.params.max_coassembly_samples
    MIN_COASSEMBLY_SAMPLES = snakemake.params.num_coassembly_samples
    MAX_RECOVERY_SAMPLES = snakemake.params.max_recovery_samples
    elusive_edges_path = snakemake.input.elusive_edges
    read_size_path = snakemake.input.read_size
    elusive_clusters_path = snakemake.output.elusive_clusters


    elusive_edges = pd.read_csv(elusive_edges_path, sep="\t", dtype={"target_ids": "string"})
    read_size = pd.read_csv(read_size_path, names = ["sample", "read_size"])

    clusters = pipeline(
        elusive_edges,
        read_size,
        MAX_COASSEMBLY_SIZE=MAX_COASSEMBLY_SIZE,
        MAX_COASSEMBLY_SAMPLES=MAX_COASSEMBLY_SAMPLES,
        MIN_COASSEMBLY_SAMPLES=MIN_COASSEMBLY_SAMPLES,
        MAX_RECOVERY_SAMPLES=MAX_RECOVERY_SAMPLES
        )
    clusters.to_csv(elusive_clusters_path, sep="\t", index=False)
