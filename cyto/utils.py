import pandas as pd
import functools as ft
import re
from networkx import MultiDiGraph
from io import StringIO
import django
from django.http import HttpResponse


def convert2cytoscapeJSON(
    G: MultiDiGraph, pagerank_scores: dict[str, float], filter_n: int
) -> dict:
    """
    load all nodes into nodes array and output top filter_n nodes to JSON
    
    arguments:
    --------
        G -- the Gene Ontology (GO) graph
        pagerank_scores -- the computed pagerank of the GO graph
        filter_n -- the top n terms to filter based on PageRank
    """

    aspect_to_color = {
        "molecular_function": "#3B185F",
        "biological_process": "#A12568",
        "cellular_component": "#539165",
    }
    g_nodes = list(G.nodes())
    top_go_terms = sorted(g_nodes, key=lambda node: pagerank_scores[node], reverse=True)
    top_go_terms = top_go_terms[:filter_n]
    G = G.subgraph(top_go_terms)
    max_pagerank = max([pagerank_scores[node] for node in G.nodes()])
    final = {}
    final["nodes"] = []
    final["edges"] = []
    for node in G.nodes():
        nx = {}
        nx["data"] = {}
        nx["data"]["id"] = node
        nx["data"]["label"] = node
        # id, name, aspect, definition
        nx["data"]["name"] = G.nodes[node]["name"]
        nx["data"]["aspect"] = G.nodes[node]["namespace"]
        nx["data"]["definition"] = G.nodes[node]["def"]
        nx["data"]["pagerank"] = pagerank_scores[node]
        nx["data"]["size"] = max(
            ((9 * ((pagerank_scores[node] / max_pagerank)) + 4) ** 2, 8)
        )
        nx["data"]["color"] = aspect_to_color[nx["data"]["aspect"]]
        final["nodes"].append(nx.copy())
    # load all edges to edges array
    for edge in G.edges():
        nx = {}
        nx["data"] = {}
        nx["data"]["id"] = edge[0] + edge[1]
        nx["data"]["source"] = edge[0]
        nx["data"]["target"] = edge[1]
        final["edges"].append(nx)
    return final


def HGNCsymbols_for_go_id(GO_ID: str, UPmerged_df: pd.DataFrame) -> list[str]:
    """return the list of genes as HGNCsymbols associated with GO_ID"""
    HGNCsymbols = (
        UPmerged_df.loc[UPmerged_df["GO_ID"] == GO_ID, "DB_Object_Symbol"]
        .unique()
        .tolist()
    )
    return HGNCsymbols


def DEGmerger(files: list[django.core.files.File]) -> pd.DataFrame:
    """merge all DEG csv/tsv files into one DataFrame"""
    dfs = []

    for file in files:
        # Check if the file is a .tsv or .csv file
        if file.name.endswith(".tsv"):
            delimiter = "\t"
        elif file.name.endswith(".csv"):
            delimiter = ","
        else:
            raise ValueError(f"Unknown file format: {file.name}")

        # Read the DEG file into a DataFrame
        file_str = file.read().decode("utf-8")
        df = pd.read_csv(StringIO(file_str), sep=delimiter)

        # Add a prefix to the column names with the filename
        df = df.add_prefix(f"{file.name[:-4]}_")
        symbol_column = [
            col
            for col in df.columns
            if col.endswith("Symbol") or col.endswith("Description")
        ][0]
        df = df.rename(columns={symbol_column: "Symbol"})
        dfs.append(df)

    # Merge all the DataFrames together into a single DataFrame
    merged_df = ft.reduce(lambda left, right: pd.merge(left, right, on="Symbol"), dfs)
    symbol_col = merged_df.pop("Symbol")
    merged_df.insert(0, "Symbol", symbol_col)
    return merged_df


def DEGassociator(
    files: list[django.core.files.File],
    ListofGOterms: list[str],
    UPmerged_df: pd.DataFrame,
) -> dict:
    """associates DEGs with GO terms in a dict"""
    # Merge all .tsv files in the specified directory
    merged_df = DEGmerger(files)

    # Create a dictionary to hold the results
    DEGsByGoTerm = {}

    for term in ListofGOterms:
        # Get the list of DEGs associated with this GO term
        goDEGs = HGNCsymbols_for_go_id(term, UPmerged_df)
        # Find DEGs that are associated with this GO term and also appear in the merged dataframe
        symbol_col = [
            col for col in merged_df.columns if col.lower().endswith("symbol")
        ][0]
        matched_genes = set(goDEGs).intersection(set(merged_df[symbol_col]))

        # If any matching genes were found, save the corresponding rows from the merged dataframe
        if matched_genes:
            DEGsByGoTerm[term] = (
                merged_df[merged_df[symbol_col].isin(matched_genes)]
                .fillna("NaN")
                .to_dict("split")
            )

    return DEGsByGoTerm


def KWNameQuery(userinput: str, graph: MultiDiGraph) -> list:
    """
    takes a user input keyword and searches all GO Terms for that keyword.
    Hits are stored in a list called keyword_nodes
    """
    name_nodes = []
    for node in graph.nodes:
        if re.search(rf"\b{userinput}\b", graph.nodes[node]["name"], re.IGNORECASE):
            name_nodes.append(node)
    return name_nodes


def KWDefQuery(userinput: str, graph: MultiDiGraph) -> list:
    """
    takes a user input keyword and searches all GO Terms definitions for that keyword.
    Hits are stored in a list called def_nodes
    """
    def_nodes = []
    for node in graph.nodes:
        if re.search(rf"\b{userinput}\b", graph.nodes[node]["def"], re.IGNORECASE):
            def_nodes.append(node)
    return def_nodes


def ConcatKeywordHits(userinput: str, graph: MultiDiGraph) -> set:
    """
    ConcatKeywordHits takes in the userinput,
    calls KWNameQuery and KWDefQuery seperately with the same input,
    then uses set to return only unique hits
    """
    return set(KWNameQuery(userinput, graph) + KWDefQuery(userinput, graph))


def fileDEGassociator(
    files: django.core.files.File,
    ListofGOterms: list[str],
    graph: MultiDiGraph,
    pagerank_scores: dict[str, float],
    UPmerged_df: pd.DataFrame,
) -> dict:
    """associate DEGs with GO terms, name, and PageRank for file output"""
    # Merge all .tsv files in the specified directory
    merged_df = DEGmerger(files)

    # Create a dictionary to hold the results
    DEGsByGoTerm = {}

    for term in ListofGOterms:
        goName = graph.nodes[term]["name"]
        # Make a one-row dataframe with GO ID and GO Name
        goNameAndTerm = pd.DataFrame(
            {"GO Term": [term], "Name": [goName], "PageRank": [pagerank_scores[term]]}
        )
        # Get the list of DEGs associated with this GO term
        goDEGs = HGNCsymbols_for_go_id(term, UPmerged_df)
        # Find DEGs that are associated with this GO term and also appear in the merged dataframe
        symbol_col = [col for col in merged_df.columns if "Symbol" in col][0]
        matched_genes = set(goDEGs).intersection(set(merged_df[symbol_col]))

        # If any matching genes were found, save the corresponding rows from the merged dataframe
        if matched_genes:
            DEGsByGoTerm[term] = merged_df[merged_df[symbol_col].isin(matched_genes)]

        # Check if the key exists in the dictionary before concatenating
        if term in DEGsByGoTerm:
            DEGsByGoTerm[term] = pd.concat([goNameAndTerm, DEGsByGoTerm[term]], axis=1)

    return DEGsByGoTerm


def output_to_csv(result: dict, file_prefix: str) -> HttpResponse:
    """output DEG GO mapped csv file"""
    # Create a Pandas DataFrame from the result dictionary
    df = pd.concat(result.values())

    # Generate the output file name with the prefix
    output_file = f"{file_prefix}_output.csv"

    # Write the DataFrame to a CSV file
    df.to_csv(output_file, index=False)

    # Print the output file name
    print(f"Results saved to {output_file}")
    with open(output_file, "rb") as deg_file:
        response = HttpResponse(deg_file.read(), content_type="text/csv")
        response["Content-Disposition"] = f'attachment; filename="{output_file}"'
        return response
