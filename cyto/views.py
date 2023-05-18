from django.http import HttpResponse, JsonResponse
from django.template import loader
import networkx as nx
import obonet
import time
import re
import os
import sys


def cyto_json(request, keyword, filter_n):
    # print("hi")
    # if request.method == 'POST':
    #     print(request.FILES)
    #     csv_files = request.FILES.getlist('de_csv_files')
    #     # Process the CSV files here
    #     print(len(csv_files))
    #     for csv_file in csv_files:
    #         print(type(csv_file))
    #         print(csv_file)
    #         print(csv_file.name)
    # print("hi")
    # print(keyword)
    # print(filter_n)
    # sys.stdout.flush()
    # Save the OBO formatted GO network from a url to a obonet graph
    url = 'http://purl.obolibrary.org/obo/go.obo'
    module_dir = os.path.dirname(__file__)  # get current directory
    file_path = os.path.join(module_dir, 'static/go.obo')
    url = file_path
    graph = obonet.read_obo(url)
    start = time.perf_counter()
    pagerank_scores = nx.pagerank(nx.DiGraph(graph))
    graph = graph.reverse(copy=False)
    print(f'{time.perf_counter()-start}')
    # From https://stackoverflow.com/a/45352658

    def convert2cytoscapeJSON(G):
        # load all nodes into nodes array
        # get top ten PageRank nodes
        aspect_to_color = {'molecular_function': '#3B185F',
                           'biological_process': '#A12568', 'cellular_component': '#539165'}
        top_go_terms = sorted(
            list(G.nodes()), key=lambda node: pagerank_scores[node], reverse=True)[:filter_n]
      #   print(top_go_terms)
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
            nx["data"]["name"] = G.nodes[node]['name']
            nx["data"]["aspect"] = G.nodes[node]['namespace']
            nx["data"]["definition"] = G.nodes[node]['def']
            nx["data"]["pagerank"] = pagerank_scores[node]
            nx["data"]["size"] = max(
                ((9 * ((pagerank_scores[node]/max_pagerank))+4)**2, 8))
            nx["data"]["color"] = aspect_to_color[nx['data']['aspect']]
          #   nx["data"]["color"] = '#%02x%02x%02x' % (180-round((pagerank_scores[node]/max_pagerank)*180), 180-round((pagerank_scores[node]/max_pagerank)*180), 255)
            final["nodes"].append(nx.copy())
        # load all edges to edges array
        for edge in G.edges():
            nx = {}
            nx["data"] = {}
            nx["data"]["id"] = edge[0]+edge[1]
            nx["data"]["source"] = edge[0]
            nx["data"]["target"] = edge[1]
            final["edges"].append(nx)
        return final

    # KWNameQuery takes a user input keyword and searches all GO Terms for that keyword.
    # Hits are stored in a list called keyword_nodes

    def KWNameQuery(userinput):
        name_nodes = []
        for node in graph.nodes:
            if re.search(rf"\b{userinput}\b", graph.nodes[node]['name'], re.IGNORECASE):
                name_nodes.append(node)
        return name_nodes

    # KWDefQuery takes a user input keyword and searches all GO Terms definitions for that keyword.
    # Hits are stored in a list called def_nodes

    def KWDefQuery(userinput):
        def_nodes = []
        for node in graph.nodes:
            if re.search(rf"\b{userinput}\b", graph.nodes[node]['def'], re.IGNORECASE):
                def_nodes.append(node)
        return def_nodes

    # ConcatKeywordHits takes in the userinput, calls KWNameQuery and KWDefQuery seperately with the same input
    # then uses set to return only unique hits

    def ConcatKeywordHits(userinput):
        return set(KWNameQuery(userinput) + KWDefQuery(userinput))
    nodes = ConcatKeywordHits(keyword)
    if len(nodes) == 0:
        return JsonResponse({})
    endConcat = time.perf_counter()
    print(f"Concat took {endConcat-start:0.4f}")
    # for node in list(nodes):
    # add the superterms to the nodes in the subgraph
    # nodes.update(nx.ancestors(graph, node))
    subgraph = graph.subgraph(nodes)
    endConcat = time.perf_counter()
    print(f"subgraph took {endConcat-start:0.4f}")
    graph_json = convert2cytoscapeJSON(subgraph)
    endConcat = time.perf_counter()
    print(f"convert2cyto took {endConcat-start:0.4f}")
    print("DONE")
    sys.stdout.flush()
    return JsonResponse(graph_json)
