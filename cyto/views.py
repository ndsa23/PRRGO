from django.http import HttpResponse, JsonResponse
import networkx as nx
import obonet
import time
import re
import os
import pandas as pd
import sys
from io import StringIO
import functools as ft
import random
import numpy as np
import urllib
from .utils import *
def mapGOTerms(request):
    # set random seeds for reproducibility
    random.seed(2023)
    np.random.seed(2023)
    if request.method == 'POST':
        # Save the OBO formatted GO network from a url to a obonet graph
        keyword = request.POST['keyword-input']
        filter_n = int(request.POST['filter-n'])
        url = 'http://purl.obolibrary.org/obo/go.obo'
        module_dir = os.path.dirname(__file__)  # get current directory
        returnFile = request.POST['isFile'] == 'True'
        graph = obonet.read_obo(url)
        # start = time.perf_counter()
        pagerank_scores = nx.pagerank(nx.DiGraph(graph))
        
        graph = graph.reverse(copy=False)

        # Check if GPI file exists
        gpi_filepath = os.path.join(module_dir, 'static/goa_human.gpi.gz')
        if not os.path.isfile(gpi_filepath):
            urllib.request.urlretrieve("http://current.geneontology.org/annotations/goa_human.gpi.gz", gpi_filepath)

        # Read the GPI file into a Pandas dataframe
        GPIdf = pd.read_csv(gpi_filepath, compression="gzip", sep='\t', comment='!', header=None,
                 names=['DB', 'DB_Object_ID', 'DB_Object_Symbol', 'DB_Object_Name', 'DB_Object_Synonym',
                        'DB_Object_Type', 'Taxon', 'Parent_Object_ID','DB_Xrefs', 'Properties'])

        # Check if GPAD file exists
        gpad_filepath = os.path.join(module_dir, 'static/goa_human.gpad.gz')
        if not os.path.isfile(gpad_filepath):
            urllib.request.urlretrieve("http://current.geneontology.org/annotations/goa_human.gpad.gz", gpad_filepath)

        # Read the GPI file into a Pandas dataframe
        GPADdf = pd.read_csv(os.path.join(module_dir, 'static/goa_human.gpad.gz'), compression="gzip", sep='\t', comment='!', header=None,
                            names=['DB', 'DB_Object_ID', 'Qualifier', 'GO_ID', 'DB:Reference', 'Evidence_Code',
                                    'WithOrFrom', 'Taxon_Intx', 'Date', 'Assigned_by', 'Annotation_Ext', 'Annotation_Prop'],
                            low_memory=False)

        # Combine annotation file information into a df merged on Uniprot IDs
        UPmerged_df = pd.merge(GPADdf, GPIdf, on='DB_Object_ID', how='inner')
        
        # Make a dictionary where these UP IDs are converted to their gene names
        # go_id_to_HGNC = UPmerged_df.groupby('GO_ID')['DB_Object_Symbol'].apply(lambda x: list(set(x))).to_dict()

        # print(f'{time.perf_counter()-start}')
        # sys.stdout.flush() 
        
        nodes = ConcatKeywordHits(keyword, graph)
        if len(nodes) == 0:
            return JsonResponse({})
        # endConcat = time.perf_counter()
        # print(f"Concat took {endConcat-start:0.4f}")
        # for node in list(nodes):
        # add the superterms to the nodes in the subgraph
        # nodes.update(nx.ancestors(graph, node))
        files = request.FILES.getlist('de-csv-files')
        if not returnFile:
            subgraph = graph.subgraph(nodes)
        # endConcat = time.perf_counter()
        # print(f"subgraph took {endConcat-start:0.4f}")
            graph_json = convert2cytoscapeJSON(subgraph, pagerank_scores, filter_n)
        # endConcat = time.perf_counter()
        # print(f"convert2cyto took {endConcat-start:0.4f}")
            print("DONE")
            mapped_degs = DEGassociator(files, nodes, UPmerged_df)
        # sys.stdout.flush()
            return JsonResponse({'graph': graph_json, 'go_data': mapped_degs})
        else:
            mapped_degs = fileDEGassociator(files, nodes, graph, pagerank_scores, UPmerged_df)
            return output_to_csv(mapped_degs, keyword)


