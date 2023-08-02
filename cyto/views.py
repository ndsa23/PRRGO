from django.http import HttpResponse, JsonResponse
from django.template import loader
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
# from .go_mapper import main

def cyto_json(request):
    # set random seeds for reproducibility
    random.seed(2023)
    np.random.seed(2023)
    if request.method == 'POST':
    # Save the OBO formatted GO network from a url to a obonet graph
        keyword = request.POST['keyword-input']
        filter_n = int(request.POST['filter-n'])
        sys.stdout.flush()
        url = 'http://purl.obolibrary.org/obo/go.obo'
        module_dir = os.path.dirname(__file__)  # get current directory
        # obo_file_path = os.path.join(module_dir, 'static/go.obo')
        # url = obo_file_path
        graph = obonet.read_obo(url)
        start = time.perf_counter()
        pagerank_scores = nx.pagerank(nx.DiGraph(graph))
        graph = graph.reverse(copy=False)
        # Read the GPAD file into a Pandas dataframe
        GPIdf = pd.read_csv(os.path.join(module_dir, 'static/goa_human.gpi.gz'), compression="gzip", sep='\t', comment='!', header=None,
                 names=['DB', 'DB_Object_ID', 'DB_Object_Symbol', 'DB_Object_Name', 'DB_Object_Synonym',
                        'DB_Object_Type', 'Taxon', 'Parent_Object_ID','DB_Xrefs', 'Properties'])

        # Read the GPI file into a Pandas dataframe
        GPADdf = pd.read_csv(os.path.join(module_dir, 'static/goa_human.gpad.gz'), compression="gzip", sep='\t', comment='!', header=None,
                            names=['DB', 'DB_Object_ID', 'Qualifier', 'GO_ID', 'DB:Reference', 'Evidence_Code',
                                    'WithOrFrom', 'Taxon_Intx', 'Date', 'Assigned_by', 'Annotation_Ext', 'Annotation_Prop'],
                            low_memory=False)

        #Combine annotation file information into a df merged on Uniprot IDs
        UPmerged_df = pd.merge(GPADdf, GPIdf, on='DB_Object_ID', how='inner')
        
        #Make a dictionary where these UP IDs are converted to their gene names
        go_id_to_HGNC = UPmerged_df.groupby('GO_ID')['DB_Object_Symbol'].apply(lambda x: list(set(x))).to_dict()

        print(f'{time.perf_counter()-start}')
        sys.stdout.flush() 

        # Based on https://stackoverflow.com/a/45352658
        def convert2cytoscapeJSON(G):
            # load all nodes into nodes array
            # get top ten PageRank nodes
            aspect_to_color = {'molecular_function': '#3B185F',
                            'biological_process': '#A12568', 'cellular_component': '#539165'}
            g_nodes = list(G.nodes())
            top_go_terms = sorted(
                g_nodes, key=lambda node: pagerank_scores[node], reverse=True)
            top_go_terms = top_go_terms[:filter_n]
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
        def HGNCsymbols_for_go_id(GO_ID):
            HGNCsymbols = UPmerged_df.loc[UPmerged_df['GO_ID'] == GO_ID, 'DB_Object_Symbol'].unique().tolist()
            return HGNCsymbols
        # accepts a list of Django Files
        def DEGmerger(files):
            dfs = []

            for file in files:
                # print(file)
                # print(file.read())
                # Check if the file is a .tsv or .csv file
                if file.name.endswith('.tsv'):
                    delimiter = '\t'
                elif file.name.endswith('.csv'):
                    delimiter = ','
                else:
                    raise ValueError(f"Unknown file format: {file.name}")

                # Read the DEG file into a DataFrame
                print(file.name)
                sys.stdout.flush()
                file_str = file.read().decode('utf-8')
                # print(file_str)
                df = pd.read_csv(StringIO(file_str), sep=delimiter)

                # Add a prefix to the column names with the filename
                df = df.add_prefix(f"{file.name[:-4]}_")
                print(df.columns)
                sys.stdout.flush()
                symbol_column = [col for col in df.columns if col.endswith("Symbol") or col.endswith("Description")][0]
                df = df.rename(columns={symbol_column: "Symbol"})
                sys.stdout.flush()
                dfs.append(df)

            # Merge all the DataFrames together into a single DataFrame
            merged_df = ft.reduce(lambda left, right: pd.merge(left, right, on='Symbol'), dfs)
            symbol_col = merged_df.pop('Symbol')
            merged_df.insert(0, 'Symbol', symbol_col)
            return merged_df

        # Function to associate DEGs with GO terms
        def DEGassociator(files, ListofGOterms):
            # Merge all .tsv files in the specified directory
            merged_df = DEGmerger(files)
            
            # Create a dictionary to hold the results
            DEGsByGoTerm = {}
            
            for term in ListofGOterms:
                goName = graph.nodes[term]['name']
                # Make a one-row dataframe with GO ID and GO Name
                goNameAndTerm = pd.DataFrame({'GO Term': [term], 'Name': [goName]})
                # Get the list of DEGs associated with this GO term
                goDEGs = HGNCsymbols_for_go_id(term)
                # Find DEGs that are associated with this GO term and also appear in the merged dataframe
                symbol_col = [col for col in merged_df.columns if col.lower().endswith("symbol")][0]
                matched_genes = set(goDEGs).intersection(set(merged_df[symbol_col]))
                
                #TODO: Double check if this is okay
                # If any matching genes were found, save the corresponding rows from the merged dataframe
                if matched_genes:
                    DEGsByGoTerm[term] = merged_df[merged_df[symbol_col].isin(matched_genes)].fillna('NaN').to_dict('split')
                    
                
            return DEGsByGoTerm
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
        files = request.FILES.getlist('de-csv-files')
        mapped_degs = DEGassociator(files, nodes)
        sys.stdout.flush()
        return JsonResponse({'graph': graph_json, 'go_data': mapped_degs})
def download(request):
    # set random seeds for reproducibility
    random.seed(2023)
    np.random.seed(2023)
    if request.method == 'POST':
        keyword = request.POST['keyword-input']
        url = 'http://purl.obolibrary.org/obo/go.obo'
        module_dir = os.path.dirname(__file__)  # get current directory
        graph = obonet.read_obo(url)
        start = time.perf_counter()
        pagerank_scores = nx.pagerank(nx.DiGraph(graph))
        graph = graph.reverse(copy=False)
        # Read the GPAD file into a Pandas dataframe
        GPIdf = pd.read_csv(os.path.join(module_dir, 'static/goa_human.gpi.gz'), compression="gzip", sep='\t', comment='!', header=None,
                 names=['DB', 'DB_Object_ID', 'DB_Object_Symbol', 'DB_Object_Name', 'DB_Object_Synonym',
                        'DB_Object_Type', 'Taxon', 'Parent_Object_ID','DB_Xrefs', 'Properties'])

        # Read the GPI file into a Pandas dataframe
        GPADdf = pd.read_csv(os.path.join(module_dir, 'static/goa_human.gpad.gz'), compression="gzip", sep='\t', comment='!', header=None,
                            names=['DB', 'DB_Object_ID', 'Qualifier', 'GO_ID', 'DB:Reference', 'Evidence_Code',
                                    'WithOrFrom', 'Taxon_Intx', 'Date', 'Assigned_by', 'Annotation_Ext', 'Annotation_Prop'],
                            low_memory=False)

        #Combine annotation file information into a df merged on Uniprot IDs
        UPmerged_df = pd.merge(GPADdf, GPIdf, on='DB_Object_ID', how='inner')
        

        print(f'{time.perf_counter()-start}')
        sys.stdout.flush() 

        def HGNCsymbols_for_go_id(GO_ID):
            HGNCsymbols = UPmerged_df.loc[UPmerged_df['GO_ID'] == GO_ID, 'DB_Object_Symbol'].unique().tolist()
            return HGNCsymbols
        # accepts a list of Django Files
        def DEGmerger(files):
            dfs = []

            for file in files:
                # Check if the file is a .tsv or .csv file
                if file.name.endswith('.tsv'):
                    delimiter = '\t'
                elif file.name.endswith('.csv'):
                    delimiter = ','
                else:
                    raise ValueError(f"Unknown file format: {file.name}")

                # Read the DEG file into a DataFrame
                file_str = file.read().decode('utf-8')
                # print(file_str)
                df = pd.read_csv(StringIO(file_str), sep=delimiter)

                # Add a prefix to the column names with the filename
                df = df.add_prefix(f"{file.name[:-4]}_")
                print(df.columns)
                sys.stdout.flush()
                symbol_column = [col for col in df.columns if col.endswith("Symbol") or col.endswith("Description")][0]
                df = df.rename(columns={symbol_column: "Symbol"})
                dfs.append(df)

            # Merge all the DataFrames together into a single DataFrame
            merged_df = ft.reduce(lambda left, right: pd.merge(left, right, on='Symbol'), dfs)
            symbol_col = merged_df.pop('Symbol')
            merged_df.insert(0, 'Symbol', symbol_col)
            return merged_df

        # Function to associate DEGs with GO terms
        def DEGassociator(files, ListofGOterms):
            # Merge all .tsv files in the specified directory
            merged_df = DEGmerger(files)
            
            # Create a dictionary to hold the results
            DEGsByGoTerm = {}
            
            for term in ListofGOterms:
                goName = graph.nodes[term]['name']
                # Make a one-row dataframe with GO ID and GO Name
                goNameAndTerm = pd.DataFrame({'GO Term': [term], 'Name': [goName], 'PageRank': [pagerank_scores[term]]})
                # Get the list of DEGs associated with this GO term
                goDEGs = HGNCsymbols_for_go_id(term)
                # Find DEGs that are associated with this GO term and also appear in the merged dataframe
                symbol_col = [col for col in merged_df.columns if 'Symbol' in col][0]
                matched_genes = set(goDEGs).intersection(set(merged_df[symbol_col]))
                
                # If any matching genes were found, save the corresponding rows from the merged dataframe
                if matched_genes:
                    DEGsByGoTerm[term] = merged_df[merged_df[symbol_col].isin(matched_genes)]
                
                # Check if the key exists in the dictionary before concatenating
                if term in DEGsByGoTerm:
                    DEGsByGoTerm[term] = pd.concat([goNameAndTerm, DEGsByGoTerm[term]], axis=1)
                
            return DEGsByGoTerm
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
        
        def output_to_csv(result: dict, file_prefix: str):
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
                response["Content-Disposition"] = f"attachment; filename=\"{output_file}\""
                return response
            
        nodes = ConcatKeywordHits(keyword)
        if len(nodes) == 0:
            return JsonResponse({})
        endConcat = time.perf_counter()
        print(f"Concat took {endConcat-start:0.4f}")

        files = request.FILES.getlist('de-csv-files')

        mapped_degs = DEGassociator(files, nodes)


        return output_to_csv(mapped_degs, keyword)
