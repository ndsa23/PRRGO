import sys
import obonet
import networkx as nx
import pandas as pd
import os

# Load the GO ontology graph
url = 'http://purl.obolibrary.org/obo/go.obo'
graph = obonet.read_obo(url)

# Convert the graph to a networkx graph
nx_graph = nx.DiGraph(graph)

# Run pagerank on the networkx graph
pagerank_scores = nx.pagerank(nx_graph)

def KWNameQuery(userinput):
    name_nodes = []
    for node in graph.nodes:
        inputWithoutDashes = userinput.replace('-', ' ')
        if (userinput.lower() in graph.nodes[node]['name'].lower() or inputWithoutDashes.lower() in graph.nodes[node]['name'].lower()):
            name_nodes.append(node)
    return (name_nodes)

# KWDefQuery takes a user input keyword and searches all GO Terms definitions for that keyword.
# Hits are stored in a list called def_nodes

def KWDefQuery(userinput):
    def_nodes = []
    inputWithoutDashes = userinput.replace('-', ' ')
    for node in graph.nodes:
        if (userinput.lower() in graph.nodes[node]['def'].lower() or inputWithoutDashes.lower() in graph.nodes[node]['def'].lower()):
            def_nodes.append(node)
    return (def_nodes)

# ConcatKeywordHits takes in the userinput, calls KWNameQuery and KWDefQuery seperately with the same input
# then uses set to return only unique hits


def ConcatKeywordHits(userinput):
    return (set(KWNameQuery(userinput) + KWDefQuery(userinput)))


# ID to Name and Name to ID fetch either ID or Name, respectively from the other. Exact matches only
id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data=True)}
name_to_id = {data['name']: id_ for id_,
              data in graph.nodes(data=True) if 'name' in data}

def Parentfinder(term):
    node = term
    for child, parent, key in graph.out_edges(node, keys=True):
        return (f'{id_to_name[child]} -> {key} -> {id_to_name[parent]}')

# Function outputs parent terms to an input term
def Childfinder(term):
    node = term
    for parent, child, key in graph.in_edges(node, keys=True):
        return (f'{id_to_name[child]} <- {key} <- {id_to_name[parent]}')

# Identifies super/subterm relationships to the input GOTERM
def SupertermIdentifier(GOTERM):
    return (sorted(id_to_name[superterm] for superterm in networkx.descendants(graph, str(GOTERM))))


def SubtermIdentifier(GOTERM):
    return (sorted(id_to_name[subterm] for subterm in networkx.ancestors(graph, str(GOTERM))))

# Identifies path to root
def AllPathsToRoot(GOTERM, ROOT_GOTERM):
    paths = networkx.all_simple_paths(
        graph,
        source=GOTERM,
        target=ROOT_GOTERM
    )
    for path in paths:
        return ' -> '.join(id_to_name[node] for node in path)

# One function to house all the relationships per GOTERM

def NetworkMapper(GOTERM):
    relationships = {
        'GOTERM': GOTERM,
        'GONAME': id_to_name[GOTERM],
        'Parents': [],
        'Children': [],
        'Superterms': [],
        'Subterms': [],
        'PathtoMF': [],
        'PathtoBP': [],
        'PathtoCC': []
    }
    relationships['Parents'] = Parentfinder(GOTERM)
    relationships['Children'] = Childfinder(GOTERM)
    relationships['Superterms'] = SupertermIdentifier(GOTERM)
    relationships['Subterms'] = SubtermIdentifier(GOTERM)
    relationships['PathtoMF'] = AllPathsToRoot(
        GOTERM, name_to_id['molecular_function'])
    relationships['PathtoBP'] = AllPathsToRoot(
        GOTERM, name_to_id['biological_process'])
    relationships['PathtoCC'] = AllPathsToRoot(
        GOTERM, name_to_id['cellular_component'])
    return (relationships)

    # Read the GPAD file into a Pandas dataframe
GPIdf = pd.read_csv('goa_human.gpi.gz', compression="gzip", sep='\t', comment='!', header=None,
                 names=['DB', 'DB_Object_ID', 'DB_Object_Symbol', 'DB_Object_Name', 'DB_Object_Synonym',
                        'DB_Object_Type', 'Taxon', 'Parent_Object_ID','DB_Xrefs', 'Properties'])

# Read the GPI file into a Pandas dataframe
GPADdf = pd.read_csv('goa_human.gpad.gz', compression="gzip", sep='\t', comment='!', header=None,
                     names=['DB', 'DB_Object_ID', 'Qualifier', 'GO_ID', 'DB:Reference', 'Evidence_Code',
                            'WithOrFrom', 'Taxon_Intx', 'Date', 'Assigned_by', 'Annotation_Ext', 'Annotation_Prop'],
                     low_memory=False)

#Combine annotation file information into a df merged on Uniprot IDs
UPmerged_df = pd.merge(GPADdf, GPIdf, on='DB_Object_ID', how='inner')

#Make a dictionary with the folowing structure k:v= GO term:[UP_ID1, UP_ID2...UP_IDN]
GO_IDtoUniprot_ID = UPmerged_df.groupby('GO_ID')['DB_Object_ID'].apply(lambda x: list(set(x))).to_dict()

#Make a dictionary where these UP IDs are converted to their gene names
GO_IDtoHGNC = UPmerged_df.groupby('GO_ID')['DB_Object_Symbol'].apply(lambda x: list(set(x))).to_dict()

def HGNCsymbols_for_go_id(GO_ID):
    HGNCsymbols = UPmerged_df.loc[UPmerged_df['GO_ID'] == GO_ID, 'DB_Object_Symbol'].unique().tolist()
    return HGNCsymbols

def DEGmerger(file_paths):
    dfs = []

    for file_path in file_paths:
        # Extract the filename from the file path
        filename = os.path.basename(file_path)

        # Check if the file is a .tsv or .csv file
        if file_path.endswith('.tsv'):
            delimiter = '\t'
        elif file_path.endswith('.csv'):
            delimiter = ','
        else:
            raise ValueError(f"Unknown file format: {file_path}")

        # Read the DEG file into a DataFrame
        df = pd.read_csv(file_path, sep=delimiter, index_col=0)

        # Add a prefix to the column names with the filename
        df = df.add_prefix(f"{filename}_")

        dfs.append(df)

    # Merge all the DataFrames together into a single DataFrame
    merged_df = pd.concat(dfs, axis=1)

    return merged_df

# Function to associate DEGs with GO terms
def DEGassociator(file_paths, ListofGOterms):
    # Merge all .tsv files in the specified directory
    merged_df = DEGmerger(file_paths)
    
    # Create a dictionary to hold the results
    DEGsByGoTerm = {}
    
    for term in ListofGOterms:
        goName = graph.nodes[term]['name']
        # Make a one-row dataframe with GO ID and GO Name
        goNameAndTerm = pd.DataFrame({'GO Term': [term], 'Name': [goName]})
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

def output_to_csv(result, output_prefix):
    # Create a Pandas DataFrame from the result dictionary
    df = pd.concat(result.values())

    # Generate the output file name with the prefix
    output_file = f"{output_prefix}_output.csv"

    # Write the DataFrame to a CSV file
    df.to_csv(output_file, index=False)

    # Print the output file name
    print(f"Results saved to {output_file}")

    # Main function
def main():
    # Read user inputs from the command line
    file_paths = sys.argv[1].split(',')
    keyword = sys.argv[2]

    # Call ConcatKeywordHits to get the list of GO terms
    go_terms = ConcatKeywordHits(keyword)

    # Call DEGassociator to get the results
    result = DEGassociator(file_paths, go_terms)

    # Output the results to a CSV file
    output_file = 'output.csv'
    output_to_csv(result, keyword)
    print(f"Results saved to {output_file}")



if __name__ == "__main__":
    main()