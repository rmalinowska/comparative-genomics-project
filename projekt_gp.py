#! /usr/bin/python

import requests
import sys, re, json, os
import copy
from Bio import SeqIO, AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo import Consensus, NewickIO
# Specify the working directory
working_dir = "C:\\Users\\roksa\\Desktop\\sem7\\GP\\projekt\\"

# You can provide either file with species names or InterPro database proteome IDs

# If you want to provide file with species names:
# Specify the name for a .txt file with species names separated by new line and set SPECIES = True
SPECIES = False
if SPECIES:
    species_filename = "fungi_v2.txt"

# If you want to provide file with proteome IDs:
# Specify the name for a .txt file with species names and proteome IDs separated by new line(example: Homo sapiens,UP000005640) and set IDS = True
IDS = True
if IDS:
    ids_file = "fungi_v2_id.txt"

# Specify the name for a fasta output file which will contain proteome sequences from all species
merged_fasta_output_filename = "all_fungi_v2.fasta"
# Specify path to mmseqs software
path_to_mmseq = "C:\\Users\\roksa\\Desktop\\sem7\\GP\\projekt\\mmseqs\\bin\\mmseqs"
# Specify minimal sequence identity for sequence clustering
IDENT = 0.4
# Specify minimal sequence coverage for sequence clustering
COVER = 0.7
# Specify coverage mode for sequence clustering (read more in mmseqs manual)
COV_MODE = 1
# Specify output filename prefix for clustering
PREF = "fungi040_v2"
# Specify the path to muscle executable
path_to_muscle = "C:\\Users\\roksa\\Desktop\\sem7\\GP\\muscle5.1.win64.exe"
# Specify tree construction method. If not provided, neighbour-joining is applied
TREE_CONSTRUCTION_METHOD = 'nj'
# Specify method for calculating distance between sequences (['blosum45', 'blosum50', 'blosum62', 'blosum80',\
# 'blosum90', 'pam250', 'pam30', 'pam70'])
DISTANCE_METHOD = 'blosum62'

# Creating directories for results
PATH_PROTEOMES = working_dir+"proteomes\\"
if not os.path.exists(PATH_PROTEOMES):
  os.makedirs(PATH_PROTEOMES)
  print("Created directory for proteomes.")

PATH_CLUSTERS = working_dir+"clusters\\"
if not os.path.exists(PATH_CLUSTERS):
  os.makedirs(PATH_CLUSTERS) 
  print("Created directory for clusters.")

PATH_TREES = working_dir+"trees\\"
if not os.path.exists(PATH_TREES):
  os.makedirs(PATH_TREES) 
  print("Created directory for trees.")

PATH_ALIGNMENTS = working_dir+"alignments\\"
if not os.path.exists(PATH_ALIGNMENTS):
  os.makedirs(PATH_ALIGNMENTS) 
  print("Created directory for alignments.")

PATH_CLUST_PARA = working_dir+"clusters_paralogs\\"
if not os.path.exists(PATH_CLUST_PARA):
    os.makedirs(PATH_CLUST_PARA)
    print("Created directory for clusters with paralogs.") 

PATH_ALIGN_PARA = working_dir+"alignments_paralogs\\"
if not os.path.exists(PATH_ALIGN_PARA):
    os.makedirs(PATH_ALIGN_PARA)
    print("Created directory for alignments with paralogs.")


os.chdir(working_dir)

def get_species(input):
    """ Returns the list of species names obtained from a file, where each species is in an individual line."""
    with open(input, "r") as inp:
        inp = inp.read().split("\n")
        return inp

def get_species_and_ids(input):
    result = {}
    with open(input, "r") as inp:
        inp = inp.read().split("\n")
        for line in inp:
            species, ID = line.strip().split(",")[:2]
            result[species] = ID
    return result

def interpro_get(path):
    """ Function that returns response from InterPro database API given path. """
    base_url = "https://www.ebi.ac.uk/interpro/api/"
    url = base_url+path
    response = requests.get(url)
    if response.status_code != 200:
        print(path)
        raise ValueError(f"Request to {url} failed with status code {response.status_code}")
    return response

def interpro_proteome_ids(res):
    """ Function that given the response from InterPro API, returns obtained accession IDs. """
    data = json.loads(res.text)
    acc = []
    for item in data["results"]:
        acc.append(item["metadata"]["accession"])
    return acc

def get_all_proteome_ids(species):
    """ Given the species list, function sends call to InterPro API and returns dictionary of all found proteome IDs for all species."""
    result = {}
    for sp in species:
        name = "+".join([x.lower() for x in sp.split(" ")])
        url = "proteome/uniprot/?search="+name
        result[sp] = interpro_proteome_ids(interpro_get(url))
    return result

def ID_per_organism(input_dict):
    output = {}
    for species in input_dict.keys():
        if species == 'Saccharomyces cerevisiae':
            ID = 'UP000002311'
        else:
            ID = input_dict[species][0]
        output[species] = ID
    return output

def download_proteomes(interpro_ids):
    """ Function runs download_proteome.py script and downloads proteomes for all given species (assuming that the first ID for the species is most likely to be the correct one.) 
    Saved proteome has a file name like [proteome_ID].fasta"""
    for ID in interpro_ids.values():
        os.system("python "+working_dir+"download_proteome.py -id " + ID+" > " + PATH_PROTEOMES+ID+".fasta")

def add_organism_ids(interpro_ids):
    """ Function that updates all proteome fasta files so each protein has the name of the organism at the beggining.
    Output files are like [proteome_ID]_withname.fasta """
    for organism, input_ID in interpro_ids.items():
        with open(PATH_PROTEOMES+input_ID+".fasta", "r") as inp:
            with open(PATH_PROTEOMES+input_ID+"_withname.fasta", "w+") as out:
                inp = inp.read(). split(">")
                for protein in inp:
                    if protein != "":
                        protein = protein.split("\n")
                        header, sequence = protein[0], protein[1:]
                        header = "_".join(organism.split(" "))+"|"+header
                        out.write(">"+header+"\n")
                        for frag in sequence:
                            if frag != '':
                                out.write(frag+'\n')

def merge_fastas(IDs):
    """ Function merges all proteomes to one file named [filename].fasta """
    with open(working_dir+merged_fasta_output_filename, "w+") as out:
        for ID in IDs.values():
            single = open(PATH_PROTEOMES+ID+"_withname.fasta").read()
            out.write(single)

def run_mmseq(path, output_prefix, identity_threshold, coverage, cov_mode):
    os.system(path+" easy-cluster "+ working_dir+merged_fasta_output_filename + " " + output_prefix + " tmp " + "--min-seq-id " + \
        str(identity_threshold) + " -c " + str(coverage) +" --cov-mode " + str(cov_mode))

def parse_mmseq(tsv_file):
    """ Based on mmseq cluster file, creates dictionary of a form: cluster_ID : [list of sequence IDs from this cluster]. """
    result = {}
    with open(working_dir+tsv_file, "r") as data:
        data = data.readlines()
        for i in range(0, len(data)):
            cluster_ID, seq_ID, _ = re.split(r"\s+", data[i])
            result.setdefault(cluster_ID, []).append(seq_ID)
    keys = copy.deepcopy(list(result.keys()))
    for k in keys:
        if len(result[k]) <= 2:
            result.pop(k)
    return result

def get_clusters_1to1(species_list, clusters):
    final_clusters = {}
    all_organisms = set(species_list)
    for k, v in clusters.items():
        organisms = set([ID.split("|")[0] for ID in v])
        if len(organisms) >= 30:
            if len(organisms.intersection(all_organisms)) == len(all_organisms):
                new_organisms = set()
                new_values = []
                for ID in v:
                    org_name = ID.split('|')[0]
                    if org_name in all_organisms:
                        if org_name not in new_organisms:
                            new_organisms.add(org_name)
                            new_values.append(ID)
                final_clusters[k] = new_values
                
    return final_clusters


def get_sequences(fasta_file):
    result = {}
    data = SeqIO.parse(fasta_file, "fasta")
    for record in data:
        result[record.id] = record.seq
    return result

def save_clusters(seq_dict, clusters_dict):
    index = 0
    for v in clusters_dict.values():
        with open(PATH_CLUSTERS+str(index)+"cluster.fasta", "w+") as f:
            for ID in v:
                f.write(">"+ID+"\n"+str(seq_dict[ID])+"\n")
            index += 1

def save_clusters_paralogs(seq_dict, clusters_dict):
    index = 0
    for v in clusters_dict.values():
        with open(PATH_CLUST_PARA+str(index)+"cluster_paralogs.fasta", "w+") as f:
            for ID in v:
                f.write(">"+ID+"\n"+str(seq_dict[ID])+"\n")
            index += 1

def align_clusters(path_muscle):
    
    index = 0
    while True:
        try:
            os.system(path_muscle+" -align "+PATH_CLUSTERS+str(index)+"cluster.fasta"+" -output "+\
                PATH_ALIGNMENTS+str(index)+"cluster_alignment.fasta")
            index += 1
        except:
            break

def align_clusters_paralogs(path_muscle):
    index = 0
    while True:
        try:
            os.system(path_muscle+" -align "+PATH_CLUST_PARA+str(index)+"cluster_paralogs.fasta"+" -output "+\
                PATH_ALIGN_PARA+str(index)+"cluster_paralogs_alignment.fasta")
            index += 1
        except:
            break

def change_ids(path):
    file_names = os.listdir(path)
    for file in file_names:
        if file[-15:] == "alignment.fasta":
            with open(file, "r") as handle:
                with open(PATH_ALIGNMENTS+"CH_"+file, "w+") as out:
                    for record in SeqIO.parse(handle, "fasta"):
                        out.write(">"+str(record.id).split("|")[0]+"\n")
                        out.write(str(record.seq)+"\n")

def construct_tree(aln_file, calc_method = 'blosum62', tree_const_method = 'nj'):
    # function that constructs single tree from alignment file given method for calculating distance 
    # and method for tree construction
    aln = AlignIO.read(PATH_ALIGNMENTS + aln_file, "fasta")
    calculator = DistanceCalculator(calc_method)
    constructor = DistanceTreeConstructor(calculator, tree_const_method)
    tree = constructor.build_tree(aln)
    return tree

def build_trees(calc_method, tree_const_method):
    trees = []
    file_names = os.listdir(PATH_ALIGNMENTS)
    for file in file_names:
        if file[-15:] == "alignment.fasta" and file[0:2] == "CH":
            f = file.open()
            tree = construct_tree(f, calc_method, tree_const_method)
            trees.append(tree)
    return trees

def write_newick(trees, output):
    writer = Phylo.NewickIO.Writer(trees)
    with open(PATH_TREES+output, "w+") as out:
        writer.write(out)

def mark_unrooted(filename):
    with open(PATH_TREES+filename, "r") as inp:
        inp = inp.read().split(";")
        with open(PATH_TREES+filename[:-4]+"_unrooted.nwk", "w+") as out:
            for tree in inp:
                tree = tree.strip()
                if len(tree) != 0:
                    out.write("[&U]"+tree+";\n")


if SPECIES:
    species = get_species(species_filename)
    proteome_ids = get_all_proteome_ids(species)
    interpro_ids = ID_per_organism(proteome_ids)
elif IDS:
    interpro_ids = get_species_and_ids(ids_file)
else:
    print("You haven't properly provided file with species names or proteom IDs. Check if you've set one of variables to True.")
    sys.exit()

#download_proteomes(interpro_ids)
#add_organism_ids(interpro_ids)
#merge_fastas(interpro_ids)
#run_mmseq(path_to_mmseq, PREF, IDENT, COVER, COV_MODE)
#clusters = parse_mmseq(PREF+"_cluster.tsv")
#sizes_of_clusters = [len(val) for val in clusters.values()]
#final_clusters = get_clusters_1to1(list(interpro_ids.keys()), clusters)
#sequences_dict = get_sequences(merged_fasta_output_filename)
#save_clusters(sequences_dict, final_clusters)
#align_clusters(path_to_muscle)
#change_ids(PATH_ALIGNMENTS)
#trees = build_trees(DISTANCE_METHOD, TREE_CONSTRUCTION_METHOD)
#write_newick(trees, "all_trees.nwk")
#mark_unrooted("all_trees.nwk")

#save_clusters_paralogs(sequences_dict, clusters)
#align_clusters_paralogs(path_to_muscle)
change_ids(PATH_ALIGN_PARA)