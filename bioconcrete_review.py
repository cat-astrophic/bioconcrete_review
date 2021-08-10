# This script performs a network + biliometric analysis of the SCOPUS bioconcrete data

# Importing required modules

import pandas as pd
import numpy as np
import networkx as nx
from matplotlib import pyplot as plt
from cpip import cpip

# Declaring the username + filepath

username = 'macary'
filepath = 'C:/Users/' + username + '/Documents/Data/bioconcrete/'

# Reading in the data

bio = pd.read_csv(filepath + 'scopus.csv')

# Creating the adjacency of the citation network :: ROW CITES COLUMN

A = np.zeros((len(bio), len(bio)))

for r in range(len(bio)): # References list
    
    for t in range(len(bio)): # Title to check for in References
        
        if str(bio.Title[t]) in str(bio.References[r]):
            
            A[r][t] = 1
            
# Creating an arc list from A

arcs = [(i+1,j+1) for i in range(len(A)) for j in range(len(A)) if A[i][j] == 1]

# Visualizing the network with netowrkx

G = nx.DiGraph() # Initialize the digraph
G.add_edges_from(arcs) # Add edges
color_map = ['red' if bio['Source title'][g] == 'Journal of Cleaner Production' else 'blue' for g in list(G.nodes)] # Create a color map for the vertices
plt.figure(figsize = (9,9)) # Create a figure with matplotlib.pyplot
nx.draw_networkx(G, with_labels = False, node_color = color_map, pos = nx.fruchterman_reingold_layout(G)) # Draw the network
plt.savefig(filepath + 'citenet.eps')

# Create overall network statistics

# Density

density = nx.density(G)

# Plot citations by year

years = bio.Year.unique()
pubs = []

for y in range(min(years),max(years)+1):
    
    bi = bio[bio.Year == y]
    pubs.append(len(bi))

pubs[-1] = np.floor(pubs[-1]*12/7) # Interpolate 2021 data

pubs2 = [sum(pubs[0:i]) for i in range(len(pubs))]
basis = [i for i in range(min(years), max(years)+1)]
labs = [i for i in range(min(years), max(years)+1) if i%2 == 1]

plt.figure()
ax = plt.axes()
plt.plot(basis, pubs, color = 'purple')
plt.fill_between(basis, pubs, color = '#33C9FF')
plt.xlabel('Year')
plt.ylabel('Publications')
plt.xticks(ticks = np.arange(min(years), max(years)+2, 2), labels = labs, rotation = 45)
plt.title('Publications on biocement/bioconcrete by year')
plt.savefig(filepath + 'pubs.eps')

plt.figure()
ax = plt.axes()
plt.plot(basis, pubs2, color = 'black')
plt.fill_between(basis, pubs2, color = '#FF7733')
plt.xlabel('Year')
plt.ylabel('Publications')
plt.xticks(ticks = np.arange(min(years), max(years)+2, 2), labels = labs, rotation = 45)
plt.title('Total publications on biocement/bioconcrete')
plt.savefig(filepath + 'pubs2.eps')

# Normalized betweenness centralities of the vertices

centralities = list(nx.betweenness_centrality(G).values()) # Betweenness centralities
centralities =  [c / max(list(nx.betweenness_centrality(G).values())) for c in centralities] # Normalized
jcp = [i for i in range(len(color_map)) if color_map[i] == 'red'] # JCP paper indicies
jcp_cty = [centralities[j] for j in jcp] # JCP centralities
c2 = [c for c in centralities if c > 0] # Nonzero centralities

# So one JCP paper is the fifth most important paper per betweenness centrality

# Get lists of countries on each paper

nats = []

for b in range(len(bio)):
    
    print(str(b+1) + ' of ' + str(len(bio)) + '.......') # Visualizing progress
    
    affs = str(bio.Affiliations[b])
    affs = affs.replace('  ', ' ')
    affs = affs.replace('United Kingdom', 'UK')
    affs = affs.replace('United States', 'US')
    affs = affs.replace('Russian Federation', 'Russia')
    flag = True
    subnats = []
    
    while flag == True:
        
        idx = affs.find(';')
        
        if idx > 0:
            
            tmp = affs[:idx]
            
        else:
            
            tmp = affs
            flag = False
            
        tmp2 = tmp
        flag2 = True
        
        while flag2 == True:
            
            idx2 = tmp2.find(' ')
            
            if idx2 > 0:
                
                tmp2 = tmp2[idx2+1:]
                
            else:
                
                flag2 = False
                
        subnats.append(tmp2)
        affs = affs[idx+2:]
        
    nats.append(subnats)

nats[205] = ['India'] # Fixing a SCOPUS affiliation misprint

n_auths = [x.count(',') + 1 for x in bio.Authors] # Number of authors on each paper
n_nats = [len(pd.Series(x).unique()) for x in nats] # Number of nations on each paper
all_nats = [x for nat in nats for x in nat] # List of all distinct nations in the data
all_nats = list(sorted(list(pd.Series(all_nats).unique())))
all_nats[0] = 'South Africa'
all_nats[22] = 'Hong Kong'
all_nats[46] = 'New Zealand'
all_nats = all_nats[:47] # drop nan
all_nats = list(sorted(all_nats))

nat_counts = [] # Papers published by nation

for n in all_nats:
    
    count = 0
    
    for b in range(len(bio)):
        
        if n in str(bio.Affiliations[b]):
            
            count += 1
            
    nat_counts.append(count)

# Updating nats to account for the above changes

for n in nats:
    
    for x in n:
        
        if x == 'Africa':
            
            n[n.index(x)] = 'South Africa'
            
        elif x == 'Kong':
            
            n[n.index(x)] = 'Hong Kong'
            
        elif x == 'Zealand':
            
            n[n.index(x)] = 'New Zealand'

M = np.zeros((len(all_nats), len(all_nats))) # Initialize the IC network

for r in all_nats:
    
    for c in all_nats:
        
        count = 0
        
        for n in nats:
            
            if r != c: # Zero diagonals
                
                if (r in n) and (c in n):
                    
                    count += 1
                    
        M[all_nats.index(r)][all_nats.index(c)] = count

# Creating an edge list from M

edges = [(i+1,j+1) for i in range(len(M)) for j in range(len(M)) if M[i][j] == 1]

# Core-periphery analysis on the IC Network

M = pd.DataFrame(M)
M.to_csv(filepath + 'M.csv', index = False)
icpath = filepath + 'M.csv'
iccore = cpip(icpath, 2.5, 4)
core_nats = [all_nats[int(i)] for i in iccore]
net_nats = [n for n in all_nats if sum(M[all_nats.index(n)]) > 0]
net_nats_ids = [all_nats.index(n) for n in net_nats]
net_nats_ids2 = [all_nats.index(n)+1 for n in net_nats]

# Visualizing the network with networkx

G2 = nx.Graph() # Initialize the digraph
G2.add_edges_from(edges) # Add edges
color_map = ['red' if str(i-1) in iccore else 'green' for i in list(G2.nodes)] # Create a color map for the vertices
shells = [[i for i in list(G2.nodes) if str(i-1) in iccore], [i for i in list(G2.nodes) if str(i-1) not in iccore]] # Define the layout options for shell_layout
plt.figure(figsize = (8,6)) # Create a figure with matplotlib.pyplot
nx.draw_networkx(G2, with_labels = True, labels = dict(zip(net_nats_ids2, net_nats)), node_color = color_map, pos = nx.shell_layout(G2, shells), font_size = 10) # Draw the network
plt.savefig(filepath + 'icnet.eps')

# Creating a bar chart for total publications by nation with bar labels

nat_pubs = []

for nat in all_nats:
    
    count = 0
    
    for b in bio.Affiliations:
        
        if nat in str(b):
            
            count += 1
            
    nat_pubs.append(count)

order = np.argsort(nat_pubs)
bars = [nat_pubs[o] for o in order]
barlabs = [all_nats[o] for o in order]

fig = plt.figure(figsize = (18,10))
ax = fig.add_axes([0,0,1,1])
ax.bar(np.arange(len(bars)), bars, color = 'b', width = 0.5)
plt.xticks(np.arange(len(bars)), barlabs, color = 'black', rotation = 45)
plt.title('Total publications including at least one author from a nation')
plt.ylabel('Publications')
plt.savefig(filepath + 'bar.eps')

# Creating a table of the most frequent keywords and their frequencies

all_keys = []

for b in bio['Author Keywords']:
    
    bb = str(b).split(';')
    bb = [x[1:] if bb.index(x) > 0 else x for x in bb]
    
    for x in bb:
        
        all_keys.append(x.lower())

all_keys = list(sorted(list(pd.Series(all_keys).unique()))) # Final list of keys

key_counts = [] # Counts for each key

for k in all_keys:
    
    count = 0
    
    for b in bio['Author Keywords']:
        
        if k in str(b).lower():
            
            count += 1
            
    key_counts.append(count)

key_order = list(np.argsort(key_counts)) # Get the order
top_key_vals = [key_counts[o] for o in key_order] # Get the values
top_keys = [all_keys[o] for o in key_order] # Get the keys
key_tab = dict(zip(top_keys, top_key_vals)) # Create a dictionary of the data to put into a dataframe
key_tab = pd.Series(key_tab, name = 'Frequency') # Create a series from the dictionary
key_tab.index.name = 'Keyword' # Rename the index so it can become a properly named column
key_tab = key_tab.reset_index() # Make the series into a dataframe
key_tab = key_tab.iloc[::-1] # Reverse the order of the dataframe so the top keys are on top
key_tab = key_tab.iloc[[x for x in range(len(key_tab)) if x != 2]] # Drop nan from key_tab (papers with no keywords)
key_tab.to_csv(filepath + 'keys.txt', index = False) # Write the dataframe to a .txt file

