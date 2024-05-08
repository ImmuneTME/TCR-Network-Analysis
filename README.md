# TCR-Network-Analysis
TCR clonal sequence similarity-based graph building is a novel approach for biomarkers' identification, a priority in immunotherapy research. 

This code in R language will allow for the creation of graphs in which each node represents a clone and each edge connects clones that differ in 1 or 2 aminoacids in their CDR3 sequence. Thus, allowing for simmilarity analysis of TCR repertoires and potentially infering functional clusters. 

For each sample it requires 2 dataframes:
1.- Clone ID as rowname and CDR3 aminoacidic sequence as variable. 
2.- Clone ID as rowname and clone's relative abundance (frequency).

When taking calculating the graph's metrics we may consider the whole sample's repertoire or only what we call the "connected" graph, meaning eliminating those clones that differ in at least 3 amino acids with all others. 
Both posibilities are calculated through this code, the latter being regarded as "con" in the variables' names. 

In the plot, node size is proportional to the clone's frequency. In order to reduce computation times (and for the sake of clarity), only the connected graph is plotted.  
