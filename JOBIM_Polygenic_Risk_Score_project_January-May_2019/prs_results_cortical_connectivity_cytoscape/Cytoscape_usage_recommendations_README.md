# Cytoscape usage recommendations

## Input file

Tab separated file:

Source_node | Target_Node | Edge_Attribute

e.g.: 
ctx_region1 | ctx_region2 | -log10(pval)

## Create network

File > Import > Network from File

Select in your table which column is source_node, target node, edge_attribute.

## Change Layout

Try different layout and find the one that suits you:

Layout > Grid layout
	   > Hierarchical layout
	   > Circular layout 
	   ...

## Compute Network Analysis

Compute networks metrics such as node degree, average_shortest_path_length, clustering_coefficient...

Tools > NetworkAnalyzer > Network Analysis > Analyse Network... 

## Change node size proportionally of node degree

Control Panel > Style > Node > Size > degree, continuous mapping

## Change edge width

Control Panel > Style > Edge > Width > -log10(pval), continuous mapping

## Color selectively certain edges 

Add another column to the inuput file, attribute ascending number to edge attribute, and same number for the edges you want to be colored with the same color, then :

Control Panel > Style > Edge > Color > color_column_edge_attribute > Discrete Mapping > Right click Mapping value generator

Note: same if you want to selectively color certain nodes: 

add column node_attribute to a new table file:
Node_name | Node_attribute

File > Import > Table from File

Control Panel > Style > Node > Color > color_column_node_attribute > Discrete Mapping > Right click Mapping value generator

## Cluster network into sub-networks

Cluster large network into sub-networks based on most connected parts.

App > clusterMaker > Community Cluster (Glay) > create new clustered network

