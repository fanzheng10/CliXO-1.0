# CliXO (Clique eXtracted Ontology) 1.0

Fan Zheng (fanzheng1101@gmail.com), Trey Ideker lab, UCSD  
Date: 08/27/2018

## About

This repository is an updated version of the CliXO (Clique eXtracted Ontology) algorithm, originally associated with the following publication:

Kramer M, Dutkowski J, Yu M, Bafna V, Ideker T. Inferring gene ontologies from pairwise similarity data. Bioinformatics, 30: i34-i42. 2014. doi: 10.1093/bioinformatics/btu282

A manuscript describing this new version is under preparation.

The input of the program is a weighted similarity network of objects (in our cases, genes), and the output is a directed acyclic graph (DAG), in which the leaf nodes are genes, and non-leaf nodes are gene sets called "terms". Term A is a descendent of another term B (i.e. in the DAG there is a path from A to B) if the gene set for A is a subset of the gene set for B. The whole DAG is called a "data-driven hierarchy" inferred from the input network.

The new version has a few key improvements over the original version, increasing the accuracy and robustness inferred from pairwise similarity (weighted network) data.

## Usage

`./clixo -i input_file -a alpha [-b beta] [-m Newman's modularity cutoff] [-z z-modularity cutoff] [-s stop_score] [-B]`

-i: input file (similarity scores); The format of this file should be three columns separated by tab. The first two columns should be two strings for the node names (using numbers may cause problem); and the third column should be a value for the edge weight.

-a: Alpha parameter; for the step size of hierarchy construction; usually, a smaller Alpha will create "deeper" hierarchies with more levels from leaves to 
 the root.
   
-b: Beta parameter; for merging overlapping terms. Two existing terms will be merged if their similarity is above the threshold defined by Beta. Usually a lower Beta will create a hierarchy with more larger terms, since small terms can be merged. However, at a lower Beta the terms look less like a clique, since the requirement for being a clique has been relaxed. Note the defintion of Beta is different from that in Kramer et al. See the manuscript for the new definition (NOT PUBLIC YET).
  
-m: Modularity parameter; calculate the contribution of each term to the Newman-Girvan's modularity in the network at the current score threshold; terms lower than this threshold will be removed from the output.
  
-z: Z-modularity parameter; another metric to remove some terms from the output. See "Miyauchi, A. & Kawase, Y. Z-Score-Based Modularity for Community Detection in Networks. PLoS One 11, e0147805 (2016)". Both -m and -z remove some terms and increase the clarity of the output hierarchy to human visualization. Note they remove different types of undesired terms; sometimes there are small but relatively isolated terms in big networks. If one believes those are important and should be kept, it is recommended to set -m to a low value or not using it at all, but use -z filter instead.
  
-s: A cutoff of similarity score, if set, the program will terminate when it reaches this point, and stop looking for more terms from scores lower than this threshold.
   
-B: if set, the program will interpret Beta (-b) as the old definition in Kramer et al.


## Differences to the old version

We observed several problems in CliXO 0.3. 

1. When we visualized the hierarchical model generated by CliXO 0.3, we always saw that biggest modules in the hierarchy forming a long chain, e.g. the top 5 biggest modules A1, A2, A3, A4, A5 form a structure in which A_{i+1} belongs to A_i. Sometimes the chain can be very long. A visual example can be found in (http://atgo.ucsd.edu/). It is suspicious that such a "Russian doll" structure is a good description of the hierarchical organization of cellular components. 
2. A robust ontology inference algorithm should be insensitive to stochastic noise in the input similarity network, i.e. slightly perturbed similarity scores should not significantly affect the resulting hierarchy. However, we found CliXO 0.3 is quite sensitive to small perturbation in input scores, sometimes the content of medium-to-big size terms are very different between two input scores with >0.99 pearson correlation.
3. Sometimes we found terms for which most of the edges weights are very low, which do not look like valid terms during visualization.

We found the major cause of this problem is due to a "missing edge inference" step. The idea of CliXO is based on finding maximal cliques in the network. However, real-world networks are noisy and there could a large number of maximal cliques that are similar to each other. In CliXO 0.3, parameter BETA was introduced to merge two maximal cliques that are highly overlapping by adding the missing edges. However, during this step, the missing edges were added to the network with edge weights equal to the edge weight threshold visited in the current iteration, which makes the future iterations see an altered network. This operation leads to the arfifacts mentioned above because:

1. The network edge density was increased every time a missing edge was inferred. Thus, a big and tight component could form within a few iterations. The iterations afterwards only incrementally added the genes that has not been covered by the big component. Thus, a "long-chain" structure was formed.
2. The missing edges that were inferred early have strong impact on the results in the future iterations. Slight perturbations of edge weights could make the order of processing edges different in early iterations, which was then augmented in future iterations. Thus the results were sensitive to slight perturbations.
3. Because missing edges became "real" edges after being inferred, it is possible that in future iterations, some missing edges were clustered together to form cliques/communities. But they are false positives since their actual edge weights are weak.

We fixed this issue and the "long-chain" structure never appeared again, and two similar inputs now generate similar output hierarchies.


## Auxillary functions (identical to CliXO_0.3): 

`extractOnt FILE THRESHOLD MIN_TERM_SIZE OUTFILE`

This script will allow the user to "peek" at the results of an in process CliXO run. It takes as input a CliXO output file (FILE), a similarity threshold above which terms will be saved and below which terms will be ignored (THRESHOLD), a minimum term size to keep (MIN_TERM_SIZE) and the name of an output file to create (OUTFILE).

`ontologyTermStats`

This utility will allow users to look at the ontology created by clixo.  There are several different types of stats that it can generate, but the most useful are the options "size" and "genes".  Option "size" will return a two column file where column one is the term identifier and column two is the number of leaf nodes / genes annotated to that term (all annotations are propagated upwards in the ontology).  Option "genes" adds a third column which is a comma separated list of all the genes annotated to the term.

# TODO: upload a new toy example
`grep -v "#" exampleOutput > exampleOutputOnt  
../ontologyTermStats exampleOutputOnt genes gene > exampleOutputOntStats  
../extractOnt exampleOutput 0.5 2 examplePeek`

## Limitations

Since maximal clique enumeration is expensive (theoretical upper bound is exponential to the number of nodes), applying CliXO on a large number of genes is computational infeasible. If the goal is to the finish the construction of the entire hierarchy,  we recommend not to exceed 1000 genes to expect the running time to be within a few hours. If either overlapping or multi-scale properties is not crucial to your research program, you should stop browsing and turn to other popular algorithms, such as Louvain clustering. 

Although there are many fast community detection algorithms, most of them partition graphs into non-overlapping parts, and in other overlapping-community compatible options only a small number of nodes on the boundaries of partitions can be assigned to multiple communities. Even fewer of them have the concept of multi-scale community detection. Thus, CliXO is unique in that it can possibly found both mutl-scale and highly-overlapping community structures. We do notice that a lot of cliques found by CliXO were unnecessary since they are later removed in post-processing steps. How to utilize such information to reduce the running time of clique finding will be in our future direction. 


## Acknowledgements

The improvements mentioned above benefited a lot from HiView (hiview.ucsd.edu), a web-based platform for visualizing hierarchical models mainly developed by Keiichiro Ono.

The author also thanks Michael Ku Yu and Anton Kratz for helpful discussion.