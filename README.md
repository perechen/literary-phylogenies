# literary-phylogenies

## Notes on the MPI cloud files

Edge lists are produced with the `network_prep.R` script.

They sit in `data/edges` folder (on the cloud). Currently, the variation comes from the percentile cutoff for determining 'meaningful' (strong, close) connections between books. 0.001 means that we take only 0.1% of the closest distances. 

- `edges_list_001.tsv` smaller network on a tight 0.001 cutoff
-  `edges_list_005.tsv` network with 0.005 of the distances considered
- `edges_001_authors.tsv` smaller network, but nodes are authors now, not individual books

The source data is Delta distances (Manhattann on z-scores) calculated for book represenations from an LDA model (200 topics). The full distance matrix is `data/distances_delta_lda.rds`
