# Chain Cluster

Cluster the chains in skempi v1 with agglomerative hierarchical clustering and blastp
## Requirement
Sklearn > version 2.1

BioPython, Skbio

Blast installed in PATH

## Version

### 1.0

```diff
# Initial version of the clustering
# Distance function is  1 - (2*Identities) / (query_seq_len + subject_seq_len)
# Use min of the distance to solve blast result inconsistency
# Use complete linkage (the distance to a intermediate cluster is the max of the distance to each item in the cluster)
# distance threshold set to 0.75, so that identity > 25% are clustered
! Depend initial label, the cluster label might be different
```

## Results

| Item | Location |
| --- | --- |
| cluster label table | ./RES/cluster_label.csv |
| skempi dataset with label | ./RES/SKP1402m.ddg_class.txt |
| Visualization of Distance Matrix | dis_matrix.png |
| Dendrogram of all chain | dendrogram.png | 