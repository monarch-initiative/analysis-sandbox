## Scripts for generate and cluster similarity and distance matrices based on phenotypic similarity

### generate-sim-matrix.py
Generates similarity matrix based on OwlSim2 (PhenoDigm)
phenotypic similarity scores.

Input is a tab delimited file in the format:
foobar:123\tHP:123|HP:345|HP:567

Output is a json formatted asymmetric 2d array

Supports batch requests using grequests; however, when using
the monarch owlsim endpoint --chunk should not be set above 3

--temp /path/to/file generates a temporary file that can be used
as a cache for subsequent runs

### cluster-dbscan.py
Cluster similarity matrix with DBSCAN and visualize with MDS

Takes the mean of distances for asymmetric matrices

Optionally perform MDS and visualize with matplotlib as well
outputting serialized (json, csv) versions of various analysis steps

Note this algorithm produces the same output as hierarchical
single linkage clustering when minSamples is set to 2
