import numpy as np
import argparse
import json
from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
"""
centers = [[1, 1], [-1, -1], [1, -1]]
X, labels_true = make_blobs(n_samples=750, centers=centers, cluster_std=0.4,
                            random_state=0)

X = StandardScaler().fit_transform(X)
"""
parser = argparse.ArgumentParser(description='description')
parser.add_argument('--input', '-i', type=str, required=True,
                    help='Location of input file')
parser.add_argument('--label', '-l', type=str, required=False,
                    help='Location of label mapping file')
args = parser.parse_args()

input_file = open(args.input, 'r')
distance_matrix = json.load(input_file)

print(distance_matrix[12][77])
print(distance_matrix[77][12])

for i in range(len(distance_matrix)):
    for k in range(len(distance_matrix[i])):
        if distance_matrix[i][k] == 0 and i != k:
            print("WARN: 0 in unexpected coordinate y:{0} x:{1}".format(i, k))


# Convert assymetric matrix to symmetric
sym_matrix = [[ ((distance_matrix[i][k] + distance_matrix[k][i]) / 2) for k in range(len(distance_matrix[i]))] for i in range(len(distance_matrix))]

highly_sim = []
for i in range(len(sym_matrix)):
    for k in range(len(sym_matrix[i])):
        if sym_matrix[i][k] <= 10 and i != k:
            highly_sim.append(sym_matrix[i][k])

print("Number of highly similar entities (90+): {0}".format(len(highly_sim)/2))

X = np.array(sym_matrix)

db = DBSCAN(eps=18, min_samples=5, metric="precomputed").fit(X)
core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
core_samples_mask[db.core_sample_indices_] = True
labels = db.labels_

print(db.labels_)

# Number of clusters in labels, ignoring noise if present.
n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
samples = sum(1 if i != -1 else 0 for i in labels)

foo = [ind for ind, val in enumerate(labels) if val == 2]
print(foo)


print('Estimated number of clusters: %d' % n_clusters_)
print(samples)
#print("Homogeneity: %0.3f" % metrics.homogeneity_score(labels_true, labels))
#print("Completeness: %0.3f" % metrics.completeness_score(labels_true, labels))
#print("V-measure: %0.3f" % metrics.v_measure_score(labels_true, labels))
#print("Adjusted Rand Index: %0.3f"
#      % metrics.adjusted_rand_score(labels_true, labels))
#print("Adjusted Mutual Information: %0.3f"
#      % metrics.adjusted_mutual_info_score(labels_true, labels))
#print("Silhouette Coefficient: %0.3f"
#      % metrics.silhouette_score(X, labels))


# Black removed and is used for noise instead.
unique_labels = set(labels)
colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))
for k, col in zip(unique_labels, colors):
    if k == -1:
        # Black used for noise.
        col = 'k'

    class_member_mask = (labels == k)

    xy = X[class_member_mask & core_samples_mask]
    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=14)

    xy = X[class_member_mask & ~core_samples_mask]
    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=6)

plt.title('Estimated number of clusters: %d' % n_clusters_)
plt.show()