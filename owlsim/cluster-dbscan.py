import numpy as np
import argparse
import json
from sklearn.cluster import DBSCAN
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from sklearn.manifold import MDS
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main():
    """
    Cluster similarity matrix with DBSCAN and visualize with MDS
    """
    parser = argparse.ArgumentParser(description='description')
    parser.add_argument('--input', '-i', type=str, required=True,
                        help='Location of input file'
                             ' that contains the sim matrix as json')
    parser.add_argument('--label', '-l', type=str, required=False,
                        help='Location of label mapping file')
    parser.add_argument('--epsilon', '-eps', type=int, required=False,
                        default=15, help='Components for MDS')
    parser.add_argument('--min_samples', '-min', type=int, required=False,
                        default=3, help='Components for MDS')
    parser.add_argument('--components', '-c', type=int, required=False,
                        default=2, help='Components for MDS')
    parser.add_argument('--visualize', '-v', action='store_true', required=False,
                        default=False, help='Plot results')
    args = parser.parse_args()

    input_file = open(args.input, 'r')
    similarity_matrix = json.load(input_file)

    distance_matrix = [
        [
            100-value for index, value in enumerate(similarity_matrix[i])
        ]
        for i in range(len(similarity_matrix))
    ]

    for i in range(len(distance_matrix)):
        for k in range(len(distance_matrix[i])):
            if distance_matrix[i][k] == 0 and i != k:
                logger.warn("WARN: 0 in unexpected coordinate y:{0} x:{1}".format(i, k))

    # Convert asymmetric matrix to symmetric
    sym_matrix = [[((distance_matrix[i][k] + distance_matrix[k][i]) / 2)
                   for k in range(len(distance_matrix[i]))] for i in range(len(distance_matrix))
    ]

    highly_sim = []
    for i in range(len(sym_matrix)):
        for k in range(len(sym_matrix[i])):
            if sym_matrix[i][k] <= 10 and i != k:
                highly_sim.append(sym_matrix[i][k])

    logger.info("Number of highly similar entities (90+): {0}".format(len(highly_sim)/2))

    X = np.array(sym_matrix)

    db = DBSCAN(eps=args.epsilon, min_samples=args.min_samples, metric="precomputed").fit(X)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_

    logger.info(db.labels_)

    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    samples = sum(1 if i != -1 else 0 for i in labels)

    logger.info('Estimated number of clusters: {0}'.format(n_clusters_))
    logger.info('Number of samples clustered {0}'.format(samples))

    if args.visualize:
        plot_dbscan_clusters(db, sym_matrix, args.components)


def plot_dbscan_clusters(db, sym_matrix, components):

    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_

    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

    # Black removed and is used for noise instead.
    unique_labels = set(labels)
    colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))

    mds = MDS(n_components=components, dissimilarity="precomputed", random_state=1)

    pos = mds.fit_transform(sym_matrix)

    if components == 3:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        for k, col in zip(unique_labels, colors):
            if k == -1:
                # Black used for noise.
                col = 'k'
                continue

            class_member_mask = (labels == k)

            xyz = pos[class_member_mask & core_samples_mask]
            ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], c=col)

            xyz = pos[class_member_mask & ~core_samples_mask]
            ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], c=col)

    if components == 2:
        for k, col in zip(unique_labels, colors):
            if k == -1:
                # Black used for noise.
                col = 'k'
                continue

            class_member_mask = (labels == k)

            xy = pos[class_member_mask & core_samples_mask]
            plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
                     markeredgecolor='k', markersize=14)

            xy = pos[class_member_mask & ~core_samples_mask]
            plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
                     markeredgecolor='k', markersize=14)

    plt.title('Estimated number of clusters: %d' % n_clusters_)
    plt.show()

if __name__ == "__main__":
    main()