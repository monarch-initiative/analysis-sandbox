import numpy as np
import argparse
import json
from sklearn.cluster import DBSCAN
from pathlib import Path
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.manifold import MDS
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main():
    """
    Cluster similarity matrix with DBSCAN and visualize with MDS

    Takes the mean of distances for asymmetric matrices

    Optionally perform MDS and visualize with matplotlib as well
    outputting serialized (json, csv) versions of various analysis
    steps

    Note this algorithm produces the same output as hierarchical
    single linkage clustering when minSamples is set to 2
    """
    parser = argparse.ArgumentParser(description='description')
    parser.add_argument('--input', '-i', type=str, required=True,
                        help='Location of input file'
                             ' that contains the sim matrix as json')
    parser.add_argument('--label', '-l', type=str, required=False,
                        help='Location of id-label mapping file')
    parser.add_argument('--epsilon', '-eps', type=int, required=False,
                        default=15, help='min distance between points in cluster')
    parser.add_argument('--min_samples', '-min', type=int, required=False,
                        default=3, help='Minimum samples to create cluster')
    parser.add_argument('--components', '-c', type=int, required=False,
                        default=2, help='Components for MDS')
    parser.add_argument('--visualize', '-v', action='store_true', required=False,
                        default=False, help='Plot results')
    parser.add_argument('--out', '-o', required=False, help='output directory')
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
                   for k in range(len(distance_matrix[i]))] for i in range(len(distance_matrix))]

    highly_sim = 0
    for i in range(len(sym_matrix)):
        for k in range(len(sym_matrix[i])):
            if sym_matrix[i][k] <= 10 and i != k:
                highly_sim += 1

    logger.info("Number of highly similar entities (90+): {0}".format(highly_sim))

    matrix = np.array(sym_matrix)

    db = DBSCAN(eps=args.epsilon, min_samples=args.min_samples, metric="precomputed").fit(matrix)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_

    logger.info(db.labels_)

    # Number of clusters in labels, ignoring noise if present.
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    samples = sum(1 if i != -1 else 0 for i in labels)

    logger.info('Number of clusters: {0}'.format(n_clusters))
    logger.info('Number of samples clustered {0}'.format(samples))

    if args.visualize:
        coordinates = convert_matrix_to_coordinates(sym_matrix, args.components)
        figure = plot_dbscan_clusters(db, coordinates, args.components)

    # Write out analysis files
    if args.out:
        """
        Make the file structure
        args.out
        |  summary.txt
        |  distance_matrix.json
        |  cluster_indices.json
        |  core_indices.json
        |
        |--clusters
        |    cluster1_ids.tsv
        |
        |--matrices
        |    cluster1_matrix.tsv
        |    cluster1_matrix.json
        |
        |--MDS
        |    2d-coordinates.json
        |    3d-coordinates.json
        """
        dir_path = Path(args.out)
        dir_path.mkdir(parents=True, exist_ok=True)
        cluster_dir = dir_path / 'clusters'
        cluster_dir.mkdir(parents=True, exist_ok=True)
        matrix_dir = dir_path / 'matrices'
        matrix_dir.mkdir(parents=True, exist_ok=True)
        mds_dir = dir_path / 'MDS'
        mds_dir.mkdir(parents=True, exist_ok=True)

        map_file = open(args.label, 'r')
        id_label_list = [line.rstrip('\n') for line in open(args.label, 'r')]
        sample_list = np.array(id_label_list)

        # Write distance matrix
        file_distance_matrix = dir_path / "distance_matrix.json"
        file_distance_fh = file_distance_matrix.open("w")
        file_distance_fh.write(json.dumps(matrix.tolist()))

        # Write core_indices
        file_core_indices = dir_path / "core_indices.json"
        file_core_fh = file_core_indices.open("w")
        file_core_fh.write(json.dumps(db.core_sample_indices_.tolist()))

        # Write cluster_indices
        file_clusters = dir_path / "clusters.json"
        file_clusters_fh = file_clusters.open("w")
        file_clusters_fh.write(json.dumps(db.labels_.tolist()))

        # Write out summary file
        file_summary = dir_path / "summary.txt"
        file_summary_fh = file_summary.open("w")
        file_summary_fh.write('Algorithm: DBSCAN\n')
        file_summary_fh.write('Epsilon: {0}\n'.format(args.epsilon))
        file_summary_fh.write('Min Samples: {0}\n'.format(args.min_samples))
        file_summary_fh.write('Number of clusters: {0}\n'.format(n_clusters))
        file_summary_fh.write('Number of samples clustered {0}\n'.format(samples))
        file_summary_fh.write('Number of input samples {0}\n'.format(len(db.labels_)))
        file_summary_fh.write("Number of highly similar entities (90+): {0}\n".format(highly_sim))
        file_summary_fh.write("Clusters: \n")

        for cluster_label in range(n_clusters):

            # For each cluster, create a tsv file with the columns
            # id, label, core
            class_member_mask = (labels == cluster_label)

            file_label = "cluster{0}_ids.tsv".format((cluster_label+1))
            cluster_label_file = cluster_dir / file_label
            cluster_label_fh = cluster_label_file.open('w')
            sample_ids = []
            sample_labels = []
            for line in sample_list[class_member_mask & core_samples_mask]:
                cluster_label_fh.write("{0}\tcore\n".format(line))
                sample_ids.append((line.split('\t')[0]))
                sample_labels.append((line.split('\t')[1]))
            for line in sample_list[class_member_mask & ~core_samples_mask]:
                cluster_label_fh.write("{0}\t\n".format(line))
                sample_ids.append((line.split('\t')[0]))
                sample_labels.append((line.split('\t')[1]))

            file_summary_fh.write("{0}\t{1}\n".format("|".join(sample_ids), "|".join(sample_labels)))

            # For each cluster create a custom distance matrix with
            # core samples first (unsorted) as tsv and json
            core_indices = np.nonzero((class_member_mask & core_samples_mask))
            non_core_indices = np.nonzero((class_member_mask & ~core_samples_mask))
            sample_indices = np.concatenate([core_indices[0], non_core_indices[0]], axis=0)
            legend = sample_list[sample_indices]
            sliced_matrix = matrix[sample_indices, :][:, sample_indices]

            # TODO merge labels with matrix

            file_json_label = "cluster{0}_matrix.json".format((cluster_label+1))
            file_json = matrix_dir / file_json_label
            file_json_fh = file_json.open('w')
            file_json_fh.write(json.dumps(sliced_matrix.tolist()))

            file_tsv_label = "cluster{0}_matrix.tsv".format((cluster_label+1))
            file_tsv = matrix_dir / file_tsv_label
            file_tsv_fh = file_tsv.open('wb')
            np.savetxt(file_tsv_fh, sliced_matrix, delimiter='\t', fmt='%2.2f')

            file_legend_label = "cluster{0}_legend.txt".format((cluster_label+1))
            file_legend = matrix_dir / file_legend_label
            file_legend_fh = file_legend.open('w')

            for line in legend:
                file_legend_fh.write("{0}\n".format(line))

        # Output position from MDS
        if args.visualize:
            file_2d_json = mds_dir / "2d-coordinates.json"
            file_2d_fh = file_2d_json.open('w')

            plot_2d_png = mds_dir / "2d-coordinates.png"
            plot_2d_fh = plot_2d_png.open('wb')

            file_3d_json = mds_dir / "3d-coordinates.json"
            file_3d_fh = file_3d_json.open('w')

            plot_3d_png = mds_dir / "3d-coordinates.png"
            plot_3d_fh = plot_3d_png.open('wb')

            if args.components == 2:
                file_2d_fh.write(json.dumps(coordinates.tolist()))
                figure.savefig(plot_2d_fh)
                three_d_coordinates = convert_matrix_to_coordinates(sym_matrix, 3)
                fig_3d = plot_dbscan_clusters(db, three_d_coordinates, 3)
                fig_3d.savefig(plot_3d_fh)
                file_3d_fh.write(json.dumps(three_d_coordinates.tolist()))

            else:
                file_3d_fh.write(json.dumps(coordinates.tolist()))
                figure.savefig(plot_2d_fh)
                two_d_coordinates = convert_matrix_to_coordinates(sym_matrix, 2)
                fig_2d = plot_dbscan_clusters(db, two_d_coordinates, 2)
                fig_2d.savefig(plot_3d_fh)
                file_2d_fh.write(json.dumps(two_d_coordinates.tolist()))


def convert_matrix_to_coordinates(sym_matrix, components):
    """
    :param sym_matrix: array, [n_samples, n_samples]
    :param components: int: 2 or 3 for MDS
    :return: Output of MDS, xy or xyz coordinates as 2d numpy array
             with shape [n_samples, components]
    """
    # Create coordinates based on multi dimensional scaling
    mds = MDS(n_components=components, dissimilarity="precomputed", random_state=1)
    coordinates = mds.fit_transform(sym_matrix)
    return coordinates


def plot_dbscan_clusters(db, coordinates, components):
    """
    :param db: dbscan object
    :param coordinates: 2d array of shape [n_samples, components]
    :param components: int: 2 or 3 for MDS
    :return: Output of MDS, xy or xyz coordinates as 2d numpy array
    """

    # Create array of Falses with the same shape as db.labels_ (list)
    # For more info on masked arrays see
    # https://docs.scipy.org/doc/numpy/reference/maskedarray.generic.html
    # although we're not using the ma module
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    # Make the core samples True
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_

    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

    # Black removed and is used for noise instead.
    unique_labels = set(labels)
    colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))
    fig = plt.figure()

    if components == 3:
        ax = fig.add_subplot(111, projection='3d')

        # Zip color array and cluster labels
        # ie [(0, array([ 0.61960787,  0.00392157,  0.25882354,  1. ])),
        #  (1, array([ 0.66189929,  0.05082661,  0.2688197 ,  1. ])), ...
        for cluster_label, color in zip(unique_labels, colors):
            if cluster_label == -1:
                # Black used for noise.
                color = 'k'
                # Remove the continue to show non clustered points
                continue

            # create a list of booleans using boolean indexing
            class_member_mask = (labels == cluster_label)

            # Plot core samples separately in case we want to
            # make them different sizes

            # Select indices that are in the cluster and also core samples
            # using boolean indexing
            xyz = coordinates[class_member_mask & core_samples_mask]

            # xyz[:, 0] = array at zero index of 2d array
            ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], c=color)

            # Slice out indices that are in the cluster and not core samples
            xyz = coordinates[class_member_mask & ~core_samples_mask]
            ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], c=color)

    if components == 2:
        for cluster_label, color in zip(unique_labels, colors):
            if cluster_label == -1:
                # Black used for noise.
                color = 'k'
                continue

            class_member_mask = (labels == cluster_label)

            xy = coordinates[class_member_mask & core_samples_mask]
            plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=color,
                     markeredgecolor='k', markersize=7)

            xy = coordinates[class_member_mask & ~core_samples_mask]
            plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=color,
                     markeredgecolor='k', markersize=7)

    plt.title('Estimated number of clusters: %d' % n_clusters_)
    plt.show()

    return fig

if __name__ == "__main__":
    main()