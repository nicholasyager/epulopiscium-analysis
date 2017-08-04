""""
density_analysis.py

Usage:
    density_analysis.py <filename> <output_file>
"""

import logging

from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from docopt import docopt
import csv

from sklearn.cluster import DBSCAN                                              
from sklearn import metrics
from scipy.spatial import ConvexHull
from matplotlib.path import Path
from sklearn.decomposition import PCA

def load_file(path):
    chromosomes = []
    with open(path, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            chromosomes.append(row)
    return np.asarray(chromosomes).astype(np.float)

def main(args):

    # Step 1: Load a file of localized chromosomes.
    path = args['<filename>']
    filename = ".".join(path.split('/')[-1].split('.')[:-1])   
    
    logger = logging.getLogger(filename)                                        
    logger.setLevel(logging.DEBUG)                                              
    ch = logging.StreamHandler()                                                
    ch.setLevel(logging.DEBUG)                                                  
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)                                                  
    logger.addHandler(ch) 
   
    logger.info("Loading chromosomes.") 
    chromosomes = load_file(path)
    logger.info("Loaded {count} chromosomes.".format(count=chromosomes.shape[0])) 
    all_chromosomes = chromosomes.copy()

    # Step 2: Perform DBSCAN to remove noise and identify unique cells.

    # Step 2a: Perform an initial scan that removes the vast majority of noise
    # This step will have a strict epsilon value but a forgiving minimum
    # cluster size.
    epsilon = 1.75
    min_points = 2 
    
    logger.info("Performing noise-removal clustering with ϵ={epsilon} ρ={min_points}.".format(
            epsilon=epsilon, 
            min_points=min_points
    ))

    db = DBSCAN(eps=epsilon, min_samples=min_points).fit(chromosomes)

    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)                   
    core_samples_mask[db.core_sample_indices_] = True                           
    labels = db.labels_

    np.delete(chromosomes, np.where(labels == -1), axis=0)

    # Step 2b: Perform a scecond DBSCAN, this time to identify large clusters.
    # This pass will have a forgiving epsilon but a strict minimum cluster
    # size
    epsilon = 2.25 
    min_points = 5

    logger.info("Performing cell differentiation clustering with ϵ={epsilon} ρ={min_points}.".format(
            epsilon=epsilon, 
            min_points=min_points
    ))


    db = DBSCAN(eps=epsilon, min_samples=min_points).fit(chromosomes)

    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)                   
    core_samples_mask[db.core_sample_indices_] = True                           
    labels = db.labels_

    largest_labels = np.where(np.bincount(labels[np.where(labels > -1)]) > 100)[0]

    cells = []

    for label in largest_labels:
        cells.append(chromosomes[np.where(labels == label)])

    #print(largest_labels)

    unique_labels = set(largest_labels)
    colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    #for k, col in zip(unique_labels, colors):
    #    if k == -1:
    #         col = 'k'
    #         continue
    #    class_member_mask = (labels == k)
    #    xy = chromosomes[class_member_mask & core_samples_mask]
    #    ax.plot(xy[:, 0], xy[:, 1], xy[:,2], 'o', markerfacecolor=col,
    #        markeredgecolor='k', markersize=6)                                  
    #    xy = chromosomes[class_member_mask & ~core_samples_mask]                    
    #    ax.plot(xy[:, 0], xy[:, 1], xy[:,2], 'o', markerfacecolor=col,         
    #        markeredgecolor='k', markersize=6)            
    #
    #plt.show()
    
    # Step 3: For each cell, find the convex hull for the points and determine 
    # the cell volume.


    output_file = open(args["<output_file>"], 'a')

    for index, cell in enumerate(cells):
        logger.info("Performing convex hull calculations for cell {index}.".format(
            index=index 
        ))

        hull = ConvexHull(cell)
        chromosome_count = 0


        for chromosome in all_chromosomes:
            new_hull = ConvexHull(np.concatenate((cell, [chromosome])))
            if np.array_equal(hull.vertices, new_hull.vertices):
                chromosome_count += 1

        #for slice_index in np.arange(np.min(cell[:,2]), np.max(cell[:,2]), 0.5):
        #    cell_slice = cell[np.where(cell[:,2] == slice_index), 0:2][0]
        #    if cell_slice.shape[0] > 2:
        #        slice_hull = ConvexHull(cell_slice)
        #        slice_hull_path = Path(cell_slice[slice_hull.vertices] )
        #
        #        remaining_chromosomes = all_chromosomes[np.where(all_chromosomes[:,2] == slice_index),:2][0]
        #        chromosome_count += np.sum(slice_hull_path.contains_points(remaining_chromosomes))
        #    else:
        #        logger.warning("Not enough chromosomes are in slice {}.".format(slice_index))

        output_file.write('"{filename}",{cell},{count},{volume},{density}\n'.format(
            filename=filename,
            cell=index,
            count=chromosome_count,
            volume=hull.volume,
            density=chromosome_count/hull.volume
        ))
        logger.info(("An estimated {count} chromosomes were found in a volume "
                    "of {volume} μm³ for a density of {density} chromosomes/μm³").format(
                        count=chromosome_count,
                        volume=hull.volume,
                        density=chromosome_count/hull.volume
                    ))

if __name__ == "__main__":
    main(docopt(__doc__))
