"""
localization.py

Usage:
    localization.py <filename>
"""

from docopt import docopt
from PIL import Image
import numpy as np
import datetime
from scipy import ndimage
import cv2

from tqdm import tqdm
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import DBSCAN
from sklearn import metrics
import random
import math

import logging

class Stack:

    def __init__(self, path):
    
        self.path = path
        self.tiff = Image.open(path)

        self.images = []

        try:
            while True:
                image = np.asarray(self.tiff, dtype=np.uint8)
                self.images.append(image)
                self.tiff.seek(self.tiff.tell()+1)
        except EOFError:
            pass

        self.images.reverse()
        self.images = np.array(self.images)
        
        self.shape = self.images.shape

class Point:

    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

        self.visited = False
        self.cluster = None

    def __str__(self):
        return "{x}, {y}, {z}".format(x=self.x, y=self.y, z=self.z)

class Dataset:

    def __init__(self, points, stack):
        self.points = points
        self.stack = stack

    def region_query(self, point, epsilon):
        """
        Return all points within the given point's epsilon radius, including 
        the given point.
        """

        neighbors = []
        for delta_z in -1 : 1

        return self.points[np.where(self.adjacency[index] <= epsilon)]

def main(args):

    # This is the path to the TIF file to be processed
    path = args['<filename>']
    filename = ".".join(path.split('/')[-1].split('.')[:-1])


    logger = logging.getLogger(filename)

    logger.setLevel(logging.DEBUG)

    # create console handler and set level to debug
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)

    # create formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    # add formatter to ch
    ch.setFormatter(formatter)

    # add ch to logger
    logger.addHandler(ch)

    logger.info("Loading {path}.".format(path=path))

    stack = Stack(path)
    logger.info("Loaded {depth} {width}x{height} images.".format(depth=stack.images.shape[0],
                                                                 width=stack.images.shape[2],
                                                                 height=stack.images.shape[1]))
    
    min_image = np.amin(stack.images, axis=0)
    max_image = np.amax(stack.images, axis=0)

    new_max_image = max_image - min_image

    # Remove background
    stack.images -= min_image

    # Normalize the images with background removed
    stack.images = (stack.images.astype('float64') * ( 255.0 / stack.images.max())).astype('uint8')

    #for image in stack.images:
    #    cv2.imshow('frame', image)
    #    cv2.waitKey(0)

    # Perform a thresholding
    #plt.hist(stack.images[np.where(stack.images > 0)].ravel(), bins = np.arange(0,255, 5))
    #plt.show()

    percentile = 0.99

    threshold = np.percentile(stack.images[np.where(stack.images > 0)], percentile*100)
    logger.info("Using a threshold of {threshold}.".format(threshold=threshold))
    logger.info("Expect a maximum of {maximum} identified pixels".format(
        maximum=percentile * stack.images.shape[0] * stack.images.shape[1] *\
                stack.images.shape[2]))

    stack.images[np.where(stack.images < threshold)] = 0
    stack.images[np.where(stack.images > threshold)] = 255

    # Binary errosion of the stack to remove noise and to differentiate between
    # chromosomes.
    stack.images = ndimage.morphology.binary_erosion(stack.images, iterations=1).astype(stack.images.dtype)
   
    # DB-SCAN clustering of light points
    # Density-based spacital clustering of applications with noise allows for h
    # clustering og points in space, grouping points that are closely packed
    # together. In this case, it will group together the associated pixels of a
    # flourecent chromosome.

    epsilon = 1 
    min_points = 1 

    data_list = [[point[0], point[1], point[2]] for point in list(zip(*np.where(stack.images >= 1)))]
    dataset = np.array(data_list)

    db = DBSCAN(eps=epsilon, min_samples=min_points).fit(dataset)

    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_

    #for layer, image in enumerate(stack.images):
    #    color_image = cv2.cvtColor(image, cv2.COLOR_GRAY2BGR)
    #    for index, chromosome in enumerate(data_list):
    #        if chromosome[0] == layer and labels[index] != -1:
    #            color_image[chromosome[1], chromosome[2]] = (0,255,0)
    #        if chromosome[0] == layer and labels[index] == -1:
    #            color_image[chromosome[1], chromosome[2]] = (255,0,0)
    #
    #    cv2.imshow("slice", color_image)
    #    cv2.waitKey(0)


    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    n_identified = len([label for label in labels if label > -1])
    n_noise = len([label for label in labels if label == -1])

    logger.info('Identified {identified} of {total}, {ratio}%'.format(
        identified=n_identified,
        total=stack.images.shape[0]*stack.images.shape[1]*stack.images.shape[2],
        ratio=n_identified/(stack.images.shape[0]*stack.images.shape[1]*stack.images.shape[2])
    ))
    logger.info('Identified {count} possible clusters in the image stack.'.\
            format(count=n_clusters_))

    # Find the com for each cluster
    clusters = {}
    for index, label in enumerate(labels):
        if label == -1:
            continue
        if label not in clusters:
            clusters[label] = []
        clusters[label].append(dataset[index,:])

    chromosomes = []

    logger.info("Writing chromosome COMs to {dest}.".format(dest=filename+'.csv'))
    with open(filename+".csv", 'w') as output_file:
        for key, cluster in clusters.items():
            chromosome = np.mean(np.array(cluster), axis=0)
            chromosomes.append(chromosome)
            output_file.write('{x},{y},{z}\n'.format(z=chromosome[0]*0.5,
                                                   y=chromosome[1]*0.1,
                                                   x=chromosome[2]*0.1))
        
    logger.info("Localization complete.")

if __name__ == "__main__":
    main(docopt(__doc__))
