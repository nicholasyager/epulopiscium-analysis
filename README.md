# Epulopiscium Analysis

**Nicholas A. Yager** and **Matthew Taylor**

This repository contains tools to analyze image stacks of flourescently-labeled
epuloposcoium chromosomes.

## Requirements:

The different tools require the following:
* Python 3
    * OpenCV 2.7
    * Numpy
    * Scipy
    * atexit
    * Scikit Learn
* R
    * rgl
    * cluster
    * plotrix
    * geometry
    * mclust
    * MASS
    * vioplot
    * ggplot2

## Python Scripts

### Localization

The localization script identifies likely flourescently-tagged chromosomes in
the image stack, clusters nearby voxels into single chromosomes, and saves the
center-of-mass for each chromosome to a CSV file.

```
localization.py

Usage:
    localization.py <filename>

Options:
    filename    The name of the TIFF image stack to analyze.
```

### Density Analysis

The density script processes a localiation file to identify real chromosomes
from spatial noise in the images, calculates the convex hull of the cell from
the chromosomes, and estimates the chromosome density of the cell.

```
density_analysis.py

Usage:
    density_analysis.py <filename> <output_file>

Options:
    filename    The name of the localized chromosome CSV file to process.
    output_file The name of the density analysis file to output.
```

## R Scripts

### Daughter Analysis
An R script to assess the differences in chromosome counts between different
offspring cells within a maternal cell. This is used to establish the relative
difference in densities between sibling cells.

### Density Analysis
The density analysis R script is used to generate figures from the CSV output
from the density analysis python script.

### Distance Plotter
An R script to generate top-down heat maps of the average distances between
chromosomes for a given stack.

### Plot Generator
A generic R script to generate top-down figures of the localized chromosomes
for all image stacks.
