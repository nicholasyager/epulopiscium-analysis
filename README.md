# Epulopiscium Analysis

**Nicholas A. Yager** and **Matthew Taylor**

This repository contains tools to analyze image stacks of flourescently-labeled
epuloposcoium chromosomes. There are three main tools for this analysis.

1. count.py: Localizes flourescent chromosomes in image stacks.
2. classifier.R: Classifies cell position, identifies independent cells, and
   calculates chromosome density.
3. density\_analysis.R: Statistical analysis of chromosome densities by stage.
4. plot\_generator.R: Script to generate figures.

## Requirements:

The different tools require the following:
* Python 3
    * OpenCV 2.7
    * Numpy
    * Scipy
    * atexit
* R
    * rgl
    * cluster
    * plotrix
    * geometry
    * mclust
    * MASS
    * vioplot
    * ggplot2

## Usage
```
usage: count.py [-h] [--batch] stackPath

Count the number of chromosomes in a stack of dark field images.

positional arguments:
  stackPath   The path to the directory containing the image stack.

optional arguments:
  -h, --help  show this help message and exit
  --batch
```

The three R scripts should be run from the directory containing the count.py
output.
