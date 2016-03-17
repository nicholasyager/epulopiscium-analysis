#!/usr/bin/python2

import atexit
import argparse
import cv2
import cv2.cv as cv
import math
import numpy as np
import os
import sys
from scipy import stats

parser = argparse.ArgumentParser(description="Count the number of chromosomes"+
                                 " in a stack of dark field images.")
parser.add_argument('stackPath', metavar='stackPath')
parser.add_argument('--batch', action="store_true",default=False)

args=parser.parse_args()

DEBUG = True

np.seterr(invalid = "raise")


if not DEBUG:
    class ExitHooks(object):
        def __init__(self):
            self.exit_code = None
            self.exception = None

        def hook(self):
            self._orig_exit = sys.exit
            sys.exit = self.exit
            sys.excepthook = self.exc_handler

        def exit(self, code=0):
            self.exit_code = code
            self._orig_exit(code)

        def exc_handler(self, exc_type, exc, *args):
            self.exception = exc

    hooks = ExitHooks()
    hooks.hook()


class Printer():
    """
    Print things to stdout on one line dynamically
    """
 
    def __init__(self,data):
 
        sys.stdout.write("\r\x1b[K"+data.__str__())
        sys.stdout.flush()


def distance(xi,yi,zi, xj, yj,zj):

    return math.sqrt( (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2  )

def distance2D(xindex, yindex, indices, xStdDev, yStdDev):

    distances = []
        
    for rowIndex, colIndex in indices:
        distance = math.sqrt( ((xindex - colIndex)**2 / xStdDev**2) +
                              ((yindex - rowIndex)**2 / yStdDev**2) )
        distances.append(distance)

    return sum(distances)/len(distances)

def umToPixels(image, coords):
    """
    Convert um measurements to coordinates.
    """
    coefs = [108.66, 108.66, 0.5]
    shape = image.shape
    newCoords = []
    for index, coord in enumerate(coords):
        newCoords.append(coord * shape[index]/ coefs)

    return(tuple(newCoords))

def main():
   
    if args.batch:
        print("Batch mode!")
        rootPath = args.stackPath
        
        stackPaths = sorted(os.listdir(rootPath))

        for stack in stackPaths:
            stackPath = rootPath +stack
       
            try:
                imagePaths = sorted(os.listdir(stackPath))
            except OSError:
                print("Skipping "+stackPath+". ");
                continue

            print(stackPath)

            #cv2.namedWindow( "Hough Circle Transform Demo", cv2.CV_WINDOW_AUTOSIZE )

            os.makedirs(stackPath+"_analysis")

            points4D = []
            values = {}

            outputFile = open(stackPath+"_analysis/points.csv","w")

            print("Processing images.")
            for index, imagePath in enumerate(imagePaths):
                Printer("\tProcessing image {0}/{1}".format(index+1,len(imagePaths)))
                orig_image = cv2.imread(stackPath+"/"+imagePath)
                image = cv2.cvtColor(orig_image,cv2.COLOR_BGR2GRAY)

                # Equalize the histogram
                

                # Threshold the image
                image = cv2.equalizeHist(image)
                ret,thresh1 = cv2.threshold(image,254,255,cv2.THRESH_TOZERO)
                kernel = np.ones((2,2),np.uint8)
                erode = cv2.erode( thresh1,kernel, iterations=1);
                equal = cv2.equalizeHist(erode)
                contours, heirarchy = cv2.findContours(equal,cv2.RETR_TREE, 
                                                       cv2.CHAIN_APPROX_SIMPLE)
                
                points = []

                shape = equal.shape

                for contour in contours:
                    M = cv2.moments(contour)
                    try:
                        cx = int(M['m10']/M['m00'])
                    except ZeroDivisionError:
                        cx = contour[0][0,0]
                    try:
                        cy = int(M['m01']/M['m00'])
                    except ZeroDivisionError:
                        cy = contour[0][0,1]
                    points.append((cx, cy))

                false_image = cv2.cvtColor(equal, cv.CV_GRAY2RGB)

                for x, y in points:
                    
                    # Add pixel to um conversion
                    z = (len(imagePaths)-index)

                    points4D.append((x,y,z))
                    values[(x,y,z)] = image[x,y]
                    

                    # Find the average cluster distance    
                    #distance = distances(colIndex, rowIndex, indices, xStdDev, yStdDev)


                    #orig_image[rowIndex,colIndex] = (0,0,distance)i
                    cv2.circle(false_image,(x,y),3, (0,255,0))

            
                #cv2.imshow("Hough Circle Transform Demo",false_image)
                #cv2.waitKey(1)
                cv2.imwrite("{0}_analysis/{1:03d}.tiff".format(stackPath,index),
                            false_image)

            #exit()

            # Perform a clustering analysis to remove duplicates.
            
            print("\nRemoving duplicate chromosomes.")
            for index, key in enumerate(values.keys()):
                x,y,z = key
                Printer("\tChecking chromosome {0}".format(index+1))
                distances = []
                continueLock = True

                for xindex in range(-3,3):
                    if not continueLock:
                        break
                    for yindex in range(-3,3):
                        if not continueLock:
                            break
                        for zindex in range(-3,3):
                            if xindex == 0 and yindex == 0 and zindex == 0:
                                continue
                            if index == 0 or index == len(values.keys()) - 1:
                                continue
                            if (x+xindex,y+yindex,z+zindex) in values:
                                try:
                                    if values[(x+xindex,y+yindex,z+zindex)] > values[(x,y,z)]:
                                        del values[(x,y,z)]
                                        points4D.remove((x,y,z))
                                        continueLock = False
                                        break
                                    else:
                                        del values[(x+xindex,y+yindex,z+zindex)]
                                        points4D.remove((x+xindex,y+yindex,z+zindex))
                                except KeyError:
                                    continueLock = False
            


            # Perform a clustering analysis to remove outliers.
            avgs = {}

            print("\nCalculating distances.")
            #medianFile = open("medians.csv","w")
            for index, coords in enumerate(points4D):
                x,y,z = coords
               
                x1 = x * 108.66 / shape[0]
                y1 = y * 108.66 / shape[1]
                z1 = z * 0.5

                Printer("\tMeasuring chromosome {0}/{1}".format(index+1,len(points4D)))
                distances = []
                for xj, yj, zj in points4D:
                
                    xj1 = xj * 108.66 / shape[0]
                    yj1 = yj * 108.66 / shape[1]
                    zj1 = zj * 0.5
    
                    if xj1 > x1 + 10 or xj1 < x1 - 10:
                        continue
                    elif yj1 > y1 + 10 or yj1 < y1 - 10:
                        continue
                    elif zj1 > z1 + 10 or zj1 < z1 - 10:
                        continue

                    if (x1,y1,z1) == (xj1,yj1,zj1):
                        continue
                    dist = distance(x1,y1,z1,xj1,yj1,zj1)
                    distances.append(dist)
                if len(distances) < 1:
                    avgs[(x,y,z)] = -1
                else:
                    if len(distances) < 10:
                        numDist = len(distances)
                    cluster = sorted(distances)[:10]
                    avgs[(x,y,z)] = stats.nanmean(np.asarray(cluster))

       
            avgslist = [ x for y, x in avgs.iteritems()]
            #print(avgslist)
            median = stats.nanmedian(np.asarray(avgslist))
            print("\tMedian: {0}".format(median))
            stddev = stats.nanstd(np.asarray(avgslist))
            print("\tStd. Dev.: {0}".format(stddev))

            """
            print("\nTrimming by distance.")
            index = 1
            removeList = []
            for key,avg in avgs.iteritems():
                x,y,z = key
                Printer("\tMeasuring chromosome {0}/{1}".format(index,len(points4D)))
                index += 1
                if avg >= median + stddev*3:
                   removeList.append((x,y,z))

            for key in removeList:
                points4D.remove(key)
                del values[key]
                del avgs[key]
            """

            # Print all of the values to a csv file
            print("\nSaving all values")
            outputFile.write('"x","y","z","value","dist"\n')
            for x, y, z in values.keys():

                try:
                    value=values[(x,y,z)]
                    avg = avgs[(x,y,z)]

                    x = x * 108.66 / shape[0]
                    y = y * 108.66 / shape[1]
                    z = z * 0.5

                    outputFile.write("{0},{1},{2},{3},{4}\n".format(x,y,z,value,avg))
                except KeyError:
                    print("Index not found!!!!1!1!!1!")
                    pass


            outputFile.close()
      
    else:

        imagePaths = sorted(os.listdir(args.stackPath))

        #cv2.namedWindow( "Hough Circle Transform Demo", cv2.CV_WINDOW_AUTOSIZE )

        os.makedirs(args.stackPath+"_analysis")

        points4D = []
        values = {}

        outputFile = open(args.stackPath+"_analysis/points.csv","w")

        print("Processing images.")
        for index, imagePath in enumerate(imagePaths):
            Printer("\tProcessing image {0}/{1}".format(index+1,len(imagePaths)))
            orig_image = cv2.imread(args.stackPath+"/"+imagePath)
            image = cv2.cvtColor(orig_image,cv2.COLOR_BGR2GRAY)

            # Equalize the histogram
            

            # Threshold the image
            image = cv2.equalizeHist(image)
            ret,thresh1 = cv2.threshold(image,254,255,cv2.THRESH_TOZERO)
            kernel = np.ones((2,2),np.uint8)
            erode = cv2.erode( thresh1,kernel, iterations=1);
            equal = cv2.equalizeHist(erode)
            contours, heirarchy = cv2.findContours(equal,cv2.RETR_TREE, 
                                                   cv2.CHAIN_APPROX_SIMPLE)
            
            points = []

            for contour in contours:
                M = cv2.moments(contour)
                try:
                    cx = int(M['m10']/M['m00'])
                except ZeroDivisionError:
                    cx = contour[0][0,0]
                try:
                    cy = int(M['m01']/M['m00'])
                except ZeroDivisionError:
                    cy = contour[0][0,1]
                points.append((cx, cy))

            false_image = cv2.cvtColor(equal, cv.CV_GRAY2RGB)

            for x, y in points:
               
                points4D.append((x,y,len(imagePaths)-index))
                values[(x,y,len(imagePaths)-index)] = image[x,y]
                

                # Find the average cluster distance    
                #distance = distances(colIndex, rowIndex, indices, xStdDev, yStdDev)


                #orig_image[rowIndex,colIndex] = (0,0,distance)i
                cv2.circle(false_image,(x,y),3, (0,255,0))

        
            #cv2.imshow("Hough Circle Transform Demo",false_image)
            #cv2.waitKey(1)
            cv2.imwrite("{0}_analysis/{1:03d}.tiff".format(args.stackPath,index),
                        false_image)

        #exit()

        # Perform a clustering analysis to remove duplicates.
        
        print("\nRemoving duplicate chromosomes.")
        for index, key in enumerate(values.keys()):
            x,y,z = key
            Printer("\tChecking chromosome {0}".format(index+1))
            distances = []
            continueLock = True

            for xindex in range(-3,3):
                if not continueLock:
                    break
                for yindex in range(-3,3):
                    if not continueLock:
                        break
                    for zindex in range(-3,3):
                        if xindex == 0 and yindex == 0 and zindex == 0:
                            continue
                        if (x+xindex,y+yindex,z+zindex) in values:
                            try:
                                if values[(x+xindex,y+yindex,z+zindex)] > values[(x,y,z)]:
                                    del values[(x,y,z)]
                                    points4D.remove((x,y,z))
                                    continueLock = False
                                    break
                                else:
                                    del values[(x+xindex,y+yindex,z+zindex)]
                                    points4D.remove((x+xindex,y+yindex,z+zindex))
                            except KeyError:
                                continueLock = False
        


        # Perform a clustering analysis to remove outliers.
        avgs = {}

        print("\nCalculating distances.")
        #medianFile = open("medians.csv","w")
        for index, coords in enumerate(points4D):
            x,y,z = coords
            x1 = x * 108.66 / shape[0]
            y1 = y * 108.66 / shape[1]
            z1 = z * 0.5

            Printer("\tMeasuring chromosome {0}/{1}".format(index+1,len(points4D)))
            distances = []
            for xj, yj, zj in points4D:
            
                xj1 = xj * 108.66 / shape[0]
                yj1 = yj * 108.66 / shape[1]
                zj1 = zj * 0.5

                if xj1 > x1 + 5 or xj1 < x1 - 5:
                    continue
                elif yj1 > y1 + 5 or yj1 < y1 - 5:
                    continue
                elif zj1 > z1 + 5 or zj1 < z1 - 5:
                    continue

                if (x1,y1,z1) == (xj1,yj1,zj1):
                    continue
                dist = distance(x1,y1,z1,xj1,yj1,zj1)
                distances.append(dist)
            if len(distances) < 1:
                avgs[(x,y,z)] = -1
            else:
                if len(distances) < 4:
                    numDist = len(distances)
                cluster = sorted(distances)[:4]
                avgs[(x,y,z)] = stats.nanmean(np.asarray(cluster))

        avgslist = [ x for y, x in avgs.iteritems()]
        #print(avgslist)
        median = stats.nanmedian(np.asarray(avgslist))
        print("\tMedian: {0}".format(median))
        stddev = stats.nanstd(np.asarray(avgslist))
        print("\tStd. Dev.: {0}".format(stddev))

        """
        print("\nTrimming by distance.")
        index = 1
        removeList = []
        for key,avg in avgs.iteritems():
            x,y,z = key
            Printer("\tMeasuring chromosome {0}/{1}".format(index,len(points4D)))
            index += 1
            if avg >= median + stddev*3:
               removeList.append((x,y,z))

        for key in removeList:
            points4D.remove(key)
            del values[key]
            del avgs[key]
        """

        # Print all of the values to a csv file
        print("\nSaving all values")
        outputFile.write('"x","y","z","value","dist"\n')
        for x, y, z in values.keys():

            try:
                value=values[(x,y,z)]
                avg = avgs[(x,y,z)]
                
                x = x * 108.66 / shape[0]
                y = y * 108.66 / shape[1]
                z = z * 0.5

                outputFile.write("{0},{1},{2},{3},{4}\n".format(x,y,z,value,avg))
            except KeyError:
                pass


        outputFile.close()

if __name__ == "__main__":
    main()
