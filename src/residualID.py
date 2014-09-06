#! /usr/bin/env python
#@Author Sean Mclaughlin

desc='''
This module analyzes the residuals left in the image. It clusters them
and attempts to ID if the object likely represents a lens. 
The function residualID can be imported.

residualID(image, c_x, c_y)
image: numpy array of the residuals
c_x, c_y: the center of the image.

Cannot currently be run as main for test cases.
'''

import numpy as np
from sklearn.cluster import DBSCAN

class Residual(object):
#object that keeps residual information in one place
    def __init__(self, pointsSet, key , image, img_cx, img_cy):
        
        self.pointsSet = pointsSet #a set of all points in tuples, for index access
        #(y,x) like how nympy does it
        self.points = np.array(list(pointsSet))

        self.key = key

        if self.key == -1:
            self.key = 'Noise'

        self.image = image
        self.imageCenter = (img_cy, img_cx)

        self.center, self.dist, self.vect, self.theta = self._findMoments()
        self.a, self.b = self.vect

        #slope from the object's center to the image center
        self.slope = (self.center[1]-imageCenter[1])/(self.center[0]-self.imageCenter[1])

    def __contains__(self, item):# for point in residual checks
        return item in self.pointsSet

    def __iter__(self):#return the tuples in an iterator
        return self.pointsSet__iter__()

    def isNoise(self):
        return self.key == 'Noise'

    def getXY(self): #return the x and y as separate arrays
        return self.points[:,1], self.points[:,0]

    def _findMoments(self):#find the moments of the distribution
        Ix = Iy = Ixx = Iyy = Ixy = I = long(0)
        
        for idx in self.pointsSet:
            y,x = idx
            val = self.image[idx]
            I+=val
            Ix+=x*val
            Iy+=y*val
            Ixx+=x*x*val
            Iyy+=y*y*val 
            Ixy+=x*y*val 
        xmean = Ix/I
        ymean = Iy/I
              
        img_y, img_x = self.image.shape
        imgCenterY, imgCenterX = img_y/2.0, img_x/2.0
        
        #the distane of the center of the light from the center of the distribution
        dist = np.sqrt((xmean-imgCenterX)**2+(ymean-imgCenterY)**2)
        
        Uxx = Ixx/I-xmean**2
        Uyy = Iyy/I-ymean**2
        Uxy = Ixy/I-xmean*ymean
        
        theta = .5*np.arctan(2.0*Uxy/(Uxx-Uyy))
        
        lambda1 = .5*(Uxx+Uyy)+.5*np.sqrt(4*Uxy**2+(Uxx-Uyy)**2)
        lambda2 = .5*(Uxx+Uyy)-.5*np.sqrt(4*Uxy**2+(Uxx-Uyy)**2)
        
        return (ymean, xmean), dist, (np.sqrt(lambda1), np.sqrt(lambda2)), theta

def makeResidualDict(image, c_x, c_y): #makes a dict of Residual objects for an image
    std = image.std()
    mean = image.mean()
    threshold = mean+std #noise v. signal threshold

    data = []
    img_y, img_x = image.shape
    for x in xrange(img_x):
        for y in xrange(img_y):
            if image[y,x]>threshold:
                data.append([y,x])

    data = np.array(data)
    
    #cluster the points
    agg = DBSCAN(eps = 2, min_samples = 6)

    agg.fit(data)
    labels = agg.labels_

    #gather the clusters into the dictionary
    nClusters = len(set(labels))
    clustRange = xrange(nCluster)
    if -1 in labels:
    #noise cluster ID'd
        clustRange = xrange(-1, nClusters-1)

    residuals = {}
    for i in clustRange:
        #gathe the poitns by cluster and create their Residual object
        cluster = set()
        d = data[labels==i]
        for idx in d:
            cluster.add(tuple(idx))

        #create the objects and put them in the dictionary
        if i == -1:
            residuals['Noise'] = Residual(cluster, 'Noise', image, c_x, c_y)
        else:
            residuals[i] = Residual(cluster, i, image, c_x, c_y)

    return residuals

def checkLens(residuals):
#return boolean if system is likely a lens
#TODO have a probability rather than a binary answer
    nObjects = len(residuals)
    if 'Noise' in residuals: #ignore noise
        nObjects -=1

    if nObjects<1:
        return False#no residuals detected

    image = residuals[0].image

    if nObjects ==1:
        residual = residuals[0]
        center_difference = 3 #allowed distance from center
        return residual.dist<=center_difference

    else: #nObjects >=2

        enoughPerp = False #enough clusters are perpendicular to the center
        closeTogether = False #objects comparable distance from the center

        imgCy, imgCx = residuals[0].c_y, residuals[0].c_x

        objs = [residuals[i] for i in xrange(nObjects)]
        slopes = [obj.slope for obj in objs]

        properlyOriented = 0
        for obj, slope in zip(objs, slopes):
            theta = -1*obj.theta#seems to come out negative. not sure why
            axesSlope = 1/np.tan(theta)
            #maybe I should be specifically checking the major axis
            properlyOriented+=1 if .1<abs(axesSlope/slope)<10 else 0

        enoughPerp = properlyOriented>1#if there is more than one properly oriented

        diffs = []
        for i, obj in enumerate(objs):
            for obj2 in objs[i+1:]:
                diffs.append(abs(obj.dist-obj2.dist))

        diffs = np.array(diffs)
        closeTogether = diffs.mean()<3 #arbitrary number

        if nObjects ==2: #special case
            slope1, slope2 = slopes[0], slopes[1]
            ratio = slope1/slope2

            #check if the objects are acorss the center from one another
            acrossOneAnother = ratio>-.2 and .1<abs(ratio)<10
            return (acrossOneAnother and enoughPerp) or (acrossOneAnother and closeTogether \
                                    or enoughPerp and closeTogether)#any 2

        return closeTogether and enoughPerp

def residualID(image, c_x, c_y):
   residuals = makeResidualDict(image, c_x, c_y) 
   return checkLens(residuals)

if __name__ == '__main__':
    print desc
