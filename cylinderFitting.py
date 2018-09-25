"""

This module fits a cylinder to data points.

Search for two sets of (three) data points.  The geometrical
center of each set determine the bottom and top coordinates.

"""
import random
import itertools
import numpy as np
from scipy.optimize import leastsq

class CylinderFitting(object):
    """
    Properties
    -----------

    """
    def __init__(self, xyz):

        self.center = center = xyz.mean(0)
        theta, phi = 0, np.pi/2
        self.W = np.array([np.cos(theta)*np.cos(phi),
                           np.sin(theta)*np.cos(phi),
                           np.sin(phi)]) # W is the axis cylinder

        params = [self.center[0], self.center[1], theta, phi,
                 (np.linalg.norm((self._project(xyz)))).mean()]     # params[4] = radius

        estParams , success = leastsq(self._cylinderFitting, params, args=(xyz,))
        self.center[0:2] = estParams[0:2]
        theta, phi = tuple(estParams[2:4])
        self.r = estParams[4]
        self.W = np.array([np.cos(theta)*np.cos(phi),np.sin(theta)*np.cos(phi),np.sin(phi)])
        self.bottom, self.top = self._computeExtremes(xyz)

    def _cylinderFitting(self,params, xyz):

        """
        Reference:

        params are variables used for computing the error function in the
        fitting procedure

        params[0] = x coordinate of the cylinder centre
        params[1] = y coordinate of the cylinder centre
        params[2] = theta, rotation angle about the z-axis
        params[3] = phi, orientation angle of the plane with normal vector W
        params[4] = r, radius of the cylinder

        xyz are the points to fit

        """
        x, y, theta, phi, r = tuple(params)
        deviation = np.zeros(xyz.shape[0])
        z = xyz[:,2].mean()
        center = np.array([x, y, z])
        W = np.array([np.cos(theta)*np.cos(phi),
                      np.sin(theta)*np.cos(phi),
                      np.sin(phi)])
        deviation = []
        for i in range(xyz.shape[0]):
            deviation.append(np.dot(self._project(xyz)[i,:],xyz[i,:]) - r**2)

        return deviation

    def _project(self,xyz):

        plane = np.identity(3) - np.dot(self.W[:,np.newaxis],self.W[np.newaxis,:])
        projection = []
        for data in range(xyz.shape[0]):
            vector = xyz[data,:] - self.center
            projection.append(np.dot(vector,plane))

        return np.asarray(projection)

    def _computeExtremes(self,xyz):

        heights = []
        for data in range(xyz.shape[0]):
            heights.append(np.inner((xyz[data,:]-self.center),self.W))
        hpoints = np.asarray(heights)
        bottom, top = self.center + min(hpoints)*self.W, self.center + max(hpoints)*self.W
        height = max(hpoints) - min(hpoints)
        self.center = (top + bottom)/2 # Recalculate the cylinder center
        bottom, top = self.center - (height/2)*self.W, self.center + (height/2)*self.W

        return bottom, top

    def _axisCylinder(self,xyz):

        p = self._project(xyz).tolist()
        V = random.sample(p, 1)/np.linalg.norm(random.sample(p, 1))
        U = np.cross(V,self.W)

        return U, V

    def vmdCommands(self):

        commands  = "set bottom {{ {} {} {} }}\n".format(*self.bottom)
        commands += "set top {{ {} {} {} }}\n".format(*self.top)
        commands += "draw material Transparent\n"
        commands += "draw color silver\n"
        commands += "draw cylinder $bottom $top radius {} resolution 100\n".format(self.r)

        return commands

    def atomsInExtremes(self, xyz, nsectors, weights=None):

        """

        Finds two combinations of points whose (possibly weighted) mean coordinates are
        the closest to the cylinder extremities. The points are split into a number of sectors
        according to their azimuth angles and the returned combinations will contain one particle
        of each sector.

        """

        delta = 2*np.pi/nsectors
        sector = [[] for i in range(nsectors)]
        U, V = self._axisCylinder(xyz)
        for i in range(xyz.shape[0]):
            xyz_u = np.inner((xyz[i,:]- self.center), U)
            xyz_v = np.inner((xyz[i,:]- self.center), V)
            angle = np.arctan2(xyz_u, xyz_v)
            isec = int((angle + np.pi)/delta)
            sector[isec].append(i)

        d0_b = d0_t = np.Inf
        for comb in itertools.product(*sector):
            masses = None if weights is None else weights[list(comb)]
            coords = np.average(xyz[list(comb),:], axis=0, weights=masses)
            d_bottom = np.linalg.norm(coords - self.bottom)
            d_top = np.linalg.norm(coords - self.top)
            if d_bottom < d0_b:
                d0_b = d_bottom
                bottom_atoms = comb
            if d_top < d0_t:
                d0_t = d_top
                top_atoms = comb

        return bottom_atoms, top_atoms

    def writeExtremesCoords(self, xyz, bottom_atoms, top_atoms, file):

        bottom_center = xyz[[int(i) for i in bottom_atoms], :].mean(0)
        top_center = xyz[[int(i) for i in top_atoms], :].mean(0)
        extremes = open(file,'w')
        extremes.write("{}\n\n".format(8))
        extremes.write("C {} {} {}\n".format(*bottom_center))
        extremes.write("C {} {} {}\n".format(*top_center))
        for i in bottom_atoms + top_atoms:
            extremes.write("C {} {} {}\n".format(*xyz[i]))
