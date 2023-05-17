"""
This is a shity first attempt at a slicer that incorperates both CFD and
generative design to allow for the rapid production of 3d printed airframes
the aim of this is to allow for easier access to the hobby through the use of 3D printers
without having to do extensive design work
"""

import numpy as np
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
class CFD():
    """This is an implementation of the Lattice Boltzmann simulation"""
    def __init__(self, Density_, Kinematic_Viscocity_, V_inf_, foil_):
        self.Density = Density
        self.mu = Kinematic_Viscocity
        self.V_inf = V_inf_
        self.foil = foil_

class Wing(Foil):
    def __init__(self, surface_, y_, chord_, sweep_, dihedral_):
        self.Sections = [Foil(surface_[i], chord_[i]) for i in range(len(chord_))]
        self.Section_angles = [[sweep_[i], dihedral_[i], y_[i]] for i in range(len(chord_))] #Sweep angle, dihedral angle, and the distance from the plater where it happens
    def generate_offsets(self):
        '''Returns an array of n, 2d arrays that describes how each layer must change'''
    def get_wing(self):
        '''Returns n, 2d arrays that descirbe the foil'''
        #Compute linear offsets


class Foil():
    def __init__(self, points, chord):
        self.Surface = self.resample(points,1000) #incoming foils are resampled
        self.Chord = chord
        self.angle = 0
    def resample(self, surface_,  samples):
        '''resamples the foil at a higher density'''
        self.x_points = surface_[:,0] #Extracts an array of all x points
        self.y_points = surface_[:,1] #Extracts an array of all y points
        self.x_points = r_[self.x_points,self.x_points[0]] # Adds the first x point to the end to form a closed curve
        self.y_points = r_[self.y_points,self.y_points[0]] # Adds the first y point to the end to form a closed curve
        self.tck , self.u = interpolate.splprep([self.x_points,self.y_points],k=5,s=0.1,per=True) #From a parameterised interpolation of the points given
        self.i_x, self.i_y = interpolate.splev(linspace(0,1,samples),self.tck) # samples n points along the interpolation
        self.resampled_foil = stack((self.i_x,self.i_y),axis=1) #Stacks all the individual samples in to 1 2d array
        return self.resampled_foil
    def cart2polar(self, point):
        '''Converts from cartesian to polar'''
        self.rho = np.sqrt(point[0]**2 + point[1]**2)
        self.phi = np.arctan2(point[1], point[0])
        return [self.rho, self.phi]
    def pol2cart(rho, phi):
        '''Converts from polar to cartesian'''
        self.x = rho * np.cos(phi)
        self.y = rho * np.sin(phi)
        return [self.x, self.y]
    def curve_model(self,x, a, b, c, d):
        '''Defines the model used during the curve fitting step'''
        return (a*x**3)+(b*x**2)+(c*x)+d
    def find_chamber_point(self, point):
        '''Find cartesian location of a point a percentage of a way along a foils chamber line using NLSR'''
        self.Coefficients, _ = curve_fit(self.curve_model, self.Surface[:,0], self.Surface[:,1], p0 = [1, 1, 1, 1]) #Computes the optimum Coefficients to fit the curve model to the foils points using NLSS
        self.a, self.b, self.c, self.d = self.Coefficients
        self.Sample_points = [[x/1000, self.curve_model(x/1000,self.a, self.b, self.c, self.d)] for x in range(0,1000)]
        self.Steps = np.diff(self.Sample_points)
        self.Step_lengths = np.hypot(self.Steps[:,0], self.Steps[:,1])
        self.Chamber_Length = np.sum(self.Step_lengths) # This is the total length of the chamber line
        self.chamber_point_length = point*self.Chamber_Length
        self.integrator = 0
        i=0
        while self.integrator<self.chamber_point_length:
            self.integrator = self.integrator +self.Step_lengths[i]
            i+=1
        self.chamber_point = [i/1000, self.curve_model(i/1000, self.a, self.b, self.c, self.d)]
        return self.chamber_point
    def set_angle(self, axis = 0.25):
        '''Rotates the foil about a given point on the chord line'''
        self.axis = self.find_chamber_point(axis) #Compute chord line
        #Rotates each point about the selected point
        self.working_foil = []
        for self.point in self.Surface:
            self.offest_point = self.point - self.axis
            self.polar_point = self.cart2polar(self.offest_point)
            self.polar_point +=[0,self.angle]
            self.new_point = self.pol2cart(self.polar_point)
            self.working_foil.append(self.new_point)
        return self.working_foil +self.axis
    def get_foil(self):
        '''returns a list of points that describe the foil'''
        self.delivered_foil = self.set_angle() #Rotate foil
        self.delivered_foil *= self.Chord      #Scale foil
        return self.delivered_foil             #return foil

class Slicer():
    def __init__(self):
        self.Layer_height = 0.1




class main():
    def __init__(self):
        self.Load_airframe_defintion()

    def Load_airframe_defintion(self):
        """XFLR5 Stores its aircraft definitions a as XML file"""

        #Waiting on response from forum
        pass
if __name__ == "__main__":
    main()
