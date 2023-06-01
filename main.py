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
from scipy import interpolate
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from tkinter import *
from tkinter import filedialog
class Material():
    def __init__(self, tensile_, shear_, tensile_strength_, shear_strength_):
        self.Tensile_Modulus = tensile_
        self.Shear_Modulus = shear_
        self.Tensile_Strength = tensile_strength_
        self.Shear_Strength = shear_strength_
class CFD():
    """This is an implementation of the Lattice Boltzmann simulation"""
    '''
    this should hopefully output a list of forces on each eadge of the foil, or the force applied to each vertex on the surface
    '''
    def __init__(self, Density_, Kinematic_Viscocity_, V_inf_, foil_):
        self.Density = Density
        self.mu = Kinematic_Viscocity
        self.V_inf = V_inf_
        self.foil = foil_

class FEA_Solver():
    def __init__(self):
        pass

class Lattice_Optimisation():
    def __init__(self):
        pass
    def Generatre_points(self, boundaries):
        '''Generates 1000 points within the foils bounding box'''
        self.surface = boundaries[0]#selecting airfoil boundary
        self.x_min = np.min(self.surface[:,0])
        self.x_max = np.max(self.surface[:,0])
        self.y_min = np.min(self.surface[:,1])
        self.y_max = np.max(self.surface[:,1])
        self.points = np. random.random((1000,2))
        self.points = self.points*np.array([self.x_max-self.x_min, self.y_max-self.y_min])
        self.points = self.points+np.array([self.x_min,self.y_min])
        '''Point in polygon to remove any points that exist in empty regions'''
        #for boundary in boundaries:

    def Number_points(self, array):
        pass
    def Do_Optimisation(self, wing):# this generates the 'infill' pattern
        '''
        in this context each sections end is a boundary to fill with the infill pattern
        '''
        pass
        #Random points
        #Points in polygon
        #point tirangularisation


class Wing():
    '''Working'''
    def __init__(self, surface_, y_, chord_, sweep_, dihedral_, twist_, spars_):
        self.Sections = [Foil(surface_[i], chord_[i], twist_[i]) for i in range(len(chord_))]
        self.Section_angles = [[sweep_[i], dihedral_[i], y_[i]] for i in range(len(chord_))] #Sweep angle, dihedral angle, and the distance from the plater where it happens
        self.y_array = y_
        self.sweep_array = sweep_
        self.dihedral_array = dihedral_
        self.twist_array = twist_
        self.Hole_Elongation = [[1,4,1],[1,4,1],[1,4,1]]
        self.Hole_Sizes = [[10,10,5],[10,10,5],[10,10,5]]
        self.Hole_Loctions = [[0.3,0.5,0.6],[0.3,0.5,0.6],[0.3,0.5,0.6]]
    def curve_model(self,x, a, b, c, d):
        '''Defines the model used during the curve fitting step'''
        return (a*x**3)+(b*x**2)+(c*x)+d
    def find_chamber_point(self, point, foil_surface):#Kidna a pain this is a duplicate
        '''Find cartesian location of a point a percentage of a way along a foils chamber line using NLSR'''
        self.Coefficients, _ = curve_fit(self.curve_model, foil_surface[:,0], foil_surface[:,1], p0 = [1, 1, 1, 1]) #Computes the optimum Coefficients to fit the curve model to the foils points using NLSS
        self.a, self.b, self.c, self.d = self.Coefficients
        self.Sample_points = np.array([[x/1000, self.curve_model(x/1000,self.a, self.b, self.c, self.d)] for x in range(int(min(foil_surface[:,0])*1000),int(max(foil_surface[:,0])*1000))])
        self.Steps = np.diff(self.Sample_points,axis=0)
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
    def generate_offsets(self):
        '''Returns an array of n, 2d arrays that describes how each layer must change'''
        self.Linear_offsets = [[self.y_array[i]*np.sin(self.sweep_array[i]),self.y_array[i]*np.sin(self.dihedral_array[i])] for i in range(len(self.y_array))]
        self.Linear_offsets = [np.sum(self.Linear_offsets[0:i],axis=0) for i in range(1,len(self.Linear_offsets)+1)]
    def generate_hole(self, radius, elongaton, centre):
        pass
    def add_spars(self, layers):
        self.locations = []
        for self.layer_index in range(0, len(layers)):
            self.boundary = layers[self.layer_index]
            for self.hole_index in range(0,len(self.Hole_Sizes)):
                self.hole_size = self.Hole_Sizes[self.layer_index, self.hole_index]
                self.hole_location = self.Hole_Loctions[self.layer_index, self.hole_index]
                self.hole_warp = self.Hole_Elongation[self.layer_index, self.hole_index]
                self.hole_defomation_factor = np.tan(self.Section_angles[self.layer_index,0:1])+1
                
                self.cart_position = self.find_chamber_point(self.hole_location, self.boundary)
                self.hole = self.generate_hole()
        self.locations = [[self.find_chamber_point(point, boundary),self.find_chamber_point(point, boundary) ]for point in self.Hole_Loctions for boundary in layers] #spar locations in each layer
        print(self.locations)

    def get_wing(self):
        '''Returns n, 2d arrays that descirbe the foil'''
        self.generate_offsets()

        '''The wing layers are each interpolation point along the wing, the slicer essentialy lofts between each of these layers'''
        self.wing_layers =  np.array([np.hstack((self.Sections[i].get_foil()+self.Linear_offsets[i],np.full((1000,1),self.y_array[i])))  for i in range(len(self.Sections))])
        self.add_spars(self.wing_layers)
class Foil():
    '''Working'''
    def __init__(self, points, chord, twist):
        self.Surface = self.resample(points,1000) #incoming foils are resampled
        self.Chord = chord
        self.angle = twist
    def resample(self, surface_,  samples):
        '''resamples the foil at a higher density'''
        self.x_points = surface_[:,0] #Extracts an array of all x points
        self.y_points = surface_[:,1] #Extracts an array of all y points
        self.tck , self.u = interpolate.splprep([self.x_points,self.y_points],k=5,s=0.1, per=True) #From a parameterised interpolation of the points given
        self.i_x, self.i_y = interpolate.splev(np.linspace(0,1,samples),self.tck) # samples n points along the interpolation
        self.resampled_foil = np.stack((self.i_x,self.i_y),axis=1) #Stacks all the individual samples in to 1 2d array
        return self.resampled_foil
    def cart2polar(self, point):
        '''Converts from cartesian to polar'''
        self.rho = np.sqrt(point[0]**2 + point[1]**2)
        self.phi = np.arctan2(point[1], point[0])
        return [self.rho, self.phi]
    def pol2cart(self, rho, phi):
        '''Converts from polar to cartesian'''
        self.x = rho * np.cos(phi)
        self.y = rho * np.sin(phi)
        return [self.x, self.y]
    def curve_model(self,x, a, b, c, d):
        '''Defines the model used during the curve fitting step'''
        return (a*x**3)+(b*x**2)+(c*x)+d
    def find_chamber_point(self, point, foil_surface):
        '''Find cartesian location of a point a percentage of a way along a foils chamber line using NLSR'''
        self.Coefficients, _ = curve_fit(self.curve_model, foil_surface[:,0], foil_surface[:,1], p0 = [1, 1, 1, 1]) #Computes the optimum Coefficients to fit the curve model to the foils points using NLSS
        self.a, self.b, self.c, self.d = self.Coefficients
        self.Sample_points = np.array([[x/1000, self.curve_model(x/1000,self.a, self.b, self.c, self.d)] for x in range(int(min(foil_surface[:,0])*1000),int(max(foil_surface[:,0])*1000))])
        self.Steps = np.diff(self.Sample_points,axis=0)
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
        self.axis = self.find_chamber_point(axis, self.Surface) #Compute chord line
        #Rotates each point about the selected point
        self.working_foil = np.empty((1,2))
        for self.point in self.Surface:
            self.offest_point = self.point - self.axis
            self.polar_point = self.cart2polar(self.offest_point)
            self.polar_point = [self.polar_point[0],self.polar_point[1]-self.angle]
            self.new_point = self.pol2cart(self.polar_point[0], self.polar_point[1])
            self.working_foil = np.vstack((self.working_foil, self.new_point))
        return self.working_foil[1:] +self.axis
    def get_foil(self):
        '''returns a list of points that describe the foil'''
        self.delivered_foil = self.set_angle() #Rotate foil
        self.delivered_foil = np.array(self.delivered_foil) * self.Chord      #Scale foil
        self.centre = self.find_chamber_point(0.3, self.delivered_foil)
        self.delivered_foil = self.delivered_foil - self.centre
        return self.delivered_foil             #return foil




class Printer():
    def __init__(self):
        self.Dimensions = []#Dimensions in mm
        self.Speed = 60#mm
        self.Filament_Diamiter = 1.75#mm
        self.Nozzle_Diamiter = 0.4
        self.Flow_Factor = ((self.Filament_Diamiter/2)**2)/((self.Nozzle_Diamiter)**2)
class Slicer():
    def __init__(self):
        self.Layer_height = 0.001
        ## TODO: Remember this will have to split the foil in to section of different heights so it fits on the y axis
    def Slice(self, model):
        '''Genrates layers that can be turned in to gcode'''
        self.layer_paths = np.empty((1,1000,3))
        for index in range(0,len(model)-1):
            self.layers = (model[index+1,0,2]-model[index,0,2])/self.Layer_height
            self.step = (model[index+1]-model[index])/self.layers
            self.paths = np.array([model[index]+self.step*layer for layer in range(0,int(self.layers)+1)])
            self.layer_paths = np.vstack((self.layer_paths,self.paths))
        return self.layer_paths

    def Genrate_Gcode(self):
        '''Turns the layer slices in to gcode commands'''

class GUI(Tk):
    def __init__(self, parent):
        Tk.__init__(self, parent)
        self.parent = parent
        '''Window setup'''
        self.title("Gen-Slice")
        self.geometry("1200x900+20+20")
        '''Setup interface'''


        '''Runs program'''
        self.mainloop()
class main():
    def __init__(self):
        self.Load_airframe_defintion()
        '''Loads foils'''
        self.foil = self.Load_foil()
        '''Build slicer instance'''
        self.slicer = Slicer()
        '''Defines Wings'''
        self.test_wing = Wing([self.foil, self.foil, self.foil],
                              [0.0, 0.2, 0.40], #y values,All array must have a leading 0 to set the reference for the bed profile
                              [0.4, 0.2, 0.15],#Chord lengths in meters
                              [0,1,0.05], #Sweep angles
                              [0,0.02,0.02], #dihedral angles
                              [0.1,0.2,0.3],
                              []) #Twist angles
        '''Test of generating wing'''
    #    ax = plt.figure().add_subplot(projection='3d')
    #    for i in self.test_wing.get_wing():
    #        ax.plot(i[:,0],i[:,1],i[:,2])
        #plt.show()
        self.Wing = self.test_wing.get_wing()


        #self.slicer.Slice(self.test_wing.get_wing())
        '''Aircraft are composed of wings which are composed of foils'''
        self.Wings = []
    def Load_foil(self):
        self.raw = open('Airfoils\s1225_norm.dat','r').read()
        self.formatted = [self.el.split(" ") for self.el in self.raw.split('\n')]
        self.formatted  = [[float(self.num) for self.num in list(filter(lambda x:x!='',self.coord))]for self.coord in self.formatted]#list(map(float,self.formatted))
        self.formatted = np.array(list(filter(lambda x:x!=[],self.formatted)))
        return self.formatted
    def Load_airframe_defintion(self):
        """XFLR5 Stores its aircraft definitions a as XML file"""

        #Waiting on response from forum
        pass
if __name__ == "__main__":
    main()
