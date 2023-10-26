"""
This is a shity first attempt at a slicer that incorperates both CFD and
generative design to allow for the rapid production of 3d printed airframes
the aim of this is to allow for easier access to the hobby through the use of 3D printers
without having to do extensive design work
"""
import serial
import time
import numpy as np
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import interpolate
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from tkinter import *
from tkinter import filedialog
from shapely.geometry import Point, Polygon
class Printer():
    def __init__(self):
        self.Dimensions = [250,250,250]#Dimensions in mm
        self.Speed = 60#mm/min
        self.Filament_Diamiter = 1.75#mm
        self.Nozzle_Diamiter = 0.4#mm
        self.Flow_Factor = ((self.Filament_Diamiter/2)**2)/((self.Nozzle_Diamiter)**2)
        self.extrusion_factor = 0.04
        self.bed_size = np.array([220,220,250])
        self.serial_port = serial.Serial('COM18',115200)
        time.sleep(2)
        self.Ready()
    def prime(self):
        self.cmds = ["G0 X2 Y2 Z0.2\r\n",
                     "G1 E7.92 X2 Y200 Z0.2\r\n",
                     "G1 E7.92 X2.1 Y2 Z0.2\r\n",
                     "G0 X2 Y2 Z10\r\n"]
        for cmd in self.cmds:
            self.serial_port.write(cmd.encode())
            self.waittillready()
        time.sleep(3)
    def home(self):
        self.serial_port.write("G28\r\n".encode())
    def retract(self):
        self.cmd = ["G1 E-25\r\n"]
        for cmd in self.cmds:
            self.serial_port(cmd.encode())
    def goto(self, pos):
        print("Going to ", pos)
        self.command = "G0 X"+str(pos[0])+" Y"+str(pos[1])+" Z"+str(pos[2])+"\r\n"
        self.serial_port.write(self.command.encode())#
        self.waittillready()
    def extrudeto(self, x,y,z,e):
        self.command = "G1 E"+str(e)+" X"+str(x)+" Y"+str(y)+" Z"+str(z)+"\r\n"
        self.serial_port.write(self.command.encode())
        self.waittillready()
    def waittillready(self):
        while True:
            if self.serial_port.read_until() == b'ok\n':
                break
    def setNozzleTemperature(self, T):
        print("Extruder heating to "+str(T)+" degrees.")
        self.cmd = "M109 S"+str(T)+"\r\n"
        self.serial_port.write(self.cmd.encode())
        self.waittillready()
        print("Extruder heated.")
    def setBedTemperature(self, T):
        print("Bed heating to "+str(T)+" degrees.")
        self.cmd = "M190 S"+str(T)+"\r\n"
        self.serial_port.write(self.cmd.encode())
        self.waittillready()
        print("Bed heated.")
    def setRelExtrusion(self):
        self.cmd = "M83\r\n"
        self.serial_port.write(self.cmd.encode())
        self.waittillready()
    def setRelMotion(self):
        self.cmd = "G91\r\n"
        self.serial_port.write(self.cmd.encode())
        self.waittillready()
    def setAbsMotion(self):
        self.cmd = "G90\r\n"
        self.serial_port.write(self.cmd.encode())
        self.waittillready()
        self.setRelExtrusion()
    def disableStepper(self):
        self.cmd = "M18\r\n"
        self.serial_port.write(self.cmd.encode())
        self.waittillready()
    def Ready(self):
        self.home()
        self.setRelExtrusion()
        self.setNozzleTemperature(260)
        self.setBedTemperature(60)
        self.prime()
        self.waittillready()
    def print(self, layers):
        '''Input format'''
        '''
        [[Sets], [of], [Paths]<- layer 1
         [Sets], [of], [Paths]]<- layer 2
        '''
        for layer in layers:#this is each vertical slice
            for element in layer:#this is each contrinious line
                self.goto(element[0])
                for index in range(1, len(element)-1):#this is is each point in the line
                    self.diff = [element[index][0]-element[index-1][0], element[index][1]-element[index-1][1]]
                    self.extrusion_distance = round(np.hypot(self.diff[0],self.diff[1])*self.extrusion_factor,5)
                    print(element[index][0],element[index][1],element[index][2],self.extrusion_distance)
                    self.extrudeto(element[index][0],element[index][1],element[index][2],self.extrusion_distance)
        self.home()
class Tools():
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

class CFD():
    """This is wrapper for XFOIL"""
    '''
    this should hopefully output a list of forces on each eadge of the foil, or the force applied to each vertex on the surface
    '''
    def Compute_Max_Loading(self):
        pass


class Material():
    def __init__(self, tensile_, shear_, tensile_strength_, shear_strength_):
        self.Tensile_Modulus = tensile_
        self.Shear_Modulus = shear_
        self.Tensile_Strength = tensile_strength_
        self.Shear_Strength = shear_strength_

class FEA_element(Material):
    def __init__(self, p1, p2, profile):
        Material.__init__(self)
        self.length = np.hypot(p2-p1)
        self.Area = np.product(np.diff(profile))
        self.Area_moment = profile[0]*(1/12)*profile[1]**3
        self.angle = np.arctan2(p2[0]-p1[0],p2[1]-p1[1])
        self.Transposition_matrix = [[np.cos(self.angle), np.sin(self.angle), 0, 0                 ,  0                 , 0],
                                     [np.sin(self.angle), np.cos(self.angle), 0, 0                 ,  0                 , 0],
                                     [0                 , 0                 , 1, 0                 ,  0                 , 0],
                                     [0                 , 0                 , 0, np.cos(self.angle), -np.sin(self.angle), 0],
                                     [0                 , 0                 , 0, np.sin(self.angle),  np.cos(self.angle), 0],
                                     [0                 , 0                 , 0, 0                 ,  0                 , 1]]

        self.Local_stiffness_matrix = [[  (self.Tensile_Modulus*self.Area)/self.length, 0                                                        , 0                                                        ,-(self.Tensile_Modulus*self.Area)/self.length, 0                                                        , 0                                                        ]
                                       [  0                                           , (12*self.Tensile_Modulus*self.Area_moment)/self.length**3, (6 *self.Tensile_Modulus*self.Area_moment)/self.length**2, 0                                           ,-(12*self.Tensile_Modulus*self.Area_moment)/self.length**3, (6 *self.Tensile_Modulus*self.Area_moment)/self.length**2],
                                       [  0                                           , (6 *self.Tensile_Modulus*self.Area_moment)/self.length**2, (4 *self.Tensile_Modulus*self.Area_moment)/self.length   , 0                                           ,-(6 *self.Tensile_Modulus*self.Area_moment)/self.length**2, (2 *self.Tensile_Modulus*self.Area_moment)/self.length   ],
                                       [ -(self.Tensile_Modulus*self.Area)/self.length, 0                                                        , 0                                                        , (self.Tensile_Modulus*self.Area)/self.length, 0                                                        , 0                                                        ],
                                       [  0                                           ,-(12*self.Tensile_Modulus*self.Area_moment)/self.length**3,-(6 *self.Tensile_Modulus*self.Area_moment)/self.length**2, 0                                           , (12*self.Tensile_Modulus*self.Area_moment)/self.length**3,-(6 *self.Tensile_Modulus*self.Area_moment)/self.length**2],
                                       [  0                                           , (6 *self.Tensile_Modulus*self.Area_moment)/self.length**2, (2 *self.Tensile_Modulus*self.Area_moment)/self.length   , 0                                           ,-(6 *self.Tensile_Modulus*self.Area_moment)/self.length**2, (4 *self.Tensile_Modulus*self.Area_moment)/self.length   ]]

        self.Global_stiffness_matrix = np.transpose(self.Transposition_matrix)*self.Local_stiffness_matrix*self.Transposition_matrix

class FEA_Solver():
    def __init__(self):
        self.Nodes = np.array([])
        self.Elements = np.array([])




class TTO_Optimisation():
    def __init__(self):
        pass
    def Generatre_mesh(self, boundaries, inclusion):
        '''Generates 1000 points within the foils bounding box'''
        self.surface = boundaries[0]#selecting airfoil boundary
        self.x_min = np.min(self.surface[:,0])
        self.x_max = np.max(self.surface[:,0])
        self.y_min = np.min(self.surface[:,1])
        self.y_max = np.max(self.surface[:,1])
        self.Points = []
        self.Num_points = 0
        while self.Num_points != 1000:
            '''Generates a point within the bounding box of the wing'''
            self.candid_point = (np.random.random((1,2))*np.array([self.x_max-self.x_min, self.y_max-self.y_min]))+np.array([self.x_min,self.y_min])
            '''Checks if the point is within the estabilished boundaries, this is quite inefficient and can be improved'''
            if all([Point(self.candid_point[0], self.candid_point[1]).within(Poly(boundary)) for boundary in boundaries]==inclusion):
                self.Num_points += 1

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

class Wing(Tools):
    '''Working'''
    def __init__(self, surface_, y_, chord_, sweep_, dihedral_, twist_, Hole_sizes_, Hole_locations_, Hole_elongations_):
        Tools.__init__(self)
        self.Sections = [Wing_section(np.array([surface_[i],surface_[i+1]]),
                                      np.array([y_[i],y_[i+1]]),
                                      np.array([chord_[i],chord_[i+1]]),
                                      np.array([twist_[i],twist_[i+1]]),
                                      np.array([sweep_[i],sweep_[i+1]]),
                                      np.array([dihedral_[i],dihedral_[i+1]]),
                                      Hole_sizes_,
                                      Hole_locations_,
                                      Hole_elongations_) for i in range(0,len(chord_)-1)]
        self.y_array = y_

        ## TODO: Make use of Spars_ fix the hole sizing, location and elongation



    def get_wing(self):
        '''Returns n, 2d arrays that descirbe the foil'''

        self.fig = plt.figure()
        self.ax = plt.axes(projection='3d')
        '''The wing layers are each interpolation point along the wing, the slicer essentialy lofts between each of these layers'''
        for self.Section in self.Sections:
            self.Section.get_wing_section()
            #self.xline = layer[:,0]
            #self.yline = layer[:,1]
            #self.zline = layer[:,2]
            #self.ax.plot3D(self.xline, self.yline, self.zline)
        #plt.show()


class Wing_section():
    def __init__(self, surfaces_, y_, chords_, twists_, sweep_, dihedral_, Hole_sizes_, Hole_locations_, Hole_elongations_):
        #Foil(surfaces_[i], chords_[i], twists_[i])
        self.Foils = np.array([Foil(surfaces_[i], chords_[i], twists_[i]) for i in range(1,len(chords_))])
        self.y_array = y_
        self.chords = chords_
        self.sweep_array = sweep_
        self.dihedral_array = dihedral_
        self.twist_array = twists_
        ## TODO: Make use of Spars_ fix the hole sizing, location and elongation
        self.Hole_Sizes = np.array(Hole_sizes_)#[10,10,5]# main spar, cable way, torque spar
        self.Hole_Elongation = np.array(Hole_elongations_)#[1,2,1]
        self.Hole_Loctions = np.array(Hole_locations_)#[0.3,0.5,0.6]
    def get_section_pressure_dat(self, section):
        self.foil_data = find_pressure_coefficients(airfoil='temp_airfoil_file_Do_Not_Delete', Reynolds = 1e6, alpha=-5,NACA=True)
        return array([self.foil_data["x"],self.foil_data["y"]])
    def generate_hole(self, radius, elongation):
        return (np.array([[radius*np.sin(angle),radius*np.cos(angle)] for angle in np.linspace(-np.pi,np.pi,100)])*elongation)
    def generate_spars(self):
        self.layer_sweep = self.sweep_array[1]
        self.layer_dihedral = self.dihedral_array[1]
        self.hole_defomation_factor = np.array([1+np.tan(self.layer_sweep),1+np.tan(self.layer_dihedral)])
        self.Holes = np.array([self.generate_hole(self.Hole_Sizes[self.index],self.hole_defomation_factor*self.Hole_Elongation[self.index]) for self.index in range(0,len(self.Hole_Sizes))])
        #working to here
        #[self.find_chamber_point(self.point, self.Foil.get_foil()) for self.point in self.Hole_Loctions]
        '''Holes are the same shape and size but in different locations'''
        print(self.chords)
        self.locations = np.array([[self.find_chamber_point(self.point, self.Foils[self.Foil_index].get_foil())*self.chords[self.Foil_index] for self.point in self.Hole_Loctions] for self.Foil_index in range(len(self.Foils)+1)]) #normalised spar locations in each layer

        self.Spars = []
        for self.index in range(self.Holes.shape[0]):
            self.offset_hole = np.array([self.point+self.locations[self.index] for self.point in self.Holes[self.index]])

            self.Spars.append(self.offset_hole)
        self.Spars = np.array(self.Spars)
        #print(self.Spars)
        for self.index in range(self.Spars.shape[0]):
            plt.plot(self.Spars[self.index][:,0],self.Spars[self.index][:,1])
        plt.show()



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
    def curve_model(self,x, a, b, c, d):
        '''Defines the model used during the curve fitting step'''
        return (a*x**3)+(b*x**2)+(c*x)+d
    def generate_offsets(self):
        '''Returns an array of n, 2d arrays that describes how each layer must change'''
        self.Linear_offsets = [[self.y_array[i]*np.sin(self.sweep_array[i]),self.y_array[i]*np.sin(self.dihedral_array[i])] for i in range(len(self.y_array))]
        self.Linear_offsets = [np.sum(self.Linear_offsets[0:i],axis=0) for i in range(1,len(self.Linear_offsets)+1)]
    def get_wing_section(self):
        '''Generates 2 profiles that can be interpolated between'''
        '''Hole Form'''
        self.Spars = self.generate_spars()
        '''Upper section'''


        '''Lower Section'''

class Foil(Tools):
    '''Working'''
    def __init__(self, points, chord, twist):
        Tools.__init__(self)
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
                              [0.0, 0.1, 0.20], #y values,All array must have a leading 0 to set the reference for the bed profile
                              [0.2, 0.2, 0.15],#Chord lengths in meters
                              [0,1,0.05], #Sweep angles
                              [0,0.02,0.02], #dihedral angles
                              [0.1,0.2,0.3], #Twist angles
                              [0.010,0.010,0.005], #Hole Size
                              [0.3,0.3,0.3],#Hole locations
                              [[1,1],[2,1],[1,1]]#Hole Elongation, vertical then horizontal streatching
                              )
        '''Test of generating wing'''
    #    ax = plt.figure().add_subplot(projection='3d')
    #    for i in self.test_wing.get_wing():
    #        ax.plot(i[:,0],i[:,1],i[:,2])
        #plt.show()
        #self.Wing = self.test_wing.get_wing()
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

    def gen_naca(self, foil_num): #genrates 4 digit naca air foils
        self.foil_num = str(foil_num)
        self.max_camber = int(self.foil_num[0])/100
        self.max_camber_pos = int(self.foil_num[1])/10
        self.thickness_ratio = int(self.foil_num[2:])/100
        self.foil_surface = []
        for x in range(100,0,-1):
            self.pp = x/100
            self.thickness = 5*self.thickness_ratio*(0.2969*sqrt(self.pp)-0.1260*self.pp-0.3516*self.pp**2+0.2843*self.pp**3-0.1015*self.pp**4)
            if self.pp<=self.max_camber_pos:
                if self.max_camber != 0:
                    self.camber_offset = (self.max_camber/self.max_camber_pos**2)*(2*self.max_camber_pos*self.pp-self.pp**2)
                    self.offset_theta = arctan((2*self.max_camber/self.max_camber_pos**2)*(self.max_camber_pos-self.pp))
                else:
                    self.camber_offset = 0
                    self.offset_theta = 0
            else:
                if self.max_camber!=0:
                    self.camber_offset = (self.max_camber/(1-self.max_camber_pos)**2)*((1-2*self.max_camber_pos)+2*self.max_camber_pos*self.pp-self.pp**2)
                    self.offset_theta = arctan((2*self.max_camber/(1-self.max_camber_pos)**2)*(self.max_camber_pos-self.pp))
                else:
                    self.camber_offset = 0
                    self.offset_theta = 0
            self.x_a = self.pp - self.thickness*sin(self.offset_theta)
            self.y_a = self.camber_offset+self.thickness*cos(self.offset_theta)
            self.foil_surface.append([self.x_a,self.y_a])
        for x in range(0,100,1):
            self.pp = x/100
            self.thickness = 5*self.thickness_ratio*(0.2969*sqrt(self.pp)-0.1260*self.pp-0.3516*self.pp**2+0.2843*self.pp**3-0.1015*self.pp**4)
            if self.pp<=self.max_camber_pos:
                if self.max_camber!=0:
                    self.camber_offset = (self.max_camber/self.max_camber_pos**2)*(2*self.max_camber_pos*self.pp-self.pp**2)
                    self.offset_theta = arctan((2*self.max_camber/self.max_camber_pos**2)*(self.max_camber_pos-self.pp))
                else:
                    self.camber_offset = 0
                    self.offset_theta = 0
            else:
                if self.max_camber!=0:
                    self.camber_offset = (self.max_camber/(1-self.max_camber_pos)**2)*((1-2*self.max_camber_pos)+2*self.max_camber_pos*self.pp-self.pp**2)
                    self.offset_theta = arctan((2*self.max_camber/(1-self.max_camber_pos)**2)*(self.max_camber_pos-self.pp))
                else:
                    self.camber_offset = 0
                    self.offset_theta = 0
            self.x_a = self.pp+self.thickness*sin(self.offset_theta)
            self.y_a = self.camber_offset-self.thickness*cos(self.offset_theta)
            self.foil_surface.append([self.x_a,self.y_a])
        self.foil_surface = array(self.foil_surface)
        return self.foil_surface

    def Load_airframe_defintion(self):
        """XFLR5 Stores its aircraft definitions a as XML file"""

        #Waiting on response from forum
        pass
if __name__ == "__main__":
    main()
