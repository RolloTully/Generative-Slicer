import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib.path as mpltPath
import numpy as np
from matplotlib.lines import Line2D
from scipy.optimize import curve_fit
from scipy import interpolate
import math
import time
class main():
    def __init__(self):
        self.grid_resolution = 5#mm
        self.mainloop()
    def resample(self, surface_,  samples):
        '''resamples the foil at a higher density'''
        self.x_points = surface_[:,0] #Extracts an array of all x points
        self.y_points = surface_[:,1] #Extracts an array of all y points
        self.tck , self.u = interpolate.splprep([self.x_points,self.y_points],k=5,s=0.1, per=True) #From a parameterised interpolation of the points given
        self.i_x, self.i_y = interpolate.splev(np.linspace(0,1,samples),self.tck) # samples n points along the interpolation
        self.resampled_foil = np.stack((self.i_x,self.i_y),axis=1) #Stacks all the individual samples in to 1 2d array
        return self.resampled_foil
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
    def gen_naca(self, foil_num): #genrates 4 digit naca air foils finished and working, heritage code.
        self.foil_num = str(foil_num)
        self.max_camber = int(self.foil_num[0])/100
        self.max_camber_pos = int(self.foil_num[1])/10
        self.thickness_ratio = int(self.foil_num[2:])/100
        self.foil_surface = []
        for x in range(500,0,-1):
            self.pp = x/500
            self.thickness = 5*self.thickness_ratio*(0.2969*np.sqrt(self.pp)-0.1260*self.pp-0.3516*self.pp**2+0.2843*self.pp**3-0.1015*self.pp**4)
            if self.pp<=self.max_camber_pos:
                if self.max_camber != 0:
                    self.camber_offset = (self.max_camber/self.max_camber_pos**2)*(2*self.max_camber_pos*self.pp-self.pp**2)
                    self.offset_theta = np.arctan((2*self.max_camber/self.max_camber_pos**2)*(self.max_camber_pos-self.pp))
                else:
                    self.camber_offset = 0
                    self.offset_theta = 0
            else:
                if self.max_camber!=0:
                    self.camber_offset = (self.max_camber/(1-self.max_camber_pos)**2)*((1-2*self.max_camber_pos)+2*self.max_camber_pos*self.pp-self.pp**2)
                    self.offset_theta = np.arctan((2*self.max_camber/(1-self.max_camber_pos)**2)*(self.max_camber_pos-self.pp))
                else:
                    self.camber_offset = 0
                    self.offset_theta = 0
            self.x_a = self.pp - self.thickness*np.sin(self.offset_theta)
            self.y_a = self.camber_offset+self.thickness*np.cos(self.offset_theta)
            self.foil_surface.append([self.x_a,self.y_a])
        for x in range(0,500,1):
            self.pp = x/500
            self.thickness = 5*self.thickness_ratio*(0.2969*np.sqrt(self.pp)-0.1260*self.pp-0.3516*self.pp**2+0.2843*self.pp**3-0.1015*self.pp**4)
            if self.pp<=self.max_camber_pos:
                if self.max_camber!=0:
                    self.camber_offset = (self.max_camber/self.max_camber_pos**2)*(2*self.max_camber_pos*self.pp-self.pp**2)
                    self.offset_theta = np.arctan((2*self.max_camber/self.max_camber_pos**2)*(self.max_camber_pos-self.pp))
                else:
                    self.camber_offset = 0
                    self.offset_theta = 0
            else:
                if self.max_camber!=0:
                    self.camber_offset = (self.max_camber/(1-self.max_camber_pos)**2)*((1-2*self.max_camber_pos)+2*self.max_camber_pos*self.pp-self.pp**2)
                    self.offset_theta = np.arctan((2*self.max_camber/(1-self.max_camber_pos)**2)*(self.max_camber_pos-self.pp))
                else:
                    self.camber_offset = 0
                    self.offset_theta = 0
            self.x_a = self.pp+self.thickness*np.sin(self.offset_theta)
            self.y_a = self.camber_offset-self.thickness*np.cos(self.offset_theta)
            self.foil_surface.append([self.x_a,self.y_a])
        self.foil_surface = np.array(self.foil_surface)
        return self.foil_surface
    def generate_isometric_triangular_grid(self, x_max, x_min, y_max, y_min, major_row_step=5):#Fucked
        self.points = np.empty((0,2))
        self.elementry_grid = np.array([[0,0],
                                        [(major_row_step/2)*np.sin(np.radians(30)),(major_row_step/2)*np.cos(np.radians(30))]])
        plt.scatter(self.elementry_grid[:,0], self.elementry_grid[:,1])
        plt.show()
        for x_index in range(0,2):#int((x_max-x_min)/major_row_step)+1):
            for y_index in range(0,2):#int((y_max-y_min)/major_row_step)+1):
                self.new_grid = self.elementry_grid+(np.array([major_row_step, major_row_step])*np.array([x_index, y_index]))+np.array([x_min, y_min])
                self.points  = np.concatenate((self.points,self.new_grid))
                print(self.new_grid)
        plt.scatter(self.points[:,0], self.points[:,1])
        plt.grid()
        plt.show()
        return self.points

    def Mesh_layer(self, bounding_polygon, holes = None, Verbose = True):
        '''Creates a convex Delaunay triangulation mesh'''
        '''
        Boundary Polygon
            a set of points that define the boundary

        holes
            a 3D array containing arrays of points that define the boundaries of the holes within the model
        '''
        if Verbose:
            plt.plot(bounding_polygon[:,0],bounding_polygon[:,1])
            plt.show()
        self.lower_x_boundary = np.min(bounding_polygon[:,0])#Finds the lower x limit of the boundary
        self.upper_x_boundary = np.max(bounding_polygon[:,0])#Finds the upper x limit of the boundary
        self.lower_y_boundary = np.min(bounding_polygon[:,1])#Finds the lower y limit of the boundary
        self.upper_y_boundary = np.max(bounding_polygon[:,1])#Finds the upper y limit of the boundary
        self.x_range = self.upper_x_boundary-self.lower_x_boundary
        self.y_range = self.upper_y_boundary-self.lower_y_boundary
        self.x_grid = np.linspace(self.lower_x_boundary, self.upper_x_boundary, int(self.x_range/self.grid_resolution))
        self.y_grid = np.linspace(self.lower_y_boundary, self.upper_y_boundary, int(self.y_range/self.grid_resolution))
        self.xv, self.yv = np.meshgrid(self.x_grid,self.y_grid) #Defines a grid that encompases the entire object
        self.mesh_points = np.stack((self.xv.flatten(),self.yv.flatten()),axis = 1)#converts the grid to a set of points
        self.mesh_points = np.concatenate((self.mesh_points, self.mesh_points+[(self.y_range/int(self.y_range/self.grid_resolution))/2,(self.x_range/int(self.x_range/self.grid_resolution))/2]),axis=0)#this overlays a second grid of points shifted by have a grid cell down and to the right, this forms a prepreating grid of triangles
        #self.mesh_points = self.generate_isometric_triangular_grid(self.upper_x_boundary,self.lower_x_boundary,self.upper_y_boundary, self.lower_y_boundary)
        if Verbose:
            plt.scatter(self.mesh_points[:,0],self.mesh_points[:,1])
            plt.gca().set_aspect('equal')
            plt.show()
        self.boundary_path = mpltPath.Path(bounding_polygon, closed = True)#Defines a path that follows outer boundary of the model
        self.contained_points = np.array([self.point  for self.point in self.mesh_points if self.boundary_path.contains_point(self.point, radius=-3)])#points that are contained within the boundary path that are not within 3mm of the boundary path
        if Verbose:
            plt.scatter(self.contained_points[:,0],self.contained_points[:,1])
            plt.gca().set_aspect('equal')
            plt.show()
        '''All remaining points now exist strictly within the boundary and not within 3mm of the 'wall' '''
        '''We must now remove points that exist within the holes and not within 3mm of the wall of the hole'''
        self.convex_points = self.contained_points #redefines array to make a loop possible; probably not needed
        if holes is not None:#Check that holes are actually given
            for hole in holes:#Go through the list of holes
                self.hole_boundary_path = mpltPath.Path(hole)#Define a path that is the boundary of the hole
                self.convex_points = np.array([self.point for self.point in self.convex_points if not self.hole_boundary_path.contains_point(self.point, radius=-3)])# We only want points that outside of the hole and not within 3mm of it
                self.convex_points = np.array([self.point for self.point in self.convex_points if not self.hole_boundary_path.contains_point(self.point, radius=3)])
        else:
            self.convex_points = self.contained_points
            #if no holes are given then you can just skip this step
        self.all_points = np.vstack((self.convex_points, bounding_polygon))#, holes[0]))#Concatinate all points for the triangularisation
        if Verbose:
            plt.scatter(self.all_points[:,0],self.all_points[:,1])
            plt.gca().set_aspect('equal')
            plt.show()
        self.all_point_length = len(self.all_points)
        if holes is not None:
            for hole in holes:
                self.all_points = np.vstack((self.all_points,hole))
        self.triangulation = tri.Triangulation(self.all_points[:,0],self.all_points[:,1])#No idea how this works, its basicly black magic, but it does the triangulation bit
        '''We know have a set of points and there connections'''
        '''We now need to filter the connections to make sure they dont cross over any holes'''
        '''this is when it becomes hell'''
        '''if a node that is contained within the holes wall has a connection that within the set of nodes within the wall of the hole but is not an immediate neighbour then this connection must be removed'''
        if holes is not None:# if holes exist
            '''this step gets us the starting index of each holes nodes'''
            self.hole_node_indexes = np.array([len(hole) for hole in holes])
            self.hole_node_indexes = np.cumsum(self.hole_node_indexes)+self.all_point_length
            self.hole_node_indexes = np.concatenate((np.array([self.all_point_length]),self.hole_node_indexes))-1
            self.connections = self.triangulation.edges
            self.filtered_connections = np.empty((0,2))
            for hole_node_index in range(0,len(holes)):
                self.index_start = self.hole_node_indexes[hole_node_index]
                self.index_stop = self.hole_node_indexes[hole_node_index+1]
                for connection in self.connections:
                    if self.index_start<=  connection[0] <= self.index_stop and self.index_start <=  connection[1] <= self.index_stop:
                        '''We have confirmed that the connection is entirly contained within the hole wall'''
                        if (connection[0]+1 == connection[1] or connection[0]-1 == connection[1]) or ((connection[0] == self.index_start or connection[0]== self.index_stop) and (connection[1] == self.index_start or connection[1] == self.index_stop)):
                            '''we have confirmed adjacencey'''
                            '''therefore its a valid connection'''
                            self.filtered_connections = np.vstack((self.filtered_connections, connection))
                        else:
                            pass
                    else:
                        self.filtered_connections = np.vstack((self.filtered_connections, connection))
        else:
            self.filtered_connections = self.triangulation.edges
        '''this next step closes all the holes, each hole has a missng connection to close it, this step inserts this connection'''
        if holes is not None:
            for i in range(0,len(holes)):
                self.index_start = self.hole_node_indexes[i]
                self.index_stop = self.hole_node_indexes[i+1]
                self.filtered_connections = np.vstack((self.filtered_connections,[self.index_start+1,self.index_stop-1]))
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
        for connection in self.filtered_connections:
            self.start_location  = self.all_points[int(connection[0])]
            self.end_location  = self.all_points[int(connection[1])]
            self.ax.plot(self.start_location[0],self.start_location[1])
            self.ax.plot(self.end_location[0],self.end_location[1])
            self.line = Line2D([self.start_location[0],self.end_location[0]],[self.start_location[1],self.end_location[1]])
            self.ax.add_line(self.line)
        plt.gca().set_aspect('equal')
        plt.show()

        #return convex_points, bounding_polygon, holes
    def mainloop(self):
        self.naca2412 = self.gen_naca(2412)
        #self.naca2412 = self.resample(self.naca2412, 400)
        self.Mesh_layer(self.naca2412*300,np.array([np.array([[ 100+20*np.cos(theta), 5+10*np.sin(theta)] for theta in np.linspace(0, 2*np.pi,100)])]))

if __name__ == "__main__":
    main()
