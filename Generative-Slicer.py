'''
what why how
is half
a section on prelim results

3 key resultsf
'''
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib.path as mpltPath
import numpy as np
from matplotlib.lines import Line2D
from scipy.optimize import curve_fit
import aeropy.xfoil_module as xf
from scipy import interpolate
import math, time, random
from tqdm import tqdm
from math import gcd, ceil
import itertools
from scipy import sparse
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point, LineString, Polygon
from tqdm import tqdm
import cProfile
import cupy as cy


class main():
    def __init__(self):
        self.grid_resolution = 7#mm
        self.Flow_Velocity = 20
        self.Angle_of_Attack = 5
        self.foil_number = '2412'
        self.DOF = 6                  #2D Truss
        self.E = 4.107e9
        self.A = 4e-4
        self.G = 0.35e5
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
        for x in range(100,0,-1):
            self.pp = x/100
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
        for x in range(0,100,1):
            self.pp = x/100
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

    def generate_isometric_triangular_grid(self, x_max, x_min, y_max, y_min, major_row_step=5):#working and fast
        self.grid_dimensions = np.array([self.grid_resolution, -2*self.grid_resolution*np.cos(60)])
        self.x_rows = (x_max-x_min)/self.grid_dimensions[0]
        self.y_rows = (y_max-y_min)/self.grid_dimensions[1]
        self.grid_1 = np.array([[self.grid_dimensions[0]*x, self.grid_dimensions[1]*y] for x in range(int(self.x_rows)+1) for y in range(int(self.y_rows)+1)])
        return np.concatenate((self.grid_1,self.grid_1+(self.grid_dimensions/2)))+np.array([x_min,y_min])


    def Mesh_layer(self, bounding_polygon, holes = None, Verbose = False, angle = 0):
        '''Creates a convex Delaunay triangulation mesh'''
        '''
        Boundary Polygon
            a set of points that define the boundary

        holes
            a 3D array containing arrays of points that define the boundaries of the holes within the model

        Comments:
        really slow, can probably be spead up, getting better iteration time is down to 0.6 seconds, which is still slow as hell but its better
        it WORKS, IT WORKS!
        i would write more comments but i think this is actually just magic at this point
        '''
        if Verbose:
            plt.plot(bounding_polygon[:,0],bounding_polygon[:,1])
            plt.show()
        '''these 4 lines are so slow, there is some fortran code you can compile to speed it up by reducing the number of function calls'''
        self.lower_x_boundary = np.min(bounding_polygon[:,0])#Finds the lower x limit of the boundary
        self.upper_x_boundary = np.max(bounding_polygon[:,0])#Finds the upper x limit of the boundary
        self.lower_y_boundary = np.min(bounding_polygon[:,1])#Finds the lower y limit of the boundary
        self.upper_y_boundary = np.max(bounding_polygon[:,1])#Finds the upper y limit of the boundary
        self.mid_x = (self.upper_x_boundary-self.lower_x_boundary)/2
        self.mid_y = (self.upper_y_boundary-self.lower_y_boundary)/2
        self.x_range = self.upper_x_boundary-self.lower_x_boundary
        self.y_range = self.upper_y_boundary-self.lower_y_boundary
        self.range = max([self.x_range, self.y_range])
        self.mid_point = np.array([self.mid_x, self.mid_y])
        self.mesh_points = self.generate_isometric_triangular_grid(self.mid_x+self.range,self.mid_x-self.range,self.mid_y+self.range,self.mid_y-self.range)


        #self.mesh_points = np.array([self.mid_point+np.dot((point - self.mid_point),np.array([[np.cos(angle), np.sin(angle)],[-np.sin(angle),np.cos(angle)]])) for point in self.mesh_points ])

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
                self.convex_points = np.array([self.point for self.point in self.convex_points if not self.hole_boundary_path.contains_point(self.point, radius=-self.grid_resolution)])# We only want points that outside of the hole and not within 3mm of it
                self.convex_points = np.array([self.point for self.point in self.convex_points if not self.hole_boundary_path.contains_point(self.point, radius=self.grid_resolution)])
        else:
            self.convex_points = self.contained_points
            #if no holes are given then you can just skip this step
        self.all_points = np.vstack((self.convex_points, bounding_polygon))#, holes[0]))#Concatinate all points for the triangularisation
        self.convex_points_indecies = np.array([0,len(self.convex_points)-1])
        self.bounding_polygon_indecies = np.array([self.convex_points_indecies[1]+1, self.convex_points_indecies[1]+len(bounding_polygon)])
        print("Convex points indecies", self.convex_points_indecies)
        print("bounding point indicices", self.bounding_polygon_indecies)
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
        '''if a node that is contained within the holes wall has a connection thats within the set of nodes within the wall of the hole but is not an immediate neighbour then this connection must be removed'''

        if holes is not None:# if holes exist
            '''this step gets us the starting index of each holes nodes'''
            self.hole_node_indecies = np.array([[0,self.all_point_length-1]])
            for hole in holes:
                self.hole_end_index = len(hole)
                self.hole_node_indecies = np.concatenate((self.hole_node_indecies,np.array([[self.hole_node_indecies[-1,1]+1,self.hole_node_indecies[-1,1]+self.hole_end_index]])))
            self.hole_node_indecies = self.hole_node_indecies[1:]
            if Verbose:
                self.fig = plt.figure()
                self.ax = self.fig.add_subplot(111)
                for index in range(0, self.convex_points_indecies[1]+1):
                    self.point = self.all_points[index]
                    #plt.scatter(self.point[0], self.point[1], 10, c = 'r')
                    #self.ax.annotate(index, (self.point[0], self.point[1]))
                for index in range(self.bounding_polygon_indecies[0], self.bounding_polygon_indecies[1]+1):
                    self.point = self.all_points[index]
                    #plt.scatter(self.point[0], self.point[1], 10, c = 'g')
                    #self.ax.annotate(index, (self.point[0], self.point[1]))

                for hole_indicies in self.hole_node_indecies:
                    for index in range(hole_indicies[0], hole_indicies[1]+1):
                        self.point = self.all_points[index]
                        #plt.scatter(self.point[0], self.point[1], 10, c = 'b')
                        #self.ax.annotate(index, (self.point[0], self.point[1]))
                #plt.gca().set_aspect('equal')
                #plt.show()
            self.connections = self.triangulation.edges.copy()
            self.invalid_connections = np.empty((0,2), np.int32)
            for hole_node_index in self.hole_node_indecies:
                self.index_start = hole_node_index[0]
                self.index_stop = hole_node_index[1]
                for index in range(hole_node_index[0], hole_node_index[1]+1):
                    '''find all connections associated with this node'''
                    self.associated_connections = np.array([connection for connection in self.connections if index in connection])#this is a horror show of inefficiency
                    '''we now test each conection'''
                    self.intra_hole_connections = np.empty((0,2), np.int32)
                    for connection in self.associated_connections:
                        if ((self.index_start<=  connection[0] <= self.index_stop) and (self.index_start <=  connection[1] <= self.index_stop)):
                            self.intra_hole_connections = np.vstack((self.intra_hole_connections, connection))

                    for connection in self.intra_hole_connections:
                        if not ((connection[0]+1 == connection[1] or connection[0]-1 == connection[1]) or ((connection[0] == self.index_start and connection[1]== self.index_stop) or (connection[1] == self.index_start and connection[0] == self.index_stop))):
                            '''the nodes are not adjacent'''
                            self.invalid_connections = np.vstack((self.invalid_connections, connection))
            self.filtered_connections = np.empty((0,2))
            for connection in self.triangulation.edges:
                if not np.any(np.all(self.invalid_connections == connection,axis=1)):
                    self.filtered_connections = np.vstack((self.filtered_connections, connection))
        else:
            self.filtered_connections = self.triangulation.edges
        self.connections = self.filtered_connections# Just tidies it up a bit.

        '''We must now do the same for the boundary of the airfoil becuase in the case that the airfoil is concave it will form external connections'''

        #self.invalid_connections = np.empty((0,2), np.int32)
        #for index in range(self.bounding_polygon_indecies[0], self.bounding_polygon_indecies[1]):
        #    '''find all connections associated with this node'''
        #    self.associated_connections = np.array([connection for connection in self.connections if index in connection])#this is a horror show of inefficiency
        #    '''we now test each conection'''
        #    self.intra_hole_connections = np.empty((0,2), np.int32)
        #    for connection in self.associated_connections:
        #        if ((self.bounding_polygon_indecies[0]<=  connection[0] <= self.bounding_polygon_indecies[1]) and (self.bounding_polygon_indecies[0] <=  connection[1] <= self.bounding_polygon_indecies[1])):
        #            self.intra_hole_connections = np.vstack((self.intra_hole_connections, connection))
        #    for connection in self.intra_hole_connections:
        #        if not ((connection[0]+1 == connection[1] or connection[0]-1 == connection[1]) or ((connection[0] == self.bounding_polygon_indecies[0] and connection[1]== self.bounding_polygon_indecies[1]) or (connection[1] == self.bounding_polygon_indecies[0] and connection[0] == self.bounding_polygon_indecies[1]))):
        #            '''the nodes are not adjacent'''
        #            self.invalid_connections = np.vstack((self.invalid_connections, connection))
        #self.filtered_connections = np.empty((0,2))
        #for connection in self.connections:
        #    '''We Now parse over each connection, if they are on the list we do not append them to the final list of connections'''
        #    if not np.any(np.all(self.invalid_connections == connection,axis=1)):
        #        self.filtered_connections = np.vstack((self.filtered_connections, connection))
        #print(len(self.filtered_connections), " Connections")

        if Verbose:
            self.fig = plt.figure()
            self.ax = self.fig.add_subplot(111)
            for connection in self.filtered_connections:
                self.start_location  = self.all_points[int(connection[0])]
                self.end_location  = self.all_points[int(connection[1])]
                self.ax.plot(self.start_location[0],self.start_location[1])
                self.ax.plot(self.end_location[0],self.end_location[1])
                self.line = Line2D([self.start_location[0],self.end_location[0]],[self.start_location[1],self.end_location[1]])
                self.ax.add_line(self.line)
            plt.title("NACA 2412 with "+str(self.grid_resolution)+"mm Lattice infill")
            plt.gca().set_aspect('equal')
            plt.show()
        self.filtered_connections = self.filtered_connections.astype(int)
        #return [self.filtered_connections, self.all_points, self.bounding_polygon_indecies, self.convex_points_indecies, self.hole_node_indecies]

    def get_pressure_distribution(self,foil):
        '''
        oh my this is so stupid, that you have to do this is maddening and took me
        so long to figure out i just when and wrote this function
        and completly went around how your meant to do it
        if i could re do my dissertation i would just fix this disaster by rewriting xfoil in python

        dont touch anything or it will 100% break
        '''
        self.foil_name = "do_not_delete_this_or_everything_breaks"
        self.output_file = open(self.foil_name,'w')
        for line in foil:
            self.output_file.write("     "+str(format(line[0],'.6f'))+"    "+str(format(line[1],'.6f'))+"\n")#
        self.output_file.close()
        self.Data = xf.find_pressure_coefficients(self.foil_name, self.Angle_of_Attack,Reynolds = 100e6, iteration=200, NACA=False)
        return self.Data

    def Generate_loading_data(self, foil, layer_thickness):
        self.loading_data = self.get_pressure_distribution(foil)
        self.x = np.asarray(self.loading_data['x'],np.float64)
        self.y = np.asarray(self.loading_data['y'],np.float64)
        self.cp = np.asarray(self.loading_data['Cp'],np.float64)
        self.pts = np.asarray(np.column_stack((self.x,self.y)),dtype = np.float64)
        self.pt_dist = np.sqrt(np.sum(np.square(np.diff(self.pts, axis = 0)),axis=1))
        self.link_area = (self.pt_dist/1000)*layer_thickness
        self.dynamic_pressure = 0.5*1.225*self.Flow_Velocity**2
        self.link_Force = self.link_area*np.mean(np.asarray(self.cp, dtype = np.float64))*self.dynamic_pressure
        self.vector_link_Force = []
        for i, pair in enumerate(np.diff(self.pts,axis=0)):
            #print(pair[1], pair[0])
            self.theta = np.arctan2(pair[1], pair[0])+np.pi/2
            self.vector_link_Force.append([self.cp[i]*np.cos(self.theta),self.cp[i]*np.sin(self.theta) ])#This is activly wrong
        self.vector_link_Force = np.array(self.vector_link_Force)
        self.node_loads = [self.vector_link_Force[0]]
        for i in range(1,len(self.link_Force)-1):
            self.node_loads.append((self.vector_link_Force[i])+(self.vector_link_Force[i+1]))
        self.node_loads.append(self.vector_link_Force[-1])
        return self.node_loads

    def Truss_Analysis(self, verbose):
        '''This thing works, dont touch it, its basicly black magic'''
        self.NE = len(self.filtered_connections)           #Number of bars
        self.d = self.all_points[self.filtered_connections[:,1],:] - self.all_points[self.filtered_connections[:,0],:]
        self.length = np.sqrt((np.square(self.d)).sum(axis=1))
        self.theta = self.d.T/self.length
        self.a = np.concatenate((-self.theta.T,self.theta.T), axis=1)
        print(self.a)

        self.Global_Stiffness = np.zeros([self.NDOF,self.NDOF])
        '''Now parsing over each element to add the mto the global stiffness matrix'''
        for index in range(self.NE):
            self.aux  = 3*self.filtered_connections[index,:]
            self.indecies = np.r_[self.aux[0]:self.aux[0]+3,self.aux[1]:self.aux[1]+3]
            self.k = ((self.E*self.A)/self.length[index])
            self.transformation_matrix = np.zeros((6,6))
            print(self.a[index].reshape((2,2)))
            print(self.a[index].reshape((2,2)).flatten())

            self.transformation_matrix[0:2,0:2] = self.a[index].reshape((2,2))
            self.transformation_matrix[3:5,3:5] = self.a[index].reshape((2,2))
            self.transformation_matrix[2,2] = 1
            self.transformation_matrix[5,5] = 1
            print(self.transformation_matrix)
            print(self.transformation_matrix.flatten())


            self.local_stiffness = np.array([[self.k,0,0,-self.k,0,0],
                                             [0,12*self.k/self.length[index]**2, 6*self.k/self.length[index],0,-12*self.k/self.length[index]**2, 6*self.k/self.length[index]],
                                             [0,6*self.k/self.length[index], 4*self.k,0,-6*self.k/self.length[index], 2*self.k],
                                             [-self.k,0,0,self.k,0,0],
                                             [0,-12*self.k/self.length[index]**2,-6*self.k/self.length[index],0,12*self.k/self.length[index],-6*self.k/self.length[index]],
                                             [0,6*self.k/self.length[index], 2*self.k,0,-6*self.k/self.length[index], 4*self.k]
                                             ])
            #with np.printoptions(precision=4, suppress=True, formatter={'float': '{:0.4f}'.format}, linewidth=100):
                #print(self.local_stiffness)
            self.ES = np.dot(np.dot(self.transformation_matrix.flatten().T,self.local_stiffness.flatten()),self.transformation_matrix[np.newaxis])
            self.Global_Stiffness[np.ix_(self.indecies,self.indecies)] = self.Global_Stiffness[np.ix_(self.indecies,self.indecies)] + self.ES
        self.supportDOF = (self.DOFCON.flatten() == 0)
        print(self.supportDOF)
        self.Kff = self.Global_Stiffness[np.ix_(self.freeDOF,self.freeDOF)]
        self.Uf = cy.linalg.solve(cy.asarray(self.Kff),cy.asarray(self.Pf)).get()
        self.U = self.DOFCON.astype(float).flatten()
        self.U[self.freeDOF] = self.Uf
        self.U[self.supportDOF] = self.Ur
        self.U = self.U.reshape(self.NN,self.DOF)
        self.Displacement_points = self.all_points+self.U
        self.new_d = self.Displacement_points[self.filtered_connections[:,1],:] - self.all_points[self.filtered_connections[:,0],:]
        self.new_length = np.sqrt((self.new_d**2).sum(axis=1))
        self.d_length = self.length-self.new_length

    def Optimise(self, verbose = True):
        '''Not all members can be removed, we must find the array of memebers that we can remove with out the k-matrix becoming singular'''
        '''Weve now found this array, we can now go throught this list and find which of there contained members carries the least stress'''
        self.current_connections = self.filtered_connections
        self.halt = False
        self.Structure_Mass = []
        self.Structure_total_displacment = []
        self.Peak_displacment = []
        self.memberlengths = self.all_points[self.filtered_connections[:,1],:] - self.all_points[self.filtered_connections[:,0],:]
        self.original_total_length = np.sum(np.sqrt((self.memberlengths**2).sum(axis=1)))
        self.frame_number = 0
        self.working_lattice = self.filtered_connections.copy()
        while not self.halt:
            self.frame_number +=1
            '''we must check that there arent any 'deadlegs' in the lattice, these cause large displacments that trips the stopping condition without contributing to its strength'''
            self.Truss_Analysis(False)


            #'''Making a video of the optimisation'''
            #self.Strain = self.displaced_s_length-self.original_s_length
            #self.min_strain = np.min(self.Strain)
            #self.max_strain = np.max(self.Strain)
            #self.strain_range = self.max_strain-self.min_strain

            #self.fig = plt.figure()
            #self.ax = self.fig.add_subplot(111)
            #self.Displacement_points = self.all_points+self.U
            #for i, connection in enumerate(self.working_lattice):
            #    self.index_strain = self.Strain[i]
            #    self.colour_encoding = (self.index_strain-self.min_strain)/self.strain_range
            #    self.colour = self.cmap(self.colour_encoding)
            #    self.start_location  = self.Displacement_points[int(connection[0])]
            #    self.end_location  = self.Displacement_points[int(connection[1])]
            #    self.ax.plot(self.start_location[0],self.start_location[1])
            #    self.ax.plot(self.end_location[0],self.end_location[1])
            #    self.line = Line2D([self.start_location[0],self.end_location[0]],[self.start_location[1],self.end_location[1]],c = self.colour)
            #    self.ax.add_line(self.line)
            #plt.title("Deformed Airfoil shape(x1000), NACA 2412, 5 Degrees, 20 ms^-1")
            #plt.gca().set_aspect('equal')
            #plt.savefig('C:\\Users\\rollo\\Documents\\GitHub\\Generative-Slicer\\Optimisation video\\'+str(self.frame_number)+'.png', dpi = 200)

            self.memberlengths = self.all_points[self.working_lattice[:,1],:] - self.all_points[self.working_lattice[:,0],:]
            self.total_length = np.sum(np.sqrt((self.memberlengths**2).sum(axis=1)))
            self.Structure_Mass.append(self.total_length/self.original_total_length)
            self.valid_changes = []
            for connection in self.working_lattice:
                self.valid = True
                if not (connection[0] or connection[1]) in range(self.bounding_polygon_indecies[0],self.bounding_polygon_indecies[1]+1):
                    for hole in self.hole_node_indecies:
                        if not (connection[0] or connection[1]) in range(hole[0],hole[1]+1):
                            pass
                        else:
                            self.valid = False
                else:
                    self.valid = False
                if self.valid:
                    self.valid_changes.append(connection)

            if self.valid_changes == []:
                self.output_lattice = self.working_lattice
                break

            self.costs = []
            self.indexes = []
            for valid_change in self.valid_changes:
                self.indexes.append(np.where(np.all(self.current_connections == valid_change, axis=1))[0])
            print("We here")
            for self.index in tqdm(self.indexes):
                try:
                    self.filtered_connections = np.delete(self.working_lattice, self.index,0)
                    self.Truss_Analysis(False)
                    self.costs.append(np.max(np.abs(self.d_length)))
                except:
                    self.costs.append(1e9)
                    pass
                #print("Removing member ", self.index," this is ",self.current_connections[self.index]," this has a cost of ",np.sum(self.d_length))
            self.index_min_cost = np.argmin(self.costs)
            print(self.costs)
            self.cost = self.costs[self.index_min_cost]
            print("minimum cost index: ", self.index_min_cost, " Cost: ", self.cost)
            if np.max(np.abs(self.d_length)) >5:
                self.output_lattice = self.working_lattice
                self.halt = True
            else:
                self.working_lattice = np.delete(self.working_lattice, self.indexes[self.index_min_cost],0)
                self.valid_changes = np.array(self.valid_changes)
                self.Structure_total_displacment.append(np.sum(np.abs(self.d_length)))
                self.Peak_displacment.append(np.max(np.abs(self.d_length)))
            '''Floppy Member removal'''
            print("and We here")
            #self.points, self.counts = np.unique(self.current_connections.flatten(), return_counts=True)
            #print(np.count_nonzero(self.counts==1))
            #while 1 in self.counts:
            #    self.points, self.counts = np.unique(self.current_connections.flatten(), return_counts=True)
            #    print(self.counts)
            #    self.extranious_memeber_end = self.points[np.argmin(self.counts)]
            #    print("Member end ", self.extranious_memeber_end)
            #    '''We must now remove any memebrs containing this point'''
            #    self.members_to_remove = np.where(self.current_connections == self.extranious_memeber_end)
            #    print("member index ", self.members_to_remove)
            #    print("Member ", self.current_connections[self.members_to_remove[0]])
            #    self.current_connections = np.delete(self.current_connections, self.members_to_remove[0],0)



            '''Live Plotting'''
            if True:
                '''Process Logging'''
                print("____________________________________________________________________________________")
                print(len(self.valid_changes), "Valid changes can be made.")
                print("Removing the ", self.index_min_cost, "Member results in the least deformation.")
                print("This is connection ",self.indexes[self.index_min_cost])
                print("The New minimum cost is ",self.cost)
                print("The Current maximum deformation is ", np.max(np.abs(self.d_length)))
                self.fig = plt.figure()
                self.ax = self.fig.add_subplot(111)
                for connection in self.working_lattice:
                    self.start_location  = self.all_points[int(connection[0])]
                    self.end_location  = self.all_points[int(connection[1])]
                    self.ax.plot(self.start_location[0],self.start_location[1])
                    self.ax.plot(self.end_location[0],self.end_location[1])
                    self.line = Line2D([self.start_location[0],self.end_location[0]],[self.start_location[1],self.end_location[1]], c='blue')
                    self.ax.add_line(self.line)
                for connection in self.valid_changes:
                    self.start_location  = self.all_points[int(connection[0])]
                    self.end_location  = self.all_points[int(connection[1])]
                    self.ax.plot(self.start_location[0],self.start_location[1])
                    self.ax.plot(self.end_location[0],self.end_location[1])
                    self.line = Line2D([self.start_location[0],self.end_location[0]],[self.start_location[1],self.end_location[1]],c = 'red')
                    self.ax.add_line(self.line)
                self.start_location = self.all_points[int(self.valid_changes[self.index_min_cost,0])]
                self.end_location = self.all_points[int(self.valid_changes[self.index_min_cost,1])]
                self.ax.plot(self.start_location[0],self.start_location[1])
                self.ax.plot(self.end_location[0],self.end_location[1])
                self.line = Line2D([self.start_location[0],self.end_location[0]],[self.start_location[1],self.end_location[1]], c='pink')
                self.ax.add_line(self.line)
                plt.title("Grid of valid changes")
                plt.gca().set_aspect('equal')
                plt.show(block=False)
                plt.pause(0.5)
                plt.close()


        '''Plotting Stuff'''
        self.fig,self.ax = plt.subplots()
        plt.grid()
        self.color = 'tab:red'
        self.ax.set_ylim(0,1.1)
        self.ax.plot(np.array(self.Structure_Mass), color = self.color, linewidth = 3.0)
        self.ax.set_xlabel("Iterations", fontsize=20)
        self.ax.set_ylabel("Volume Fraction (%)", color = self.color, fontsize=20)
        self.ax.tick_params(axis='y',labelcolor = self.color)

        self.ax2 = self.ax.twinx()  # instantiate a second axes that shares the same x-axis
        self.color = 'tab:blue'
        self.ax2.set_ylabel('Total Structual Displacment(mm)', color = self.color, fontsize=20)  # we already handled the x-label with ax1
        self.ax2.plot(self.Structure_total_displacment, color = self.color, linewidth = 3.0)
        self.ax2.yaxis.tick_right()
        self.ax2.yaxis.set_label_position("right")
        self.ax2.tick_params(axis='y',labelcolor = self.color)

        self.ax3 = self.ax.twinx()
        self.ax3.spines.right.set_position(("outward",60))
        self.make_patch_spines_invisible(self.ax3)
        self.ax3.spines["right"].set_visible(True)
        self.ax3.semilogy(self.Peak_displacment, color = "Green", linewidth = 3.0)
        self.ax3.set_ylabel('Maximum Structual Displacment(mm)', color = 'Green', fontsize=20)
        self.ax3.tick_params(axis='y',labelcolor = 'Green')
        plt.title("Optmisation of NACA 2412 for 5 degree AoA at Re = 1e5")
        plt.show()


        self.Truss_Analysis(False)
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
        self.Displacement_points = self.all_points+self.U*100
        self.original_memberlengths = self.all_points[self.filtered_connections[:,1],:] - self.all_points[self.filtered_connections[:,0],:]
        self.original_s_length = np.sum(np.sqrt((self.original_memberlengths**2).sum(axis=1)))
        self.displaced_memberlengths = self.Displacement_points[self.filtered_connections[:,1],:] - self.all_points[self.filtered_connections[:,0],:]
        self.displaced_s_length = np.sqrt((self.displaced_memberlengths**2).sum(axis=1))
        self.Strain = self.displaced_s_length-self.original_s_length
        self.min_strain = np.min(self.Strain)
        self.max_strain = np.max(self.Strain)
        self.strain_range = self.max_strain-self.min_strain
        self.cmap = plt.cm.get_cmap('turbo')
        for i, point in enumerate(self.Displacement_points):
            print(i, point)
            self.ax.annotate(str(i), (int(point[0]), int(point[1])))
        #plt.scatter(self.point[0], self.point[1], 10, c = 'r')
        #self.ax.annotate(index, (self.point[0], self.point[1]))
        for i, connection in enumerate(self.filtered_connections):
            self.index_strain = self.Strain[i]
            self.colour_encoding = (self.index_strain-self.min_strain)/self.strain_range
            self.colour = self.cmap(self.colour_encoding)
            self.start_location  = self.Displacement_points[int(connection[0])]
            self.end_location  = self.Displacement_points[int(connection[1])]
            self.ax.plot(self.start_location[0],self.start_location[1])
            self.ax.plot(self.end_location[0],self.end_location[1])
            self.line = Line2D([self.start_location[0],self.end_location[0]],[self.start_location[1],self.end_location[1]],c = self.colour)
            self.ax.add_line(self.line)
        #for connection in self.filtered_connections:
        #    self.start_location  = self.all_points[int(connection[0])]
        #    self.end_location  = self.all_points[int(connection[1])]
        #    self.ax.plot(self.start_location[0],self.start_location[1])
        #    self.ax.plot(self.end_location[0],self.end_location[1])
        #    self.line = Line2D([self.start_location[0],self.end_location[0]],[self.start_location[1],self.end_location[1]], c='blue')
        #    self.ax.add_line(self.line)
        plt.title("Deformed Airfoil shape(x100), NACA 2412, 5 Degrees, 20 ms^-1")
        plt.gca().set_aspect('equal')
        plt.show()

    def make_patch_spines_invisible(self, ax):
        ax.set_frame_on(True)
        ax.patch.set_visible(False)
        for sp in ax.spines.values():
            sp.set_visible(False)

    def mainloop(self):
        self.foil = self.gen_naca(self.foil_number)
        self.third_chord = self.find_chamber_point(0.3, self.foil)
        self.three_quater_chord = self.find_chamber_point(0.6, self.foil)
        self.Layer = self.Mesh_layer(self.foil*300,np.array([
                                                            np.array([[ self.third_chord[0]*300+10*np.cos(theta), self.third_chord[1]*300+10*np.sin(theta)] for theta in np.linspace(0, 2*np.pi,20)])[:-1],
                                                            np.array([[ self.three_quater_chord[0]*300+5*np.cos(theta), self.three_quater_chord[1]*300+5*np.sin(theta)] for theta in np.linspace(0, 2*np.pi,20)])[:-1]]))
        self.all_points = self.all_points.astype(np.float64)
        self.forces = self.Generate_loading_data(self.foil, 0.2)
        '''Pre-processing as much as possible'''
        self.DOFCON = np.ones_like(self.all_points*3).astype(int)
        self.Ur = []
        self.Forces = np.zeros_like(self.all_points)
        for index in range(self.bounding_polygon_indecies[0], self.bounding_polygon_indecies[1]+1):
            self.Forces[index,:] = self.forces[index-self.bounding_polygon_indecies[1]]
        for indexs in self.hole_node_indecies:
            for i in range(indexs[0], indexs[1]+1):
                self.DOFCON[i,:] = 0
                self.Ur.append(0)
                self.Ur.append(0)
                self.Ur.append(0)
                self.Ur.append(0)
                self.Ur.append(0)
                self.Ur.append(0)

        self.freeDOF = self.DOFCON.flatten().nonzero()[0]
        self.Pf = self.Forces.flatten()[self.freeDOF]

        self.NN = len(self.all_points)          #Number of nodes
        self.NDOF = self.DOF*self.NN            #Total number of degree of freedom
        self.Truss_Analysis(True)
        self.Displacement_points = self.all_points+self.U*1000
        self.original_memberlengths = self.all_points[self.filtered_connections[:,1],:] - self.all_points[self.filtered_connections[:,0],:]
        self.original_s_length = np.sum(np.sqrt((self.original_memberlengths**2).sum(axis=1)))
        self.displaced_memberlengths = self.Displacement_points[self.filtered_connections[:,1],:] - self.all_points[self.filtered_connections[:,0],:]
        self.displaced_s_length = np.sqrt((self.displaced_memberlengths**2).sum(axis=1))
        self.Strain = self.displaced_s_length-self.original_s_length
        self.min_strain = np.min(self.Strain)
        self.max_strain = np.max(self.Strain)
        self.strain_range = self.max_strain-self.min_strain
        self.cmap = plt.cm.get_cmap('turbo')

        self.Optimise(False)
if __name__ == "__main__":
    #cProfile.run('main()')
    main()
