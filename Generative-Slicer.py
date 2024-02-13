'''
what why how
is half
a section on prelim results

3 key results
'''
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib.path as mpltPath
import numpy as np
from matplotlib.lines import Line2D
from scipy.optimize import curve_fit
import aeropy.xfoil_module as xf
from scipy import interpolate
import math, time, random, numba
from tqdm import tqdm
from math import gcd, ceil
import itertools
from scipy import sparse
import numpy as np
import cvxpy as cvx
import matplotlib.pyplot as plt
from shapely.geometry import Point, LineString, Polygon
class Optimisation():
    def calcB(Nd, Cn, dof):
        self.m, self.n1, self.n2 = len(Cn), Cn[:,0].astype(int), Cn[:,1].astype(int)
        self.l, self.X, self.Y = Cn[:,2], Nd[self.n2,0]-Nd[self.n1,0], Nd[self.n2,1]-Nd[self.n1,1]
        self.d0, self.d1, self.d2, self.d3 = dof[self.n1*2], dof[self.n1*2+1], dof[self.n2*2], dof[self.n2*2+1]
        self.s = np.concatenate((-self.X/self.l * self.d0, -self.Y/self.l * self.d1, self.X/self.l * self.d2, self.Y/self.l * self.d3))
        self.r = np.concatenate((self.n1*2, self.n1*2+1, self.n2*2, self.n2*2+1))
        self.c = np.concatenate((np.arange(self.m), np.arange(self.m), np.arange(self.m), np.arange(self.m)))
        return sparse.coo_matrix((self.s, (self.r, self.c)), shape = (len(Nd)*2, self.m))
    #Solve linear programming problem
    def solveLP(Nd, Cn, f, dof, st, sc, jc):
        self.l = [col[2] + jc for col in Cn]
        self.B = calcB(Nd, Cn, dof)
        self.a = cvx.Variable(len(Cn))
        self.obj = cvx.Minimize(np.transpose(self.l) * self.a)
        self.q, self.eqn, self.cons= [], [], [self.a>=0]
        for self.k, self.fk in enumerate(f):
            self.q.append(cvx.Variable(len(Cn)))
            self.eqn.append(self.B * self.q[k] == self.fk * dof)
            self.cons.extend([self.eqn[k], self.q[k] >= -sc * self.a, self.q[k] <= st * self.a])
        self.prob = cvx.Problem(self.obj, self.cons)
        self.vol = prob.solve()
        self.q = [np.array(self.qi.value).flatten() for self.qi in self.q]
        self.a = np.array(self.a.value).flatten()
        self.u = [-np.array(eqnk.dual_value).flatten() for self.eqnk in self.eqn]
        return self.vol, self.a, self.q, self.u
    #Check dual violation
    def stopViolation(Nd, PML, dof, st, sc, u, jc):
        self.lst = np.where(PML[:,3]==False)[0]
        self.Cn = PML[self.lst]
        self.l = Cn[:,2] + jc
        self.B = calcB(Nd, self.Cn, dof).tocsc()
        self.y = np.zeros(len(Cn))
        for self.uk in u:
            self.yk = np.multiply(self.B.transpose().dot(self.uk) / self.l, np.array([[st], [-sc]]))
            self.y += np.amax(self.yk, axis=0)
        self.vioCn = np.where(y>1.0001)[0]
        self.vioSort = np.flipud(np.argsort(self.y[self.vioCn]))
        self.num = ceil(min(len(self.vioSort), 0.05*max( [len(Cn)*0.05, len(self.vioSort)])))
        for i in range(self.num):
            PML[self.lst[self.vioCn[self.vioSort[i]]]][3] = True
        return num == 0
    def trussopt(self, Links, Nd, Boundary, Convex, Ground, Forces, st, sc, jc):
        # [self.filtered_connections, self.all_points,
        #print(self.hole_node_indecies)
        #print(self.bounding_polygon_indecies)
        #print(self.convex_points_indecies)
        print(Boundary, Convex, Ground, Forces, st, sc, jc)
        self.dof, self.f, self.PML = np.ones((len(Nd),2)), np.zeros((len(Nd),2)), []
        '''Ground Constraints'''
        for hole in Ground:
            for index in range(hole[0], hole[1]):
                self.dof[index,:] = [0, 0]
        '''Loads constraints'''
        for index in range(Boundary[0], Boundary[1]):
            self.f[i,:] = Forces[index-Boundary[1]]
        '''Ground Structure'''

        for i, j in Links:#Parses over link start and end indecies
            self.dx, self.dy = abs(Nd[i][0] - Nd[j][0]), abs(Nd[i][1] - Nd[j][1])
            if gcd(int(self.dx), int(self.dy)) == 1 or jc != 0:
                    self.seg = [] if convex else LineString([Nd[i], Nd[j]])
                    if convex or poly.contains(self.seg) or poly.boundary.contains(self.seg):
                        self.PML.append( [i, j, np.sqrt(self.dx**2 + self.dy**2), False] )
        self.PML, self.dof = np.array(self.PML), np.array(self.dof).flatten()
        self.f = [self.f[i:i+len(Nd)*2] for i in range(0, len(self.f), len(Nd)*2)]
        print('Nodes: %d Members: %d' % (len(Nd), len(self.PML)))
        for pm in [p for p in self.PML if p[2] <= 1.42]:
            pm[3] = True
        #Start the ’member adding’ loop
        for itr in range(1, 100):
            self.Cn = self.PML[self.PML[:,3] == True]
            self.vol, self.a, self.q, self.u = solveLP(Nd, self.Cn, self.f, self.dof, st, sc, jc)
            print("Itr: %d, vol: %f, mems: %d" % (itr, self.vol, len(self.Cn)))
            plotTruss(Nd, Cn, a, q, max(a) * 1e-3, "Itr:" + str(itr))
            if stopViolation(Nd, PML, dof, st, sc, u, jc): break
        print("Volume: %f" % (vol))
        plotTruss(Nd, Cn, a, q, max(a) * 1e-3, "Finished", False)

class main():
    def __init__(self):
        self.grid_resolution = 5#mm
        self.Flow_Velocity = 10
        self.optimisation = Optimisation()
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

    #@numba.jit(nopython = False) #this is literaly my crack
    def Mesh_layer(self, bounding_polygon, holes = None, Verbose = True):
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
        self.mesh_points = self.generate_isometric_triangular_grid(self.upper_x_boundary,self.lower_x_boundary,self.upper_y_boundary, self.lower_y_boundary)
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
                    plt.scatter(self.point[0], self.point[1], 10, c = 'r')
                    #self.ax.annotate(index, (self.point[0], self.point[1]))
                for index in range(self.bounding_polygon_indecies[0], self.bounding_polygon_indecies[1]+1):
                    self.point = self.all_points[index]
                    plt.scatter(self.point[0], self.point[1], 10, c = 'g')
                    #self.ax.annotate(index, (self.point[0], self.point[1]))

                for hole_indicies in self.hole_node_indecies:
                    for index in range(hole_indicies[0], hole_indicies[1]+1):
                        self.point = self.all_points[index]
                        plt.scatter(self.point[0], self.point[1], 10, c = 'b')
                        #self.ax.annotate(index, (self.point[0], self.point[1]))
                plt.gca().set_aspect('equal')
                plt.show()
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
        self.invalid_connections = np.empty((0,2), np.int32)
        for index in range(self.bounding_polygon_indecies[0], self.bounding_polygon_indecies[1]):
            '''find all connections associated with this node'''
            self.associated_connections = np.array([connection for connection in self.connections if index in connection])#this is a horror show of inefficiency
            '''we now test each conection'''
            self.intra_hole_connections = np.empty((0,2), np.int32)
            for connection in self.associated_connections:
                if ((self.bounding_polygon_indecies[0]<=  connection[0] <= self.bounding_polygon_indecies[1]) and (self.bounding_polygon_indecies[0] <=  connection[1] <= self.bounding_polygon_indecies[1])):
                    self.intra_hole_connections = np.vstack((self.intra_hole_connections, connection))
            for connection in self.intra_hole_connections:
                if not ((connection[0]+1 == connection[1] or connection[0]-1 == connection[1]) or ((connection[0] == self.bounding_polygon_indecies[0] and connection[1]== self.bounding_polygon_indecies[1]) or (connection[1] == self.bounding_polygon_indecies[0] and connection[0] == self.bounding_polygon_indecies[1]))):
                    '''the nodes are not adjacent'''
                    self.invalid_connections = np.vstack((self.invalid_connections, connection))
        self.filtered_connections = np.empty((0,2))
        for connection in self.connections:
            '''We Now parse over each connection, if they are on the list we do not append them to the final list of connections'''
            if not np.any(np.all(self.invalid_connections == connection,axis=1)):
                self.filtered_connections = np.vstack((self.filtered_connections, connection))
        print(len(self.filtered_connections), " Connections")
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
            plt.gca().set_aspect('equal')
            plt.show()
        return [self.filtered_connections, self.all_points, self.bounding_polygon_indecies, self.convex_points_indecies, self.hole_node_indecies]

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
        self.Data = xf.find_pressure_coefficients(self.foil_name, 0., iteration=200, NACA=False)
        return self.Data

        #return convex_points, bounding_polygon, holes
    def Generate_loading_data(self, foil, layer_thickness):
        self.loading_data = self.get_pressure_distribution(foil)
        self.x = np.asarray(self.loading_data['x'],np.float16)
        self.y = np.asarray(self.loading_data['y'],np.float16)
        self.cp = np.asarray(self.loading_data['Cp'],np.float16)
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
        self.ax.plot(self.x,self.cp)
        plt.gca().set_aspect('equal')
        plt.show()
        print(self.cp)

        self.pts = np.asarray(np.column_stack((self.x,self.y)),dtype = np.float32)
        self.pt_dist = np.sqrt(np.sum(np.square(np.diff(self.pts, axis = 0)),axis=1))
        self.link_area = (self.pt_dist/1000)*layer_thickness
        self.dynamic_pressure = 0.5*1.225*self.Flow_Velocity**2
        self.link_Force = self.link_area*np.mean(np.asarray(self.cp, dtype = np.float32))*self.dynamic_pressure
        self.vector_link_Force = []
        plt.plot(self.link_Force)
        plt.show()
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
        for i, pair in enumerate(np.diff(self.pts,axis=0)):
            print(pair[1], pair[0])
            #input()
            self.theta = np.arctan2(pair[1], pair[0])+np.pi/2
            #print(self.theta)
            #if self.theta<0:
            #    print("Case1")
            #    self.theta += np.pi/2+((15*np.pi)/360)
            #elif self.theta>0:
            #    print("Case2")
            #    self.theta -= np.pi/2+((15*np.pi)/360)
            #print("Theta ", self.theta)
            #print([self.link_Force[i]*np.sin(self.theta),self.link_Force[i]*np.cos(self.theta)])
            self.vector_link_Force.append([self.cp[i]*np.cos(self.theta),self.cp[i]*np.sin(self.theta) ])
            self.ax.scatter(i,self.theta+np.pi/2)
            #self.vector_link_Force.append([5*np.sin(self.theta),5*np.cos(self.theta)])
        #plt.gca().set_aspect('equal')
        plt.grid()
        plt.show()
        self.vector_link_Force = np.array(self.vector_link_Force)
        self.node_loads = [self.vector_link_Force[0]]
        for i in range(1,len(self.link_Force)-1):
            self.node_loads.append((self.vector_link_Force[i])+(self.vector_link_Force[i+1]))
        self.node_loads.append(self.vector_link_Force[-1])
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
        self.ax.plot(self.pts[:,0],self.pts[:,1])
        print(self.node_loads)
        for i, point in enumerate(foil):
            try:
                print(i)
                print((self.node_loads[i][0], self.node_loads[i][1]))
                self.ax.arrow(self.x[i], self.y[i],self.node_loads[i][0],self.node_loads[i][1])
                #self.ax.annotate("", xy=(point[0], point[1]), xytext=(point[0]+self.node_loads[i][1]*1000, point[1]+self.node_loads[i][0]*1000),arrowprops=dict(arrowstyle="->"))
                self.ax.plot(self.x[i],self.y[i])
            except IndexError:
                pass
        plt.gca().set_aspect('equal')
        plt.show()
        print(self.node_loads)



    def mainloop(self):
        self.start_time = time.time()
        self.foil_number = '9918'
        self.foil = self.gen_naca(self.foil_number)
        self.third_chord = self.find_chamber_point(0.3, self.foil)
        self.three_quater_chord = self.find_chamber_point(0.6, self.foil)
        self.Layer = self.Mesh_layer(self.foil*300,np.array([
                                                            np.array([[ self.third_chord[0]*300+10*np.cos(theta), self.third_chord[1]*300+10*np.sin(theta)] for theta in np.linspace(0, 2*np.pi,20)])[:-1],
                                                            np.array([[ self.three_quater_chord[0]*300+5*np.cos(theta), self.three_quater_chord[1]*300+5*np.sin(theta)] for theta in np.linspace(0, 2*np.pi,20)])[:-1]]))
        self.forces = self.Generate_loading_data(self.foil, 0.2)
        #self.optimisation.trussopt(self.filtered_connections, self.all_points, self.bounding_polygon_indecies, self.convex_points_indecies, self.hole_node_indecies, self.forces, st = 1, sc =1, jc = 1)

        print(time.time()-self.start_time)
if __name__ == "__main__":
    main()
