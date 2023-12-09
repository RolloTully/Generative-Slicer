import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib.path as mpltPath
import numpy as np
from matplotlib.lines import Line2D
import math
import time
grid_resolution = 20#mm at most one vertex per mm
def generate_random_points(N, x_range, y_range):
    # Generate random x and y coordinates within the given ranges
    x_points = np.random.uniform(x_range[0], x_range[1], N)
    y_points = np.random.uniform(y_range[0], y_range[1], N)

    # Combine x and y coordinates to form random points
    random_points = np.column_stack((x_points, y_points))

    return random_points

def generate_isometric_triangular_grid(rows, columns, size):
    points = []
    for r in range(rows):
        for c in range(columns):
            x = c * size * math.sqrt(3)
            y = r * size * 1.5
            if c % 2 != 0:
                y += size * 0.5
            points.append([x, y])
    return np.array(points)

def Mesh_layer(bounding_polygon, holes = None):
    '''Creates a convex Delaunay triangulation mesh'''
    '''
    Boundary Polygon
        a set of points that define the boundary

    holes
        a 3D array containing arrays of points that define the boundaries of the holes within the model
    '''
    #plt.scatter(bounding_polygon[:,0],bounding_polygon[:,1])
    #plt.show()
    lower_x_boundary = np.min(bounding_polygon[:,0])#Finds the lower x limit of the boundary
    upper_x_boundary = np.max(bounding_polygon[:,0])#Finds the upper x limit of the boundary
    lower_y_boundary = np.min(bounding_polygon[:,1])#Finds the lower y limit of the boundary
    upper_y_boundary = np.max(bounding_polygon[:,1])#Finds the upper y limit of the boundary
    x_range = upper_x_boundary-lower_x_boundary
    y_range = upper_y_boundary-lower_y_boundary
    x_grid = np.linspace(lower_x_boundary, upper_x_boundary, int(x_range/grid_resolution))
    y_grid = np.linspace(lower_y_boundary, upper_y_boundary, int(y_range/grid_resolution))
    xv, yv = np.meshgrid(x_grid,y_grid) #Defines a grid that encompases the entire object
    mesh_points = np.stack((xv.flatten(),yv.flatten()),axis = 1)#converts the grid to a set of points
    #mesh_points = np.concatenate((mesh_points, mesh_points+[(y_range/int(y_range/grid_resolution))/2,(x_range/int(x_range/grid_resolution))/2]),axis=0)#this overlays a second grid of points shifted by have a grid cell down and to the right, this forms a prepreating grid of triangles
    #mesh_points = generate_isometric_triangular_grid(100,100,5)
    #mesh_points = generate_random_points(1000, [lower_x_boundary,upper_x_boundary],[lower_y_boundary,upper_y_boundary])
    plt.scatter(mesh_points[:,0],mesh_points[:,1])
    plt.gca().set_aspect('equal')
    #plt.show()
    boundary_path = mpltPath.Path(bounding_polygon)#Defines a path that follows outer boundary of the model
    contained_points = np.array([point  for point in mesh_points if boundary_path.contains_point(point, radius=-3)])#points that are contained within the boundary path that are not within 3mm of the boundary path
    '''All remaining points now exist strictly within the boundary and not within 3mm of the 'wall' '''
    '''We must now remove points that exist within the holes and not within 3mm of the wall of the hole'''
    convex_points = contained_points #redefines array to make a loop possible
    if holes is not None:#Check that holes are actually given
        for hole in holes:#Go through the list of holes
            hole_boundary_path = mpltPath.Path(hole)#Define a path that is the boundary of the hole
            convex_points = np.array([point for point in convex_points if not hole_boundary_path.contains_point(point, radius=-3)])# We only want points that outside of the hole and not within 3mm of it
            convex_points = np.array([point for point in convex_points if not hole_boundary_path.contains_point(point, radius=3)])
    else:
        pass
    #plt.scatter(convex_points[:,0],convex_points[:,1])
    #plt.plot(holes[0,:,0],holes[0,:,1])
    #plt.show()
        #if no holes are given then you can just skip this step
    all_points = np.vstack((convex_points, bounding_polygon))#, holes[0]))#Concatinate all points for the triangularisation
    all_point_length = len(all_points)
    if holes is not None:
        for hole in holes:
            all_points = np.vstack((all_points,hole))
    triangulation = tri.Triangulation(all_points[:,0],all_points[:,1])#No idea how this works, its basicly black magic, but it does the triangulation bit
    #infill_point_index = #all thrse get altered
    #bounding_polygon_index = #allthese have forcres applied to them
    #holes_index = #all these indexes are grounded
    '''We know have a set of points and there connections'''
    '''We now need to filter the connections to make sure they dont cross over any holes'''
    '''this is when it becomes hell'''
    '''if a node that is contained within the holes wall has a connection that within the set of nodes within the wall of the hole but is not an immediate neighbour then this connection must be removed'''
    if holes is not None:# if holes exist
        '''this step gets us the starting index of each holes nodes'''
        hole_node_indexes = np.array([len(hole) for hole in holes])
        hole_node_indexes = np.cumsum(hole_node_indexes)+all_point_length
        hole_node_indexes = np.concatenate((np.array([all_point_length]),hole_node_indexes))
        connections = triangulation.edges
        filtered_connections = np.empty((0,2))
        for hole_node_index in range(0,len(holes)):
            index_start = hole_node_indexes[hole_node_index]
            index_stop = hole_node_indexes[hole_node_index+1]
            for connection in connections:
                if index_start<=  connection[0] <=index_stop and index_start <=  connection[1] <= index_stop:
                    '''We have confirmed that the connection is entirly contained within the hole wall'''
                    if (connection[0]+1 == connection[1] or connection[0]-1 == connection[1]) or ((connection[0] == index_start or connection[0]==index_stop) and (connection[1] == index_start or connection[1]== index_stop)):
                        '''we have confirmed adjacencey'''
                        '''therefore its a valid connection'''
                        filtered_connections = np.vstack((filtered_connections,connection))
                    else:
                        pass
                else:
                    filtered_connections = np.vstack((filtered_connections,connection))
    '''this next step closes all the holes, each hole has a missng connection to close it, this step inserts this connection'''
    for i in range(0,len(holes)):
        index_start = hole_node_indexes[hole_node_index]
        index_stop = hole_node_indexes[hole_node_index+1]
        filtered_connections = np.vstack((filtered_connections,[index_start+1,index_stop-1]))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for connection in filtered_connections:
        start_location  = all_points[int(connection[0])]
        end_location  = all_points[int(connection[1])]
        ax.plot(start_location[0],start_location[1])
        ax.plot(end_location[0],end_location[1])
        line = Line2D([start_location[0],end_location[0]],[start_location[1],end_location[1]])
        ax.add_line(line)
    plt.gca().set_aspect('equal')
    plt.show()

    #return convex_points, bounding_polygon, holes

    plt.triplot(triangulation, 'bo-', lw=0.5)
    plt.scatter(convex_points[:,0],convex_points[:,1],s = 2)
    plt.plot(bounding_polygon[:,0],bounding_polygon[:,1], 'r')
    plt.plot(holes[0,:,0],holes[0,:,1])
    plt.gca().set_aspect('equal')
    #plt.show()
start = time.process_time()
Mesh_layer(np.array([[200*np.cos(theta),200*np.sin(theta)] for theta in np.linspace(0, 2*np.pi,100)]),
           np.array([np.array([[ 50-50*np.cos(theta), 50-50*np.sin(theta)] for theta in np.linspace(0, 2*np.pi,100)]),
                     np.array([[ -50-50*np.cos(theta), -50-50*np.sin(theta)] for theta in np.linspace(0, 2*np.pi,100)])]))
print("Time to execute: ",time.process_time()-start)
