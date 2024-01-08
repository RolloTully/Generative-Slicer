import matplotlib.pyplot as plt
import numpy as np

#making an equaliater grid using 2 overlaying square grids
link_length = 5
grid_dimensions = np.array([link_length, -2*link_length*np.cos(60)])
grid_1 = np.array([[grid_dimensions[0]*x, grid_dimensions[1]*y] for x in range(10) for y in range(10)])
grid_2 = grid_1+(grid_dimensions/2)
grid = np.concatenate((grid_1,grid_2))
