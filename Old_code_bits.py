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
