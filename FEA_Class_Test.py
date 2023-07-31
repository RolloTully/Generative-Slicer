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
from pprint import pprint

class Material():
    def __init__(self, tensile_, shear_, tensile_strength_, shear_strength_):
        self.Tensile_Modulus = tensile_
        self.Shear_Modulus = shear_
        self.Tensile_Strength = tensile_strength_
        self.Shear_Strength = shear_strength_

class FEA_element(Material):
    def __init__(self, p1, p2, node1, node2, profile):
        Material.__init__(self, 50e6, 25.74e6,3.5e9,30e6)
        self.Node1 = node1
        self.Node2 = node2
        self.length = np.hypot(p2[0]-p1[0],p2[1]-p1[1])
        self.Area = np.product(np.diff(profile))
        self.Area_moment = profile[0]*(1/12)*profile[1]**3
        self.angle = np.arctan2(p2[0]-p1[0],p2[1]-p1[1])
        self.Transposition_matrix = [[np.cos(self.angle),np.sin(self.angle),0,0],
                                     [-np.sin(self.angle),np.cos(self.angle),0,0],
                                     [0,0,np.cos(self.angle),np.sin(self.angle)],
                                     [0,0,-np.sin(self.angle),np.cos(self.angle)]
                                     ]
        '''
        [[np.cos(self.angle), np.sin(self.angle), 0, 0                 ,  0                 , 0],
                                     [np.sin(self.angle), np.cos(self.angle), 0, 0                 ,  0                 , 0],
                                     [0                 , 0                 , 1, 0                 ,  0                 , 0],
                                     [0                 , 0                 , 0, np.cos(self.angle), -np.sin(self.angle), 0],
                                     [0                 , 0                 , 0, np.sin(self.angle),  np.cos(self.angle), 0],
                                     [0                 , 0                 , 0, 0                 ,  0                 , 1]]
        '''
        self.Local_stiffness_matrix = [[(self.Tensile_Modulus*self.Area)/self.length,0, -(self.Tensile_Modulus*self.Area)/self.length,0],
                                       [0,0,0,0],
                                       [-(self.Tensile_Modulus*self.Area)/self.length,0, (self.Tensile_Modulus*self.Area)/self.length,0],
                                       [0,0,0,0]]
        '''
        self.Local_stiffness_matrix = [[  (self.Tensile_Modulus*self.Area)/self.length, 0                                                        , 0                                                        ,-(self.Tensile_Modulus*self.Area)/self.length, 0                                                        , 0                                                        ],
                                       [  0                                           , (12*self.Tensile_Modulus*self.Area_moment)/self.length**3, (6 *self.Tensile_Modulus*self.Area_moment)/self.length**2, 0                                           ,-(12*self.Tensile_Modulus*self.Area_moment)/self.length**3, (6 *self.Tensile_Modulus*self.Area_moment)/self.length**2],
                                       [  0                                           , (6 *self.Tensile_Modulus*self.Area_moment)/self.length**2, (4 *self.Tensile_Modulus*self.Area_moment)/self.length   , 0                                           ,-(6 *self.Tensile_Modulus*self.Area_moment)/self.length**2, (2 *self.Tensile_Modulus*self.Area_moment)/self.length   ],
                                       [ -(self.Tensile_Modulus*self.Area)/self.length, 0                                                        , 0                                                        , (self.Tensile_Modulus*self.Area)/self.length, 0                                                        , 0                                                        ],
                                       [  0                                           ,-(12*self.Tensile_Modulus*self.Area_moment)/self.length**3,-(6 *self.Tensile_Modulus*self.Area_moment)/self.length**2, 0                                           , (12*self.Tensile_Modulus*self.Area_moment)/self.length**3,-(6 *self.Tensile_Modulus*self.Area_moment)/self.length**2],
                                       [  0                                           , (6 *self.Tensile_Modulus*self.Area_moment)/self.length**2, (2 *self.Tensile_Modulus*self.Area_moment)/self.length   , 0                                           ,-(6 *self.Tensile_Modulus*self.Area_moment)/self.length**2, (4 *self.Tensile_Modulus*self.Area_moment)/self.length   ]]
        '''
        self.Global_stiffness_matrix = np.matmul(np.transpose(self.Transposition_matrix),np.matmul(self.Local_stiffness_matrix,self.Transposition_matrix))

class FEA_Solver():
    def __init__(self,nodes_, elements_, forces_, D_C_):
        self.Element_array = []
        self.Nodes = np.array([[0,0],#0
                               [1,0],#1
                               [0,1],#2
                               [1,1]#3
                               ])# corners of a square, the position of each joint
        self.Elements = np.array([[0,1],
                                  [0,2],
                                  [2,3],
                                  [1,3]
                                  ]) #How the joints are linked
        self.Loads = np.array([[0,0],
                               [0,0],
                               [-20,20],
                               [20,20]
                                ])
        self.Dirichlet_Conditions = [0,1]#which nodes are fixed

        for self.pair in self.Elements:
            self.Element_array.append(FEA_element(self.Nodes[self.pair[0]],
                                                  self.Nodes[self.pair[1]],
                                                  self.pair[0],
                                                  self.pair[1],[0.005,0.002]))

        #print(np.array([self.Element.Local_stiffness_matrix for self.Element in self.Element_array]))
        print(np.array([self.Element.Global_stiffness_matrix for self.Element in self.Element_array]))
        #print(self.Convert_To_Solvable(self.Elements))
        self.Assemble_GSM()
        print(self.Solve())
        #pprint(self.GSM)
    def Solve(self):
        return self.Eliminate(np.linalg.inv(self.GSM))*self.Eliminate(elf.self.Convert_To_Solvable(self.Loads))
    def Assemble_GSM(self):
        self.GSM = np.zeros((len(self.Nodes)*4,len(self.Nodes)*4))
        for self.Element in self.Element_array:
            self.node1 = self.Element.Node1
            self.node2 = self.Element.Node2
            self.Stiffness = self.Element.Global_stiffness_matrix
            print(self.Stiffness)
            print(self.GSM[self.node1*4:(self.node1*4)+4,self.node2*4:(self.node2*4)+4])
            self.GSM[self.node1*4:(self.node1*4)+4,self.node2*4:(self.node2*4)+4] = self.Stiffness
        print(self.GSM)
        print(np.linalg.det(self.GSM.T))
    def Convert_To_Solvable(self,array):
        return array.flatten()[np.newaxis].T
    def Eliminate(self,array):
        for row in self.Dirichlet_Conditions:
            array = np.delete(array,row,0)
            array = np.delete(array,row,1)
        return array



FEA_Solver([[0,0],#0
            [1,0],#1
            [0,1],#2
            [1,1]#3
            ],
            [[0,1],
            [0,2],
            [2,3],
            [1,3]
            ],
            [[0,0],
            [0,0],
            [-20,20],
            [20,20]
            ],
            [0,1])
