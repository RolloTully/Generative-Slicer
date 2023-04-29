"""
This is a shity first attempt at a slicer that incorperates both CFD and
generative design to allow for the rapid production of 3d printed airframes
the aim of this is to allow for easier access to the hobby through the use of 3D printers
without having to do extensive design work
"""

import numpy as np
import xml.etree.ElementTree as ET
class CFD():
    """This is an implementation of the Lattice Boltzmann simulation"""
    def __init__(self, Density_, Kinematic_Viscocity_, V_inf_, foil_):
        self.Density = Density
        self.mu = Kinematic_Viscocity
        self.V_inf = V_inf_
        self.foil = foil_

class Wing(Foil):
    def __init__(self, root_, tip_, Sweep_, Dihedral_, ):
        Foil.__init__()
        self.Foils =


class Foil():
    def __init__(self, points):




class main():
    def __init__(self):
        self.Load_airframe_defintion()

    def Load_airframe_defintion(self):
        """XFLR5 Stores its airfract definitions a as XML file"""

        #Waiting on response from forum
        pass
if __name__ == "__main__":
    main()
