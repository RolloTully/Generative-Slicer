import numpy as np
import cv2
class Foil():
    def __init__(self, address_, chord_):
        self.foil_address = address_
        self.Chord = chord_
        self.Foil = self.Generate_Foil()
        self.thickness = self.Find_thickness()
        print(self.thickness)
    def Find_thickness(self):
        """Finds the thickness of the current foil"""
        return (np.max(self.Foil[:,1])-np.min(self.Foil[:,1]))*self.Chord

    def format_dat(self,data):
        # TODO: Needs Improvment for greater compatibility
        self.data = data.split("\n")
        self.formatted = [self.el.split(" ") for self.el in self.data]
        self.formatted  = [[float(self.num) for self.num in list(filter(lambda x:x!='',self.coord))]for self.coord in self.formatted]#list(map(float,self.formatted))
        self.formatted = list(filter(lambda x:x!=[],self.formatted))
        return self.formatted
    def Load_foil(self):
        """Load the Airfoil for testing"""
        self.raw = open(self.foil_address,'r').read()
        self.foil_dat = np.array(self.format_dat(self.raw))[2:-2]
        return self.foil_dat
    def Generate_Foil(self):
        return self.Load_foil()*self.Chord
    def Bounding_box(self):
        return np.array([np.max(self.Foil[:,0])-np.min(self.Foil[:,0]),
                         np.max(self.Foil[:,1])-np.min(self.Foil[:,1])])

class Mesh_Element():
    def __init__(self):
        self.Flow_State = np.array([0,0])
        self.Density = 1.
class Source_Element():
    """Defines the Inlet face"""
    """This has a fixed flow condition and will not vary"""
    def __init__(self):
        self.Flow_State = np.array([1,0])
class Surface_Element():
    """Defines the surface of the test element"""
    def __init__(self):

class Sink_Element():
    """Defines the outlet face"""
    def __init__(self):

class LBM():
    def __init__(self):
        self.Foil_address = "C:/Users/rollo/Documents/GitHub/Generative-Slicer/Airfoils/s1225_norm.dat"
        self.foil = Foil(self.Foil_address,0.2)
        self.element_size = 0.001 #Meters, so 1mm
        self.mainloop()

    def Generate_Foil_Mask(self):
        """Generates a logical mesh of the airfoil"""
        self.Foil_Dimensions = self.foil.Bounding_box()/self.element_size
        self.Mesh_Dimensions = np.int32(self.Foil_Dimensions*np.array([4,8]))
        self.foil_mask = np.zeros(tuple(self.Mesh_Dimensions[::-1]))#Creates a 0 initialised array
        self.Foil_copy = self.foil.Foil*(1/self.element_size)
        self.offset = np.array([self.Mesh_Dimensions/2])-self.foil.Bounding_box()/(2*self.element_size)
        self.Foil_copy = self.Foil_copy+self.offset
        self.Foil_copy = np.int32([self.Foil_copy])
        self.foil_mask = cv2.fillPoly(self.foil_mask,self.Foil_copy,color=255)
        self.foil_mask = (self.foil_mask>0)
        #self.foil_mask = (self.foil_mask>125)*1
        cv2.imshow("Foil", self.foil_mask)
        cv2.waitKey(2000)
    def Generate_Lattice(self):
        """Creates the LBM lattice"""
        pass

    def Streaming_step(self):
        """Defines the LBM streaming step"""
        pass
    def Collision_step(self):
        """Defines the LBM collision step"""
        pass


    def mainloop(self):
        self.Generate_Foil_Mask()
        pass

if __name__ == "__main__":
    LBM()
