import numpy as np
class LBM():
    def __init__(self, foil):
        self.Foil_address = "C:\Users\rollo\Documents\GitHub\Generative-Slicer\Airfoils\s1225_norm.dat"
        self.foil = self.Load_foil()
        self.mainloop()

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
    def Generate_Foil_Mask(self):
        """Generates a logical mesh of the airfoil"""
    def Generate_Lattice(self):
        """Creates the LBM lattice"""
    def Streaming_step(self):
        """Defines the LBM streaming step"""
    def Collision_step(self):
        """Defines the LBM collision step"""


    def mainloop():

if __name__ == "__main__":
    LBM()
