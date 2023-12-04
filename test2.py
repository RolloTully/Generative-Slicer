from aeropy.xfoil_module import find_pressure_coefficients, find_coefficients, create_input, file_name
import matplotlib.pyplot as plt
from numpy import *
class Main():
    def __init__(self):
        self.mainloop()
    def generate_naca(self, foil_num, closed = False): #genrates 4 digit naca air foils
            self.foil_num = str(foil_num)
            print(self.foil_num)
            self.max_camber = int(self.foil_num[0])/100
            self.max_camber_pos = int(self.foil_num[1])/10
            self.thickness_ratio = int(self.foil_num[2:])/100
            self.foil_surface = []
            for x in range(100,0,-1):
                self.pp = x/100
                self.thickness = 5*self.thickness_ratio*(0.2969*sqrt(self.pp)-0.1260*self.pp-0.3516*self.pp**2+0.2843*self.pp**3-0.1015*self.pp**4)
                if self.pp<=self.max_camber_pos:
                    if self.max_camber != 0:
                        self.camber_offset = (self.max_camber/self.max_camber_pos**2)*(2*self.max_camber_pos*self.pp-self.pp**2)
                        self.offset_theta = arctan((2*self.max_camber/self.max_camber_pos**2)*(self.max_camber_pos-self.pp))
                    else:
                        self.camber_offset = 0
                        self.offset_theta = 0
                else:
                    if self.max_camber!=0:
                        self.camber_offset = (self.max_camber/(1-self.max_camber_pos)**2)*((1-2*self.max_camber_pos)+2*self.max_camber_pos*self.pp-self.pp**2)
                        self.offset_theta = arctan((2*self.max_camber/(1-self.max_camber_pos)**2)*(self.max_camber_pos-self.pp))
                    else:
                        self.camber_offset = 0
                        self.offset_theta = 0
                self.x_a = self.pp - self.thickness*sin(self.offset_theta)
                self.y_a = self.camber_offset+self.thickness*cos(self.offset_theta)
                self.foil_surface.append([self.x_a,self.y_a])
            for x in range(0,100,1):
                self.pp = x/100
                self.thickness = 5*self.thickness_ratio*(0.2969*sqrt(self.pp)-0.1260*self.pp-0.3516*self.pp**2+0.2843*self.pp**3-0.1015*self.pp**4)
                if self.pp<=self.max_camber_pos:
                    if self.max_camber!=0:
                        self.camber_offset = (self.max_camber/self.max_camber_pos**2)*(2*self.max_camber_pos*self.pp-self.pp**2)
                        self.offset_theta = arctan((2*self.max_camber/self.max_camber_pos**2)*(self.max_camber_pos-self.pp))
                    else:
                        self.camber_offset = 0
                        self.offset_theta = 0
                else:
                    if self.max_camber!=0:
                        self.camber_offset = (self.max_camber/(1-self.max_camber_pos)**2)*((1-2*self.max_camber_pos)+2*self.max_camber_pos*self.pp-self.pp**2)
                        self.offset_theta = arctan((2*self.max_camber/(1-self.max_camber_pos)**2)*(self.max_camber_pos-self.pp))
                    else:
                        self.camber_offset = 0
                        self.offset_theta = 0
                self.x_a = self.pp+self.thickness*sin(self.offset_theta)
                self.y_a = self.camber_offset-self.thickness*cos(self.offset_theta)
                self.foil_surface.append([self.x_a,self.y_a])
            if closed:
                self.foil_surface.append(self.foil_surface[0])
            self.foil_surface = array(self.foil_surface)
            return self.foil_surface

    def convert_to_xfoil(self, foil):
        self.file = open("test_airfoil","w")
        for point in foil:
            self.row = "    "+str(round(point[0],5))+"   "+str(round(point[1],5))+"\n"
            self.file.write(self.row)
        self.file.close
    def Compute_Max_Loading(self,input_foil=[], alphas=[-15,15], Reynolds = 1e6, naca = True):
        if naca:#This will be a string that descirbes its shape so can be passed directly
            pass
        else:#This will be a 2D list of that describes the normalised surface of the airfoil
            pass

        self.Loading_series = []
        for self.angle in range(alphas[0],alphas[1]+1):
            self.Loading_series.append(find_pressure_coefficients(airfoil='naca0012', Reynolds = Reynolds, alpha=self.angle, NACA=naca))
        self.Loading_data = [data["Cp"] for data in self.Loading_series]
        self.Foil = [self.Loading_series[0]["x"],self.Loading_series[0]["y"]]

    def mainloop(self):
        #self.Compute_Max_Loading()
        self.foil = self.generate_naca(2412)
        self.convert_to_xfoil(self.foil)
        self.alfa = -5
        self.reynold = 1e6
        self.foil_data = find_pressure_coefficients(airfoil="naca2412", Reynolds = self.reynold, alpha=self.alfa, NACA=True)
        #print(find_coefficients(airfoil='naca2412',alpha=-5))
        print(self.foil_data.keys())
        #print(foil_data["x"])
        #print(foil_data["y"])
        #print("\nCp-----")
        #print(len(foil_data["Cp"]))
        #plt.plot(self.foil_data["x"],self.foil_data["Cp"])
        #plt.plot(self.foil_data["x"][:80],self.foil_data["y"][:80])
        #plt.plot(self.foil_data["x"][80:],self.foil_data["y"][80:])
        #plt.show()


if __name__ == "__main__":
    Main()
