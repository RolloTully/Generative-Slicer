import numpy as np
import serial
import time
class Printer():
    def __init__(self):
        self.Dimensions = [250,250,250]#Dimensions in mm
        self.Speed = 60#mm/min
        self.Filament_Diamiter = 1.75#mm
        self.Nozzle_Diamiter = 0.4#mm
        self.Flow_Factor = ((self.Filament_Diamiter/2)**2)/((self.Nozzle_Diamiter)**2)
        self.extrusion_factor = 0.04
        self.bed_size = np.array([220,220,250])
        self.serial_port = serial.Serial('COM18',115200)
        time.sleep(2)
        self.Ready()
    def prime(self):
        self.cmds = ["G0 X2 Y2 Z0.2\r\n",
                     "G1 E7.92 X2 Y200 Z0.2\r\n",
                     "G1 E7.92 X2.1 Y2 Z0.2\r\n",
                     "G0 X2 Y2 Z10\r\n"]
        for cmd in self.cmds:
            self.serial_port.write(cmd.encode())
            self.waittillready()
        time.sleep(3)
    def home(self):
        self.serial_port.write("G28\r\n".encode())
    def retract(self):
        self.cmd = ["G1 E-25\r\n"]
        for cmd in self.cmds:
            self.serial_port(cmd.encode())
    def goto(self, pos):
        print("Going to ", pos)
        self.command = "G0 X"+str(pos[0])+" Y"+str(pos[1])+" Z"+str(pos[2])+"\r\n"
        self.serial_port.write(self.command.encode())#
        self.waittillready()
    def extrudeto(self, x,y,z,e):
        self.command = "G1 E"+str(e)+" X"+str(x)+" Y"+str(y)+" Z"+str(z)+"\r\n"
        self.serial_port.write(self.command.encode())
        self.waittillready()
    def waittillready(self):
        while True:
            if self.serial_port.read_until() == b'ok\n':
                break
    def setNozzleTemperature(self, T):
        print("Extruder heating to "+str(T)+" degrees.")
        self.cmd = "M109 S"+str(T)+"\r\n"
        self.serial_port.write(self.cmd.encode())
        self.waittillready()
        print("Extruder heated.")
    def setBedTemperature(self, T):
        print("Bed heating to "+str(T)+" degrees.")
        self.cmd = "M190 S"+str(T)+"\r\n"
        self.serial_port.write(self.cmd.encode())
        self.waittillready()
        print("Bed heated.")
    def setRelExtrusion(self):
        self.cmd = "M83\r\n"
        self.serial_port.write(self.cmd.encode())
        self.waittillready()
    def setRelMotion(self):
        self.cmd = "G91\r\n"
        self.serial_port.write(self.cmd.encode())
        self.waittillready()
    def setAbsMotion(self):
        self.cmd = "G90\r\n"
        self.serial_port.write(self.cmd.encode())
        self.waittillready()
        self.setRelExtrusion()
    def Ready(self):
        self.home()
        self.setRelExtrusion()
        self.setNozzleTemperature(260)
        self.setBedTemperature(60)
        self.prime()
        self.waittillready()
    def print(self, layers):
        for layer in layers:#this is each vertical slice
            for element in layer:#this is each contrinious line
                self.goto(element[0])
                for index in range(1, len(element)-1):#this is is each point in the line
                    self.diff = [element[index][0]-element[index-1][0], element[index][1]-element[index-1][1]]
                    self.extrusion_distance = round(np.hypot(self.diff[0],self.diff[1])*self.extrusion_factor,5)
                    print(element[index][0],element[index][1],element[index][2],self.extrusion_distance)
                    self.extrudeto(element[index][0],element[index][1],element[index][2],self.extrusion_distance)


printer = Printer()
'''Generate print path'''
angles = np.linspace(-np.pi,np.pi,500)
path = []
for z in range(2,200,2):
    c1 = []
    c2 = []
    c3 = []
    for angle in angles:
        c1.append([round(30*np.sin(angle)+110,3),round(30*np.cos(angle)+110,3),z/10])
    for angle in angles:
        c2.append([round(50*np.sin(angle)+110,3),round(50*np.cos(angle)+110,3),z/10])
    for angle in angles:
        c3.append([round(70*np.sin(angle)+110,3),round(70*np.cos(angle)+110,3),z/10])
    path.append([c1,c2,c3])
printer.print(path)
