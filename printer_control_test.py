import serial
import time
import numpy as np
def Load_foil():
    raw = open('Airfoils\e591_norm.dat','r').read()
    formatted = [el.split(" ") for el in raw.split('\n')]
    formatted  = [[float(num) for num in list(filter(lambda x:x!='',coord))]for coord in formatted]#list(map(float,self.formatted))
    formatted = np.array(list(filter(lambda x:x!=[],formatted)))
    return formatted
def home():
    serial_port.write("G28\r\n".encode())
def goto(pos):
    command = "G0 X"+str(pos[0])+" Y"+str(pos[1])+" Z"+str(pos[2])+"\r\n"
    serial_port.write(command.encode())#
def extrudeto(x,y,z,e):
    command = "G1 E"+str(e)+" X"+str(x)+" Y"+str(y)+" Z"+str(z)+"\r\n"
    serial_port.write(command.encode())

def waittillready():
    while True:
        if serial_port.read_until() == b'ok\n':
            break
def setNozzleTemperature(T):
    cmd = "M109 S"+str(T)+"\r\n"
    serial_port.write(cmd.encode())
def setBedTemperature(T):
    cmd = "M190 S"+str(T)+"\r\n"
    serial_port.write(cmd.encode())
def setRelExtrusion():
    cmd = "M83\r\n"
    serial_port.write(cmd.encode())
extrusion_factor = 0.03
serial_port = serial.Serial('COM18',115200)
time.sleep(2)
home()
setRelExtrusion()
setNozzleTemperature(260)
waittillready()
setBedTemperature(60)
waittillready()

foil = np.array((Load_foil()*200)+[10,100])
foil = np.hstack((foil, np.full_like(foil, 0.2)))
print(foil)
'''
angles = np.linspace(-np.pi,np.pi,500)
positions = []
for angle in angles:
    #print(angle)
    y = 70*np.sin(angle)+100
    x = 70*np.cos(angle)+100
    z = 0.2
    positions.append([x,y,z])
'''

goto(foil[0])
for z in range(2,100,3):
    for i in range(1,len(foil)-1):
        current_pos = foil[i-1]
        new_pos = foil[i]
        diff = new_pos-current_pos
        extrusion_distance = np.hypot(diff[0],diff[1])*extrusion_factor
        print(new_pos[0],new_pos[1],new_pos[2],extrusion_distance)
        extrudeto(new_pos[0],new_pos[1],z/10,extrusion_distance)
        waittillready()


goto([0,220,30])
serial_port.close()
