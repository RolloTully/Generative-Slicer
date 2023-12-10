import aeropy.xfoil_module as xf
from aeropy.aero_module import Reynolds
import matplotlib.pyplot as plt
from aeropy.geometry.airfoil import CST, create_x
foil_name = "test_input_foil"
def make_aeropy_file(foil):
    output_file = open(foil_name,'w')
    for line in foil:
        output_file.write("     "+str(format(line[0],'.6f'))+"    "+str(format(line[1],'.6f'))+"\n")
    output_file.close()

foil_array = []
file = open("e591_norm.dat","r")
for line in file:
    array = line.strip().split()
    foil_array.append([float(array[0]),float(array[1])])
make_aeropy_file(foil_array)


#print("Input created")
#data = xf.find_pressure_coefficients("naca2412", 5., Reynolds=16e5, iteration=10, NACA=True)
#print(data)
Data = xf.find_pressure_coefficients(foil_name, 0., iteration=100, NACA=False)
x = Data['x']
y = Data['y']
cp = Data['Cp']
for i in range(0,len(x)):
    plt.scatter(float(x[i]),float(y[i]))
    print(x[i],y[i])
    plt.scatter(float(x[i]),float(cp[i]))
plt.axis('equal')
plt.show()
