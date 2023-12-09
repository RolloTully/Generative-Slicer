import aeropy.xfoil_module as xf
from aeropy.aero_module import Reynolds
import matplotlib.pyplot as plt
foil_name = "test_input_foil"

x=[]
y=[]
file = open("e591_norm.dat","r")
for line in file:
    array = line.strip().split()
    x.append(float(array[0]))
    y.append(float(array[1]))
yu = y[:31]#elements up to the 29th element
yl = y[30:-1]#elements after the 30th element
xu = x[:31]
xl = x[30:-1]
print(x)
print(y)
plt.plot(xu,yu,'r')
plt.plot(xl,yl,'g')
plt.show()

xf.create_input(x,yu,yl,foil_name,different_x_upper_lower=True)

#print("Input created")
#data = xf.find_pressure_coefficients("naca2412", 5., Reynolds=16e5, iteration=10, NACA=True)
#print(data)
