import aeropy.xfoil_module as xf
from aeropy.aero_module import Reynolds
foil_name = "test_input_foil"

x=[]
y=[]
file = open("e591_norm.dat","r")
for line in file:
    array = line.strip().split()
    x.append(float(array[0]))
    y.append(float(array[1]))
x=x
y=y
print(y[:30])
print(y[29:])

print(len(x))
try:
    xf.create_input(x,y[30:],y[:30],foil_name,different_x_upper_lower=False)
except IndexError:
    pass
data = xf.find_pressure_coefficients(foil_name, 5., Reynolds=1600, iteration=10, NACA=False)
print(data)
