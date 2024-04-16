import numpy as np#
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
mesh = [883,650,611,566, 539, 533]
time = [8025,3917, 3174,2298,1126, 1040]
def curve_model(x, a, b, c):
    '''Defines the model used during the curve fitting step'''
    return (a*x**2)+(b*x)+c

Coefficients, _ = curve_fit(curve_model, mesh, time, p0 = [1, 1, 1])
a, b, c = Coefficients
print(Coefficients)
fit_curve = [curve_model(x, a, b, c) for x in range(533,884)]
fig, ax = plt.subplots()
plt.grid()
ax.plot(mesh, time, label = "Measured Data")
ax.plot(range(533,884), fit_curve, label = "Fit line: "+str(a)[0:5]+"$x^2$+"+str(b)[0:5]+"x"+str(c)[0:5])
plt.legend()
ax.set_xlabel("Mesh Elements")
ax.set_ylabel("Optimsation Time [Seconds]")
plt.title("Optimisation time")
#plt.x_label("Grid Resolution")
#plt.y_label("Optimisation Time")
plt.show()
