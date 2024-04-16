import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt

angles = np.array([0,20,40,60])
resolutions = np.array([533,539,566,611,650,883])
#resolutions = np.array([18,16,14,12,10,8])
angles, resolutions = np.meshgrid(angles, resolutions)
optimised_mass = np.array([[0.5984027952668227, 0.6162953483753347, 0.6256867792678744, 0.6037040339219235],
                           [0.5998473176771821, 0.5986107033992938, 0.6155441235193964, 0.6026964883128804],
                           [0.5706846138923475, 0.5767564800990483, 0.5935401435256663, 0.5942101728699716],
                           [0.5379086179012356, 0.5486451024730493, 0.5322901436912285, 0.5337808000662],
                           [0.5177331075553816, 0.5073158995248792, 0.5045517104346265, 0.5094293845457155],
                           [0.4438062145913483, 0.4546148802973241, 0.4435170772595853, 0.4295840860546163]
                            ])
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf = ax.plot_surface(resolutions, angles, optimised_mass, cmap=cm.coolwarm, linewidth=0, antialiased=False)
ax.set_xlabel("Mesh Resolution[mm]")
ax.set_ylabel("Infill Offset angle")
ax.set_zlabel("Optimised Mass Ratio")
plt.title("Hyperparameter Variance")
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()
