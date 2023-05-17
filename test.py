import numpy as np
import cv2
import cProfile
import re
cProfile.run('re.compile("foo|bar")')
# Simulation parameters
Nx          = 400    # resolution x-dir
Ny          = 100    # resolution y-dir
rho0        = 50    # average density
tau         = 0.6    # collision timescale
Nt          = 40000   # number of timesteps

# Lattice speeds / weights
NL = 9
idxs = np.arange(NL)
cxs = np.array([0, 0, 1, 1, 1, 0,-1,-1,-1])
cys = np.array([0, 1, 1, 0,-1,-1,-1, 0, 1])
weights = np.array([4/9,1/9,1/36,1/9,1/36,1/9,1/36,1/9,1/36]) # sums to 1
X, Y = np.meshgrid(range(Nx), range(Ny))

# Initial Conditions - flow to the right with some perturbations
F = np.ones((Ny,Nx,NL)) + 0.1*np.random.randn(Ny,Nx,NL)
F[:,:,3] += 2 * (1+0.2*np.cos(2*np.pi*X/Nx*4))
rho = np.sum(F,2)
for i in idxs:
  F[:,:,i] *= rho0 / rho

# Cylinder boundary
cylinder = (X - Nx/4)**2 + (Y - Ny/2)**2 < (Ny/4)**2

print(cylinder)
# Simulation Main Loop
iter=0
while True:
    iter=iter+1
    # Drift
    for i, cx, cy in zip(idxs, cxs, cys):
      F[:,:,i] = np.roll(F[:,:,i], cx, axis=1)
      F[:,:,i] = np.roll(F[:,:,i], cy, axis=0)

    # Set reflective boundaries
    bndryF = F[cylinder,:]
    bndryF = bndryF[:,[0,5,6,7,8,1,2,3,4]]

    # Calculate fluid variables
    rho = np.sum(F,2)
    ux  = np.sum(F*cxs,2) / rho
    uy  = np.sum(F*cys,2) / rho
    #velocity_field *= 255.0/velocity_field.max()
    #rho2 = rho
    print(np.int8(rho))
    rho2 = np.int8(rho)*cylinder
    rho2 = np.interp(rho2, (rho2.min(), rho2.max()), (10, 250))
    
    cv2.imshow("Density", np.int8(rho2)*cylinder)
    cv2.waitKey(1)
    print("Saving image")
    print("C:/Users/rollo/Documents/GitHub/Generative-Slicer/Images/"+str(iter)+".JPG")
    cv2.imwrite("C:/Users/rollo/Documents/GitHub/Generative-Slicer/Images/"+str(iter)+".JPG", np.int8(rho))
    # Apply Collision
    Feq = np.zeros(F.shape)
    for i, cx, cy, w in zip(idxs, cxs, cys, weights):
      Feq[:,:,i] = rho*w* (1 + 3*(cx*ux+cy*uy) + 9*(cx*ux+cy*uy)**2/2 - 3*(ux**2+uy**2)/2)

    F += -(1.0/tau) * (F - Feq)

    # Apply boundary
    F[cylinder,:] = bndryF
