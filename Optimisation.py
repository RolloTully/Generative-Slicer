import numpy as np
import matplotlib.pyplot as plt
E = 4.107e9
p = 1240
s_lim = 13.8e9

#Assign geometry variable
def AssignGeometry(nodes, Var,d,d1):
    NCi = Var[d1:d]          #Node Conditions
    print(NCi)
    #Assign geometry variable
    nodes[1,0] = NCi[0]      #Node 2 (Axis X)
    nodes[5,0] = NCi[0]      #Node 6 (Axis X)
    nodes[2,0] = NCi[1]      #Node 3 (Axis X)
    nodes[6,0] = NCi[1]      #Node 7 (Axis X)
    nodes[1,1] = NCi[2]      #Node 2 (Axis Y)
    nodes[2,1] = NCi[3]      #Node 3 (Axis Y)
    nodes[3,1] = NCi[4]      #Node 4 (Axis Y)
    nodes[5,1] = NCi[5]      #Node 6 (Axis Y)
    nodes[6,1] = NCi[6]      #Node 7 (Axis Y)
    nodes[7,1] = NCi[7]      #Node 8 (Axis Y)
    return nodes

#%% Truss structure analysis
def TrussAnalysis(nodes, bars, Var, d, d1, DOFCON, P, Ur, Ai):
    NN = len(nodes)          #Number of nodes
    NE = len(bars)           #Number of bars
    DOF = 2                  #2D Truss
    NDOF = DOF*NN            #Total number of degree of freedom

    A = np.copy(Var[0:d1])
    for i in range(NE):
        A[i] = Ai[Var[i].astype(int)]
    nodes = AssignGeometry(nodes, Var, d, d1)

    #Structural analysis
    d = nodes[bars[:,1],:] - nodes[bars[:,0],:]
    L = np.sqrt((d**2).sum(axis=1))
    angle = d.T/L
    a = np.concatenate((-angle.T,angle.T), axis=1)
    K = np.zeros([NDOF,NDOF])         #Stiffness matrix
    for k in range(NE):
        #DOFs
        aux = 2*bars[k,:]
        index = np.r_[aux[0]:aux[0]+2, aux[1]:aux[1]+2]
        #Stiffness matrix
        ES = np.dot(a[k][np.newaxis].T*E*A[k],a[k][np.newaxis])/L[k]
        K[np.ix_(index,index)] = K[np.ix_(index,index)] + ES
    freeDOF = DOFCON.flatten().nonzero()[0]
    print(DOFCON)
    supportDOF = (DOFCON.flatten() == 0).nonzero()[0]
    Kff = K[np.ix_(freeDOF,freeDOF)]
    Pf = P.flatten()[freeDOF]
    print(Kff, Pf)
    Uf = np.linalg.solve(Kff,Pf)
    U = DOFCON.astype(float).flatten()
    print(supportDOF)
    U[freeDOF] = Uf
    U[supportDOF] = Ur
    U = U.reshape(NN,DOF)                                    #Displacement (in)
    u = np.concatenate((U[bars[:,0]],U[bars[:,1]]),axis=1)
    N = E*A[:]/L[:]*(a[:]*u[:]).sum(axis=1)                     #Force (kips)
    S = N/A
    Mass = (p*A*L).sum()
    return S, Mass

def limitV(V,d, Vlim):
    for j in range(d):
        if V[j] > Vlim[j,1]:
            V[j] = Vlim[j,1]
        if V[j] < Vlim[j,0]:
            V[j] = Vlim[j,0]
    return V

def FindNearest(array,value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


class Particle:
    def __init__(self,nodes, bars, ps, d, d1, A_index, Xlim, Vlim, MaxIt, DOFCON, P, Ur, Ai):
        self.Nodes = nodes
        self.Bars = bars
        self.ps = p
        self.d = d
        self.d1 = d1
        self.Xlim  = Xlim
        self.Vlim = Vlim
        self.MaxIt = 500
        self.ps = 9
        self.c1 = 2
        self.c2 = 2
        self.DOFCON = DOFCON
        self.P = P
        self.Ur = Ur
        self.Ai = Ai
        self.w = 0.9-((0.9-0.4)/self.MaxIt)*np.linspace(0,self.MaxIt,self.MaxIt)
        self.A_index = A_index
        self.position = np.zeros([ps,d])
        self.velocity = np.zeros([ps,d])
        self.cost = np.zeros(ps)
        self.stress = np.zeros([ps,d1])
        for i in range(ps):
            for j in range(d):
                if j < d1:
                    self.position[i,j] = np.random.choice(A_index[A_index > 15])
                else:
                    self.position[i,j] = np.random.uniform(Xlim[j,0],Xlim[j,1])
                self.velocity[i,j] = np.random.uniform(Vlim[j,0],Vlim[j,1])
            self.stress[i], self.cost[i] = TrussAnalysis(self.Nodes, self.Bars, self.position[i],d,d1,DOFCON, P, Ur, Ai)
        print(self.stress)
        self.pbest = np.copy(self.position)
        self.pbest_cost = np.copy(self.cost)
        self.index = np.argmin(self.pbest_cost)
        self.gbest = self.pbest[self.index]
        self.gbest_cost = self.pbest_cost[self.index]
        self.BestCost = np.zeros(MaxIt)
        self.BestPosition = np.zeros([MaxIt,d])
    def Evaluate(self):
        for it in range(self.MaxIt):
            for i in range(self.ps):
                self.velocity[i] = (self.w[it]*self.velocity[i]
                                   + self.c1*np.random.rand(self.d)*(self.pbest[i] - self.position[i])
                                   + self.c2*np.random.rand(self.d)*(self.gbest - self.position[i]))
                self.velocity[i] = limitV(self.velocity[i],self.d,self.Vlim)
                self.position[i] = self.position[i] + self.velocity[i]
                for p in range(self.d1):
                    self.position[i,p] = FindNearest(self.A_index,self.position[i,p])
                self.stress[i], self.cost[i] = TrussAnalysis(self.Nodes, self.Bars, self.position[i],self.d,self.d1,self.DOFCON, self.P, self.Ur, self.Ai)
                self.C_total = 0
                for cd in range(self.d1):
                    if np.abs(self.stress[i,cd]) > s_lim:
                        self.C1 = np.abs((self.stress[i,cd] - s_lim)/s_lim)
                    else:
                        self.C1 = 0
                    self.C_total = self.C_total + self.C1
                phi = (1 + self.C_total)
                self.cost[i] = self.cost[i]*phi
                if self.cost[i] < self.pbest_cost[i]:
                    self.pbest[i] = self.position[i]
                    self.pbest_cost[i] = self.cost[i]
                    if self.pbest_cost[i] < self.gbest_cost:
                        self.gbest = self.pbest[i]
                        self.gbest_cost = self.pbest_cost[i]
            self.BestCost[it] = self.gbest_cost
            self.BestPosition[it] = self.gbest
        print(self.Nodes)
        return self.Nodes
        plt.scatter(self.Nodes[:,0],self.Nodes[:,1],c='red')
        plt.show()
class Optimisation():
    def __init__(self, Nodes, Bars):#, Forces, Grounds):
        self.MaxIt = 1000
        self.ps = 100
        self.c1 = 2
        self.c2 = 2
        self.w = 0.9-((0.9-0.4)/self.MaxIt)*np.linspace(0,self.MaxIt,self.MaxIt)
        self.nodes = np.array(Nodes).astype(np.float32)
        self.bars = np.array(Bars)
        self.Ur = [0, 0, 0, 0]
        self.DOFCON = np.ones_like(self.nodes).astype(int)
        self.DOFCON[0,:] = 0
        self.DOFCON[4,:] = 0

        self.P = np.zeros_like(self.nodes)
        self.P[4,1] = -20

        #for i in range(Grounds[0], Grounds[1]):
        #    self.DOFCON[i,:] = 0

        self.Ai = [0.111, 0.141, 0.174, 0.220, 0.270, 0.287, 0.347, 0.440, 0.539, 0.954,
              1.081, 1.174, 1.333, 1.488, 1.764, 2.142, 2.697, 2.800, 3.131, 3.565,
              3.813, 4.805, 5.952, 6.572, 7.192, 8.525, 9.300, 10.850, 13.330, 14.290,
              17.170, 19.180]
        self.A_index = np.linspace(0,len(self.Ai)-1,len(self.Ai)).astype(int)
        '''How Much each node can move. for convex point this is huge, for boundary points this must be small'''
        self.NC = []                                   #Node Conditions
        self.NC.append([220, 260])                     #X3 [min, max]
        self.NC.append([100, 140])                     #X2 [min, max]
        self.NC.append([100, 140])                     #Y2 [min, max]
        self.NC.append([100, 140])                     #Y3 [min, max]
        self.NC.append([50, 90])                       #Y4 [min, max]
        self.NC.append([-20, 20])                      #Y6 [min, max]
        self.NC.append([-20, 20])                      #Y7 [min, max]
        self.NC.append([20, 60])                       #Y8 [min, max]
        self.NC = np.array(self.NC)
        plt.scatter(self.nodes[:,0],self.nodes[:,1])
        plt.show()
        self.d1 = len(self.bars)
        self.d2 = len(self.NC)
        self.d = self.d1 + self.d2

        self.Xlim = [min(self.A_index), max(self.A_index)]*np.ones([len(self.bars),2])
        self.Xlim = np.concatenate([self.Xlim,self.NC],axis=0)
        self.Vlim = np.zeros([self.d,2])
        self.Vlim[:,1] = (self.Xlim[:,1] - self.Xlim[:,0])*0.2
        self.Vlim[:,0] = -self.Vlim[:,1]

        self.MaxIt = 500
        self.ps = 9
        self.c1 = 2
        self.c2 = 2
        self.w = 0.9-((0.9-0.4)/self.MaxIt)*np.linspace(0,self.MaxIt,self.MaxIt)

        a = Particle(self.nodes, self.bars, self.ps, self.d, self.d1, self.A_index,self.Xlim, self.Vlim,self.MaxIt,self.DOFCON, self.P, self.Ur, self.Ai)
        a.Evaluate()

nodes = []
bars = []

nodes.append([0,120])
nodes.append([120,120])
nodes.append([240,120])
nodes.append([360,120])
nodes.append([0,0])
nodes.append([120,0])
nodes.append([240,0])
nodes.append([360,0])
nodes.append([440,0])

bars.append([0,1])
bars.append([1,2])
bars.append([2,3])
bars.append([4,5])
bars.append([5,6])
bars.append([6,7])
bars.append([5,1])
bars.append([6,2])
bars.append([7,3])
bars.append([0,5])
bars.append([4,1])
bars.append([1,6])
bars.append([5,2])
bars.append([2,7])
bars.append([6,3])
bars.append([7,8])
bars.append([8,6])
Optimisation(nodes, bars)
