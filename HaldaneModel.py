import numpy as np
import matplotlib.pyplot as plt
from math import cos,sin,sqrt,pi

def Ham_Haldane(ky,M,J,J2,phi,Nx):
    M1 = 2*J2*(cos(phi)*cos(sqrt(3)*ky)-sin(phi)*sin(sqrt(3)*ky)) + M
    M2 = 2*J2*(cos(phi)*cos(sqrt(3)*ky)-sin(phi)*sin(sqrt(3)*ky)) - M
    Jp = 2*J*cos(sqrt(3)*ky/2)
    J2p = 2*J2*cos(phi-sqrt(3)*ky/2)

    d = []
    for i in range(2*Nx):
        if i % 2 == 0:
            d.append(M1)
        else:
            d.append(M2)
    ds = []
    for i in range(2*Nx-1):
        if i % 2 == 0:
            ds.append(J)
        else:
            ds.append(Jp)
    dss = []
    for i in range(2*Nx-2):
        dss.append(J2p)

    H = np.diag(d,0)+np.diag(ds,1)+np.diag(ds,-1)+np.diag(dss,-2)+np.diag(dss,2)

    return H

if __name__ == "__main__":
    Nx = 12
    Ny = 40
    J = 1
    J2 = 0
    phi = pi/2
    M = 1

    ky_list = np.linspace(0,1,Ny+1,endpoint=True)*2*pi/sqrt(3)
    Energy = np.empty([2*Nx,Ny+1])
    for i in range(Ny+1):
        ky = ky_list[i]
        H = Ham_Haldane(ky,M,J,J2,phi,Nx)
        Energy[:,i],States = np.linalg.eigh(H)
    
    fig = plt.figure()
    for i in range(2*Nx):
        plt.plot(ky_list,Energy[i])
    plt.xlim(0,2*pi/sqrt(3))
    plt.xlabel('ky')
    plt.ylabel('Energy')
    fig.savefig('EnergyHaldane.pdf')

    
