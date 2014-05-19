#def clearall():
#    """clear all globals"""
#    for uniquevar in [var for var in globals().copy() if var[0] != "_" and var != 'clearall']:
#        del globals()[uniquevar]
#clearall()

from numpy import arange, array, pi, zeros, empty, copy, sin, cos, nonzero, ravel, append
from matplotlib import pyplot as plt
#import time

phi = arange(0.,2*pi,2*pi/1000)
phis = [0.05*pi, 0.51*pi]

mu0 = 4*pi*1E-7
N_Ni = 9.14E28
N_Gd = 3.02E28
muB = 9.274E-24
M = array([0.7, 2.9])/mu0 #A/m 
K = [0.5E4, 1E4] #J/m**3
theta = pi/6
Js = [1., 0.2]
Jex = -1.
exNi = -Js[0]*mu0*M[0]**2
exGd = -Js[1]*mu0*M[1]**2
exi = -Jex*mu0*M[0]*M[1]

# Define the thickness of the magnetic stack
Ni = [0]
Gd = [1]
layers = 8*Ni + 16*Gd + 8*Ni

# magnetic_energy calculates the total energy (Zeeman, anisotropy, exchange and 
# dipolar of a single domain particle. Energy is per unit volume.
# 
# return the calculated energy and the angle at which the energy is minimu
def magnetic_energy(phi,H,Ku,theta,Ms,Ms1,phi1,J1,Ms2,phi2,J2):
    """
    magnetic_energy(phi, H, Ku, theta, Ms[i], Ms[i+1]...)
    """
    Zeeman = -mu0*Ms*H*cos(phi)
    anisotropy = Ku*sin(phi-theta)**2
    exchange1 = -J1*mu0*Ms1*Ms*cos(phi-phi1)
    exchange2 = -J2*mu0*Ms2*Ms*cos(phi-phi2)
    #dipolar = 0.
    energy = Zeeman + anisotropy + exchange1 + exchange2 #+ dipolar
    Emin = min(energy)
    phi_ind = ravel(nonzero(energy==Emin))
    phi_min = phi[phi_ind]
    return energy, phi_min
    
def create_stack(n,layers,Ms,Ku,J,phi_start):
    Ms_vector = []
    Ku_vector = []
    J_vector = []
    phi_vector = []
    for i in range(0,n,1):
        Ms_vector = append(Ms_vector,Ms[layers[i]])
        Ku_vector = append(Ku_vector,Ku[layers[i]])
        J_vector = append(J_vector,J[layers[i]])
        phi_vector = append(phi_vector,phi_start[layers[i]])
    Ms_vector = append(append(0,Ms_vector),0)
    Ku_vector = append(append(0,Ku_vector),0)
    J_vector = append(append(0,J_vector),0)
    phi_vector = append(append(0,phi_vector),0)
    return Ms_vector, Ku_vector, J_vector, phi_vector

parameters = create_stack(len(layers),layers,M,K,Js,phis)
Ms = parameters[0]
Ku = parameters[1]
J = parameters[2]
phi0 = parameters[3]
n = len(layers)
nGd = len(nonzero(layers)[0])
nNi = n - nGd
J[nNi/2+1] = Jex
J[nGd+nNi/2+1] = Jex

H = arange(-0.5,0.5,.01)/mu0
Mtot_Ni = empty(len(H))
Mtot_Gd = empty(len(H))
Mtot_Ni_bw = empty(len(H))
Mtot_Gd_bw = empty(len(H))
phi_plot = [0]*len(H)
e = [0]*len(H)
phi_check = empty(len(phi0))

iterations = 0;
for k in range(0,len(H),1):
    delta = abs(zeros(len(phi0))-phi0)
    while delta[2:-2].max()>0.001:
        phi_check=copy(phi0)
        for i in range(1,len(Ms)-1,1):
            layer_energy = magnetic_energy(phi,H[k],Ku[i],theta,Ms[i],\
                Ms[i-1],phi0[i-1],J[i],Ms[i+1],phi0[i+1],J[i+1])
            phi0[i] = layer_energy[1][0]
        delta = abs(phi_check - phi0)
        e[k] = layer_energy[0][:]
        iterations+=1
    Mtot_Ni[k] = 5*(sum(Ms[1:nNi/2]*cos(phi0[1:nNi/2])) + \
        sum(Ms[nNi/2+nGd+1:-1]*cos(phi0[nNi/2+nGd+1:-1])))
    Mtot_Gd[k] = sum(Ms[nNi/2+1:nNi/2+nGd]*cos(phi0[nNi/2+1:nNi/2+nGd]))
    Mtot_Ni_bw[k] = -5*((sum(Ms[1:nNi/2]*cos(phi0[1:nNi/2])) + \
        sum(Ms[nNi/2+nGd+1:-1]*cos(phi0[nNi/2+nGd+1:-1]))))
    Mtot_Gd_bw[k] = -sum(Ms[nNi/2+1:nNi/2+nGd]*cos(phi0[nNi/2+1:nNi/2+nGd]))
    phi_plot[k] = phi0

print 'Ground state found after ' + str(iterations) + ' iterations'

# Plot the hysteresis loop for each layer     
plt.figure(1)
plt.plot(mu0*H*1E3,array(Mtot_Ni)/(muB*N_Ni)/nNi,'g')
plt.plot(mu0*H*1E3,array(Mtot_Gd)/(muB*N_Gd)/nGd,'b')
plt.plot(-mu0*H*1E3,array(Mtot_Ni_bw)/(muB*N_Ni)/nNi,'g')
plt.plot(-mu0*H*1E3,array(Mtot_Gd_bw)/(muB*N_Gd)/nGd,'b')
plt.xlabel('$\mu0$ H (mT)')
plt.ylabel('magnetic moment ($\mu_B$/atom)')
plt.xlim([min(1E3*mu0*H), max(1E3*mu0*H)])
plt.grid(True)

# Plot the configuration of each layer as arrows, whose length is proportional
# to the magnetization. 
plt.figure(2)
for i in range(0,len(Ms),1):
    plt.arrow(0.5, 0.01+0.015*i, 0.1*mu0*Ms[i]*cos(phi_plot[4][i]), \
        0.1*mu0*Ms[i]*sin(phi_plot[4][i]),
        fc="k", ec="k",head_width=0.02, head_length=0.03 )
plt.axis('off')
        
plt.show()