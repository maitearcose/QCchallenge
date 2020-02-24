#!/usr/bin/env python
# coding: utf-8


import numpy as np




def u(alpha, beta, gamma, delta):
        rx_beta = np.array([[np.cos(beta/2), -1j*np.sin(beta/2)],
             [-1j*np.sin(beta/2),np.cos(beta/2) ]])
        ry_gamma = np.array([[np.cos(gamma/2), -np.sin(gamma/2)],
             [np.sin(gamma/2), np.cos(gamma/2)]])
        rx_delta = np.array([[np.cos(delta/2), -1j*np.sin(delta/2)],
             [-1j*np.sin(delta/2),np.cos(delta/2) ]])
        miau= np.matmul( rx_beta, ry_gamma )
        U= np.power( np.e, 1j*alpha)*np.matmul(miau,rx_delta)
        return U
        
            
        



print(u(0,1,0,-1))



h = np.array([[1/np.sqrt(2), 1/np.sqrt(2)],
            [1/np.sqrt(2), -1/np.sqrt(2)]])




def f(alpha,beta, gamma, delta):
    p= u(alpha, beta, gamma, delta)-h
    return(np.linalg.norm(p))




epsilon = 0.0000001



def grad (alpha,beta, gamma, delta):
    a= f(alpha,beta, gamma, delta)-f(alpha + epsilon ,beta, gamma, delta)
    b= f(alpha,beta, gamma, delta)-f(alpha ,beta + epsilon, gamma, delta)
    c= f(alpha,beta, gamma, delta)-f(alpha ,beta, gamma + epsilon, delta)
    d= f(alpha, beta, gamma, delta)-f(alpha ,beta, gamma, delta + epsilon)
    return(1/epsilon*np.array([a,b,c,d]))



print(grad(1,1,1,1))




wow= np.array([1.0,1.0,1.0,1.0])
print(f(*wow))
c = 0.001
wow= wow +  c*grad(*wow)
print(f(*wow))
wow= wow +  c*grad(*wow)
print(f(*wow))

for i in range(10000):
    wow =  wow +  c*grad(*wow)
print(f(*wow))





print(u(*wow))



print(h)




print(wow)



alpha = wow[0]
beta = wow[1]
gamma = wow[2]
delta = wow[3]
import qiskit as ma
simulator = Aer.get_backend('qasm_simulator')




from qiskit.tools.visualization import plot_bloch_multivector


circuit= ma.QuantumCircuit(2,2)

circuit.rx(delta,0)
circuit.ry(gamma,0)
circuit.rx(beta,0)
circuit.cx(0,1)
circuit.measure([0,1], [0,1])
circuit.measure([0,1], [0,1])
job = ma.execute(circuit, simulator, shots=100000)
result = job.result()
counts = result.get_counts(circuit)
print("\nTotal count for 00 and 11 are:",counts)
circuit.draw()




#to get - instead of + you do an x gate at the beginning on the first qubit. 







