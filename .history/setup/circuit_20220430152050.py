import numpy as np
from qiskit import QuantumCircuit

class DTC_Circ:

    def __init__(self,
                 N=None,
                 J=None,
                 t=None,
                 dt=None,
                 lamb=None,
                 omega=None,
                 h=None
                 ):
        self.parameters = {
            "N": N,
            "J": J,
            "t": t,
            "dt": dt,
            "lamb": lamb,
            "omega": omega,
            "h": h,
        }
        self.dtc_circ = QuantumCircuit(N, N, name='N=' + str(N) +
                                       ',J=' + str(round(J,2)) +
                                       ',h=' + str(round(h,2)) +
                                       ',lamb=' + str(round(lamb,2)) +
                                       ',omega=' + str(round(omega,2)) +
                                       ',dt=' + str(dt) +
                                       ',step=' + str(int(t/dt))
                                       )
        self.N = N
        self.J = J
        self.lamb = lamb
        self.omega = omega
        self.h = h
        self.t = t
        self.dt = dt

    def prepare(self):
        for i in range(self.N):
            self.dtc_circ.h(i)
            self.dtc_circ.s(i)

    def measure(self):
        for i in range(self.N):
            self.dtc_circ.sdg(i)
            self.dtc_circ.h(i)
            self.dtc_circ.measure(i, i)

    def hamiltonian(self, t):

        def coeff(h, t, omega):
            return -h*np.cos(omega*t/2)**2

        # H0 terms with sigma^y and sigma^z
        for i in range(self.N):
            self.dtc_circ.ry(2*self.lamb*self.dt, i)
            self.dtc_circ.rz(2*self.lamb*self.dt, i)

        self.dtc_circ.barrier()
        # H1 terms with sigma^z
        for i in range(0, self.N-1, 2):
            self.dtc_circ.cx(i, i+1)
            self.dtc_circ.rz(-2*self.J*self.dt, i+1)
            self.dtc_circ.cx(i, i+1)

        for i in range(1, self.N-1, 2):
            self.dtc_circ.cx(i, i+1)
            self.dtc_circ.rz(-2*self.J*self.dt, i+1)
            self.dtc_circ.cx(i, i+1)

        self.dtc_circ.barrier()
        # H2 time dependent terms with sigma^x
        for i in range(self.N):
            self.dtc_circ.rx(2*coeff(self.h, t , self.omega)*self.dt, i)

    def generate_circuit(self):
        self.prepare()
        for i in range(int(self.t//self.dt)):
            self.hamiltonian(i*self.dt)
        self.measure()
        
def prepare(circuit,N):
    for i in range(N):
        circuit.h(i)
        circuit.s(i)

def measure(circuit,N):
    for i in range(N):
        circuit.sdg(i)
        circuit.h(i)
        circuit.measure(i,i)

def hamiltonian(circuit,N,t,dt,lamb,J,h,omega):
    
    def coeff(t,h,omega):
        return -h*np.cos(omega*t/2)**2

    #H0 terms with sigma^y and sigma^z
    for i in range(N):
        circuit.ry(2*lamb*dt,i)
        circuit.rz(2*lamb*dt,i)

    circuit.barrier()
    #H1 terms with sigma^z
    for i in range(0,N-1,2):
        circuit.cx(i,i+1)
        circuit.rz(-2*J*dt,i+1)
        circuit.cx(i,i+1)

    for i in range(1,N-1,2):
        circuit.cx(i,i+1)
        circuit.rz(-2*J*dt,i+1)
        circuit.cx(i,i+1)

    circuit.barrier()
    #H2 time dependent terms with sigma^x
    for i in range(N):
        circuit.rx(2*coeff(t,h,omega)*dt,i)