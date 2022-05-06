import numpy as np
from qiskit import *
from qiskit import Aer, transpile
from qiskit.tools.jupyter import *
from matplotlib import *
from pylatexenc import *
import itertools
from mitiq import *
from mitiq import zne
from typing import List
from tqdm import *

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

def convert_site(site,N):
    # this function perform the site index conversion as follows:
    # 01234, then convert, e.g. 0 to 4, 1 to 3, etc.
    if N%2 == 0:
        return int(N/2-site+1)
    elif N%2 == 1:
        M = (N-1)/2
        return int(M + (M-site))

def spin_combinations(j,N):
    # j is the site index which counts from the l.h.s. in the physical dimension
    # we convert j to ibm_site which counts from the r.h.s.
    combi = list(itertools.product(['0', '1'], repeat=N-1))
    bits = []
    for i in range(2**(N-1)):
        bits.append(''.join(combi[i]))
    spin_up = []
    ibm_site = convert_site(j,N)
    for i in range(2**(N-1)):
        spin_up.append(bits[i][:ibm_site] + '0' +  bits[i][ibm_site:]) # spin up is '0'. previously '1'
    spin_down = []
    for i in range(2**(N-1)):
        spin_down.append(bits[i][:ibm_site] + '1' +  bits[i][ibm_site:]) # spin down is '1'. previously '0'
    return spin_up,spin_down

def expectation_value(N,site,qobj,shot_num,circuit):
    spin_up,spin_down = spin_combinations(site,N)
    spin_up_amp = []
    spin_down_amp = []
    for i in range(2**(N-1)):
        spin_up_amp.append(qobj.get_counts(circuit).get(spin_up[i]))
        spin_down_amp.append(qobj.get_counts(circuit).get(spin_down[i]))
    spin_down_amp = [i/shot_num for i in spin_down_amp if i]
    spin_up_amp = [i/shot_num for i in spin_up_amp if i]
    return sum(spin_up_amp) - sum(spin_down_amp)

def run_zne(IsSimulator,simulator,backend,shot_num,N,qc_list):

    def batch_run_zne(circuits: List[QuantumCircuit]) -> List[float]:

        def executor(circuit: QuantumCircuit,backend_name: str = "aer_simulator", shots: int = shot_num) -> float:
            if IsSimulator == 'Yes':
                qobj = execute(circuit, simulator,seed_simulator=4, shots = shot_num).result()
            elif IsSimulator == 'No':
                qobj = execute(circuit, backend, shots = shot_num).result()
            
            return expectation_value(N,1,qobj,shot_num,circuit)

        return [zne.execute_with_zne(circuits[j],executor) for j in tqdm(range(len(circuits)))]
    
    return Executor(batch_run_zne, max_batch_size=75).evaluate(qc_list,force_run_all = True)

def run_no_zne(IsSimulator,simulator,backend,shot_num,N,qc_list):

    def batch_run_no_zne(circuits: List[QuantumCircuit]) -> List[float]:

        def executor(circuit: QuantumCircuit,backend_name: str = "aer_simulator", shots: int = shot_num) -> float:
            if IsSimulator == 'Yes':
                qobj = execute(circuit, simulator,seed_simulator=4, shots = shot_num).result()
            elif IsSimulator == 'No':
                qobj = execute(circuit, backend, shots = shot_num).result()
            
            return expectation_value(N,1,qobj,shot_num,circuit)

        return [executor(circuits[j]) for j in tqdm(range(len(circuits)))]
    
    return Executor(batch_run_no_zne, max_batch_size=75).evaluate(qc_list,force_run_all = True)