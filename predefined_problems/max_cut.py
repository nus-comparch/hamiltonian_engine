import os 
os.path.sys.path.append('../hamiltonian_engine/')
from hamiltonian import phase_hamiltonian as phs_ham
from hamiltonian import mixer_hamiltonian as mix_ham
from expectation_value import expectation_value as ex_v

from scipy import optimize as opt
from scipy.optimize import Bounds
from qiskit import *
from qiskit.visualization import plot_histogram
import numpy as np
import networkx as nx

# TODO: Change class to MAX-CUT then add different sub-class types

class max_cut:
    objective_function = '~x_1 & x_2 | ~x_2 & x_1' # Try to state the actual functions 
    objective_functionDi = '~x_1 & x_2'
    variables = ['x_1','x_2']
    graph = None
    shots = 0
    

    def __init__(self, graph:nx.Graph, shots=1024, directed=False):
        self.shots = shots
        self.graph = graph
        self.directed = directed

        if self.directed == False:
            self.phse_ham = phs_ham(self.objective_function,self.variables)
            self.expectation = ex_v(self.objective_function, self.variables, is_graph=True)
        else:
            self.phse_ham = phs_ham(self.objective_functionDi, self.variables)
            self.expectation = ex_v(self.objective_functionDi, self.variables, is_graph=True)

        self.mx_ham = mix_ham(self.phse_ham)

        # generate Phase Hamiltonian
        self.phse_ham.Hamify(boolean=True)


    def generate_quantumCircuit(self, graph:nx.Graph, hyperparams:list):
        if graph == None:
            raise ValueError('Missing Argument: {} for "graph:nx.Graph"'.format(graph))
        else:
            self.circuit = self.phse_ham.perEdgeMap(graph,hyperparams[0],True,True)

            self.circuit += self.mx_ham.generalXMixer(hyperparams[1],0,True,True,graph)

            return self.circuit.draw(output='mpl')
    
    def run_circuit(self, graph:nx.Graph, shots=1024):
        # Add backend for actual quantum chip
        backend      = Aer.get_backend("qasm_simulator")
        print('backend setup: Complete running circuit')

        simulate     = execute(self.circuit, backend=backend, shots=shots)
        results = simulate.result()

        print('Simulation: Complete!')

        print("Expectation Value : {}".format(self.expectation.get_expectationValue(results,shots,graph)))

        return results

    def MAX_CUT(self, hyperparameters:list):
        self.generate_quantumCircuit(self.graph, hyperparameters)

        backend = Aer.get_backend("qasm_simulator")
        simulate     = execute(self.circuit, backend=backend, shots=self.shots)
        results = simulate.result()
        res_maxcut = self.expectation.get_expectationValue(results,self.shots,self.graph)

        return -1 * res_maxcut


    def run_QAOA(self, init_hyperparams:list, method:str):
        #define the bounds for the hyperparameters
        bounds = [[0, 2*np.pi], [0, np.pi]]
        cons = []
        for factor in range(len(bounds)):
            lower, upper = bounds[factor]
            l = {'type': 'ineq',
                'fun': lambda x, lb=lower, i=factor: x[i] - lb}
            u = {'type': 'ineq',
                'fun': lambda x, ub=upper, i=factor: ub - x[i]}
            cons.append(l)
            cons.append(u)
        
        res = opt.minimize(self.MAX_CUT, init_hyperparams,constraints=cons, tol= 1e-3, method=method)

        print(res)

        return res.x        

