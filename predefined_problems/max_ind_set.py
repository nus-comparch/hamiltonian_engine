import os 
os.path.sys.path.append('../hamiltonian_engine/')
from hamiltonian import phase_hamiltonian as phs_ham
from hamiltonian import mixer_hamiltonian as mix_ham
from expectation_value import expectation_value as ex_v

from scipy import optimize as opt
from scipy.optimize import Bounds
from qiskit import *
import numpy as np
import networkx as nx


class max_ind_set:
    objective_function = 'x_1' # Try to state the actual functions 
    variables = ['x_1']
    graph = None
    shots = 0
    qubit_map = None

    def __init__(self, p:int, graph:nx.Graph, shots=1024):
        self.shots = shots
        self.graph = graph
        self.p = p


        self.phse_ham = phs_ham(self.objective_function,self.variables)
        self.expectation = ex_v(self.objective_function, self.variables, is_graph=True)

        self.mx_ham = mix_ham()

        # generate Phase Hamiltonian
        self.phse_ham.Hamify(boolean=True)


    def generate_quantumCircuit(self, graph:nx.Graph, hyperparams:list):
        if graph == None:
            raise ValueError('Missing Argument: {} for "graph:nx.Graph"'.format(graph))
        else:
            assert len(hyperparams) == 2*self.p
            
            l = len(hyperparams)
            gammas = hyperparams[:l//2]
            betas  = hyperparams[l//2:]

            self.phse_ham.perEdgeMap(gammas, self.p, graph, True, True)

            self.qubit_map = self.phse_ham.get_qubitMap()
            self.expectation.use_qubitMap(self.qubit_map)

            self.mx_ham.controlledXMixer(betas, self.p, graph, True, True)

            self.circuit = self.phse_ham / self.mx_ham

            return self.circuit.draw(output='mpl')
    
    def run_circuit(self, shots=1024):
        # Add backend for actual quantum chip
        backend      = Aer.get_backend("qasm_simulator")
        print('backend setup: Complete running circuit')

        simulate     = execute(self.circuit, backend=backend, shots=shots)
        results = simulate.result()

        print('Simulation: Complete!')

        print("Expectation Value : {}".format(self.expectation.get_expectationValue(results,shots,self.graph)))

        return results

    def MAX_IND(self, hyperparameters:list):
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
        
        res = opt.minimize(self.MAX_IND, init_hyperparams,constraints=cons, tol= 1e-3, method=method)

        print(res)

        return res.x 









    