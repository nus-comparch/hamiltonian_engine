import networkx as nx
import numpy as np
from qiskit.visualization import plot_histogram
from qiskit import *
from scipy.optimize import Bounds
from scipy import optimize as opt
from expectation_value import expectation_value as ex_v
from hamiltonian import mixer_hamiltonian as mix_ham
from hamiltonian import phase_hamiltonian as phs_ham
import os
os.path.sys.path.append('../hamiltonian_engine/')


class skeletor:
    objective_function = None  # Try to state the actual functions
    variables = None
    graph = None
    shots = 0

    def __init__(self, p: int, obj_fun: str, variables: list, boolean: bool, shots=1024, graph: nx.Graph = None):
        self.shots = shots
        self.graph = graph
        self.p = p
        self.objective_function = obj_fun
        self.variables = variables

        self.phse_ham = phs_ham(self.objective_function, self.variables)

        # generate Phase Hamiltonian
        self.phse_ham.Hamify(boolean=boolean)

        if graph != None:
            self.expectation = ex_v(
                self.objective_function, self.variables, is_graph=True)
        else:
            self.expectation = ex_v(
                self.objective_function, self.variables, is_graph=False)

        self.mx_ham = mix_ham()

    def get_objFun(self):
        return self.phse_ham.get_objFun()

    def get_pHamil(self):
        return self.phse_ham.get_pHamil()

    def generate_quantumCircuit(self, hyperparams: list):
        assert len(hyperparams) == 2*self.p

        l = len(hyperparams)
        gammas = hyperparams[:l//2]
        betas = hyperparams[l//2:]

        if self.graph != None:
            self.phse_ham.perEdgeMap(gammas, self.p, self.graph, True, True)
        else:
            self.phse_ham.perQubitMap(gammas, self.p, True, True)

        phse_map = self.phse_ham.qubit_map

        self.expectation.use_qubitMap(phse_map)

        #self.mx_ham.controlledXMixer(betas, self.p, self.graph, inverse=True, measure=True)
        self.mx_ham.generalXMixer(betas, self.p, phse_map, True)

        self.circuit = self.phse_ham / self.mx_ham

        return self.circuit.draw(output='mpl')

    def run_circuit(self, shots=1024):
        # Add backend for actual quantum chip
        backend = Aer.get_backend("qasm_simulator")
        print('backend setup: Complete running circuit')

        simulate = execute(self.circuit, backend=backend, shots=shots)
        results = simulate.result()

        print('Simulation: Complete!')

        print("Expectation Value : {}".format(
            self.expectation.get_expectationValue(results, shots, self.graph)))

        return results

    def run_skeletor(self, hyperparameters: list):
        self.generate_quantumCircuit(hyperparameters)

        backend = Aer.get_backend("qasm_simulator")
        simulate = execute(self.circuit, backend=backend, shots=self.shots)
        results = simulate.result()
        res_maxcut = self.expectation.get_expectationValue(
            results, self.shots, self.graph)

        return -1 * res_maxcut
    def run_QAOA(self, init_hyperparams: list, method: str):
        # define the bounds for the hyperparameters
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

        res = opt.minimize(self.run_skeletor, init_hyperparams,
                           constraints=cons, tol=1e-3, method=method)

        print(res)

        return res.x
