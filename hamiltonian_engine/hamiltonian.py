from sympy import *
from sympy.abc import I
from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit
import networkx as nx
import numpy as np

import os
os.path.sys.path.append('../hamiltonian_engine/')
from utils import circuit_builder as cir_build

class hamiltonian:

    # Symbols for converting objective function into Pauli expressions
    X = {}
    Y = {}
    Z = {}
    constants = []
    symbolic_map = {}

    # Meta data for circuit generations
    quanCir_list = []
    qubit_map = None

    quantum_circuit = []
    an = QuantumRegister(1, 'ancilla')

    # Equations of phase hamiltonians
    obj_exp = ""
    Hamil_exp = ""
    full_hamiltonian = None


    def __init__(self, expr_str, var):
        if expr_str != None:

            self.__checkVariables(expr_str, var)

            # Initialize expressions and variables
            self.expression_str = expr_str
            self.variables = var

            self.obj_exp = sympify(self.expression_str)

            # Number of pauli objects are limited to the number of variables/qubits being used
            i = 0
            for sym in self.obj_exp.free_symbols:
                if self.variables.index(str(sym)) >= 0:
                    i = self.variables.index(str(sym))
                    ind = self.variables[i].find('_')
                    subscrpt = self.variables[i][ind + 1]      
                    self.Z[sym] = symbols('Z_{}'.format(subscrpt))
                    self.X[sym] = symbols('X_{}'.format(subscrpt))
                    self.Y[sym] = symbols('Y_{}'.format(subscrpt))
                    self.symbolic_map[Not(sym)] = 0.5 * (I + self.Z[sym])
                    self.symbolic_map[sym] = 0.5*(I - self.Z[sym])
                    #i += 1
                else:
                    self.symbolic_map[sym] = sym * I
                    self.constants.append(sym)
        else:
            self.variables = var

    def __add__(self, other):
        if isinstance(other, phase_hamiltonian) & isinstance(self, phase_hamiltonian):
            temp_func = str(self.obj_exp + other.obj_exp - (self.obj_exp * other.obj_exp))
            temp_var  = self.variables + list(set(other.variables) - set(self.variables))  

            return phase_hamiltonian(temp_func, temp_var) 

        elif (isinstance(other, float) | isinstance(other, int)) & isinstance(self, phase_hamiltonian):
            temp_func = str(float(other) + self.obj_exp)

            return phase_hamiltonian(temp_func, self.variables)
        

    def __mul__(self, other):
        if isinstance(other, phase_hamiltonian) & isinstance(self, phase_hamiltonian):
            temp_func = str(self.obj_exp * other.obj_exp)
            temp_var  = self.variables + list(set(other.variables) - set(self.variables))

            return phase_hamiltonian(temp_func, temp_var)

        elif (isinstance(other, float) | isinstance(other, int)) & isinstance(self, phase_hamiltonian):
            temp_func = str(float(other) * self.obj_exp)

            return phase_hamiltonian(temp_func, self.variables)
    
    def __invert__(self):
        return phase_hamiltonian(str(1 - (self.obj_exp)), self.variables)

    def __truediv__(self, other):
        if isinstance(self, hamiltonian) & isinstance(other, hamiltonian):
            assert len(self.quantum_circuit) == len(other.quantum_circuit)
            cirq = None
            for i in range(len(self.quantum_circuit)):
                if i == 0 :
                    cirq = self.quantum_circuit[i] 
                    temp = other.quantum_circuit[i]
                    if temp.has_register(self.an):
                        cirq.add_register(self.an)
                    cirq = cirq + temp
                else:
                    temp1 = self.quantum_circuit[i] 
                    temp2 = other.quantum_circuit[i]
                    if temp.has_register(self.an):
                        temp1.add_register(self.an)
                    cirq = cirq + temp1 + temp2
        return cirq

    # Need to check more cases so as to prevent future errors
    def __checkVariables(self, expression: str, variables: list):
        for v in variables:
            if(v in expression):
                pass
            else:
                raise ValueError(
                    'Variables Mismatch! Unable to find {} in the Objective Function: {}'.format(v, expression))

    def get_objFun(self):
        return self.obj_exp

    def get_pHamil(self):
        if self.full_hamiltonian == None:
            return self.Hamil_exp
        else:
            return self.full_hamiltonian

    def get_qclist(self):
        return self.quanCir_list

    def get_exprstr(self):
        return self.expression_str

    def get_variables(self):
        return self.variables

    def get_qubitMap(self):
        return self.qubit_map

    def get_quantumCircuitAslist(self):
        return self.quantum_circuit
    
    def get_quantumCircuit(self):
        for i in range(len(self.quantum_circuit)):
            if i  == 0:
                c = self.quantum_circuit[i]
            else:
                temp = self.quantum_circuit[i]
                c  = c + temp

        return c


class phase_hamiltonian(hamiltonian):
    def __init__(self, expr_str:str, var):
        super().__init__(expr_str, var)

    def Hamify(self, pwr_args=True, boolean=False):
        if boolean == True:
            self.Hamil_exp = simplify_logic(self.obj_exp)
            self.Hamil_exp = self.Hamil_exp.replace(Or, Add)
            self.Hamil_exp = self.Hamil_exp.replace(And, Mul)

        else:
            self.Hamil_exp = self.obj_exp

        # maps x -> 1/2(I - Zi) or ~x -> 1/2(I + Zi)
        self.Hamil_exp = self.Hamil_exp.xreplace(self.symbolic_map)

        self.Hamil_exp = expand(self.Hamil_exp)

        self.Hamil_exp = self.Hamil_exp.replace(
            lambda expr: expr.is_Pow and (expr.count(I) > 0), lambda expr: expr.base**1)

        # # Remove all Identity matrices and global phases
        for sym in self.obj_exp.free_symbols:
            self.Hamil_exp = self.Hamil_exp.subs(self.Z[sym]* I, self.Z[sym])
        # coeff = self.Hamil_exp.as_coefficients_dict()

        # # Reduce variables with >= power(1) to power(1)
        if pwr_args == True:
            self.Hamil_exp = self.Hamil_exp.replace(
                lambda expr: expr.is_Pow, lambda expr: expr.base**1)

        self.Hamil_exp = self.Hamil_exp.subs([(c, 0) for c in self.constants])

        # Convert to expression into a sympy poly to get list of monomial expression to build the QC
        # However for simplicity we will still reduce expressions with power > 1 to 1
        if pwr_args == False:
            temp = self.Hamil_exp.replace(
                lambda expr: expr.is_Pow, lambda expr: expr.base**1)
            self.quanCir_list = Poly(temp).monoms()
        else:
            temp = self.Hamil_exp.subs(I , 1)
            
            coeff = temp.as_coefficients_dict()
            gbl_phse = coeff.get(1)

            if gbl_phse != None:
                temp = temp - gbl_phse
            
            self.quanCir_list = Poly(temp).monoms()

    # multidimensional qubit mapping for problems involving k-dits
    def perDitMap(self, gamma, p, k_dits, graph: nx.Graph, sub_expr={'i': 'i'}, barrier=False, initial_Hadamard=False):
        assert p == len(gamma)

        self.quantum_circuit = []

        self.qubit_map = cir_build.circuit_builder.map_qubits(self.variables, k_dits, graph)

        no_qubits = len(graph.nodes) * k_dits

        for i in range(p):

            cir = QuantumCircuit(no_qubits)

            if i == 0 and initial_Hadamard == True:
                for j in range(no_qubits):
                    cir.h(j)
                cir.barrier() 

            for e in graph.edges:
                cir += cir_build.circuit_builder.generate_Ditcircuit(self.quanCir_list, gamma[i], self.qubit_map, e, sub_expr, k_dits, no_qubits)
            
            if barrier == True:
                cir.barrier()

            self.quantum_circuit.append(cir)

    # Map the qubits directly to each variable
    def perQubitMap(self, gamma:list, p, barrier=False, initial_Hadamard=False):
        assert p == len(gamma)

        self.quantum_circuit = []

        self.qubit_map = cir_build.circuit_builder.map_qubits(self.variables, 0)

        no_qubits = len(self.qubit_map.values())

        for i in range(p):
            cir = QuantumCircuit(no_qubits)

            if i == 0 and initial_Hadamard == True:
                for j in range(no_qubits):
                    cir.h(j)

            cir += cir_build.circuit_builder.generate_Zcircuit(self.quanCir_list, gamma[i], self.qubit_map)       

            if barrier == True:
                cir.barrier()

            self.quantum_circuit.append(cir)

    # Only for 2 variable Expressions since each edge is an interaction between 2 vertices(qubits)
    def perEdgeMap(self, gamma:list, p:int, graph:nx.Graph,  barrier=False, initial_Hadamard=False):
        assert p == len(gamma) 

        self.__add_defaultWeights(graph)

        self.full_hamiltonian = 0
        self.quantum_circuit = []

        self.qubit_map = cir_build.circuit_builder.map_qubits(self.variables, 0, graph)

        no_qubits = len(self.qubit_map.values())

        for i in range(p):
            cir = QuantumCircuit(no_qubits)

            if i == 0 and initial_Hadamard == True:
                for j in range(no_qubits):
                    cir.h(j)

            if len(self.variables) == 2:
                for e in graph.edges:
                    if i == 0:
                        temp = self.Hamil_exp
                        l = 0
                        for sym in self.Hamil_exp.free_symbols:
                            if not (sym == I):
                                temp = temp.subs(sym, symbols('Z_{}'.format(e[l])))
                                l = (l + 1) % 2
                            
                        self.full_hamiltonian += graph.get_edge_data(e[0],e[1])["weight"]* temp

                    cir += cir_build.circuit_builder.generate_Zcircuit(self.quanCir_list, gamma[i], self.qubit_map, e)
            else:
                if i == 0:
                    for v in graph.nodes:
                        temp = self.Hamil_exp
                        for sym in self.Hamil_exp.free_symbols:
                            if not (sym == I):
                                temp = temp.subs(sym, symbols('Z_{}'.format(v)))
                        self.full_hamiltonian += temp

                cir += cir_build.circuit_builder.generate_Zcircuit(self.quanCir_list, gamma[i], self.qubit_map, edge=(-1,-1))

            if barrier == True:
                cir.barrier()

            self.quantum_circuit.append(cir)

    def __add_defaultWeights(self, graph:nx.Graph, weights=1):
        for u,v in graph.edges:
            if not bool(graph.get_edge_data(u,v)):
                graph.add_edge(u, v, weight = weights)



class mixer_hamiltonian(hamiltonian):
    def __init__(self, var=None, expr_str=None):
        super().__init__(expr_str, var)

    # Mixer Hamiltonian for qubits to have dynamicity between {0,1}
    def generalXMixer(self, betas:list, p:int, qubit_map:dict, measure=False):

        self.quantum_circuit = []

        for i in range(p):
            cir = QuantumCircuit()

            cir += cir_build.circuit_builder.generate_RXcircuit(qubit_map, betas[i])
                
            if measure == True and i == p - 1:
                cir.measure_all()

            self.quantum_circuit.append(cir)


    def controlledXMixer(self, beta:list, p:int, graph: nx.Graph, inverse:bool= False, measure=False):
        # allow for ancillary qubits so that controlled rotations can be performed.
        # Include the controlled X-not version by adding X gates to each side of the controlled qubit line.
        self.quantum_circuit = []

        for i in range(p):
            cir = QuantumCircuit(len(graph.nodes))

            # Get all the q-regs
            quantum_regs = cir.qregs[0]

            for n in graph.nodes:
                bfs = dict(nx.traversal.bfs_successors(graph, n, depth_limit=1))
                for source in bfs:
                    if len(bfs[source]) > 0:
                        control_bits = list(quantum_regs[n] for n in bfs[source])
                        
                        if inverse == True:
                            cir.x(control_bits)

                        cir.mcrx(beta[i], control_bits, quantum_regs[int(source)])
                        
                        if inverse == True:
                            cir.x(control_bits)
                                         
            if measure == True and i == p - 1:
                classical_regs = ClassicalRegister(len(graph.nodes))
                cir.add_register(classical_regs)
                cir.measure(quantum_regs, classical_regs)

            self.quantum_circuit.append(cir)


    # def single_XYMixer(self, xy: str, beta: float, qubit_1: int, qubit_2: int, ancillary_qubit: int = None):
    #     if qubit_1 == qubit_2:
    #         raise ValueError(
    #             "Error: qubit_1 and qubit_2 cannot have the same int values.")
    #     else:
    #         if xy == "xx":
    #             self.mixer_circuit.h(qubit_1)
    #             self.mixer_circuit.h(qubit_2)
    #             self.mixer_circuit.cx(qubit_1, qubit_2)

    #             if ancillary_qubit == None:
    #                 self.mixer_circuit.rz(beta, qubit_2)
    #             else:
    #                 self.mixer_circuit.crz(beta, ancillary_qubit, qubit_2)

    #             self.mixer_circuit.cx(qubit_1, qubit_2)
    #             self.mixer_circuit.h(qubit_1)
    #             self.mixer_circuit.h(qubit_2)

    #         elif xy == 'yy':
    #             x_rot = np.pi / 2
    #             self.mixer_circuit.rx(x_rot, qubit_1)
    #             self.mixer_circuit.rx(x_rot, qubit_2)
    #             self.mixer_circuit.cx(qubit_1, qubit_2)

    #             if ancillary_qubit == None:
    #                 self.mixer_circuit.rz(beta, qubit_2)
    #             else:
    #                 self.mixer_circuit.crz(beta, ancillary_qubit, qubit_2)

    #             self.mixer_circuit.cx(qubit_1, qubit_2)
    #             self.mixer_circuit.rx(x_rot, qubit_1)
    #             self.mixer_circuit.rx(x_rot, qubit_2)

    #         else:
    #             raise ValueError(
    #                 "Incorrect xy string; it can only be either 'xx' or 'yy'.")

    #     return self.mixer_circuit
