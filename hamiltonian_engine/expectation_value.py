from sympy import *
from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit
import networkx as nx
import numpy as np

class expectation_value:
    
    obj_exp = ""
    variables = []
    is_graph = False
    is_boolean = False

    def __init__(self, objective_function:str, variables:list , is_graph:bool):
        self.obj_exp = sympify(objective_function)
        self.variables = variables
        self.is_graph = is_graph
        #self.qubit_map = qubit_map
        self.is_boolean = self.obj_exp.is_Boolean

    # If no weights have been set all the edges in the graphs without a weight defaults to 1
    def __add_defaultWeights(self, graph:nx.Graph, weights=1):
        for u,v in graph.edges:
            if not bool(graph.get_edge_data(u,v)):
                graph.add_edge(u, v, weight = weights)

    def use_qubitMap(self, qubit_map:dict):
        self.qubit_map = qubit_map

    def get_expectationValue(self, ciruit_results, shots:int, graph:nx.Graph=None): 
        
        counts = ciruit_results.get_counts()

        expectation_value = 0
        
        #Assumption: the coefficients are already in the obj function thus no weights are required
        #            weights are assumed to used for graphs thus adding weights to the edges that 
        #            are connected by 2 vertices.

        if self.is_graph == False:
            for bitstr in counts:
                temp_func = self.obj_exp

                for sym in self.obj_exp.free_symbols:
                    temp_func = temp_func.subs(sym, bitstr[self.qubit_map[str(sym)]])      

                expectation_value = expectation_value + (counts[bitstr] / shots) * temp_func
        else:
            if graph != None:
                self.__add_defaultWeights(graph, 1)
                if len(self.variables) > 1:
                    edges = graph.edges
                    
                    for bitstr in counts:
                        reverse_bitstr = bitstr[::-1]
                        graph_value = 0

                        for e in edges:
                            temp_func1 = self.obj_exp
                            temp_func1 = temp_func1.subs(self.variables[0], reverse_bitstr[e[0]]) 
                            temp_func1 = temp_func1.subs(self.variables[1], reverse_bitstr[e[1]])

                            if self.is_boolean:
                                if temp_func1:
                                    temp_func = 1
                                else:
                                    temp_func = 0
                                graph_value = graph_value + graph.get_edge_data(e[0],e[1])["weight"] * temp_func    
                            else:
                                graph_value = graph_value + graph.get_edge_data(e[0],e[1])["weight"] * temp_func1
                            
                        
                        expectation_value = expectation_value + counts[bitstr] * graph_value      
                else:
                    vertices = graph.nodes
                    
                    for bitstr in counts:
                        reverse_bitstr = bitstr[::-1]
                        graph_value = 0

                        for v in vertices:
                            temp_func1 = self.obj_exp
                            temp_func1 = temp_func1.subs(self.variables[0], reverse_bitstr[int(self.qubit_map[v])])

                            if self.is_boolean:
                                if temp_func1:
                                    temp_func = 1
                                else:
                                    temp_func = 0
                                graph_value = graph_value + temp_func    
                            else:
                                graph_value = graph_value + temp_func1
                            
                        expectation_value = expectation_value + counts[bitstr] * graph_value   
            else:
                raise ValueError('Missing Argument: {} for "graph:nx.Graph"'.format(graph))

        return float(expectation_value / shots)





