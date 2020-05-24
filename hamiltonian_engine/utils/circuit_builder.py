from qiskit import *
import networkx as nx


class circuit_builder:

    @staticmethod
    def map_qubits(variables:list, k_dit:int ,graph:nx.graph=None):
        map = {}
        if graph == None:
            i = 0
            for v in variables:
               map[v] =  i
               i += 1

        elif graph != None and k_dit == 0:
            for n in graph.nodes:
                map[n] = n

        else:
            no_of_qubits = len(graph.nodes) * k_dit
            i = 0
            for v in graph.nodes:
                qbits = tuple()
                for j in range(k_dit):
                    qbits = qbits + (i,)
                    i += 1
                map[v] = qbits
            assert i == no_of_qubits

        return map
        


    @staticmethod
    def generate_Zcircuit(quanCir_list:list, gamma, qubit_map:dict, edge=None):
        qcircuit = QuantumCircuit(len(qubit_map.values()))

        for sub_cir in quanCir_list:
            if edge == None:
                if sum(sub_cir) > 1:
                    indices = [i for i, x in enumerate(sub_cir) if x == 1]
                    for i in range(len(indices)):
                        if i == len(indices) - 1:
                            qcircuit.rz(gamma, indices[i])
                        else:
                            qcircuit.cx(indices[i], indices[i + 1])

                    for i in range(len(indices)-1, 0, -1):
                        qcircuit.cx(indices[i-1], indices[i])
                else:
                    qcircuit.rz(gamma, sub_cir.index(1))
            else:
                if len(sub_cir) > 1:
                    if sum(sub_cir) > 1:
                        qcircuit.cx(edge[0], edge[1])
                        qcircuit.rz(gamma, edge[1])
                        qcircuit.cx(edge[0], edge[1])
                    else:
                        qcircuit.rz(gamma, edge[sub_cir.index(1)])
                else:
                    for v in qubit_map:
                        qcircuit.rz(gamma, int(qubit_map[v]))
        
        return qcircuit

    @staticmethod
    def generate_Ditcircuit(quanCir_list:list, gamma, qubit_map:dict, edge, sub_expr, k_dits, no_qubits):
        qcircuit = QuantumCircuit(no_qubits)

        for sub_cir in quanCir_list:
            for i in range(k_dits):
                    nxt_i = eval(sub_expr['i']) % k_dits

                    if sum(sub_cir) > 1:
                        qcircuit.cx(qubit_map[edge[0]][i], qubit_map[edge[1]][nxt_i])
                        qcircuit.rz(gamma, qubit_map[edge[1]][nxt_i])
                        qcircuit.cx(qubit_map[edge[0]][i], qubit_map[edge[1]][nxt_i])
                    else:
                        qcircuit.rz(gamma, edge[sub_cir.index(1)])

        return qcircuit

    @staticmethod
    def generate_RXcircuit(qubit_map:dict, gamma):
        value = list(qubit_map.values())[0]
        if isinstance(value, tuple):
            l = 0
            for v in qubit_map.values():
                l += len(v)
            cir = QuantumCircuit(l)

            for i in range(l):
                cir.rx(gamma, i)

            return cir
        else:
            l = len(qubit_map)

            cir = QuantumCircuit(l)

            for i in range(l):
                cir.rx(gamma, i)

            return cir

    # @staticmethod
    # def generate_Toffoli(q)