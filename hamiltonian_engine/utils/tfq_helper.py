from hamiltonian import hamiltonian, phase_hamiltonian
import networkx as nx
import numpy as np
import cirq
import qiskit
import tensorflow as tf
import tensorflow_quantum as tfq
import sympy
import os
os.path.sys.path.append('../hamiltonian_engine/')


def _qubit_mapping(qubit):
    q = str(qubit).split("_")
    return cirq.GridQubit(0, q[len(q)-1])


def _get_qubitNo(qubit):
    qubit_arr = str(qubit).split("_")
    value = int(qubit_arr[len(qubit_arr)-1])
    return value


class tfq_helper:

    @staticmethod
    def to_tfq_circuit(hamil: type([hamiltonian, qiskit.QuantumCircuit])):
        if isinstance(hamil, hamiltonian):
            qasm_cir = hamil.get_quantumCircuit().qasm()
        elif isinstance(hamil, qiskit.QuantumCircuit):
            qasm_cir = hamil.qasm()

        from cirq.contrib.qasm_import import circuit_from_qasm

        cirq_cir = circuit_from_qasm(qasm_cir)
        cirq_cir_gridded = cirq.Circuit([cirq.ops.Moment(operation.transform_qubits(
            _qubit_mapping)for operation in moment.operations)for moment in cirq_cir.moments])

        if isinstance(hamil, phase_hamiltonian):
            expanded_PauliSum = hamil.get_pHamil()

            Z_list = list(expanded_PauliSum.free_symbols)
            Z_list.remove(sympy.abc.I)

            Z_list.sort(key=_get_qubitNo)

            _cir = cirq.GridQubit.rect(1, len(Z_list))

            cirq_zList = []

            for q in _cir:
                cirq_zList.append(cirq.Z(q))

            conversion_dict = dict(zip(Z_list, cirq_zList))
            conversion_dict[sympy.abc.I] = 1

            pauli_sum = 0
            for exp in expanded_PauliSum.args:
                pauli_exp = 1
                for var in exp.args:
                    if var.is_Number:
                        number = float(var)
                        if number < 0:
                            number = number * -1
                            pauli_exp *= number
                    else:
                        pauli_exp *= conversion_dict[var]

                pauli_sum += cirq.PauliString(pauli_exp)

            return (cirq_cir_gridded, pauli_sum)
        else:
            return(cirq_cir_gridded, None)

    @staticmethod
    def exponentiate_parameters(param: list, cirq_circuit: cirq.Circuit):
        qaoa_parameters = sympy.symbols(param)
        p = 0
        cir_exponent = cirq.Circuit()

        while p != len(qaoa_parameters):
            new_cir = cirq_circuit.copy()
            for i in range(len(new_cir.moments)):
                if new_cir.moments[i].operations[0].gate == cirq.rz(1):
                    new_cir.moments[i] = new_cir.moments[i].__pow__(
                        qaoa_parameters[p])
                elif new_cir.moments[i].operations[0].gate == cirq.rx(1):
                    new_cir.moments[i] = new_cir.moments[i].__pow__(
                        qaoa_parameters[p + 1])
            p += 2
            cir_exponent += new_cir

        return cir_exponent

