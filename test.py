# import the predefined library
import os 
os.path.abspath(os.curdir)
os.path.sys.path.append('predefined_problems/')
import skeletor as sk
import numpy as np
# Create a graph for the problem its is also adviseable to instantiate the weight of the graphs else it 
# be defaulted to 1
import networkx as nx 
import numpy as np
import matplotlib.pyplot as plt 
from   matplotlib import cm
from   matplotlib.ticker import LinearLocator, FormatStrFormatter

n     = 5
V     = np.arange(0,n,1)
E     =[(0,1,6.0),(2,0,7.0),(1,2,8.0),(3,2,9.0),(3,4,10.0),(4,2,4.0)] 

G     = nx.Graph()
G.add_nodes_from(V)
G.add_weighted_edges_from(E)

print("Finished generating  graph.......")

obj1 = '~x_1 & x_2 | ~x_2 & x_1'
v1 = ['x_1', 'x_2']

hamiltonian = sk.skeletor(60, obj1, v1, True, graph=G)
print("Instantiated Skeletor .......")

hyperparams = [np.pi, np.pi, np.pi, np.pi/2, np.pi/3,
      np.pi, np.pi, np.pi, np.pi/2, np.pi/3,
       np.pi, np.pi, np.pi, np.pi/2, np.pi/3,
      np.pi, np.pi, np.pi, np.pi/2, np.pi/3,
      np.pi, np.pi, np.pi, np.pi/2, np.pi/3,
       np.pi, np.pi, np.pi, np.pi/2, np.pi/3,
       np.pi, np.pi, np.pi, np.pi/2, np.pi/3,
       np.pi, np.pi, np.pi, np.pi/2, np.pi/3,
      np.pi, np.pi, np.pi, np.pi/2, np.pi/3,
       np.pi, np.pi, np.pi, np.pi/2, np.pi/3,
       np.pi, np.pi, np.pi, np.pi/2, np.pi/3,
       np.pi, np.pi, np.pi, np.pi/2, np.pi/3,
       np.pi, np.pi, np.pi, np.pi/2, np.pi/3,
       np.pi, np.pi, np.pi, np.pi/2, np.pi/3,
       np.pi, np.pi, np.pi, np.pi/2, np.pi/3,
       np.pi, np.pi, np.pi, np.pi/2, np.pi/3,
       np.pi, np.pi, np.pi, np.pi/2, np.pi/3,
      np.pi, np.pi, np.pi, np.pi/2, np.pi/3,
       np.pi, np.pi, np.pi, np.pi/2, np.pi/3,
       np.pi, np.pi, np.pi, np.pi/2, np.pi/3,
       np.pi, np.pi, np.pi, np.pi/2, np.pi/3,
       np.pi, np.pi, np.pi, np.pi/2, np.pi/3,
       np.pi, np.pi, np.pi, np.pi/2, np.pi/3,
       np.pi, np.pi, np.pi, np.pi/2, np.pi/3]

hamiltonian.generate_quantumCircuit(hyperparams=hyperparams)
print("finshied circuit map.......")
opt_param = hamiltonian.run_QAOA(hyperparams,'COBYLA')
print("finshied obtaining hyperparams....... running circuit once more......")
hamiltonian.generate_quantumCircuit(opt_param)
res = hamiltonian.run_circuit()

print(res.get_counts())
