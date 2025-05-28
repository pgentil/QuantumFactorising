
# !pip install -U qiskit-aer
# !pip install -U pylatexenc
# !pip install -U qiskit-ibm-runtime
# !pip install -U qiskit
import sympy
import numpy as np
from numpy import pi
from qiskit import QuantumCircuit, QuantumRegister, transpile, ClassicalRegister
from qiskit_aer import Aer
from qiskit.visualization import plot_histogram
from qiskit.circuit.library import UnitaryGate
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2 as Sampler, IBMBackend
from factorising_functions import *

def orderFindingCircuit(N:int, a:int, backend, show_circuit: bool = False, transpileCirc: bool = False, optimization: int = 2) -> QuantumCircuit:
  n = int(np.ceil(np.log2(N)))
  estimation = QuantumRegister(2 * n, name = "estimation")
  state = QuantumRegister(n, name="state")
  creg = ClassicalRegister(2 * n, name="measurements")
  circ = QuantumCircuit(state, estimation, creg, name = f"Period Finding SHOR N{N} a{a}")

  circ.h(estimation[:])
  circ.x(state[0])
  for i in range(2*n):
    circ.append(AExpXModNControlled(a, i, N),  state[:] + [estimation[i]])
  circ.barrier();
  circ.append(iqft(2*n), estimation[:])
  circ.measure(estimation, creg)

  if (show_circuit):
    circ.draw(output="mpl", reverse_bits=True)
  if(transpileCirc):
    circ = transpile(circ, backend, optimization_level=optimization)
  return circ

def orderFinding(N:int, a:int, backend: IBMBackend, shot_treshold: int = 80, show_circuit: bool = False, show_results: bool = True, simulation: bool = True) -> int:
  print(f"N: {N}, a: {a}")
  n = int(np.ceil(np.log2(N)))
  estimation = QuantumRegister(2 * n, name = "estimation")
  state = QuantumRegister(n, name="state")
  creg = ClassicalRegister(2 * n, name="measurements")
  circ = QuantumCircuit(state, estimation, creg, name = f"Period Finding SHOR N{N} a{a}")

  circ.h(estimation[:])
  circ.x(state[0])
  for i in range(2*n):
    circ.append(AExpXModNControlled(a, i, N),  state[:] + [estimation[i]])
  circ.barrier();
  circ.append(qft(2*n).inverse(), estimation[:])
  circ.measure(estimation, creg)

  if (show_circuit):
    circ.draw(output="mpl", reverse_bits=True)

  if (simulation):
    backend = Aer.get_backend('aer_simulator')
    tc = transpile(circ, backend, optimization_level=3)
    job = backend.run(tc, shots=1000)
    result= job.result()
    counts = result.get_counts()
  else:
    sampler = Sampler(backend)
    tc = transpile(circ, backend, optimization_level=3)
    job = sampler.run([tc], shots=1000)
    print(f"Job ID: {job.job_id()}")
    result= job.result()
    counts = result[0].data.measurements.get_counts()
    metrics = job.metrics()
    print(f"Total execution time (from metrics): {metrics['usage']}")

  if (show_results):
        plot_histogram(counts, title = "Shor's period finding")


  possibilities = []
  for key in counts:
    if (counts[key] > shot_treshold):
      contFrac = continuedFractions(key, 15)
      retq = []
      convergent(contFrac, len(contFrac) - 1, retq, [])
      aux = [i for i in retq if i < N]
      aux = max(aux)
      possibilities.append(aux)
  order = max(possibilities)
  return order

import random
import time

def shorAlgorithm(N: int, excluded: list, backend: IBMBackend, simulation : bool = True) -> list:
  random.seed(time.time())
  factors = []
  found = False
  if not sympy.isprime(N): ##paso I
    order = 1
    while(not found):
      back = False
      a = random.randint(2, N - 1) ## paso II
      print(f"a: {a}")
      gcd = sympy.gcd(a, N)
      print(f"gcd: {gcd}")
      if (gcd > 1): ##paso III
        if (gcd in excluded):
          back = True
        else:
          factors.append(gcd)
          found = True
      if (not back and not found):
        order = orderFinding(N, a, backend, simulation=simulation, shot_treshold=30) ## paso IV
        print(f"order: {order}")
      while (not back and not found):
        if (order % 2 == 0):
          x = pow(a, order // 2, N) #paso VI
          x1 = (x + 1) % N
          x2 = (x - 1) % N
          print(f"x1: {x1}, x2: {x2}")
          if (x1 == 0):
            back = True #paso VIIa
          elif (x2 == 0):
            print("order / 2")
            order = order // 2 #paso VIIb
          else: #paso VIII
            p = sympy.gcd(x1, N)
            q = sympy.gcd(x2, N)
            if (p != 1 and p not in excluded):
              factors.append(p)
            if (q != 1 and q not in excluded):
              factors.append(q)
            if (len(factors) > 0):
              found = True
        else: #paso V
          back = True
  return factors
