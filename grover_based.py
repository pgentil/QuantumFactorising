# -*- coding: utf-8 -*-
"""## Factorising with Grover's search algorithm

This is the article where this algorithm was taken from
ref:
* Whitlock, S., & Kieu, T. D. (2023). Quantum Factoring Algorithm using Grover Search. arXiv preprint arXiv:2312.10054 [physics.gen-ph]
https://arxiv.org/pdf/2312.10054

We are going to use Grover's search to factorize a biprime N = pq; such that p, q are prime and 3 < p < N, 3 < q < N

"""
import sympy
import numpy as np
from numpy import pi
from qiskit.circuit.library import UnitaryGate, PhaseGate, CPhaseGate, ZGate
from qiskit.quantum_info import Statevector
from qiskit import QuantumCircuit, QuantumRegister, transpile, ClassicalRegister
from qiskit_aer import Aer
from qiskit.visualization import plot_histogram
from qiskit_ibm_runtime import QiskitRuntimeService
from factorising_functions import *

"""### Algorithm's Implementation"""
def groverQuantumCircuit(N:int, s: int, backend, showCircuit: bool = False, transpileCirc: bool = False, optimization: int = 2, dis: int=0) -> QuantumCircuit:
  ## number of bits necessary to store M
  bS = int(1.5 - 0.5 * (N % 6))
  M = int((N - bS)/6 - 1)

  nn = int(np.ceil(np.log2(N)))
  n = int(np.ceil(np.log2(M)))
  print(f"n (M bits): {n}, nn (N bits): {nn}")
  nxaux = int(np.floor(nn/2 - 2 - dis))
  nyaux = int(np.ceil(nn/2 - 2 + dis))
  nx = max(nxaux, 1)
  ny = max(nyaux, 1)
  print(f"nx: {nx}, ny: {ny}")
  nz = nx + ny + 3
  K = int(np.floor((pi/4)*pow(2, (nx + ny)/2)))

  ## initialize registers and circuit
  X = QuantumRegister(nx, name="x")
  Y = QuantumRegister(ny, name="y")
  Z = QuantumRegister(nz, name="z")
  cXY = ClassicalRegister(nx + ny, name="cXY")
  circ = QuantumCircuit(X, Y, Z, cXY, name=f"FactoriseGroverN{N}")

  ## initialize X and Y to uniform state
  for i in range(nx + ny):
    circ.h(i)

  circ.append(numberToCircuit(M, nz), Z[:])
  ## starting qft
  circ.append(qft(nz, swaps=False), Z[:])
  for i in range(K):
    ## ORACLE
    circ.append(phiQMA(nx, ny, nz, -6, -s*bS, - s, 1), X[:] + Y[:] + Z[:])
    circ.append(diff(nz), Z[:])
    circ.append(phiQMA(nx, ny, nz, 6, s*bS, s, -1), X[:] + Y[:] + Z[:])
    ## XY DIFFUSER
    circ.append(diff(nx + ny), X[:] + Y[:])
  ## ending qft inverse
  circ.append(qft(nz, swaps=False).inverse(), Z[:])
  circ.measure(X[:] + Y[:], cXY[:])
  if (showCircuit):
    circ.draw(output="mpl", reverse_bits = True)
  if (transpileCirc):
    circ = transpile(circ, backend, optimization_level=optimization)
  return circ

"""If we have some hints about p and q we can take this circuit even more efficient by reducing the number of qbits we are going to need."""

def factoriseQuantumCircuit(M: int, dis: int, bS: int, s: int, N:int, backend, simulation: bool = True) -> list:
  """
    Quantum part of the Factorising algorithm using Grover's search.
    We establish a relation between M, x and y such that choosing the correct x
    and y makes the NQMA (QMA inverse) operator plus the inverse QFT make the Z
    register go to |0>, state that is then marked by the oracle. Then we revert
    the operation and apply the diffuser operator to the X and Y registers to
    amplify the probability of choosing the correct x and y integer we are looking
    for.
    WARNING - This circuit does not give the factors of N, it just outputs an x
    and a y that can then be used to compute those factors.

    It is important to know that:
      bS  = 1.5 - 0.5*(N mod 6)
      M   = (N - bS)/6 - 1
      dis = distance between numbers of bits for X register and number of bits
            for Y register (in most cases it should be 0)
      N = number to be factorized
      s = +/- 1
  """
  assert(abs(bS) == 1 and abs(s) == 1)
  ## number of bits necessary to store M
  nn = int(np.ceil(np.log2(N)))
  n = int(np.ceil(np.log2(M)))
  print(f"n (M bits): {n}, nn (N bits): {nn}")
  nxaux = int(np.floor(nn/2 - 2 - dis))
  nyaux = int(np.ceil(nn/2 - 2 + dis))
  nx = max(nxaux, 1)
  ny = max(nyaux, 1)
  print(f"nx: {nx}, ny: {ny}")
  nz = nx + ny + 3
  K = int(np.floor((pi/4)*pow(2, (nx + ny)/2)))

  ## initialize registers and circuit
  X = QuantumRegister(nx, name="x")
  Y = QuantumRegister(ny, name="y")
  Z = QuantumRegister(nz, name="z")
  cXY = ClassicalRegister(nx + ny, name="cXY")
  circ = QuantumCircuit(X, Y, Z, cXY, name=f"FactoriseGroverN{N}")

  ## initialize X and Y to uniform state
  for i in range(nx + ny):
    circ.h(i)

  circ.append(numberToCircuit(M, nz), Z[:])
  ## starting qft
  circ.append(qft(nz, swaps=False), Z[:])
  for i in range(K):
    ## ORACLE
    circ.append(phiQMA(nx, ny, nz, -6, -s*bS, - s, 1), X[:] + Y[:] + Z[:])
    circ.append(diff(nz), Z[:])
    circ.append(phiQMA(nx, ny, nz, 6, s*bS, s, -1), X[:] + Y[:] + Z[:])
    ## XY DIFFUSER
    circ.append(diff(nx + ny), X[:] + Y[:])
  ## ending qft inverse
  circ.append(qft(nz, swaps=False).inverse(), Z[:])
  circ.measure(X[:] + Y[:], cXY[:])


  ## We run the circuit
  if (simulation):
    backend = Aer.get_backend("aer_simulator")
    tc = transpile(circ, backend, optimization_level=3)
    job = backend.run(tc, shots=1000)
    # metrics = job.metrics()
    # print(f"Total execution time (from metrics): {metrics['usage']}")
    result = job.result()
    counts = result.get_counts()
  else:
    sampler = Sampler(backend)
    tc = transpile(circ, backend, optimization_level=3)
    job = sampler.run([tc], shots=1000)
    print(f"Job ID: {job.job_id()}")
    result= job.result()
    counts = result[0].data.cXY.get_counts()
    metrics = job.metrics()
    print(f"Total execution time (from metrics): {metrics['usage']}")
  plot_histogram(counts)
  maxKey = ""
  maxValue = -1
  for i in counts:
    if counts[i] > maxValue:
      maxKey = i
      maxValue = counts[i]
  print(maxKey)
  resY = maxKey[0:ny]
  resX = maxKey[ny: ny + nx]
  print(f"bin x: {resX},bin y: {resY}")
  if (len(resY) > 0):
    resY = int(resY, 2)
  else: resY = 0
  if (len(resX) > 0):
    resX = int(resX, 2)
  else: resX = 0
  print(f"x: {resX}, y: {resY}")

  return [resX, resY]

def factoriseWithGrover(N: int, backend, dis: int = 0, simulation: bool = True) -> list:
  """
  Only works with biprime numbers N = pq, s.t. p and q are primes other than 2
  or 3. The extended version of this circuit is yet to be implemented.
  Distance parameter should be let at 0. The optimized version of this circuit is
  yet to be implemented.

  Parameters:
    N - biprime number to be factorised
    dis - distance between X and Y registers' number of bits, if it is known that
    N 's factors do not have the same bitlength.
  """
  fact = []
  ## We check that N is not a multiple of 2 or 3
  if N % 2 == 0:
    fact.append(2)
    N = N // 2
  if N % 3 == 0:
    fact.append(3)
    N = N // 3
  ## ----------
  if (len(fact) == 0):
    end = False
    ## We initialize our arguments for the algorithm
    S = int(1.5 - 0.5*(N % 6))
    print(f"S = {S}")
    M = (N - S)//6 - 1
    print(f"M = {M}")
    s = 1
    print(f"s = {s}")
    while(not end):
      ## We enter the quantum part of the algorithm
      aux = factoriseQuantumCircuit(M, dis, S, s, N, backend, simulation=simulation)

      ## We construct the prime numbers
      p = (aux[0]) * 6 + s
      q = (aux[1]) * 6 + s*S
      print(f"p: {p}, q: {q}")
      ## Check if they are correct
      if (p*q != N and s == 1):
        s = -1 ## We try with the other value for s
        print("s = -1")
      elif (p*q == N): ## We have our primes
        fact.append(p)
        fact.append(q)
        end = True
      else:
        end = True
    pass
  return fact

