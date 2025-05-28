import sympy
import numpy as np
from numpy import pi
from qiskit import QuantumCircuit, QuantumRegister, transpile, ClassicalRegister
from qiskit_aer import Aer
from qiskit.visualization import plot_histogram
from qiskit.circuit.library import UnitaryGate
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2 as Sampler, IBMBackend

"""

### Functions, gates and subcircuits declaration

We are going to define some function that will allow us to implement Shor's Algorithm. 
We define: The Quantum Fourier Transform (QFT), the Approximate Quantum Fourier Transform (AQFT), the quantum function of $a^{2^x} \ mod \ N$,
the quantum controlled function of $a^{2^x} \ mod \ N$ plus two functions to do fast exponentiation (both normal and modular). 
We define a function to convert a binary number (string) to decimal integer, a function to compute the continued fraction expresion of a real number,
and a function that computes the convergents of a continued fraction given a limit k. This last one is partiicularly important as it does not return the actual convergents,
but 2 lists containing the numerator and denominator of each one of them. The information that we are interested in is the denominator, because by runnning Shor's algorithm
- specifically the order-finding circuit -  we are going to get a fraction, which we will express as a continued fraction, and then see all its convergents,
taking the denominator of the convergent with the biggest limit (k) such that the denominator is smaller than the number we are trying to factorize (N). 
This denominator will be the order/period we are looking for.
\
_**Remark:**_ A function has been declared to create the controlled version of the (a^{2^x} mod N) gate by using a unitary matrix. 
The .control() function form Qiskit takes too long to create the gate otherwise. In case of wanting to modify the implementation of either of both $a^{2^x} \ mod \ N$ gates
(controlled or normal) do not hesitate to do so. The same kind of problem applies to the .inverse() method in Qiskit,
 so to increase the efficiency of the circuit's construction it would be advisable to define manually the inverse QFT or AQFT.
"""

def numberToCircuit(n: int, nbits: int) -> QuantumCircuit:
  """
  Auxiliary function to hardcode a number to a register of nbits
  """
  circ = QuantumCircuit(nbits)
  if (nbits >= int(np.ceil(np.log2(n)))):
    binary = bin(n)
    binary = binary[2:len(binary)]
    binary = binary[::-1]
    for i in range(len(binary)):
      if binary[i] == '1':
        circ.x(i)
  return circ


def qft(n: int, swaps: bool = True) -> QuantumCircuit:
  circuit = QuantumCircuit(n, name="QFT")
  PIx2 = 2 * pi
  for i in range(n - 1, -1, -1):
    circuit.h(i)
    for j in range(0, i):
      circuit.cp((PIx2 / (2**(1 + i - j))), j, i)
  circuit.barrier()
  upper = n//2
  if (swaps):
    for i in range(upper):
      if (i != n - i - 1):
        circuit.swap(i, n - i - 1)
  return circuit

def qftC(n: int, swaps: bool = True, constant: int = 1) -> QuantumCircuit:
  circuit = QuantumCircuit(n, name="QFT")
  PIx2 = 2 * pi
  for i in range(n - 1, -1, -1):
    circuit.h(i)
    for j in range(0, i):
      circuit.cp((PIx2 / (2**(1 + i - j))), j, i)
    circuit.p((constant - 1)) * pi
  circuit.barrier()
  upper = n//2
  if (swaps):
    for i in range(upper):
      if (i != n - i - 1):
        circuit.swap(i, n - i - 1)
  return circuit

def iqft(n: int, swaps: bool = True) -> QuantumCircuit:
  circuit = QuantumCircuit(n, name="IQFT")
  upper = n//2
  if (swaps):
    for i in range(upper):
      if (i != n - i - 1):
        circuit.swap(i, n - i - 1)
  PIx2 = 2 * pi
  for i in range(0, n):
    for j in range(0, i):
      circuit.cp(-(PIx2 / (2**(1 + i - j))), j, i)
    circuit.h(i)
  circuit.barrier()


  return circuit

def aqft(n: int, max_rot: int, swaps: bool = True) -> QuantumCircuit:
  circuit = QuantumCircuit(n, name="AQFT")
  PIx2 = 2 * pi
  for i in range(n - 1, -1, -1):
    circuit.h(i)
    for j in range(i - 1, i - min(i, max_rot) - 1, -1): #hacemos el mínimo entre la cantidad de rotaciones que debe hacer el qbit y las rotaciones máximas que se pueden hacer
      circuit.cp((PIx2 / (2**(1 + i - j))), j, i)
  circuit.barrier()
  upper = n//2
  if (swaps):
    for i in range(upper):
      if (i != n - i - 1):
        circuit.swap(i, n - i - 1)
  return circuit

def fexp(val, power) -> int:
    result = pow(val, power//2)
    result = result * result

    if power % 2 != 0:
        result = result * val
    return result

def FastModularExponentiation(b, k, m):
    return pow(b, pow(2, k), m)

def put_one(row: int, col: int, ax: int, N: int) -> int:
  res = 0
  AUX = col * ax
  if (col >= N):
    if (row == col):
      res = 1
  elif (((AUX) % N) == row):
    res = 1
  return res

def AExpXModN(a: int, x: int, N: int) -> QuantumCircuit:
  n = int(np.ceil(np.log2(N)))
  circ = QuantumCircuit(n, name=f"{a}^(2^{x}) % {N}")
  ax = FastModularExponentiation(a, x, N)

  matrix = [[put_one(i, j, ax, N) for j in range(fexp(2, n))] for i in range (fexp(2, n))]
  gate = UnitaryGate(matrix)
  circ.append(gate, [i for i in range(n)])

  return circ

def printMatrix(matrix):
  for i in matrix:
    for j in i:
      print(j, end=" ")
    print("")

def AExpXModNControlled(a: int, x: int, N: int) -> QuantumCircuit:
  """
  Eng:
    Implementation of the controlled version of the quantum gate a^(2^x) mod N.
    The control bit will be the most significant bit of the circuit.

    Parameters:
      a - number to be powered
      x - number to be used for exponent 2^x
      N - number to be used for modulo operation as divisorç

    Returns: Quantum circuit that envelopes the behaviour of a controlled
    a^(2^x) mod N quantum gate.

  Esp:
    Implementación de la puerta a^(2^x) mod N controlada. Se debe
    tener en cuenta que el qbit de control será el más significativo,
    ya que hace que la creación de la matriz sea más fácil.

    Parametros:
      a - número a exponenciar
      x - exponente (en realidad el exponente será 2^x)
      N - módulo

    Salida: Circuito cúantico que contiene el comportamiento de la
    puerta a^(2^x) mod N controlada
  """
  n = int(np.ceil(np.log2(N)))
  circ = QuantumCircuit(n + 1, name=f"{a}^(2^{x}) % {N} CONTROLLED")
  ax = FastModularExponentiation(a, x, N)
  pnm1 = fexp(2, n)
  pn = fexp(2, n + 1)
  matrix = [
      [put_one(i - pnm1 + 1, j - pnm1 + 1, ax, N)
      if (i - pnm1 + 1 >= 0 and j - pnm1 + 1 >= 0)
      else (1 if i == j else 0) for j in range(pn)]
      for i in range(pn)]
  gate = UnitaryGate(matrix)
  circ.append(gate, [i for i in range(n + 1)])

  return circ

def AExpXModNControlled(a: int, x: int, N: int) -> QuantumCircuit:
  """
  Eng:
    Implementation of the controlled version of the quantum gate a^(2^x) mod N.
    The control bit will be the most significant bit of the circuit.

    Parameters:
      a - number to be powered
      x - number to be used for exponent 2^x
      N - number to be used for modulo operation as divisorç

    Returns: Quantum circuit that envelopes the behaviour of a controlled
    a^(2^x) mod N quantum gate.

  Esp:
    Implementación de la puerta a^(2^x) mod N controlada. Se debe
    tener en cuenta que el qbit de control será el más significativo,
    ya que hace que la creación de la matriz sea más fácil.

    Parametros:
      a - número a exponenciar
      x - exponente (en realidad el exponente será 2^x)
      N - módulo

    Salida: Circuito cúantico que contiene el comportamiento de la
    puerta a^(2^x) mod N controlada
  """
  n = int(np.ceil(np.log2(N)))
  circ = QuantumCircuit(n + 1, name=f"{a}^(2^{x}) % {N} CONTROLLED")
  ax = FastModularExponentiation(a, x, N)
  pnm1 = fexp(2, n)
  pn = fexp(2, n + 1)
  matrix = [
      [put_one(i - pnm1 + 1, j - pnm1 + 1, ax, N)
      if (i - pnm1 + 1 >= 0 and j - pnm1 + 1 >= 0)
      else (1 if i == j else 0) for j in range(pn)]
      for i in range(pn)]
  gate = UnitaryGate(matrix)
  circ.append(gate, [i for i in range(n + 1)])

  return circ

circ = AExpXModNControlled(2, 0, 15)
#display(circ.decompose().draw(output="mpl", reverse_bits=True))

def controlledAddAndScale2(na: int, c: int) -> QuantumCircuit:
  """
  Controlled add and scale operator:
  For |d> (control register of length 1), |a> in register A, |b> in register B and a constant C this operator will
  return |a + b*C>in register A and |b> in register B if |d> == |1>

  Input - control register -> Less significant bit
          A -> Less significant bits following the control register
          B -> Most significant bits

  Parameters:
    na - number of bits in A
    nb - number of bits in B
    c - constant C
  Control bit is set to be the LSB of the circuit
  """
  control = QuantumRegister(1, name="control")
  A = QuantumRegister(na, name="a")
  circ = QuantumCircuit(control, A, name="MUL")
  for i in range(0, na):
    circ.append(customRGate(1, c * (i+1)).control(1), [0, i])
  return circ

def fromBinToDecimal(binn: str) ->int:
  dec = int(binn, 2)
  dec = dec
  return dec


def continuedFractions(dec: str, limit: int) -> list:
  ret = []
  den = pow(2, len(dec))
  num = fromBinToDecimal(dec)
  ret.append(0)
  i = 0
  while num != 0 and i < limit:
    num, den = den, num
    aux = num // den
    num = num - aux * den
    ret.append(aux)
    i = i + 1
  return ret

#Precondición: len(contFrac) > 0 y len(contFrac) > k >= 0 y ret = [] y k <= len(contFrac)
def convergent(contFrac: list, k: int, retq: list, retp: list):
  p = -1
  q = -1
  if (k == 0):
    p = contFrac[0]
    q = 1
  elif (k == 1):
    convergent(contFrac, k - 1, retq, retp)
    p = contFrac[0] * contFrac[1] + 1
    q = contFrac[1]
  else:
    convergent(contFrac, k - 1, retq, retp)
    p = retp[k-1] * contFrac[k] + retp[k-2]
    q = retq[k-1] * contFrac[k] + retq[k-2]
  retp.append(p)
  retq.append(q)




"""
### Grover's approach algorithm's auxiliary functions
#### Functions
"""



def customRGate(i: int, c: int) -> PhaseGate:
  """
  Wrapper for a phase rotation gate using the PhaseGate from
  Qiskit. The custom R gate applies a phase shift to a qbit when
  it is |1> (e^(c*pow / 2^(i))). This gate will allow us to perform
  addition and multiplications using i and c as parameters.

  """
  angle = c*pi/pow(2, i)
  gate = PhaseGate(angle, label=f"R{i}({c})")
  return gate


def customCRGate(i: int, c:int) ->CPhaseGate:
  """
  Wrapper for a controlled phase rotation gate using the CPhaseGate from
  Qiskit. The custom R gate applies a phase shift to a qbit when
  it is |1> (e^(c*pow / 2^(i))) and control bit is set to |1>.
  This gate will allow us to perform addition and multiplications using i and
  c as parameters.
  """
  angle = c*pi/pow(2, i)
  gate = CPhaseGate(angle, label=f"CR{i}({c})")
  return gate


def addAndScale(na: int, nb: int, c: int) -> QuantumCircuit:
  """
  Add and scale operator:
  For |a> in register A, |b> in register B, and a constant C this operator will
  return |a + b*C>in register A and |b> in register B

  Input - A -> Less significant bits
          B -> Most significant bits

  Parameters:
    na - number of bits in A
    nb - number of bits in B
    c - constant C
  """
  A = QuantumRegister(na, name="a")
  B = QuantumRegister(nb, name="b")
  circ = QuantumCircuit(A, B, name="ADD")
  if (na >= nb):
    for i in range(nb):
      for j in range(i, na):
        circ.append(customCRGate(j - i, c), [na + i, j])
  return circ


def controlledAddAndScale(na: int, nb:int, c: int) -> QuantumCircuit:
  """
  Controlled add and scale operator:
  For |d> (control register of length 1), |a> in register A, |b> in register B and a constant C this operator will
  return |a + b*C>in register A and |b> in register B if |d> == |1>

  Input - control register -> Less significant bit
          A -> Less significant bits following the control register
          B -> Most significant bits

  Parameters:
    na - number of bits in A
    nb - number of bits in B
    c - constant C
  Control bit is set to be the LSB of the circuit
  """
  control = QuantumRegister(1, name="control")
  A = QuantumRegister(na, name="a")
  B = QuantumRegister(nb, name="b")
  circ = QuantumCircuit(control, A, B, name="CADD")
  if (na >= nb):
    for i in range(nb):
      for j in range(i, na):
        circ.append(customRGate(j - i, c).control(2), [0, na + i + 1, j + 1])
  return circ


def phiQMA(nx: int, ny: int, nz: int, a:int, b:int, c: int, d: int) -> QuantumCircuit:
  """
  QUANTUM MULTIPLY ADD
  This circuit computes the calculation z + axy + bx + cy + d mod 2^nz
  and stores the result in register Z.
  """
  X = QuantumRegister(nx, name="x")
  Y = QuantumRegister(ny, name="y")
  Z = QuantumRegister(nz, name="z")
  circ = QuantumCircuit(X, Y, Z, name="phiQMA")
  for i in range(nx):
    circ.append(controlledAddAndScale(nz, ny, a*pow(2, i)), [i] +  Z[:] + Y[:]) # controlledAddAndScale
  circ.append(addAndScale(nz, nx, b), Z[:] + X[:])
  circ.append(addAndScale(nz, ny, c), Z[:] + Y[:])
  for i in range(nz):
    circ.append(customRGate(i, d), [nx + ny + i])
  return circ

def diff(n : int) -> QuantumCircuit:
  """
    Diffuser operator for Grover's search algorithm.
    In this case we will use it for the result registers and the
    "auxiliary" register (Z) as the oracle can be implemented with this
    operator.
  """
  circ = QuantumCircuit(n)
  for i in range(n):
    circ.ry(pi/2, i)
  circ.append(ZGate().control(n - 1), [i for i in range(1, n)] + [0])
  for i in range(n):
    circ.ry(-pi/2, i)
  return circ

def diffZ(nz: int) -> QuantumCircuit:
  """
  Other way (less eficient) of implementing the oracle for the factorising algorithm.
  """
  circ = QuantumCircuit(nz)
  circ.append(qft(nz, swaps=False).inverse(), [i for i in range(nz)])
  for i in range(nz):
    circ.x(i)
  circ.append(ZGate().control(nz - 1), [i for i in range(1, nz)] + [0])
  for i in range(nz):
    circ.x(i)
  circ.append(qft(nz, swaps=False), [i for i in range(nz)])
  return circ
