import random
from blueqat import Circuit, BlueqatGlobalSetting
import blueqat_gk_backend

print('GHZ State')
print('Circuit().h[0].cx[0, 1].cx[0, 2].m[:].run(shots=100)')
print('numpy')
print(Circuit().h[0].cx[0, 1].cx[0, 2].m[:].run(shots=100))
print('Gottesman-Knill')
print(Circuit().h[0].cx[0, 1].cx[0, 2].m[:].run_with_gk(shots=100))

print('2 qubit grover (Expected: |11>)')
print('Circuit().h[:].h[1].cx[0, 1].h[1].h[:].x[:].h[1].cx[0, 1].h[1].x[:].h[:].m[:].run(shots=100)')
print('numpy')
print(Circuit().h[:].h[1].cx[0, 1].h[1].h[:].x[:].h[1].cx[0, 1].h[1].x[:].h[:].m[:].run(shots=100))
print('Gottesman-Knill')
print(Circuit().h[:].h[1].cx[0, 1].h[1].h[:].x[:].h[1].cx[0, 1].h[1].x[:].h[:].m[:].run_with_gk(shots=100))

def make_circuit(n):
    s = 'Circuit().h[:]'
    c = Circuit().h[:]
    for i in range(n - 1):
        c.cx[i, i + 1]
        s += f'.cx[{i}, {i + 1}]'
    c.y[:]
    s += f'.y[:]'
    for i in range(n - 1):
        c.cx[i + 1, i]
        s += f'.cx[{i + 1}, {i}]'
    c.h[:].m[:]
    s += '.h[:].m[:]'
    return c, s

print('Circuit (4 qubits)')
c, s = make_circuit(4)
print(s)
print('numpy')
print(c.run(shots=100))
print('Gottesman-Knill')
print(c.run_with_gk(shots=100))

print('Circuit (10 qubits)')
c, s = make_circuit(10)
print(s)
print('numpy')
print(c.run(shots=100))
print('Gottesman-Knill')
print(c.run_with_gk(shots=100))

print('Circuit (100 qubits)')
c, s = make_circuit(100)
print(s)
print('Gottesman-Knill')
print(c.run_with_gk(shots=100))
