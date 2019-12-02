# Blueqat Gottesman-Knill backend
Blueqat backend for Clifford gates (X, Y, Z, H, S, S†, CNOT, CZ) and measurement

## Usage

```py
from blueqat import Circuit
import blueqat_gk_backend

print(Circuit().x[0].cx[0, 1].ccx[0, 1, 2].m[:].run_with_gk())
```
