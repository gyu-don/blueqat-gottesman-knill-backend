"""The Blueqat backend by Gottesman-Knill Theorem.

    This module provides Blueqat backend for X, Y, Z, S, S†, CX gates and measurement.
    "gk" backend is faster and less memory than general quantum computer simulator.
    However, available gates are not enough to make quantum algorithm.
"""

from collections import Counter
from dataclasses import dataclass
from enum import Enum
import random
from typing import List, Tuple
import numpy as np
from blueqat import BlueqatGlobalSetting
from blueqat.backends.backendbase import Backend

gk_table = {
        # (gate, pauli): (next pauli, swap sign? (+1: keep, -1: swap))
        ('x', 'X'): ('X', +1),
        ('x', 'Y'): ('Y', -1),
        ('x', 'Z'): ('Z', -1),
        ('y', 'X'): ('X', -1),
        ('y', 'Y'): ('Y', +1),
        ('y', 'Z'): ('Z', -1),
        ('z', 'X'): ('X', -1),
        ('z', 'Y'): ('Y', -1),
        ('z', 'Z'): ('Z', +1),
        ('h', 'X'): ('Z', +1),
        ('h', 'Y'): ('Y', -1),
        ('h', 'Z'): ('X', +1),
        ('s', 'X'): ('Y', +1),
        ('s', 'Y'): ('X', -1),
        ('s', 'Z'): ('Z', +1),
        ('sdg', 'X'): ('Y', -1),
        ('sdg', 'Y'): ('X', +1),
        ('sdg', 'Z'): ('Z', +1),
}

_mul = {
    'I': {
        'I': ('I', 0),
        'X': ('X', 0),
        'Y': ('Y', 0),
        'Z': ('Z', 0),
    },
    'X': {
        'I': ('X', 0),
        'X': ('I', 0),
        'Y': ('Z', -1),
        'Z': ('Y', +1),
    },
    'Y': {
        'I': ('Y', 0),
        'X': ('Z', +1),
        'Y': ('I', 0),
        'Z': ('X', -1),
    },
    'Z': {
        'I': ('Z', 0),
        'X': ('Y', -1),
        'Y': ('X', +1),
        'Z': ('I', 0),
    },
}

@dataclass
class Stabilizer:
    pauli: List[str]
    sign: int

    def __mul_pauli_and_phase(self, other):
        pauli = [_mul[lc][rc][0] for lc, rc in zip(self.pauli, other.pauli)]
        sign = self.sign * other.sign
        phase = sum(_mul[lc][rc][1] for lc, rc in zip(self.pauli, other.pauli))
        assert phase in (-2, 0, 2)
        if phase != 0:
            sign = -sign
        return pauli, sign

    def __mul__(self, other: 'Stabilizer') -> 'Stabilizer':
        return Stabilizer(*self.__mul_pauli_and_phase(other))

    def __imul__(self, other: 'Stabilizer'):
        pauli, sign = self.__mul_pauli_and_phase(other)
        self.pauli = pauli
        self.sign = sign


class _GKReturnType(Enum):
    SHOTS = 'shots'
    STABILIZERS = 'stabilizers'


@dataclass
class _GottesmanKnillBackendContext:
    stabilizers: List[Stabilizer]
    measured: List[bool]
    returns: _GKReturnType


class GottesmanKnillBackend(Backend):
    """Backend by Gottesman-Knill Theorem

    This backend supports only Clifford gates and measurement.
    """
    @staticmethod
    def _make_z_stabilizer(n, n_qubits):
        s = []
        for _ in range(n):
            s.append('I')
        s.append('Z')
        for _ in range(n_qubits - n - 1):
            s.append('I')
        return Stabilizer(s, 1)

    def _preprocess_run(self, gates, n_qubits, args, kwargs):
        stabilizers = [self._make_z_stabilizer(i, n_qubits) for i in range(n_qubits)]
        rettype = kwargs['returns']
        return gates, _GottesmanKnillBackendContext(stabilizers, [False] * n_qubits, rettype)

    def _postprocess_run(self, ctx):
        if ctx.returns == _GKReturnType.SHOTS:
            return ''.join('1' if b else '0' for b in ctx.measured)
        return ctx.stabilizers

    def run(self, gates, n_qubits, *args, **kwargs):
        returns = kwargs.get('returns', 'shots')

        if returns == 'shots':
            rettype = _GKReturnType.SHOTS
            shots = kwargs.get('shots', 100)
        elif returns == 'stabilizers':
            rettype = _GKReturnType.STABILIZERS
            shots = kwargs.get('shots', 1)
        else:
            raise ValueError('Unknown return type.')
        kwargs['returns'] = rettype

        counter = Counter()
        stabilizers = []
        for _ in range(shots):
            result = self._run(gates, n_qubits, args, kwargs)
            if rettype == _GKReturnType.SHOTS:
                counter[result] += 1
            else:
                stabilizers.append(result)

        if rettype == _GKReturnType.SHOTS:
            return counter
        else:
            return stabilizers


    def _single_qubit_gate(self, gate, ctx):
        sts = ctx.stabilizers
        for idx in gate.target_iter(len(sts)):
            for st in sts:
                if st.pauli[idx] != 'I':
                    pauli, sign = gk_table[gate.lowername, st.pauli[idx]]
                    st.pauli[idx] = pauli
                    st.sign *= sign
        return ctx

    gate_x = _single_qubit_gate
    gate_y = _single_qubit_gate
    gate_z = _single_qubit_gate
    gate_h = _single_qubit_gate
    gate_s = _single_qubit_gate
    gate_sdg = _single_qubit_gate

    def gate_cx(self, gate, ctx):
        sts = ctx.stabilizers
        for (c, t) in gate.control_target_iter(len(sts)):
            for st in sts:
                phase = 0
                cpauli = st.pauli[c]
                tpauli = st.pauli[t]
                if cpauli in 'XY':
                    st.pauli[t], dphase = _mul['X'][tpauli]
                    phase += dphase
                if tpauli in 'YZ':
                    st.pauli[c], dphase = _mul[cpauli]['Z']
                    phase += dphase
                assert phase in (-2, 0, 2)
                if phase != 0:
                    st.sign = -st.sign
        return ctx


    @staticmethod
    def _meas_all_comm(idx, stabilizers):
        n_qubits = len(stabilizers)
        op = GottesmanKnillBackend._make_z_stabilizer(idx, n_qubits).pauli

        def check(sts):
            try:
                idx = [st.pauli for st in sts].index(op)
            except ValueError:
                return None
            return sts[idx]

        def trans(sts, xs):
            if xs:
                for st in xs[1:]:
                    st *= xs[0]
                sts.remove(xs[0])

        chk = check(stabilizers)
        if chk:
            return chk.sign == -1

        sts = stabilizers[:]
        for i in range(n_qubits):
            if i == idx:
                continue
            xs = [st for st in sts if st.pauli[i] == 'X']
            trans(sts, xs)
            ys = [st for st in sts if st.pauli[i] == 'Y']
            trans(sts, ys)
            zs = [st for st in sts if st.pauli[i] == 'Z']
            trans(sts, zs)
        # Now, sts is I..IZI..I shaped.
        chk = check(sts)
        if not chk:
            raise ValueError(f'Backend is buggy. idx: {idx}, stabilizers: {stabilizers}')
        return chk.sign == -1

    @staticmethod
    def _meas_partial_noncomm(idx, noncomms):
        n_qubits = len(noncomms[0].pauli)
        for st in noncomms[1:]:
            st *= noncomms[0]
            
        noncomms[0].pauli = GottesmanKnillBackend._make_z_stabilizer(idx, n_qubits).pauli
        r = random.random()
        if r < 0.5:
            noncomms[0].sign = +1
            return False
        else:
            noncomms[0].sign = -1
            return True

    def gate_measure(self, gate, ctx):
        sts = ctx.stabilizers
        n_qubits = len(ctx.stabilizers)

        for idx in gate.target_iter(len(sts)):
            noncommutatives = [st for st in sts if st.pauli[idx] in 'XY']
            if noncommutatives:
                meas = self._meas_partial_noncomm(idx, noncommutatives)
            else:
                meas = self._meas_all_comm(idx, sts)
            ctx.measured[idx] = meas
        return ctx


BlueqatGlobalSetting.register_backend('gk', GottesmanKnillBackend)
