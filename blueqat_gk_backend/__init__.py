"""The Blueqat backend by Gottesman-Knill Theorem.

    This module provides Blueqat backend for X, Y, Z, S, S†, CX gates and measurement.
    "gk" backend is faster and less memory than general quantum computer simulator.
    However, available gates are not enough to make quantum algorithm.
"""

from collections import Counter
from dataclasses import dataclass
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
        ('s', 'Y'): ('Z', -1),
        ('s', 'Z'): ('Z', +1),
        ('sdg', 'X'): ('Y', -1),
        ('sdg', 'Y'): ('Z', +1),
        ('sdg', 'Z'): ('Z', +1),
}

gk_mult_x = {
        'I': 'X',
        'X': 'I',
        'Y': 'Z',
        'Z': 'Y',
}

gk_mult_y = {
        'I': 'Y',
        'X': 'Z',
        'Y': 'I',
        'Z': 'X',
}

gk_mult_z = {
        'I': 'Z',
        'X': 'Y',
        'Y': 'X',
        'Z': 'I',
}

gk_mult_i = {
        'I': 'I',
        'X': 'X',
        'Y': 'Y',
        'Z': 'Z',
}

gk_mult_dict = {
        'I': gk_mult_i,
        'X': gk_mult_x,
        'Y': gk_mult_y,
        'Z': gk_mult_z,
}

@dataclass
class _Stabilizer:
    pauli: List[str]
    sign: int

@dataclass
class _GottesmanKnillBackendContext:
    stabilizers: List[_Stabilizer]
    measured: List[bool]

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
        return _Stabilizer(s, 1)

    def _preprocess_run(self, gates, n_qubits, args, kwargs):
        stabilizers = [self._make_z_stabilizer(i, n_qubits) for i in range(n_qubits)]
        return gates, _GottesmanKnillBackendContext(stabilizers, [False] * n_qubits)

    def _postprocess_run(self, ctx):
        return ''.join('1' if b else '0' for b in ctx.measured)

    def run(self, gates, n_qubits, *args, **kwargs):
        shots = kwargs.get('shots', 100)
        counter = Counter()
        for i in range(shots):
            result = self._run(gates, n_qubits, args, kwargs)
            counter[result] += 1
        return counter

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
                if st.pauli[c] in 'XY':
                    st.pauli[t] = gk_mult_x[st.pauli[t]]
                if st.pauli[t] in 'YZ':
                    st.pauli[c] = gk_mult_z[st.pauli[c]]
        return ctx

    @staticmethod
    def _mult_pauli(lhs, rhs):
        return[gk_mult_dict[l][r] for l, r in zip(lhs, rhs)]

    @staticmethod
    def _meas_all_comm(idx, stabilizers):
        n_qubits = len(stabilizers)
        op = GottesmanKnillBackend._make_z_stabilizer(idx, n_qubits).pauli
        mul = GottesmanKnillBackend._mult_pauli

        def check(sts):
            try:
                idx = [st.pauli for st in sts].index(op)
            except ValueError:
                return None
            return sts[idx]

        def trans(sts, xs):
            if xs:
                for st in xs[1:]:
                    st.pauli = mul(st.pauli, xs[0].pauli)
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
            st.pauli[idx] = gk_mult_z[st.pauli[idx]]
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
