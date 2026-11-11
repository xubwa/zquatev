import numpy as np
from zquatev import eigh

def test_eigh():
    np.random.seed(2)
    n = 6
    a, b = np.random.rand(2,n,n) + 1j*np.random.rand(2,n,n)
    a = a + a.conj().T
    b = b - b.T
    mat = np.empty((2,2,n,n), dtype=np.complex128)
    mat[0,0] = a
    mat[1,1] = a.conj()
    mat[0,1] = b
    mat[1,0] = -b.conj()
    mat = mat.transpose(0,2,1,3).reshape(n*2,n*2)

    e1, v1 = np.linalg.eigh(mat)

    e, v = eigh(mat, iop=1)
    assert abs(e1 - e).max() < 1e-14
    assert abs(mat.dot(v) - v*e).max() < 1e-14

    e, v = eigh(mat)
    assert abs(mat.dot(v) - v*e).max() < 1e-14


def test_eigh_blocked():
    # n = 300 > the panel size (nb = 20 in zquatev.cc): the blocked
    # Householder path, which n = 6 does not reach
    rng = np.random.default_rng(5)
    n = 300
    a, b = rng.random((2,n,n)) + 1j*rng.random((2,n,n))
    a = a + a.conj().T
    b = b - b.T
    mat = np.block([[a, b], [-b.conj(), a.conj()]])
    e1 = np.linalg.eigvalsh(mat)
    e, v = eigh(mat)
    scale = abs(e1).max()
    # Kramers pairs: the partners of the first n eigenpairs are the last n
    assert abs(e[:n] - e[n:]).max() < 1e-12 * scale
    assert abs(np.sort(e) - e1).max() < 1e-12 * scale
    assert abs(mat.dot(v) - v*e).max() < 1e-11 * scale
    assert abs(v.conj().T.dot(v) - np.eye(2*n)).max() < 1e-12
