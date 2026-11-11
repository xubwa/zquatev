# zquatev

Python interface to ZQUATEV by Toru Shiozaki
([qsimulate-open/zquatev](https://github.com/qsimulate-open/zquatev)), an
eigensolver for Hermitian quaternionic (Kramers-symmetric) matrices.  It
diagonalizes

    H = [ A  -B* ]      A Hermitian, B antisymmetric (B^T = -B),
        [ B   A* ]

with a quaternionic Householder reduction and returns eigenvectors in
time-reversal pairs: if (u, v) is an eigenvector, so is (-v*, u*), with the
same eigenvalue.  Such matrices appear, e.g., as Fock matrices of
Kramers-restricted relativistic (two- and four-component) SCF in a
Kramers-paired basis.

Reference: T. Shiozaki, Mol. Phys. 115, 5-12 (2017),
DOI: [10.1080/00268976.2016.1158423](https://doi.org/10.1080/00268976.2016.1158423),
[arXiv:1512.08934](https://arxiv.org/abs/1512.08934).

## Install

    pip install zquatev

Wheels are provided for Linux (x86_64, aarch64, with OpenBLAS) and macOS
(x86_64, arm64, with Accelerate).  The library is loaded with ctypes, so one
wheel serves all Python 3 versions.

### Building from source

Requirements: a C++ compiler, CMake (>= 3.10) and BLAS/LAPACK.

    pip install --no-binary zquatev zquatev                 # from the sdist on PyPI
    pip install git+https://github.com/sunqm/zquatev        # from the repository

BLAS/LAPACK are found by CMake (FindBLAS/FindLAPACK).  Options are passed to
CMake through the environment variable `CMAKE_ARGS`, e.g.

    CMAKE_ARGS="-DBLA_VENDOR=OpenBLAS" pip install --no-binary zquatev zquatev
    CMAKE_ARGS="-DBLA_VENDOR=Intel10_64lp_seq" pip install ...   # MKL, sequential
    CMAKE_ARGS="-DBLAS_LIBRARIES=/path/to/libmkl_rt.so" pip install ...

On macOS the default is Apple Accelerate (`-DBLA_VENDOR=Apple`).

The shared library alone (for use from C/C++ or other ctypes code) is built
with CMake:

    cmake -B build -DCMAKE_INSTALL_PREFIX=/path/to/prefix
    cmake --build build && cmake --install build      # -> prefix/lib/libzquatev.so

## Usage

    import numpy as np
    from zquatev import eigh, geigh

    n = 100
    a = np.random.rand(n, n) + 1j*np.random.rand(n, n)
    b = np.random.rand(n, n) + 1j*np.random.rand(n, n)
    a = a + a.conj().T
    b = b - b.T
    h = np.block([[a, -b.conj()], [b, a.conj()]])

    e, c = eigh(h)        # H C = C E

* `eigh(h, iop=0)`: `h` is (2n, 2n); only its left half `h[:, :n]` is read.
  `e[:n]` are the n distinct eigenvalues in ascending order and `e[n:] = e[:n]`;
  `c[:, :n]` are the eigenvectors and `c[:, n:]` their Kramers partners,
  `c[:n, n:] = -c[n:, :n].conj()` and `c[n:, n:] = c[:n, :n].conj()`.
  With `iop=1` the pairs are interleaved instead: `e[0::2] = e[1::2]`, the
  partner of column 2k is column 2k+1.
* `geigh(f, s)`: the generalized problem F C = S C E (F and S both of the
  form above, S positive definite); eigenpairs interleaved as for `iop=1`,
  C^+ S C = 1.

## License

BSD 2-Clause; see [LICENSE](LICENSE).
