# Algebraic Complex Numbers
A library implementing precise representation of complex numbers in a form of $(4+1)$ tuple.
A complex number $c$ is represented as $(n_1, n_2, m_1, m_2, k) \in \mathbb{Z}^5$
where the value of $c$ is given by
$$c = (\frac{1}{\sqrt{2}})^k n_1 + \sqrt{2} n_2 + m_1 i + m_2 i \sqrt{2}.$$
The library is developed to be used in formal-method tools for quantum computing
such as simulators and verifiers.


## Building
Dependencies (for fedora):
- A C++ compiler
- `cmake, make`: build tools
- `catch2-devel`: For tests

To build the project, use the following command sequence:
```
mkdir build
cd build
cmake ../src
make -j$(nproc)
```

After the build process finishes, the library is placed `build/libalgebraic_complex_numbers.a`.
