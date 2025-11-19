# Algebraic Complex Numbers
A library implementing precise representation of complex numbers in a form of $(4+1)$ tuple.
A complex number $c$ is represented as $(n_0, n_1, n_2, n_3, k) \in \mathbb{Z}^5$
where the value of $c$ is given by

$$c = (\frac{1}{\sqrt{2}})^k \big(n_0 e^{0 \cdot \frac{2 \pi i}{4}} + n_1 e^{1 \cdot \frac{2 \pi i}{4}} + n_2 e^{2 \cdot \frac{2 \pi i}{4}} + n_3 e^{2 \cdot \frac{2 \pi i}{4}} \big).$$

We are currently working on extending the library to support representation with arbitrary division of the unit half-circle, i.e.,

$$c = (\frac{1}{\sqrt{2}})^k \big(n_0 e^{0 \cdot \frac{2 \pi i}{N}} + \cdots + n_{(N-1)} e^{(N-1) \cdot \frac{2 \pi i}{N}} \big)$$

where $N$ is some power of two.

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
