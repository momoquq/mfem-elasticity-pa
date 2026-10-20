# MFEM Elasticity PA

[简体中文](README.zh-CN.md)

This repository contains the implementation accompanying the paper **“Shifting the Sweet Spot: High-Performance Matrix-Free Method for High-Order Elasticity.”**

**Publication:** *ACM Transactions on Architecture and Code Optimization (TACO)*, 23(4), Article 139 (2026).

**Paper:** [ACM Digital Library](https://dl.acm.org/doi/full/10.1145/3844725)

**Authors:** Dali Chang, Chong Zhang, Kaiqi Zhang, Mingguan Yang, Huiyuan Li, and Weiqiang Kong.

## Overview

The code implements optimized Partial Assembly (PA) for three-dimensional linear elasticity on affine tensor-product hexahedral meshes. It combines sum factorization, Voigt notation, macro-kernel fusion, and slice-wise loop organization, with a GMG-PCG driver based on [MFEM v4.8](https://github.com/mfem/mfem/tree/v4.8).

Across `p ∈ {1, 2, 4, 8}`, the [paper](https://dl.acm.org/doi/full/10.1145/3844725) reports `7–83×` kernel speedup and `3.6–16.8×` end-to-end speedup on AMD EPYC 7713. The operator-throughput sweet spot shifts to `p ≥ 6`, with the same trend observed on Huawei Kunpeng 920 (ARMv8.2). These figures refer to the experiments reported in the paper.

## Requirements

- An MPI-enabled MFEM 4.8 with HYPRE 2.31.0 and [METIS 5.1.0](https://github.com/mfem/tpls).

## Build dependencies (AMD Zen 3)

Run the following in a separate working directory with an MPI toolchain, Make, and CMake available. This example uses the AMD compiler flags from the experiments; adjust `CMAKE_CXX_FLAGS` for other CPUs.

### 1. HYPRE 2.31.0

```bash
wget https://github.com/hypre-space/hypre/archive/v2.31.0.tar.gz
tar -xzf v2.31.0.tar.gz
(
  cd hypre-2.31.0/src
  ./configure --prefix="$PWD/hypre" --with-MPI
  make -j$(nproc)
  make install
)
```

### 2. METIS 5.1.0

```bash
git clone https://github.com/mfem/tpls.git
tar -xzf tpls/metis-5.1.0.tar.gz
(
  cd metis-5.1.0
  make config shared=1 prefix="$PWD"
  make -j$(nproc)
  make install
)
```

### 3. MFEM 4.8

```bash
git clone --branch v4.8 https://github.com/mfem/mfem.git
(
  cd mfem
  mkdir build && cd build
  cmake .. -DCMAKE_INSTALL_PREFIX=../mfem \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_FLAGS="-O3 -march=znver3 -mavx2 -mfma -mprefetchwt1 -DNDEBUG -fopenmp-simd" \
    -DMFEM_USE_MPI=YES \
    -DMFEM_USE_SIMD=YES \
    -DHYPRE_DIR=../../hypre-2.31.0/src/hypre \
    -DMETIS_DIR=../../metis-5.1.0
  make -j$(nproc)
  make install
)
```

MFEM is installed in `mfem/mfem/`. Use its absolute path as `MFEM_INSTALL_DIR` below. For other configurations, see the [MFEM 4.8 build instructions](https://github.com/mfem/mfem/blob/v4.8/INSTALL).

## Build this project

From the repository root:

```bash
make -j MFEM_INSTALL_DIR=/path/to/mfem-install
```

Select a build profile with `PLATFORM`:

| Option | Target | Output directory |
|---|---|---|
| `PLATFORM=generic` | Default compiler target | `build/generic/` |
| `PLATFORM=amd` | AMD Zen 3 | `build/amd/` |
| `PLATFORM=arm` | ARMv8.2-A with SVE | `build/arm/` |

For example:

```bash
make -j MFEM_INSTALL_DIR=/path/to/mfem-install PLATFORM=amd
```

## Run

Use `beam-hex.mesh` from the [MFEM 4.8 source tree](https://github.com/mfem/mfem/blob/v4.8/data/beam-hex.mesh). From the repository root, run the optimized PA example with two MPI processes:

```bash
mpirun -np 2 ./build/generic/myex2 \
  -m /path/to/mfem/data/beam-hex.mesh \
  -r 0 -gr 0 -or 2 -pa -op -no-vis
```

`-m` points to the MFEM mesh data. Replace `-op` with `-no-op` to select MFEM's elasticity integrator.

| Option | Meaning |
|---|---|
| `-m <file>` | Input mesh |
| `-o <p>` | Coarse polynomial degree, default `1` |
| `-r <n>` | Parallel mesh refinements before building the multigrid hierarchy |
| `-gr <n>` | Geometric refinement levels in the multigrid hierarchy |
| `-or <p>` | Finest polynomial degree in the multigrid hierarchy |
| `-pa` | Enable Partial Assembly on the fine multigrid levels |
| `-op` / `-no-op` | Select the optimized / MFEM elasticity integrator |
| `-no-vis` | Disable GLVis visualization |

## Clean and rebuild

From the repository root, clean the AMD build and the generated mesh and solution files:

```bash
make clean PLATFORM=amd
```

Clean all build profiles and rebuild:

```bash
make distclean
make -j MFEM_INSTALL_DIR=/path/to/mfem-install PLATFORM=amd
```

For a custom `BUILD_DIR`, use the same path in the clean and build commands.

## Source layout

| Path | Content |
|---|---|
| `src/kernels/` | Elasticity action and diagonal kernels |
| `src/integ/` | MFEM integrator |
| `myex2.cpp` | GMG-PCG elasticity driver |
| `makefile` | Library and driver build |

## Citation

If our work helps you in your research, please consider citing: qwq

```bibtex
@article{10.1145/3844725,
  author = {Chang, Dali and Zhang, Chong and Zhang, Kaiqi and Yang, Mingguan and Li, Huiyuan and Kong, Weiqiang},
  title = {Shifting the Sweet Spot: High-Performance Matrix-Free Method for High-Order Elasticity},
  journal = {ACM Trans. Archit. Code Optim.},
  year = {2026},
  volume = {23},
  number = {4},
  articleno = {139},
  numpages = {23},
  doi = {10.1145/3844725},
  url = {https://doi.org/10.1145/3844725}
}
```

## License

This repository is licensed under the [BSD 3-Clause License](LICENSE).
