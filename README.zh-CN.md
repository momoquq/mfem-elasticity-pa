# MFEM Elasticity PA

[English](README.md)

本仓库提供论文 **“Shifting the Sweet Spot: High-Performance Matrix-Free Method for High-Order Elasticity”** 的配套实现。

**发表信息：** *ACM Transactions on Architecture and Code Optimization (TACO)*，23(4)，文章编号 139（2026）。

**论文：** [ACM Digital Library](https://dl.acm.org/doi/full/10.1145/3844725)

**作者：** Dali Chang、Chong Zhang、Kaiqi Zhang、Mingguan Yang、Huiyuan Li 和 Weiqiang Kong。

## 项目简介

代码面向 affine tensor-product 六面体网格上的三维线弹性问题，实现优化的 Partial Assembly（PA）。核心方法包括 sum factorization、Voigt notation、macro-kernel fusion 和 slice-wise 循环组织，并提供基于 [MFEM v4.8](https://github.com/mfem/mfem/tree/v4.8) 的 GMG-PCG 驱动。

在 `p ∈ {1, 2, 4, 8}` 范围内，[论文](https://dl.acm.org/doi/full/10.1145/3844725)在 AMD EPYC 7713 上报告了 `7–83×` kernel 加速和 `3.6–16.8×` 端到端加速。operator-throughput sweet spot 移动到 `p ≥ 6`，并在 Huawei Kunpeng 920（ARMv8.2）上观察到相同趋势。这些数字来自论文中的实验。

## 依赖

- 启用 MPI 并链接 HYPRE 2.31.0 和 [METIS 5.1.0](https://github.com/mfem/tpls) 的 MFEM 4.8。

<details><summary>构建依赖（AMD Zen 3）</summary>

在已配置 MPI 工具链、Make 和 CMake 的环境中，选择单独的工作目录执行以下命令。示例采用实验中的 AMD 编译选项；其他 CPU 请相应调整 `CMAKE_CXX_FLAGS`。

```bash
# HYPRE
wget https://github.com/hypre-space/hypre/archive/v2.31.0.tar.gz
tar -xzf v2.31.0.tar.gz
(
  cd hypre-2.31.0/src
  ./configure --prefix="$PWD/hypre" --with-MPI
  make -j$(nproc)
  make install
)

# METIS
git clone https://github.com/mfem/tpls.git
tar -xzf tpls/metis-5.1.0.tar.gz
(
  cd metis-5.1.0
  make config shared=1 prefix="$PWD"
  make -j$(nproc)
  make install
)

# MFEM
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

MFEM 安装在 `mfem/mfem/` 中。下文的 `MFEM_INSTALL_DIR` 使用该目录的绝对路径。其他配置见 [MFEM 4.8 构建说明](https://github.com/mfem/mfem/blob/v4.8/INSTALL)。

</details>

## 构建

在仓库根目录执行：

```bash
make -j MFEM_INSTALL_DIR=/path/to/mfem-install
```

通过 `PLATFORM` 选择构建配置：

| 选项 | 目标平台 | 输出目录 |
|---|---|---|
| `PLATFORM=generic` | 编译器默认目标 | `build/generic/` |
| `PLATFORM=amd` | AMD Zen 3 | `build/amd/` |
| `PLATFORM=arm` | ARMv8.2-A，支持 SVE | `build/arm/` |

例如：

```bash
make -j MFEM_INSTALL_DIR=/path/to/mfem-install PLATFORM=amd
```

## 运行

使用 [MFEM 4.8 源码目录](https://github.com/mfem/mfem/blob/v4.8/data/beam-hex.mesh)中的 `beam-hex.mesh`。在仓库根目录，以两个 MPI 进程运行优化 PA 示例：

```bash
mpirun -np 2 ./build/generic/myex2 \
  -m /path/to/mfem/data/beam-hex.mesh \
  -r 0 -gr 0 -or 2 -pa -op -no-vis
```

`-m` 指向 MFEM 网格数据。把 `-op` 改为 `-no-op` 即可选择 MFEM 的线弹性 integrator。

| 选项 | 含义 |
|---|---|
| `-m <file>` | 输入网格 |
| `-o <p>` | 粗层多项式阶数，默认 `1` |
| `-r <n>` | 建立多重网格层次前的并行网格加密次数 |
| `-gr <n>` | 多重网格层次中的几何加密层数 |
| `-or <p>` | 多重网格层次中最细层的多项式阶数 |
| `-pa` | 在多重网格细层启用 Partial Assembly |
| `-op` / `-no-op` | 选择优化的 / MFEM 的线弹性 integrator |
| `-no-vis` | 关闭 GLVis 可视化 |

## 源码结构

| 路径 | 内容 |
|---|---|
| `src/kernels/` | 线弹性 action 与对角组装 kernel |
| `src/integ/` | MFEM integrator |
| `myex2.cpp` | GMG-PCG 线弹性驱动 |
| `makefile` | 库与驱动的构建入口 |

## 引用

如果我们的工作对你的研究有所帮助，请考虑引用： qwq

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

## 许可证

本仓库采用 [BSD 3-Clause License](LICENSE)。
