# MFEM Elasticity PA

[English](README.md)

本仓库将用于发布论文 **“Shifting the Sweet Spot: High-Performance Matrix-Free Method for High-Order Elasticity”** 的官方开源实现。

**出版状态：** 已被 *ACM Transactions on Architecture and Code Optimization (TACO)* 接收。

**预印本：** [arXiv:2601.08374](https://arxiv.org/abs/2601.08374)

**作者：** Dali Chang、Chong Zhang、Kaiqi Zhang、Mingguan Yang、Huiyuan Li 和 Weiqiang Kong。

## 项目简介

[MFEM](https://github.com/mfem/mfem) 原生 [linear-elasticity Partial Assembly path](https://github.com/mfem/mfem/blob/v4.8/fem/bilininteg.hpp#L3148) 的 element operator 使用 `O((p+1)^6)` contraction，在论文的 baseline measurements 中，CPU operator-throughput sweet spot 位于 `p ≈ 2` 附近。论文面向 affine tensor-product hexahedral meshes，结合 sum factorization、Voigt notation、macro-kernel fusion 和 slice-wise loop organization 弥补这一实现差距，并在 GMG-PCG solver workflow 中评测优化后的 operator。

在 `p ∈ {1, 2, 4, 8}` 范围内，论文在 AMD EPYC 7713（x86-64）上报告了 `7–83×` kernel speedup 和 `3.6–16.8×` end-to-end speedup。operator-throughput sweet spot 移动到 `p ≥ 6`，并在 Huawei Kunpeng 920（ARMv8.2）上复现了相同趋势。论文还提供了 per-stage ablation 和 hardware-counter characterization。

## 代码发布

> **源代码正在整理和正确性验证中，完成后将在本仓库公开。**

首个版本预计包含：

- optimized elasticity PA kernel；
- 带 automatic fallback 的 MFEM-compatible integrator；
- GMG-PCG driver 和 beam-hex 示例；
- portable build 和 validation script。

以上是论文报告的实验结果。目前尚未发布源码，因此暂时无法从本仓库复现这些结果。你可以 Watch 本仓库以获取后续发布通知。

## 相关 upstream 工作

- [MFEM PR #5438](https://github.com/mfem/mfem/pull/5438) — 面向 MFEM `ElasticityIntegrator` 的 correctness fixes、compact PA data 和 tensor-product optimization 工作。这份近期的 MFEM PR 方向与我们的工作较为接近，本项目后续计划向 MFEM 贡献互补的 high-order CPU 路径。

## 引用

论文已被 *ACM Transactions on Architecture and Code Optimization (TACO)* 接收，正式 DOI 分配后将更新本条目。在此之前请引用 arXiv 预印本：

```bibtex
@misc{chang2026shiftingsweetspothighperformance,
      title={Shifting the Sweet Spot: High-Performance Matrix-Free Method for High-Order Elasticity},
      author={Dali Chang and Chong Zhang and Kaiqi Zhang and Mingguan Yang and Huiyuan Li and Weiqiang Kong},
      year={2026},
      eprint={2601.08374},
      archivePrefix={arXiv},
      primaryClass={cs.DC},
      url={https://arxiv.org/abs/2601.08374},
}
```

## 许可证

本仓库采用 [BSD 3-Clause License](LICENSE)。
