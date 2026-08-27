# MFEM Elasticity PA

[English](README.md)

本仓库将用于发布论文 **“Shifting the Sweet Spot: High-Performance Matrix-Free
Method for High-Order Elasticity”** 的官方开源实现。

**作者：** Dali Chang、Chong Zhang、Kaiqi Zhang、Mingguan Yang、Huiyuan Li 和
Weiqiang Kong。

## 项目简介

论文面向 MFEM 中的高阶三维线弹性问题，研究高性能 matrix-free Partial
Assembly operator。实现结合了 sum factorization、Voigt notation、
macro-kernel fusion 和 slice-wise loop organization，并在 GMG-PCG solver
workflow 中进行评测。

## 代码发布

> **maintained source code 正在整理和正确性验证中，完成后将在本仓库公开。**

首个版本预计包含：

- optimized elasticity PA kernel；
- 带 automatic fallback 的 MFEM-compatible integrator；
- GMG-PCG driver 和 beam-hex 示例；
- portable build 和 validation script。

目前尚未发布源码，因此暂时无法从本仓库复现论文的数值和性能结果。你可以 Watch
本仓库以获取后续发布通知。

本项目是独立的论文软件 companion。后续开发会与相关 MFEM upstream 工作协调，
包括 [MFEM PR #5438](https://github.com/mfem/mfem/pull/5438)，但不暗示 MFEM
endorsement。

## 引用

目前尚无稳定的 DOI 或 arXiv 标识。论文公开后，这里会补充正式引用信息。软件引用
元数据见 [`CITATION.cff`](CITATION.cff)。

## 许可证

本仓库采用 [BSD 3-Clause License](LICENSE)。
