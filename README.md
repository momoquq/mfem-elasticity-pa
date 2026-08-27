# MFEM Elasticity PA

[简体中文](README.zh-CN.md)

This repository will host the official open-source implementation accompanying
the paper **“Shifting the Sweet Spot: High-Performance Matrix-Free Method for
High-Order Elasticity.”**

**Authors:** Dali Chang, Chong Zhang, Kaiqi Zhang, Mingguan Yang, Huiyuan Li,
and Weiqiang Kong.

## About

The paper develops a high-performance matrix-free Partial Assembly operator for
high-order 3D linear elasticity in MFEM. It combines sum factorization, Voigt
notation, macro-kernel fusion, and slice-wise loop organization, and evaluates
the operator in a GMG-PCG solver workflow.

## Code release

> **The maintained source code is being prepared and will be released in this
> repository after correctness validation.**

The first release is expected to include:

- the optimized elasticity PA kernel;
- an MFEM-compatible integrator with automatic fallback;
- the GMG-PCG driver and beam-hex example;
- a portable build and validation script.

No source release is available yet, so the paper's numerical and performance
results cannot currently be reproduced from this repository. Watch this
repository for release updates.

This is an independent paper companion. Development will be coordinated with
relevant MFEM upstream work, including
[MFEM PR #5438](https://github.com/mfem/mfem/pull/5438), without implying MFEM
endorsement.

## Citation

A stable DOI or arXiv identifier is not yet available. The paper citation will
be added here when it becomes public. Software citation metadata is available
in [`CITATION.cff`](CITATION.cff).

## License

This repository is licensed under the [BSD 3-Clause License](LICENSE).
