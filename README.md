# MFEM Elasticity PA

[简体中文](README.zh-CN.md)

This repository will host the official open-source implementation accompanying the paper **“Shifting the Sweet Spot: High-Performance Matrix-Free Method for High-Order Elasticity.”**

**Publication status:** Accepted for publication in *ACM Transactions on Architecture and Code Optimization (TACO)*.

**Preprint:** [arXiv:2601.08374](https://arxiv.org/abs/2601.08374)

**Authors:** Dali Chang, Chong Zhang, Kaiqi Zhang, Mingguan Yang, Huiyuan Li, and Weiqiang Kong.

## Paper overview

[MFEM](https://github.com/mfem/mfem)'s native [linear-elasticity Partial Assembly path](https://github.com/mfem/mfem/blob/v4.8/fem/bilininteg.hpp#L3148) applies an `O((p+1)^6)` contraction in the element operator, leaving the CPU operator-throughput sweet spot near `p ≈ 2` in the paper's baseline measurements. The paper addresses this gap for affine tensor-product hexahedral meshes by combining sum factorization, Voigt notation, macro-kernel fusion, and slice-wise loop organization, and evaluates the resulting operator in a GMG-PCG solver workflow.

Across `p ∈ {1, 2, 4, 8}`, the paper reports `7–83×` kernel speedup and `3.6–16.8×` end-to-end speedup on AMD EPYC 7713 (x86-64). The operator-throughput sweet spot shifts to `p ≥ 6`, and the same trend is reproduced on Huawei Kunpeng 920 (ARMv8.2). The paper also includes per-stage ablation and hardware-counter characterization.

## Code release

> **The source code is being prepared and will be released in this repository after correctness validation.**

The first release is expected to include:

- the optimized elasticity PA kernel;
- an MFEM-compatible integrator with automatic fallback;
- the GMG-PCG driver and beam-hex example;
- a portable build and validation script.

The results above are reported in the paper. No source release is available yet, so they cannot currently be reproduced from this repository. Watch this repository for release updates.

## Related upstream work

- [MFEM PR #5438](https://github.com/mfem/mfem/pull/5438) — correctness fixes, compact PA data, and tensor-product optimization for MFEM's `ElasticityIntegrator`. This recent MFEM PR is close to our work in direction, and we plan to contribute a complementary high-order CPU path to MFEM.

## Citation

The paper has been accepted for publication in *ACM Transactions on Architecture and Code Optimization (TACO)*. This entry will be updated with the publisher DOI once it is assigned. Until then, please cite the arXiv preprint:

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

## License

This repository is licensed under the [BSD 3-Clause License](LICENSE).
