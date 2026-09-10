# Helper consolidation — September 11, 2026

Ten helper files added during the mathematical repairs have been merged into
their primary class files:

| Primary file | Consolidated parts |
| --- | --- |
| [Maths.cs](../sources/Core/Maths.cs) | Arithmetic |
| [Matrice.cs](../sources/Core/Matrice.cs) | Statistics, Interpolation |
| [Special.cs](../sources/Core/Special.cs) | Numerics, Elementary, Series, ErrorFunctions, Integrals, Bessel, DistributionKernels |

The original WaveletPacket and Weyl–Heisenberg partial files are retained at the
user's request. The [source layout rules](../AGENTS.md) require future helpers to
stay in the primary file and include meaningful English XML documentation.

XML comments now describe 104 helper methods: their numerical purpose,
parameters, result, and relevant domains or branch conventions. This includes
the moved helpers, DistributionNumerics, the local mean-filter helper, and
trapezoidal moment calculation. Existing explanations of delicate numerical
steps are retained near their implementations.

Executable source tokens, declarations, and constants match the previous
implementations after accounting for the file moves. The generated XML
documentation was checked for all 104 helper methods.

The complete Windows Release audit contains **14,707 cases: 14,524 passed,
183 failed, none skipped**. Comparing every test ID against the B05 snapshot
shows no changed outcomes, added tests, or removed tests. B01–B05 remain closed;
all 183 failing cases retain their B06–B12 assignments.

Coverage remains **27,806 / 33,393 lines (83.26%)** and
**10,529 / 13,764 branches (76.49%)**. The source inventory now contains
**427 C# files** after removing the ten companion files.

- [Test-ID comparison](audit-consolidation/comparison.json).
- [Run summary and source digest](audit-consolidation/summary.json).
- [Current source inventory](audit-consolidation/source-inventory.md).
- [Current failures and block assignments](audit-consolidation/repair-blocks.json).
- [Old-to-current source paths](audit-consolidation/source-relocations.json).

Historical audit JSON/CSV files retain the paths and hashes from their original
runs. Markdown source links point to the current implementations.

The verification command was:

```powershell
./tools/Run-MathAudit.ps1 -NoRestore -ResultsDirectory artifacts/math-audit/consolidation/verified
```

It exits with status 1 because the 183 known failing expectations remain enabled.
