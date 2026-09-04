# Provenance

This directory was prepared from the exact formulas printed in
`k2_higgs_zero_reconnection_v189.tex` for the charge-two dyonic-instanton
sharp nodal theorem and its equal-weight isosceles auxiliary calculation.

The central proof script reconstructs the nodal polynomial from the physical
JNR coefficient identity rather than hard-coding the final `F_nod` polynomial.
All resultants, factorizations, derivative checks and Bernstein coefficients
are recomputed from exact rational/algebraic input using SymPy.

The auxiliary Sturm script reproduces the displayed scaled Sturm chain and
endpoint signs for the equal-weight isosceles reconnection quintic.  That
Sturm calculation is not logically required by the sharp `Xi_J >= 4` proof;
it is archived because it is another exact symbolic calculation printed in
the manuscript.
