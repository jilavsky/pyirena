# Size Distribution Fitting Methods

pyIrena implements three independent inversion methods for recovering a particle
size distribution P(r) from small-angle scattering (SAS) data.  All three share
the same G-matrix framework — the forward model is always:

```
I(Q) = G · x + background
```

where `x` is the bin amplitude vector (related to P(r) by bin width), and each
row of G contains the form-factor integral for one (Q, r) pair.  The three
methods differ only in how they choose x from the (severely underdetermined)
inverse problem.

---

## 1. Maximum Entropy (MaxEnt)

### Concept

MaxEnt finds the **smoothest possible** distribution that is **just consistent**
with the data.  Two competing objectives are balanced at every iteration:

| Objective | Direction |
|-----------|-----------|
| Maximise entropy S | Push x toward a flat, featureless distribution |
| Achieve χ² = M | Push x until the fit matches the data within noise |

**Entropy** here is defined as:

```
S = -Σⱼ (xⱼ/Z) ln(xⱼ/sky)
```

where Z = Σ xⱼ and `sky` is a small positive reference value (the prior).
S is maximised by the flattest possible distribution, so MaxEnt produces the
least-structured solution consistent with the data — it introduces no peaks or
features that the data do not require.

**Chi-squared target** is χ² = M (number of data points), meaning on average
each measured point is reproduced within 1σ.

### Algorithm

Each iteration proceeds as follows:

1. Compute the **entropy gradient** `sgrad` — which direction in the N-dimensional
   bin space increases entropy (pushes toward the flat prior).
2. Compute the **chi-squared gradient** `cgrad` — which direction decreases χ²
   (moves the model toward the data).
3. Build two orthonormal search directions ξ₀, ξ₁ from these gradients
   (Gram-Schmidt); reduce the N-dimensional problem to a 2×2 linear system.
4. Solve for step sizes β₀, β₁ that simultaneously drive χ² toward M and
   keep entropy high (Lagrange balance: λ = (M − χ²)/χ²).
5. Apply the step with backtracking (halve if χ² blows up); enforce x > 0.
6. Declare convergence when |χ² − M| < tolerance **and** the two gradients
   are nearly collinear — the KKT condition for the constrained optimum.

### Parameters

| Parameter | Default | Meaning |
|-----------|---------|---------|
| **Sky background** | `1e-6` | Prior value — all bins start here; sets the entropy reference and positivity floor. Should be much smaller than the expected volume fraction (typically 10⁻⁴ to 10⁻¹). The result is insensitive to the exact value as long as it is negligibly small. |
| **Max iterations** | `1000` | Safety cap. Well-posed problems converge in 20–200 iterations. If routinely hitting the limit, check whether the background subtraction has left negative corrected intensities, or whether the contrast value is correct. |

> **Note — Stability (internal constant, not exposed):**
> The convergence tolerance is `tol = 0.01 × √(2M)`.  This is the
> standard Skilling-Bryan constant used in all MEMSYS-class implementations
> and was a hardcoded internal parameter in the original Igor Pro IR1R code.
> It does not need user adjustment.

### When to use MaxEnt

- When the data are **well resolved** and errors are reliable.
- When you want the **most conservative interpretation** — no artificial peaks.
- For **unimodal** distributions it is particularly well-behaved.
- For **bimodal** distributions with a large size ratio MaxEnt may merge
  the two peaks; try Regularization or TNNLS in that case.

---

## 2. Tikhonov Regularization

### Concept

Regularization replaces the ill-posed inversion with a penalised least-squares
problem:

```
minimise  ||W^½ (I − G·x)||²  +  α ||L·x||²   subject to  x ≥ 0
```

- **Data term** `||W^½ (I − G·x)||²` — weighted chi-squared (W = diag(1/err²)).
- **Smoothness penalty** `α ||L·x||²` — penalises roughness via the
  2nd-order finite-difference matrix L (rows `[1, −2, 1]`), which approximates
  the second derivative of x(r).
- **α** (the regularisation parameter) controls the trade-off: large α → very
  smooth; small α → follows data closely.

The algorithm finds α by **binary search** in log₁₀ space so that the
final χ² equals M (same target as MaxEnt).  If the target is not achievable
(e.g. because background-corrected data contain negative values), it falls back
to the L-curve elbow (the α that minimises χ² without imposing the target).

Non-negativity (x ≥ 0) is enforced exactly via
`scipy.optimize.nnls` (Lawson-Hanson algorithm) on the augmented system
`[W^½·G ; √α·L]`.  This is more reliable than unconstrained least-squares
followed by clamping.

### Parameters

| Parameter | Default | Meaning |
|-----------|---------|---------|
| **χ² tolerance** | `1.0` | Binary search stops when `|χ² − M| < evalue × √(2M)`. Larger values accept a looser fit; `1.0` is the standard (same scale as MaxEnt). |
| **Min ratio** | `1e-4` | After convergence, bins below `min_ratio × max(x)` are raised to this floor. Suppresses numerical noise in near-zero tails without affecting the main peaks. |
| **Max iterations** | 100 (internal binary search) | The binary search converges in <50 steps for any physical dataset; this is not exposed to the user. |

### When to use Regularization

- The most **generally reliable** method for a wide range of sample types.
- Good default choice when you are unsure which method to use.
- Handles **bimodal** distributions better than MaxEnt in most cases.
- If background subtraction leaves some negative corrected intensities,
  the L-curve fallback ensures a useful result is still returned.

---

## 3. TNNLS / IPG (Truncated Non-Negative Least Squares)

### Concept

TNNLS (also called Interior-Point Gradient, IPG) solves the plain
non-negative least squares problem:

```
minimise  ||A·x − b||²   subject to  x ≥ 0
```

where A = G/err and b = I/err (error-weighted).  There is **no explicit
smoothness penalty** — the only constraint is non-negativity.

The method uses a diagonally preconditioned gradient descent that stays
strictly inside the positive orthant:

1. Compute gradient: `grad = AᵀAx − Aᵀb`
2. Diagonal preconditioner: `D = x / (AᵀA·x)` — scales the step so that
   negative-gradient components cannot overshoot zero.
3. Search direction: `P = −D · grad`
4. Exact line search (optimal step for quadratic objective).
5. Step limiting: if any component of x would go negative, the step is
   reduced by the `approach_param` safety factor.
6. Convergence when reduced χ² ≤ 1 (i.e. χ² ≤ M).

### Parameters

| Parameter | Default | Meaning |
|-----------|---------|---------|
| **Approach param** | `0.95` | Safety factor (< 1) applied when a step would drive any bin to zero. `0.95` means "go 95 % of the way to the zero boundary". Reducing toward 0.5 makes the algorithm more conservative but slower. |
| **Max iterations** | `1000` | Safety cap. TNNLS typically converges in 50–500 iterations. If routinely hitting the limit, the data may be too noisy or the model range may be mismatched to the data. |

### When to use TNNLS

- When you expect **sharp, narrow peaks** — TNNLS imposes no smoothness,
  so it can resolve closely spaced populations that MaxEnt or Regularization
  would merge.
- For **exploratory analysis** — the unregularised result reveals the natural
  resolution limit of the data.
- Be aware that TNNLS results can be **noisy** (many small spurious peaks)
  when data quality is limited; always cross-check with Regularization.

---

## 4. McSAS (Monte Carlo Size Analysis)

### Concept

McSAS (Bressler et al., *J. Appl. Cryst.* **48**, 2015) approaches the inversion
problem differently from the three deterministic methods above.  Instead of
solving a linear system for the entire distribution at once, it builds the
model from **N_c discrete particle contributions**, each with an independently
optimised radius.

The scattered intensity from N_c contributions with a common scale factor A is:

```
I(Q) ≈ A × Σⱼ g(Q, rⱼ)
```

where `g(Q, rⱼ)` is the form-factor column for radius `rⱼ` (the same quantity
as one column of the G-matrix used by the other methods).

The scale factor A is updated at each step by exact weighted least squares:

```
A = Σᵢ [I(Qᵢ) g_sum(Qᵢ) / σᵢ²]  /  Σᵢ [g_sum(Qᵢ)² / σᵢ²]
```

so that it always minimises the weighted χ² for the current set of radii.

### Algorithm

Each **repetition** (there are `N_repetitions` independent ones) proceeds as:

1. Draw `N_c` radii log-uniformly from [R_min, R_max].
2. Compute the total G-matrix sum `g_sum = Σⱼ g(Q, rⱼ)`.
3. Find the optimal scale factor A by the formula above.
4. Compute χ².
5. **Replacement loop** (up to `max_iter` steps):
   a. Pick a random contribution index j.
   b. Draw a new candidate radius r_new.
   c. Compute the trial `g_sum_trial = g_sum − g(Q, rⱼ) + g(Q, r_new)`.
   d. Compute A_trial and χ²_trial.
   e. **Accept** if χ²_trial ≤ χ² (Metropolis criterion at T=0).
   f. Stop early if χ²/M ≤ `convergence` (default 1.0).
6. Histogram the accepted radii {rⱼ} onto the output r_grid:
   `x_raw[k] = A × count[k]`   (volume fraction in bin k).

After all repetitions the per-bin mean and **standard deviation** of x_raw
are computed.  The standard deviation gives a per-bin uncertainty that reflects
how reproducibly the Monte Carlo reconstruction recovers the distribution —
it is shown as ±1σ error bars on the P(r) plot.

### Parameters

| Parameter | Default | Meaning |
|-----------|---------|---------|
| **N contrib** | `200` | Number of discrete particle contributions per repetition. More contributions give a smoother, more stable result but slow convergence. For unimodal distributions 50–200 is usually sufficient; for complex multi-modal systems increase to 500–1000. |
| **N reps** | `10` | Number of independent repetitions. Error bars on P(r) are estimated from the spread across repetitions. Fewer than 5 gives unreliable uncertainty estimates; 10–20 is the practical range. |
| **Convergence** | `1.0` | Stop the replacement loop when χ²/M ≤ this value (1.0 means on average one σ per point). Values > 1 stop earlier (rougher fit); values < 1 require a tighter fit (may never converge for noisy data). |
| **Max iter** | `100000` | Safety cap on the number of replacement attempts per repetition. Typical convergence requires 10–50 × N_c steps; 100 000 is generous for N_c = 200. |

### When to use McSAS

- When you need a **model-independent uncertainty estimate** on P(r) — the
  ±1σ error bars directly reflect information content of the data.
- For **exploratory analysis** of data with uncertain noise levels, where the
  deterministic methods' χ² targets may be poorly calibrated.
- When the number of features in the distribution is not known a priori —
  McSAS imposes no smoothness, so it can resolve features that MaxEnt or
  Regularization would suppress.
- **Caution**: the result depends on `N_c` (contribution count).  If N_c is
  too small the histogram is sparsely populated and spurious gaps appear.
  Always check that the P(r) shape is stable as N_c is increased.

---

## Comparison Summary

| | MaxEnt | Regularization | TNNLS | McSAS |
|---|---|---|---|---|
| **Regularisation** | Entropy (information-theoretic) | Tikhonov smoothness penalty | None (positivity only) | None (positivity only) |
| **χ² target** | Hard (χ² = M) | Hard, with L-curve fallback | Soft (χ² ≤ M) | Soft (χ²/M ≤ convergence) |
| **Smoothness** | Implicit (maximum flatness) | Explicit (2nd-derivative penalty) | None | None |
| **Sharpness** | Low | Medium | High | High |
| **Uncertainty** | None | None | None | Yes (±1σ from repetitions) |
| **Robustness** | Good, sensitive to negative data | Very good, fallback available | Good | Good |
| **Speed** | Medium (20–200 iters) | Fast (<50 binary-search steps) | Medium (50–500 iters) | Slow (MC loops × repetitions) |
| **Best for** | Conservative/unimodal | General use | Sharp/exploratory | Uncertainty quantification |

All four methods produce physically equivalent results when the data are
high quality and the size range is well chosen.  Significant disagreement
between methods indicates that the problem is underdetermined — the data do not
uniquely constrain the distribution.

---

## Shared Parameters (all methods)

### Radius grid

| Parameter | Default | Meaning |
|-----------|---------|---------|
| **R min** | 10 Å | Smallest particle radius modelled. Should be consistent with the maximum Q in your data: `R_min ≈ π/Q_max`. |
| **R max** | 1000 Å | Largest particle radius modelled. Should be consistent with the minimum Q: `R_max ≈ π/Q_min`. |
| **N bins** | 200 | Number of radius bins. More bins give finer resolution but make the problem more underdetermined. 50–200 is typical. |
| **Log spacing** | Yes | Logarithmically spaced bins give equal resolution per decade — strongly recommended for SAS data spanning multiple decades of r. |

### Background

| Parameter | Default | Meaning |
|-----------|---------|---------|
| **Background** | 0 | Flat incoherent background [cm⁻¹] subtracted before fitting. |
| **B (power-law amplitude)** | 0 | Amplitude of the low-Q power-law contribution B·Q⁻ᴾ. |
| **P (power-law exponent)** | 4 | Exponent of the low-Q power law (4 = Porod). |

### Other

| Parameter | Default | Meaning |
|-----------|---------|---------|
| **Contrast (Δρ)²** | 1.0 | Scattering length density contrast in units of 10²⁰ cm⁻⁴. Scales the absolute amplitude of P(r): volume fraction = ∫P(r)dr / contrast. Getting this right is essential for a physically meaningful volume fraction. |
| **Error scale** | 1.0 | Multiplies all measurement errors before fitting. Values > 1 make the fit trust the data less (more regularisation); < 1 does the opposite. Use only if your error bars are known to be systematically under- or over-estimated. |

---

## References

- J. Skilling & R.K. Bryan, *MNRAS* **211**, 111 (1984) — MaxEnt algorithm.
- P.C. Hansen, *Rank-Deficient and Discrete Ill-Posed Problems*, SIAM (1998) — Tikhonov regularisation, L-curve.
- J. Ilavsky & P.R. Jemian, *J. Appl. Cryst.* **42**, 347 (2009) — Irena package (Igor Pro original).
- C.L. Lawson & R.J. Hanson, *Solving Least Squares Problems*, Prentice-Hall (1974) — NNLS.
- I. Bressler, B.R. Pauw & A.F. Thünemann, *J. Appl. Cryst.* **48**, 962 (2015) — McSAS Monte Carlo method.
