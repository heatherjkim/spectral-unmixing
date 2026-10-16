# spectral-unmixing
UV-Vis spectral unmixing using non-negative least squares.

# NNLS unmixing: `BSA_PPy_unmixing.m`

It turns raw UV-Vis spectra of Franz-cell receptor samples into **per-sample BSA concentrations (µg/mL)**. It separates BSA absorbance from PPy (polypyrrole) absorbance and from background drift with a non-negative least-squares (NNLS) fit.

It writes the per-sample concentrations to `BSA_PPy_unmixing_results.csv`.

---

## Why NNLS

BSA is quantified at 278 nm, but PPy leaching from the hydrogel also absorbs in the 250–310 nm UV range. A single-wavelength reading would therefore overestimate BSA. Instead, the whole 250–310 nm spectrum of each sample is fit as:

```
A(λ) = a_BSA · S_BSA(λ)  +  a_PPy · S_PPy(λ)  +  b0·P0(λ) + b1·P1(λ) + b2·P2(λ)
       └── BSA shape ──┘    └── PPy shape ──┘    └──── smooth baseline ────┘
```

- `S_BSA` and `S_PPy` are measured reference spectra (the "basis").
- `a_BSA, a_PPy ≥ 0`: neither analyte can be present in a negative amount.
- `b0, b1, b2` are free. The baseline can be positive or negative, which absorbs scattering, cuvette and buffer offsets.

Because `S_BSA` is scaled to 1 at 278 nm, the fitted `a_BSA` **is** the BSA-only absorbance at 278 nm. The standard BSA calibration line then converts it to µg/mL.

---

## Quick start

```matlab
cd('<folder containing the script>')   % optional: paths don't depend on the working folder
BSA_Drug_Delivery_NNLS_Cam_Ready
```

From VS Code: open the file and press **F5** (or run **MATLAB: Run File**). It takes under a minute.

**Output file.** Each run overwrites `BSA_PPy_unmixing_results.csv` in this folder (path can be changed using `outCsv` at the end of the script). It has a header row and one row per included sample:

```
file,geom,t_min,isPass,Vr_mL,C_BSA_ugmL
T0_B3C.csv,C,0,0,9.850898,-2.1199809779715002
```

`C_BSA_ugmL` is written at full precision (`%.17g`), so downstream analysis reproduces the values exactly. Nothing else is written to disk.

---

## Requirements

| Requirement | Notes |
|---|---|
| MATLAB (developed on R2024b) | R2019a or later (`readcell`, `readmatrix`, tables) |
| **Optimization Toolbox** | `lsqlin` (primary) and `quadprog` (backup) solve the constrained fit. Without the toolbox, the script stops with an error. |
| `Data/` folder next to the script | All paths are relative to the script's own folder (`scriptDir`, `dataDir` at the top of the settings). The script reads `Data/BSA_nml_cc_combined_final/`, `Data/BSA_PPy-PBS/`, `Data/BSA_PPy-PBS-02/` and `Data/Delivery_Test_Samples/`, and writes its CSV next to itself. Move the folder anywhere and it still runs. |

---

## Inputs

### Spectrum CSVs

There's no header. Column 1 is the wavelength (nm) and columns 2 onward are absorbance. Multiple absorbance columns are averaged. Rows are sorted, de-duplicated, and rows that aren't finite numbers are dropped.

### Reference and calibration folders

| Folder | Role | Baseline subtracted |
|---|---|---|
| `BSA_nml_cc_combined_final/` | All spectra are averaged into the **BSA basis** `S_BSA` | None |
| `BSA_PPy-PBS/` | PPy standards S3–S6. Used for the **PPy basis** and **PPy calibration**. | `baseline-pbs.csv` |
| `BSA_PPy-PBS-02/` | PPy standards S7-S9. Used for the PPy basis and calibration. | `baseline-pbs2.csv` |

The PPy basis `S_PPy` is the average of the two folders' mean spectra (after baseline subtraction). Note that the basis uses **every** CSV in these folders except the baselines, so it includes S7 and S9. The S-number exclusions below apply only to the calibration.

### PPy calibration standards (`calib.ratioList`)

| Standard | PPy ratio (v/v) | Used |
|---|---|---|
| S3 | 1 | ✔ |
| S4 | 0.556 | ✔ |
| S5 | 0.25 | ✔ |
| S6 | 0.125 | ✔ |
| S8 | 0.0625 | ✔ |
| S1, S2, S7, S9 | — | ✘ (`calib.excludeExact`; S7 and S9 excluded for high linearity) |

### Delivery samples

The script reads eight batch folders under `Delivery_Test_Samples/` (listed in `folders.samples`). Each file name encodes the sample:

```
T<min>_B<batch><GEOM>[_Pass][_DIL[n]].csv        e.g.  T180_B4H2_Pass_DIL.csv
```

| Part | Meaning |
|---|---|
| `T<min>` | Sampling time |
| `B<n><GEOM>` | Batch + geometry (`R`, `C`, `H2`, `H3`, `I`) |
| `_Pass` | PASS condition. Files without it are CIP. |
| `_DIL[n]` | Diluted before measurement. The factor comes from the metadata, not the name. |

These files are skipped: `baseline*` files, and the 7 stems in `EXCLUDE_FILES`. The latter were removed to keep 3 replicates per group (see `Delivery_Test_Samples/README.rtf`).

### Metadata: `Delivery_Test_Samples/metadata_map_final.csv`

No header. Each row is `stem, dilution factor, receptor type`:

```
T30_B1R_DIL,3,Rec
T0_B3C,1,Cir
```

Receptor type `Rec` means 10.255033 mL and `Cir` means 9.850898 mL (`VR_TYPES_mL`). Matching ignores case. If a sample has no DF or no receptor type, it is skipped and logged.

---

## Algorithm

### 1. Build the basis
1. Average the BSA references into `S_BSA`, and the baseline-subtracted PPy references into `S_PPy`.
2. Interpolate both onto a common wavelength grid (the overlapping range, at the finest step) and crop it to **250–310 nm** (`use_nm_window`).
3. Divide `S_BSA` by its value at **278 nm** (`targetNm`). If that value is ≤ 0, the script stops with an error.
4. Build a quadratic baseline `P = [1, λ, λ²]`, then **orthogonalise** it against `[S_BSA, S_PPy]`:
   ```
   P_orth = P − Q (QᵀQ)⁻¹ QᵀP,   Q = [S_BSA, S_PPy]
   ```
   This removes from the baseline any part that looks like BSA or PPy, so the free baseline terms can't absorb analyte signal.

### 2. Calibrate PPy amplitude → ratio
1. Baseline-subtract each kept standard and fit it with `[S_PPy, P_orth]` only (`nnls_one_analyte`). The fit alternates between `lsqnonneg` for the PPy amplitude and ordinary least squares for the baseline, for 40 iterations.
2. Fit a straight line from fitted amplitude to known ratio: `C_PPy = α·a_PPy + β`.
3. Set **PPy bounds** from the standards: amplitude and ratio ranges, widened by ±10 %. These are used only for flagging (step 3.2).

**Current calibration** (2026-10-13 run): α = 0.76666, β = −0.48856, R² = 0.994, RMSE = 0.026, n = 5.

### 3. Fit each sample
For every sample CSV (sorted by geometry, then time):

| # | Action | Outcome if it fails |
|---|---|---|
| 1 | Read the raw A(278) | — |
| 2 | Run a **pre-DF fit** on the undiluted-scale spectrum to get `aP_raw` and `C_PPy_raw`, then flag `ppy_preDF_high/low` if outside the PPy bounds | Flag only. The sample is kept. |
| 3 | **Screen:** raw A(278) > **0.4 AU** (top of BSA's linear range) | Sample **excluded** |
| 4 | Look up DF and receptor type in the metadata | Sample **skipped** |
| 5 | Multiply the spectrum by DF | — |
| 6 | **Fit** `[S_BSA, S_PPy, P_orth]` with `a_BSA, a_PPy ≥ 0` using `lsqlin`, falling back to `quadprog` | Script **stops with an error** |
| 7 | Compute RMSE and R² of the fit | — |
| 8 | `C_BSA = (a_BSA − 0.001112) / 0.000524533` µg/mL; `C_PPy = max(α·a_PPy + β, 0)` | — |

In the solver, the constraint is written as `−a ≤ 0` on only the first two coefficients (`Aineq = −eye(2, 5)`), which leaves the baseline coefficients free.

---

## Outputs (console)

| Block | Contents |
|---|---|
| `PPy calibration file scan` | KEPT/SKIP for each calibration file, with the reason |
| `PPy amplitude calibration` | α, β, R², RMSE, n, and the files used (printed twice) |
| `METADATA keys with NO matching file` | Metadata rows with no matching CSV. They are usually the excluded files, which is expected. |
| `[STRICT] Excluding …` / `EXCLUDED by pre-DF A(278)` | Samples removed by the 0.4 AU screen |
| `Ingest summary` | Number of candidate files after the exclusion list |
| `Results (key columns)` | Full per-sample table (see below) |
| `Metadata audit` | DF and receptor type used for each sample, and the metadata key that matched |
| **`Wrote N samples to …`** | Confirms `NNLS_Cam_Ready_results.csv` was written (`file, geom, t_min, isPass, Vr_mL, C_BSA_ugmL`) |

### Per-sample table columns

| Column | Meaning |
|---|---|
| `file, geom, t_min, isPass` | Sample identity, parsed from the file name |
| `DF, Vr_type, Vr_mL` | From the metadata |
| `A278_raw`, `A278_corr` | A(278) as measured, and after multiplying by DF |
| `a_BSA, a_PPy, b0, b1, b2` | Fit coefficients (DF-corrected spectrum) |
| `A_BSA_278` | BSA-only absorbance at 278 nm (= `a_BSA`) |
| `C_BSA_ugmL` | **BSA concentration**: the value written to the CSV |
| `C_PPy` | PPy ratio (v/v), clipped at 0 |
| `rmse, r2` | Fit quality |
| `aP_raw, C_PPy_raw` | PPy from the pre-DF fit |
| `ppy_preDF_high/low` | PPy out-of-range flags |

### Current run (2026-10-13)

- 188 candidate files after the exclusion list. **8 excluded** by the A278 screen (7 are H3, 1 is R). That leaves **180 samples** in `NNLS_Cam_Ready_results.csv`.
- BSA ranges from −2.1 to 916 µg/mL.

---

## Settings

| Variable | Value | Meaning |
|---|---|---|
| `use_nm_window` | `[250 310]` | Fit window (nm). Also used for the calibration (`calib.use_nm_window`). |
| `targetNm` | `278` | BSA quantitation wavelength / basis normalisation |
| `m_BSA`, `b_BSA` | `0.000524533`, `0.001112` | BSA calibration line: A = m·C + b (C in µg/mL) |
| `A278_preDF_max` | `0.4` | Raw A(278) exclusion threshold (AU) |
| `VR_TYPES_mL` | `rec` 10.255033, `cir` 9.850898 | Receptor volumes |
| `EXCLUDE_FILES` | 7 stems | Samples removed to keep replicates balanced |
| `calib.excludeExact` | baselines, S1, S2, S7, S9 | Calibration files that are not used |
| `calib.ratioList` | S3, S4, S5, S6, S8 | Known PPy ratios |
| `clip_neg_PPy` | `true` | Clip `C_PPy` below 0 to 0 |
| `metaMapFile` | `…/metadata_map_final.csv` | DF and receptor-type map |

---

## Strict mode

The script never guesses paths, baselines or dilution factors.

| Situation | Result |
|---|---|
| Metadata file, PPy baseline file, or all sample CSVs missing | **Error** |
| BSA basis ≤ 0 at 278 nm | **Error** |
| Fewer than 2 usable PPy standards | **Error** |
| Receptor type not `rec`/`cir` | **Error** |
| `lsqlin` and `quadprog` both fail | **Error** (check the Optimization Toolbox) |
| DF or receptor type missing for a sample | Sample skipped and logged |
| Raw A(278) > 0.4 AU | Sample excluded and listed |
| A sample folder has no usable CSVs | Warning, then the script continues |

---

## Caveats

- **Negative BSA values.** `a_BSA` is constrained to be ≥ 0, but the calibration intercept `b_BSA` is 0.001112. Any sample with `a_BSA < b_BSA` therefore gets a small negative concentration (as low as −2.1 µg/mL). In the current run, 26 of 180 samples are negative, mostly early time points with near-zero BSA. These values are **not** clipped; they are written to the CSV as they are.
- **Fit quality isn't used to filter.** `rmse` and `r2` are reported, but no cutoff is applied. Check them by eye in the results table.
- **PPy flags are informational.** Samples outside the calibrated PPy range are still fit and kept.
- **Extrapolation.** Every interpolation uses `'extrap'`. If a sample spectrum doesn't cover the full 250–310 nm range, the edge values are extrapolated without a warning.
- **The basis and the calibration use different PPy file sets** (see [Inputs](#reference-and-calibration-folders)).
- **No validation in this script.** I do have validation against samples of known concentration, which lives in `BSA_nnls_validate*.m` and the `NNLS_val/` data.
