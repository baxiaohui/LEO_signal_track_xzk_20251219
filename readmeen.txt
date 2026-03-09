# LEO Tracking and Demodulation Simulation Code (BPSK + CSK Time-Division Structure)

This repository provides a receiver **tracking and demodulation** simulation codebase for **LEO navigation augmentation signals** in which the BPSK and CSK components are time-division multiplexed. The code supports multiple tracking strategies (PLL/FLL/KF/EKF/UKF, etc.) and can output tracking-error statistics, CSK demodulation BER statistics, and batch simulation summary tables. It is mainly intended for paper reproducibility and comparative evaluation.

> Note: This code is developed for research reproducibility. The default settings follow the simulation configurations used in the paper. For engineering or real-time deployment, further refactoring is required according to the target platform and interfaces.

---

## 1. Overview

### 1.1 Supported tracking/demodulation modes (`mode`)

Different algorithm branches are selected via `mode` (the following list illustrates the intended meaning; please refer to in-code comments for the exact implementation):

* `mode = 0`: PLL (phase-locked loop)
* `mode = 1`: FLL (frequency-locked loop, including HDC-FLL-related logic)
* `mode = 2`: LKF (linear Kalman filter)
* `mode = 3`: EKF (extended Kalman filter, with data-bit wipeoff)
* `mode = 4`: EKF (without data-bit wipeoff; **marked as not available / requires further derivation in the current version**)
* `mode = 5`: CSK-only demodulation (no tracking adjustment; mainly for observation/statistics)
* `mode = 6`: UKF (unscented Kalman filter, with data-bit wipeoff)
* `mode = 8`: Extrapolated/scheduled LKF (example: 1 ms prediction + 2 ms correction)
* `mode = 9`: Extrapolated/scheduled EKF (example: 1 ms prediction + 2 ms correction, with data-bit wipeoff)

### 1.2 Metrics and outputs

* Tracking errors: carrier-phase error, Doppler error, C/N0 estimation, etc.
* Demodulation errors: coherent/non-coherent CSK BER and (optional) BER after LDPC decoding
* Saved results: log `.txt`, summary `.xlsx`, and per-run workspace `.mat` snapshots

---

## 2. Requirements

* MATLAB: R2021a or later is recommended
* OS: Windows is more convenient (due to Excel I/O via `xlswrite/xlsread`)

  * On Linux/macOS: consider replacing `xlswrite` with `writetable` / `writecell`
* Toolboxes (depending on your code usage): Signal Processing / Communications (optional)

---

## 3. Code and data dependencies (must be on the MATLAB path)

### 3.1 Required scripts/functions (examples)

* `CAgenerate.m`: generate local PRN codes (low-rate/high-rate)
* `BIThightrans.m`: high-rate data-bit preprocessing/mapping
* `EMSdecoderforcsk.m`: CSK decoding (e.g., EMS decoding)
* `countBERandSER.m`: BER/SER statistics
* `Trackplotandcount.m`: tracking-error statistics and plotting/output

### 3.2 Required data files (example)

* `LDPC200_100high.mat`: LDPC-coded bit sequences (e.g., variable `resultbin1`)

---

## 4. Quick start (minimal reproducible run)

1. Save the main script as `main.m` (or your preferred filename), and make sure all functions/data files listed in Section 3 are accessible via the MATLAB path.
2. Configure key parameters in the script (see the next section), e.g.:

   * `CN0` (dB-Hz)
   * `mode`
   * `AcqSatNUM` (PRN)
   * `LocalDop` (initial Doppler)
3. Run:

   * Execute `main` in the MATLAB command window, or click the Run button.
4. After the run, check the outputs in Section 6.

---

## 5. Key parameters (recommended starting point for modifications)

* `CN0`: carrier-to-noise density ratio (dB-Hz); affects signal amplitude, demodulation enabling, and statistics
* `mode`: algorithm selection (see Section 1.1)
* `AcqSatNUM`: target satellite PRN
* `LocalDop`: initial Doppler (recommended to be close to the true/assumed value)
* `Freq_sample`: sampling frequency (Hz)
* `Fc`: intermediate frequency (Hz)
* `Freq_code`: code rate (Hz)
* `runtime`: run time (s); in some versions it may be computed based on the visible arc
* `anger`: geometry-related angle parameter (used to construct pass duration/Doppler profiles)

---

## 6. Outputs (generated after running)

* `LKF.txt` (example filename): run log (e.g., per-update records such as `loopcnt`, carrier frequency, discriminator/innovation, etc.)
* `simulation_results.xlsx`: simulation summary table (automatically appended)

  * Typical columns: CN0, mode, BER (coherent/non-coherent/LDPC-coherent), carrier-phase RMSE, Doppler RMSE, C/N0 mean, etc.
* `workspace_data_<CN0>_<mode>.mat`: per-run workspace snapshot for reproducibility and plotting

---

## 7. Reproducibility tips (for paper comparisons)

* It is recommended to follow the paper setting and perform sweeps over `mode × CN0`.
* After each run:

  1. Read `simulation_results.xlsx` and generate summary plots
  2. If detailed inspection is needed, load the corresponding `workspace_data_<CN0>_<mode>.mat`
* For repeatability:

  * Fix the random seed (e.g., `rng(1)`)
  * Keep parameters consistent (loop bandwidth, coherent integration time, non-coherent accumulation count, etc.)

---

## 8. FAQ

**Q1: The run reports missing `LDPC200_100high` or `resultbin1`.**
A: Please make sure `LDPC200_100high.mat` is in the working directory or on the MATLAB path, and that it contains the variable `resultbin1`.

**Q2: `xlswrite/xlsread` is not available or fails on non-Windows systems.**
A: Replace them with `writecell/writetable`, or run on Windows to generate the Excel output and post-process on other systems.

**Q3: `mode=4` is marked unavailable or behaves abnormally.**
A: This branch is an experimental EKF variant without data-bit wipeoff (including a data-bit state). It requires further derivation and validation before being used for paper conclusions.

---

## 9. Citation and disclaimer

If you use this repository in your paper/report, please cite the source in the acknowledgements or the code-availability statement (author/repository link/version or commit hash).

* Author: Xuzhekai
* Purpose: research reproducibility and algorithm comparison

---

## 10. Contact

For help with reproducing paper settings, interpreting output metrics, or adding comparison scripts, please open an issue in this repository or contact the author.
