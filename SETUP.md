# Environment Setup

## Why Python 3.11?

Python 3.11 is a stable LTS-era release that is widely available on both macOS and Windows via Anaconda/conda-forge. It avoids the floating-point and numerical-library differences observed between Python 3.9 (original Mac environment) and Python 3.13 (Windows), while still being fully compatible with all pinned package versions listed below. Python 3.12+ introduced changes to internal floating-point handling in some libraries that can shift peak-detection results; 3.11 provides a stable, reproducible middle ground.

**Both Windows and Mac users must use this exact environment to obtain identical numerical results.**

---

## Creating the Environment

```bash
conda env create -f environment_frozen.yml
```

This will create a conda environment named `tfnt-analysis` with all packages pinned to the exact versions from the reference Mac run.

## Activating the Environment

```bash
conda activate tfnt-analysis
```

Run this command before launching any analysis scripts or notebooks.

## Verifying the Environment

```bash
conda activate tfnt-analysis
python -c "import pandas, numpy, scipy, sklearn; print('OK')"
```

## Notes

- All package versions are pinned with `=` (exact match) to ensure bit-for-bit reproducible numerical results across operating systems.
- `conda-forge` is the primary channel, as it provides pre-built binaries for all listed packages on both macOS and Windows.
- Do **not** update individual packages without coordinating with both machines — even a minor version bump (e.g. numpy 2.0.2 → 2.1.0) can shift floating-point results in feature detection.
