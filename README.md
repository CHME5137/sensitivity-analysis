# Sensitivity analysis

Notebooks for CHME 5137 (2026): sensitivity analysis with SALib, after Lecture 15.

| Notebook | What it is |
|---|---|
| `cstr-sensitivity.ipynb` | The sensitivity homework: local, Morris and Sobol analyses of a cooled CSTR running A ⇌ B |

Each notebook has a paired `.py` file holding the same cells as plain text (see
**Jupytext** below).

## Getting started

1. **Fork** this repository to your own GitHub account (button at top right).
2. **Clone your fork**, not this one:
   ```bash
   git clone git@github.com:YOUR-USERNAME/sensitivity-analysis.git   # SSH, e.g. on Explorer
   # or
   git clone https://github.com/YOUR-USERNAME/sensitivity-analysis.git
   ```
   On Explorer, clone it into `/courses/CHME5137.202710/students/$USER/`.
3. **Install SALib 1.6** into the environment your kernel uses, if you haven't already.
   On Explorer, in the course venv, from a terminal on a compute node:
   ```bash
   source ~/venvs/chme5137/bin/activate
   pip install --upgrade SALib ipywidgets
   ```
   With Anaconda on a laptop: `conda install -c conda-forge salib ipywidgets`. Or, either
   way, run `%pip install --upgrade SALib` in a notebook cell, then restart the kernel.
4. Open `cstr-sensitivity.ipynb` and work through it, using your own kernel (`CHME5137` on
   Explorer). Commit and push as you go:
   ```bash
   git add cstr-sensitivity.ipynb cstr-sensitivity.py
   git commit -m "Task 1: base case and coolant temperature sweep"
   git push
   ```
   This is your own fork, so committing on `main` is fine. Don't open a pull request back
   to `CHME5137/sensitivity-analysis`: that would publish your answers to everyone.

## Jupytext

`jupytext.toml` pairs each notebook with a `.py` script in the *percent* format, as in the
`differential-equations` and `parameter-estimation` repositories. Saving the notebook updates
the `.py` as well. If the two get out of step, `jupytext --sync cstr-sensitivity.ipynb` brings
them back in line.

---

*Written with Claude Code (Claude Opus 5.5) and checked by Prof. West.*
