# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.5
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Sensitivity analysis of a cooled CSTR
#
# **CHME 5137, homework after Lecture 15 (October 30, 2026).** In class we ran SALib on
# someone else's model, the electric-car battery. Here you run the whole workflow yourself,
# on a reactor.
#
# *Written for 2026 by Claude Code (Claude Opus 5.5) and checked by Prof. West.*
#
# ---
#
# A jacketed CSTR runs a liquid-phase, reversible, exothermic reaction
# $$\mathrm{A} \rightleftharpoons \mathrm{B}, \qquad r = k_f C_A - k_r C_B,$$
# first order both ways. A cooling jacket holds it near the temperature you want. Two things
# compete:
#
# - **Kinetics.** $k_f = A\exp(-E_a/RT)$ rises steeply with temperature, so a hotter reactor
#   gets closer to equilibrium in the time the liquid spends in it.
# - **Equilibrium.** The reaction is exothermic, so the equilibrium constant
#   $K = k_f/k_r = \exp(\Delta S/R - \Delta H/RT)$ *falls* as the temperature rises, and so
#   does the conversion you could ever reach.
#
# So there's a best temperature, and how close you get to it depends on the energy balance:
# the heat the reaction releases, the feed temperature, and how hard the jacket cools.
#
# The model is written for you, in the next cell. **Read it before you start**: every task
# depends on knowing what goes in and what comes out. It has 12 inputs. Some are properties
# of the chemistry that you only know roughly (rate parameters, thermodynamics, heat-transfer
# coefficient). Others are design and operating choices (residence time, feed and coolant
# temperatures, jacket area, feed concentration).
#
# The coolant is water at atmospheric pressure. The solvent boils well above any temperature
# the reactor reaches, so the reactor itself never boils.
#
# ### How to work
#
# - **Use SALib 1.6** (the newest). The first cell checks. Installation is explained at the top
#   of the Lecture 15 notebook.
# - **Follow the Lecture 15 notebook.** Every SALib call you need is there: the `problem`
#   dictionary, `sample`, running the model, `analyze`, the plots.
# - **Use `seed=2026`** wherever SALib takes a seed, so that your numbers can be compared with
#   other people's.
# - **Commit as you go.** After each numbered task works, Restart & Run All, then commit with
#   a message that says what you did, and push. **Commit the `.ipynb` and its paired `.py`
#   together.**
# - **Answer the questions in words**, in the markdown cells marked *Answer*. A sentence or
#   two each.

# %%
from importlib.metadata import version
import time
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['figure.figsize'] = [5, 3]
print('SALib', version('SALib'))

# %% [markdown]
# ## The model

# %%
R = 8.314  # gas constant, J/mol/K


def cstr_steady_state(log10_A=9.1, Ea=80e3, dH=-60e3, dS=-146.0,
                      tau=600.0, T_feed=310.0, T_cool=365.0,
                      U=500.0, area=6.0, C_feed=3000.0,
                      rho_cp=1.7e6, W_stir=1000.0):
    """
    Steady state of a jacketed CSTR running the reversible, exothermic,
    first-order liquid-phase reaction A <=> B.

    Arguments, all in SI units:
        log10_A  log10 of the forward pre-exponential factor (A in 1/s)
        Ea       forward activation energy, J/mol
        dH       enthalpy of reaction, J/mol (negative: exothermic)
        dS       entropy of reaction, J/mol/K
        tau      residence time V/q, s
        T_feed   feed temperature, K
        T_cool   coolant temperature in the jacket, K (the jacket is well mixed)
        U        overall heat-transfer coefficient, W/m2/K
        area     heat-transfer area of the jacket, m2
        C_feed   concentration of A in the feed, mol/m3 (no B in the feed)
        rho_cp   volumetric heat capacity of the liquid, J/m3/K
        W_stir   power the stirrer puts into the liquid, W

    Returns:
        X        conversion of A, between 0 and 1
        T        reactor temperature, K
    """
    V = 1.0          # reactor volume, m3
    q = V / tau      # volumetric flow rate, m3/s

    # We don't know the reactor temperature yet, so try every temperature
    # from 250 K to 700 K, in steps of 0.1 K, all at once (NumPy arrays).
    T = np.linspace(250.0, 700.0, 4501)

    # Kinetics: the forward rate constant, from Arrhenius.
    k_f = 10**log10_A * np.exp(-Ea / (R * T))
    # Thermodynamics: K = k_f / k_r, and -RT ln K = dG = dH - T dS.
    K = np.exp(dS / R - dH / (R * T))
    k_r = k_f / K

    # Mole balance on A at steady state, in = out + reacted:
    #     q C_feed = q C_A + V (k_f C_A - k_r C_B),
    # with C_A = C_feed (1 - X) and C_B = C_feed X. Solved for X:
    X = k_f * tau / (1 + (k_f + k_r) * tau)

    # Energy balance at steady state, in watts. Heat comes from the reaction and
    # the stirrer; it goes into warming the feed up to T, and into the jacket.
    heat_released = -dH * q * C_feed * X + W_stir
    heat_removed = q * rho_cp * (T - T_feed) + U * area * (T - T_cool)
    imbalance = heat_released - heat_removed

    # A steady state is a temperature where the imbalance is zero. There can be
    # more than one (ignition and extinction). Take the lowest: it's the one the
    # reactor settles at when it's started up cold.
    crossings = np.nonzero(np.diff(np.sign(imbalance)))[0]
    i = crossings[0]
    # Interpolate between the grid points either side of that crossing.
    fraction = imbalance[i] / (imbalance[i] - imbalance[i + 1])
    T_ss = T[i] + fraction * (T[i + 1] - T[i])
    X_ss = X[i] + fraction * (X[i + 1] - X[i])
    return X_ss, T_ss


# %% [markdown]
# The inputs, their base-case values, and the ranges to explore:
#
# | Input | Base | Range | Why that range |
# |---|---|---|---|
# | `log10_A` | 9.1 | 8.8 – 9.4 | *task 1(c)* |
# | `Ea` (kJ/mol) | 80 | 78 – 82 | Fitted from rate data, ± 2 kJ/mol |
# | `dH` (kJ/mol) | −60 | −65 – −55 | Calorimetry, ± 5 kJ/mol |
# | `dS` (J/mol/K) | −146 | −153 – −139 | Estimated by group additivity |
# | `tau` (s) | 600 | 300 – 1000 | The flow rate the pumps can deliver |
# | `T_feed` (K) | 310 | 290 – 330 | The feed tank is outdoors: winter to summer |
# | `T_cool` (K) | 365 | 320 – 373 | *task 1(c)* |
# | `U` (W/m²/K) | 500 | 300 – 700 | From a correlation, ± 40% |
# | `area` (m²) | 6 | 4 – 8 | A design choice: the jacket size |
# | `C_feed` (mol/m³) | 3000 | 2500 – 3500 | The feed specification |
# | `rho_cp` (MJ/m³/K) | 1.7 | 1.6 – 1.8 | The solvent's density × heat capacity |
# | `W_stir` (kW) | 1 | 0.5 – 2 | The stirrer motor |
#
# In the code, everything is in SI units (J/mol, J/m³/K, W).

# %%
problem = {
    'num_vars': 12,
    'names': ['log10_A', 'Ea', 'dH', 'dS', 'tau', 'T_feed', 'T_cool',
              'U', 'area', 'C_feed', 'rho_cp', 'W_stir'],
    'bounds': [[8.8, 9.4],         # log10_A: see task 1(c)
               [78e3, 82e3],       # Ea, J/mol
               [-65e3, -55e3],     # dH, J/mol
               [-153.0, -139.0],   # dS, J/mol/K
               [300.0, 1000.0],    # tau, s
               [290.0, 330.0],     # T_feed, K
               [320.0, 373.0],     # T_cool, K: see task 1(c)
               [300.0, 700.0],     # U, W/m2/K
               [4.0, 8.0],         # area, m2
               [2500.0, 3500.0],   # C_feed, mol/m3
               [1.6e6, 1.8e6],     # rho_cp, J/m3/K
               [500.0, 2000.0]],   # W_stir, W
}
base = dict(log10_A=9.1, Ea=80e3, dH=-60e3, dS=-146.0, tau=600.0, T_feed=310.0,
            T_cool=365.0, U=500.0, area=6.0, C_feed=3000.0, rho_cp=1.7e6, W_stir=1000.0)

# %% [markdown]
# ## Task 1. The base case, one input at a time
#
# **(a)** Run the model at the base case, and print the conversion and the reactor
# temperature. (`cstr_steady_state(**base)` passes the dictionary as keyword arguments.)

# %%

# %% [markdown]
# **(b)** Vary the coolant temperature alone, from 320 to 400 K, with everything else at
# its base value. Plot the conversion, and the reactor temperature, against `T_cool`. Mark
# 373 K. Where is the best coolant temperature? Why does the conversion rise, and then
# fall?

# %%

# %% [markdown]
# **Answer:**

# %% [markdown]
# **(c)** Two ranges in the table have no reason given: `log10_A` and `T_cool`. Write one
# for each. Why might the upper limit on `T_cool` be a hard limit, not just an uncertainty?

# %% [markdown]
# **Answer:**

# %% [markdown]
# ## Task 2. Running the model on a sample, and local sensitivities
#
# **(a)** The vehicle-to-grid model in class took arrays, so `model(*sample.T)` ran every
# row at once. This one doesn't: it takes one value of each input, and returns one $X$ and
# one $T$. (Why not? Look at what it does with `T`.)
#
# Write a function `run_model(sample)` that takes a SALib sample, with one row per run and one
# column per input in the order of `problem['names']`, and returns two arrays, $X$ and $T$, with
# one value per row. A loop is fine. Test it on the first five rows of a Sobol sample, and time
# it on 1000 rows. Then use your timing to estimate how long 10 000 and 100 000 runs would
# take.

# %%
from SALib.sample import sobol as sobol_sample

def run_model(sample):
    """X and T for each row of a SALib sample."""
    


# your test and timing here


# %% [markdown]
# **Answer:**

# %% [markdown]
# **(b) Local sensitivity.** At the base case, compute the normalised local sensitivity
# $$S_i^\mathrm{local} = \frac{\partial \ln X}{\partial \ln p_i}
# \approx \frac{\ln X(p_i(1+\delta)) - \ln X(p_i(1-\delta))}{2\delta}$$
# of the conversion to each input $p_i$, by a central finite difference with $\delta = 0.01$
# (a 1% change). It's the percentage change in $X$ for a 1% change in the input. Print a table,
# largest magnitude first. By this measure, which inputs matter most, and which hardly at all?
# What does it say about `T_cool`?
#
# *(For `log10_A`, a 1% change in log₁₀ A is a factor of about 1.23 in A. That's fine; just
# note it.)*

# %%

# %% [markdown]
# **Answer:**

# %% [markdown]
# ## Task 3. Morris screening
#
# Run a Morris analysis of the conversion, as in class: 100 trajectories, 4 levels,
# `seed=2026`. How many model runs is that? Print μ, μ\* and σ for each input, and make the
# two Morris plots (`horizontal_bar_plot` and `covariance_plot`).
#
# Which inputs come out on top, and which could you probably ignore? Which have σ larger
# than μ\*, and what does that tell you?

# %%
from SALib.sample import morris as morris_sample
from SALib.analyze import morris as morris_analyze
from SALib.plotting import morris as morris_plot

# your code here


# %% [markdown]
# **Answer:**

# %% [markdown]
# ## Task 4. Sobol indices
#
# **(a)** Compute the first-order, total and second-order Sobol indices of the conversion, with
# $N = 1024$ (a power of 2), `calc_second_order=True` and `seed=2026`. How many model runs is
# that? Print $S_1$, $S_T$ and the 95% confidence interval on $S_T$ for each input, largest
# $S_T$ first.

# %%
from SALib.analyze import sobol as sobol_analyze

# your code here


# %% [markdown]
# **(b)** Answer in words:
#
# 1. Rank the inputs by $S_T$. Does it agree with the Morris ranking from task 3? Morris used
#    far fewer runs: was it good enough to tell you which inputs matter?
# 2. Which input has the biggest gap between $S_T$ and $S_1$? What does that gap mean?
# 3. Compare with your local sensitivities from task 2(b). Pick the input where local and
#    global disagree most, and explain why. (Look back at your plot from task 1(b).)
# 4. Find the largest second-order index, `Si['S2']`, and its interval, `Si['S2_conf']`. Is it
#    clearly bigger than zero? Does that pair interacting make physical sense?

# %%

# %% [markdown]
# **Answer:**

# %% [markdown]
# ## Task 5. Fixing the unimportant inputs
#
# The point of factor fixing is a simpler model that gives the same answers. Choose the inputs
# with $S_T < 0.03$, and fix them at their base values. Then **check** it:
#
# 1. Take your Sobol sample from task 4, set the fixed columns to their base values, and run
#    the model again on the same rows. (Copy the array first: `reduced = sobol_X.copy()`.)
# 2. Plot the reduced model's conversion against the full model's, one point per row, with
#    the line $y = x$.
# 3. What fraction of the variance in $X$ does the reduced model keep? Use
#    $1 - \sum (X_\mathrm{full} - X_\mathrm{reduced})^2 / \sum (X_\mathrm{full} - \bar X_\mathrm{full})^2$.
#    Compare it with $1 - \sum S_T$ of the inputs you fixed.
# 4. What's the largest difference for any single row? Is the reduced model good enough to
#    use?

# %%

# %% [markdown]
# **Answer:**

# %% [markdown]
# ## Task 6. What should you measure better? It depends how you run the reactor
#
# You've designed the reactor, so now the design and operating choices are fixed at their base
# values. What's left uncertain is the chemistry and the heat transfer: `log10_A`, `Ea`, `dH`,
# `dS`, `U` and `rho_cp`. You have money for one set of measurements.
#
# **(a)** Make a new `problem` with just those six inputs (same ranges). For three coolant
# temperatures, 325 K, 350 K and 373 K, compute their total Sobol indices for the conversion
# ($N = 1024$, first and total order only, `seed=2026`), with every other input at its base
# value. Print a table: one row per input, one column per coolant temperature.
#
# *(Hint: write a small function that takes a row of six values and the coolant temperature,
# builds the full set of keyword arguments from `base`, and calls `cstr_steady_state`.)*

# %%

# %% [markdown]
# **(b)** Which property would you spend the money measuring if you run the reactor with cold
# coolant? With coolant near boiling? Explain the change using kinetics and equilibrium. Is there
# anything you'd never bother measuring?

# %% [markdown]
# **Answer:**

# %% [markdown]
# **(c) Think about your answer.** From task 1, the best coolant temperature is only a few
# kelvin below the boiling point of water. Would you run the reactor there? Give at least two
# reasons for caution, including something the model doesn't know about. How much conversion
# would it cost you to back off by 10 K? (Your plot from task 1(b) has the numbers.)

# %% [markdown]
# **Answer:**

# %% [markdown]
# ---
#
# ### Stretch (optional)
#
# - **Scatter plots.** Plot $X$ against each input for the Sobol sample, as in class. Which
#   inputs show a clear trend? Can you see `T_cool`'s optimum?
# - **Reactor temperature.** Repeat task 4 for the reactor temperature (`sobol_T`) instead of
#   the conversion. Which inputs control the temperature, and how is that different?
# - **Convergence.** Repeat task 4(a) with $N$ = 128, 256 and 512. How much do the $S_T$ values
#   and their intervals change? What $N$ would you have been happy with?
# - **More than one steady state.** The model takes the lowest steady state when there's more than
#   one. For how many rows of your Sobol sample is there more than one? Where in the input space
#   are they? How does that help explain the largest difference you found in task 5?
