## Introduction
This [Julia](https://julialang.org) and [Maple](https://maplesoft.com/products/Maple) code accompanies the research article "*Numerical evaluation of the transition probability of the simple birth-and-death process*" by Pessia A. and Tang J. (2019).
A pre-print version of the article can be found on [arXiv](https://arxiv.org/abs/1909.10765).

## Setup

First of all, you obviously need a working Julia installation.
Download and install the latest stable version from the official Julia [website](https://julialang.org/downloads).

You now need to install package *SimpleBirthDeathProcess* and third-party packages *KahanSummation*, *SpecialFunctions*, and *Plots*.
Use the following commands in a Julia environment to install all required packages:

- Enter the *Pkg* mode by pressing  `]`
- Issue the command `add KahanSummation SpecialFunctions Plots`
- Issue the command `add https://github.com/albertopessia/SimpleBirthDeathProcess.jl`
- Press the *Backspace* key to return to the normal Julia prompt

## How to use the code

To generate numerically accurate values in arbitrary precision you need a working installation of *Maple* 2018.
Simply issue the command `read "log_trans_prob.mpl"`.

Figures 1, 3, and 4 can be recreated by copy/pasting the corresponding source code.
Published figures were edited and finalized with [Inkscape](https://inkscape.org).

Values in Tables 1-4 can be obtained by using code found in `Tables.jl`.
Note that this particular operation requires very long computation times and should only be tried in parallel on cluster servers.
