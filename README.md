# Simplicity of confinement in SU(3) Yang-Mills theory — Monte Carlo Code

# Program for simulation of lattice Yang-Mills theories

## Contents
- [Simplicity of confinement in SU(3) Yang-Mills theory — Monte Carlo Code](#simplicity-of-confinement-in-su3-yang-mills-theory--monte-carlo-code)
- [Program for simulation of lattice Yang-Mills theories](#program-for-simulation-of-lattice-yang-mills-theories)
  - [Contents](#contents)
  - [Installation](#installation)
  - [Action](#action)
  - [Configuration parameters](#configuration-parameters)
  - [Input file and some conventions](#input-file-and-some-conventions)
  - [A note on topological charge](#a-note-on-topological-charge)
  - [Abelian projection for SU(3) lattice gauge theory](#abelian-projection-for-su3-lattice-gauge-theory)

## Installation

This code base is designed to operate with a linux or compatible distribution. Detailed installation instructions are included in [`INSTALL`](INSTALL).

To configure the program, use

    ./configure

See [configuration options](#configuration-parameters) below.

To build the program, use

    make

and to clean the build, accordingly, use

    make clean

To clean the entire program distribution for reconfiguration, use

    make distclean


## Action
The action used is the standard Wilson action, which for the $N$ colors can be generically written as the sum of

$$-\beta \frac{1}{N} \sum_{\text{plaquettes}} Re(Tr(\text{plaquette}))$$

and possibly a theta term. If `--enable-imaginary-theta` is used during configuration an imaginary theta term is used of the form

$$-\theta_{im} Q$$

where the simplest discretisation with definite parity of the topological
charge $Q$ is used (see later [A note on the topological charge](#a-note-on-topological-charge) for more details).

In `yang_mills_tracedef`, a trace deformation is also added to the action of the form

$$\sum_{\text{sites on a spatial slice}} \sum_{i=0}^{N/2} h_i|Tr(\text{polyakov}^{i+1})|^2$$

The update is performed by means of heatbath and overrelaxation updates (and in some case also Metropolis updates) implemented a la Cabibbo-Marinari.

In `yang_mills_local_fundadj`, the action is

$$-\beta \frac{1}{N} \sum_{\text{plaquettes}} Re(Tr(\text{plaquette})) 
-\beta_{\text{adj}} \frac{1}{N^2 -1} \sum_{\text{plaquettes}} Tr_{\text{adj}}(\text{plaquette})$$

and $Tr_\text{adj}=|Tr_\text{fund}|^2-1$.

In `yang_mills_higgs`, the action is 

$$-\beta \frac{1}{N} \sum_{\text{plaquettes}} Re(Tr(\text{plaquette})) 
-N_{\text{higgs}} \beta_\text{higgs} \sum_{r, \mu} \sum_{i=1}^{N_\text{higgs}} Re([H_r^{(i)}]^{\dag} U_{r, \mu} H_{r+\mu}^{(i)} )$$

where the sum on $\mu$ is just on positive orientations and Higgs fields are normalized according to

$$ \sum_{i=1}^{N_\text{higgs}} | H_r^{(i)} |^2 = 1 $$

for every site $r$.


## Configuration parameters

The choice of the gauge group happens at configuration time through the macros `Gauge_Group` and `N_c`. Two classes of gauge groups are implemented, the `SoN` and
`SuN` classes (with `N_c=1` for `SuN` corresponding to the $U(1)$ gauge theory). As an example, to simulate the $SU(5)$ gauge theory the program has to be configured using

    ./configure Gauge_Group=SuN N_c=5

Relevant configure options are

- `--enable-use-openmp`     to enable the use of openmp
- `--enable-use-theta`      to enable the use of an imaginary theta term

and the following macro are available

- `Gauge_Group`  the gauge group (`SuN` or `SoN`, default `SuN`)
- `N_c`          the number of colors (default 2)
- `Num_levels`   the number of levels in multilevel (default 1)
- `Num_threads`  number of threads to be used in OpenMP (default 1)
- `ST_dim`       spacetime dimensionality (default 4)
- `N_higgs`      the number of Higgs fields (default 1)


After the configuration, the compilation is performed as usual by

    make

Defining the macro `DEBUG` several sanity checks are activated, while the macros `OPT_MULTIHIT` and `OPT_MULTILEVEL` can be used to optimize the number of multihit steps and the number of updates to be used in the multilevel algorithm.

## Input file and some conventions

A template input file `template_input.in` is created when calling executables without input file and everything following `#` (up to carriage return) in the input file is interpreted as a comment.

Some conventions used in the code are the following:

- time is direction 0 and the numbers following "size" in the imput are in the order "time space1 space2 ..."

- when computing Polyakov loop correlators, Polyakov loops are separated along direction 1 and, when computing flux tube profiles, the transverse direction is direction 2. 

- Level 0 of the multilevel is the one corresponding to the largest timeslice.

## A note on topological charge

Topology for the general group: in the continuum the topological charge is defined in 4 space-time dimensions by

$$Q=\frac{1}{64 \pi^2} \int \epsilon_{\mu\nu\rho\sigma} 
   F_{\mu\nu}^a F_{\rho\sigma}^a d^4 x$$

with the normalization $Tr(T_iT_j)=K \delta_{ij}$ such that the longest root of
the representation is normalized to 1. See [C. W. Bernard, N. H. Christ, A. H. Guth, E. J. Weinberg Phys. Rev. D 16, 2967 (1977)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.16.2967).

For $SU(N)$, this base is that of the generalized Gell-Mann matrices divided by
two and in this representation $Tr(T^aT^b)=(1/2)\delta^{ab}$, so that 

$$F_{\mu\nu}^aF_{\rho\sigma}^a=2Tr(F_{\mu\nu}F_{\rho\sigma})$$

and thus 

$$Q=\frac{1}{32 \pi^2} \int \epsilon_{\mu\nu\rho\sigma} 
Tr(F_{\mu\nu} F_{\rho\sigma}) d^4 x$$

The discretisation used on the lattice is obtained by using the clover form of the discretised field-strength, i.e. the clover $Q_{\mu\nu}$ is given by

$$Q_{\mu\nu}  \sim 4 + 4 i a^2 F_{\mu\nu}$$

and thus

$$F_{\mu\nu} = \frac{1}{8i}( Q_{\mu\nu} - Q_{\mu\nu}^{\dag} )$$

By using this expression we get

$$Tr(F_{\mu\nu}F_{\rho\sigma}) = -\frac{1}{2^5} ReTr(Q_{\mu\nu}[Q_{\rho\sigma}-Q_{\rho\sigma}^{\dag}])$$

We now note that of the 24 terms
$\epsilon_{\mu\nu\rho\sigma}F_{\mu\nu}F_{\rho\sigma}$ only 3 are independent:

- single exchange on the first term   $F_{\mu\nu}F_{\rho\sigma} \to F_{\nu\mu}F_{\rho\sigma}$
- single exchange on the second term  $F_{\mu\nu}F_{\rho\sigma} \to F_{\mu\nu}F_{\sigma\rho}$
- exchange $F_{\mu\nu}F_{\rho\sigma} \to F_{\rho\sigma}F_{\mu\nu}$ 
- remaining permutations of indices: e.g 0123, 0213, 0312

thus $$\frac{1}{32} \epsilon_{\mu\nu\rho\sigma}Tr(F_{\mu\nu}F_{\rho\sigma}) \\ =-\frac{1}{32} \cdot \frac{1}{2^5} \cdot 2^3 \cdot [ \text{sum on independent permutations of } \epsilon_{\mu\nu\rho\sigma} ReTr(Q_{\mu\nu}[Q_{\rho\sigma}-Q_{\rho\sigma}^{\dag}]) ] \\ =-\frac{1}{128} [ \text{sum on independent permutations of } \epsilon_{\mu\nu\rho\sigma} ReTr(Q_{\mu\nu}[Q_{\rho\sigma}-Q_{\rho\sigma}^{\dag}]) ]$$

This is the expression used in the code.

## Abelian projection for SU(3) lattice gauge theory

A small section of the codebase has been modified to facilitate integration with an analysis pipeline specified by the repository [Simplicity of confinement in SU(3) Yang-Mills — Analysis Code](https://github.com/xavier-crean/su3_mon):
- [`src/yang_mills_local.c`](src/yang_mills_local.c) has been modified to include the function [`perform_measure_u1subg()`](src/yang_mills_local.c#L23) which calls a function [`U1_extract()`](src/yang_mills_local.c#L56).
- [`lib/gauge_conf_meas.c`](lib/gauge_conf_meas.c) has been to modified so that [`U1_extract()`](lib/gauge_conf_meas.c#L1459) prints the two $U(1)$-sector configurations to respective binary data files.
- respective header files have been modified accordingly.
---

**N.B.** This README file has been modified and expanded from an original plaintext version $-$ see [`AUTHORS`](AUTHORS).
