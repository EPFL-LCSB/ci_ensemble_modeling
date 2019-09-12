Statistical inference in ensemble modeling of cellular metabolism

By Tuure Hameri, Marc-Olivier Boldi*, Vassily Hatzimanikatis*

*corresponding authors: marc-olivier.boldi@unil.ch (MOB), vassily.hatzimanikatis@epfl.ch (VH)

This code is an early release. You will need Matlab 2016a to run it. Other versions of Matlab should work but have not been tested.

You will need to have Git LFS in order to properly download some binary files:

git clone https://github.com/EPFL-LCSB/ci_ensemble_modeling.git /path/to/ci_ensemble_modeling

cd /path/to/ci_ensemble_modeling

git lfs install

git lfs pull

The main run file is run_CI_EM_paper.m and can be run to reproduce figures presented in the paper. 

Data (flux control coefficients) used to demonstrate construction of confidence intervals in this paper are contained in /path/to/ci_ensemble_modeling/rawData. This data was obtained from the publication:

Kinetic models of metabolism that consider alternative steady-state solutions of intracellular fluxes and concentrations

Metabolic Engineering

Volume 52, March 2019, Pages 29-41

Tuure Hameri, Georgios Fengos, Meric Ataman, Ljubisa Miskovic, Vassily Hatzimanikatis

https://doi.org/10.1016/j.ymben.2018.10.005.

The software in this repository is put under an APACHE-2.0 licensing scheme - please see the LICENSE file for more details.