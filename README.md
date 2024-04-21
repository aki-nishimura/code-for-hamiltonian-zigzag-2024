
# Code for "Zigzag path connects two Monte Carlo samplers: Hamiltonian counterpart to a piecewise deterministic Markov process" by Nishimura et al. 2024

This repository contains instruction and scripts to reproduce the results in the article "Zigzag path connects two Monte Carlo samplers: Hamiltonian counterpart to a piecewise deterministic Markov process" by Nishimura et al. 2024.


## Setting up BEAST

The simulation results of Section 4 uses the implementation of Hamiltonian zigzag in the [BEAST](https://beast.community/) software package. 
On Unix systems, you can use the commands below to install BEAST.
The `ant` command requires Apache Ant, which should be available through a package manger such as Homebrew on macOS.

```
git clone https://github.com/beast-dev/beast-mcmc.git
cd beast-mcmc
git checkout hmc_develop
ant
```

Under `beast-mcmc/build/dist/`, you should now find the `beast.jar` application file for running BEAST on a model as specified by an input XML file.


## Reproducing numerical results using BEAST

Having built BEAST, you can run the software on the compound symmetric targets of Sec 4.3 and on the phylogenetic probit posterior of Sec 4.4 by inputting corresponding XML files. 
The command for this is the form

```
java -jar ${path_to_beast_jar} -seed ${seed} ${path_to_input_xml}
```

where `path_to_beast_jar` and `path_to_input_xml` are variables specifying paths to the `beast.jar` and XML file and `seed` specifying the random seed for the Markov chain. 
For example, if you cloned this and `beast-mcmc` repositories to your home directory, you would run zigzag NUTS on the phylogenetic probit posterior with command

```
java -jar ~/beast-mcmc/build/dist/beast.jar -seed 111 ~/code-for-hamiltonian-zigzag-2024/ZNUTS_phylo_application.xml
```

The corresponding BEAST run generates a text file `ZNUTS_phylo_application_samples.log` which contains zigzag NUTS samples from the target.
(The log file is updated dynamically as the sampling proceeds;
given the high-dimensionality of the target, the sampling takes hours to complete.) 
Similarly, to run zigzag HMC with relative integration time (or "step size") of 0.707 on the 256-dimensional compound symmetric target with correlation 0.9, you would use a command

```
java -jar ~/beast-mcmc/build/dist/beast.jar -seed 111 ~/code-for-hamiltonian-zigzag-2024/cs_example/HZZ_cs_rho0.9_d256_ss0.707.xml
```

which generate a text file `HZZ_cs_rho0.9_d256_ss0.707_samples.log` which contains corresponding zigzag HMC samples.


## Reproducing supplemental numerical results using hdtg

The `supplement` folder provides R scripts to reproduce the numerical results in Supplement Section S6, S7, S11, and S12.
These scripts use a Hamiltonian and Markovian zigzag implementations in the R package [hdtg](https://github.com/suchard-group/hdtg).
You can install the package from an R session via

```
remotes::install_github("https://github.com/suchard-group/hdtg", build = FALSE)
```
<!-- See https://github.com/suchard-group/hdtg/issues/10 for why `build = FALSE` is needed -->

The results for Section S7 and S12 also require the parameters of the phylogenetic probit posterior.
These are provided by the files `precision_matrix.csv`, `orthant_indicator.txt`, and `mean.txt` available on Zenodo at http://doi.org/10.5281/zenodo.4679720.
Place these files under the `supplement/data` folder within this Git repo.
The folder also contains `marginal_variances.txt` and `principal_component.txt`, precomputed from the precision matrix to save time and used by the simulation scripts.

The results of Section S6 can be reproduced by running the script `supplement/compare_zigzag_against_rejection_sampler.R`.

To reproduce the results of Section S7 ("Comparison of zigzag HMC/NUTS with existing samplers"), run 
```
Rscript supplement/compare_zigzag_against_harmonic_hmc.R \
    $sampler_type $sampler_seed $harmonic_integ_time
```
where the variable `sampler_type` takes a value of either `"zigzag"` (zigzag NUTS) or `"harmonic"` (harmonic HMC of Pakman and Paninski, 2014).
The `sampler_seed` takes any integer; `1`, `2`, `3`, `4`, and `5` are used for the results in the article.
The `harmonic_integ_time` is used only when the sampler type is set to `"harmonic"` and takes either `"short"`, `"medium"`, and `"long"`, corresponding to the integration time of `[pi / 8, pi / 4]`, `[3 pi / 16, 3 pi / 8]`, and `[pi / 4, pi / 2]`.
For the result on the method of Botev (2017), run
```
Rscript supplement/benchmark_botev_method.R
```

To reproduce the results of Section S11 ("Two zigzags' performance on rotated versions of compound symmetric targets"), run
```
Rscript supplement/benchmark_zigzag_on_rotated_cs.R \
    $target_seed $sampler_seed $rho $n_dim $stepsize_multiplier $sampler_type 
```
where the variable `target_seed` specifies the seed for random rotation of the compound symmetric target; 
`"null"` means no rotation and `111`, `222`, are `333` used for the results in the article.
The `rho` specifies the correlation parameter and takes either `0.9` and `0.99`.
The `sampler_type` takes either `"markovian"` (Markovian zigzag), `"nuts"` (zigzag NUTS), and `"hmc"` (zigzag HMC).
The `stepsize_multiplier` is ignored when the sampler is set to Markovian zigzag but otherwise scales the relative integration time;
it takes either `0.1` or `1` <!-- labeled as "shorter" and "longer" in corresponding output files --> for zigzag NUTS
and either `0.5`, `0.7071`, `1`, `1.414`, or `2` <!-- labeled as "shortest", "shorter", "medium", "longer", and "longest" in corresponding output files --> for zigzag HMC.


To reproduce the results of S12 ("Preconditioning target under Hamiltonian zigzag"), run
```
Rscript supplement/compare_zigzag_with_and_without_precond.R \
    $precond_type $sampler_seed
```
where the variable `precond_type` takes a value of either `"none"` (no preconditioning), `"cond_var"` (diagonal preconditioning by conditional variance), and `"marg_var"` (diagonal preconditioning by marginal variance).

Effective sample sizes and algorithms' run times are saved to files named `*_summary_seed_${sampler_seed}.rds` and MCMC samples in `*_samples_seed_${sampler_seed}.rds`.
The default location for these outputs is `supplement/results/` but you can change it to any other location by editing the path at the top of `supplement/helper.R`;
you can also change the path for the input phylogenetic data there.
