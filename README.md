

### This code is associated with the paper from Tian et al., "Allosteric mechanism of the V. vulnificus adenine riboswitch resolved by four-dimensional chemical mapping". eLife, 2017. http://dx.doi.org/10.7554/eLife.29602


# M2R Simulations

*RFAM simulations of mutate-map-rescue quartets and helix frequency estimates*


## Installation

1. Need `Biers` (and maybe `HiTRACE` for some util scripts)

> Which also requires local `RNAstructure` (likely `5.6`) running. In particular, we need `partition` (default) or `allsub`.

2. Clone repo and add subfolders to path in MATLAB


## Examples

#### 1. _in silico_ RFAM M2R quartet simulations

See `examples/simulate_RFAM.m`. Initial digested RFAM sequences are provided in `examples/simulate_RFAM.mat`. Previous run result of this job (~330 RFAM families) are in dropbox:

- RFAM raw results [0811](https://www.dropbox.com/s/yitfgaovuhi1tjx/add_20160811_RFAM_bpp.mat?dl=0)
- RFAM summary [0831](https://www.dropbox.com/s/bqynunh1j6itwlh/add_20160831_RFAM_helix_partial.mat?dl=0) (loaded from 0811 and saved in separate)

#### 2. MLE (maximum likelihood estimates) summary of RFAM results

See `examples/calc_MLE.m`. Need to load from RFAM summary from step 1. This step produces the estimator of helix frequency ratio vs. in silico bpp; both based on the above simulated RFAM data.

- MLE result [0716](https://www.dropbox.com/s/i57bsxbajapxk2s/add_20170716_RFAM_MLE_95.mat?dl=0) (used by M2R scores for experiments)

#### 3. Helix frequency of experimental M2R data based on MLE

See `examples/estimate_M2R.m`. Related MATALAB sessions are:

- MLE result (see above)
- add M2R + LM2R exp. quartets [0811](https://www.dropbox.com/s/urdya0o6a4vcl3o/add_20160811_a140T_score_cleanup_more.mat?dl=0) and [1014](https://www.dropbox.com/s/ylf7yk2lc0b3b7x/add_20161014_3Dlock_quartet.mat?dl=0)


## Notes:

> All add riboswitch related MATLAB scripts and sessions are available in [dropbox](https://www.dropbox.com/sh/7fir8bn7hgh72kp/AADDUC_0aZdsascQKKKWr1Hca?dl=0); in case digging out data and commands are needed.
