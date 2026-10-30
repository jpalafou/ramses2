The two namelist blocks `&CORE_PARAMS` and `&CORE_ACCRETION_PARAMS` control the core particle module. Cores are sink-like particles that form in dense gas clumps found by the PHEW clump finder, grow by flux accretion, can merge, and turn their mass into star particles. Cores are switched on at runtime with `core=.true.`; there is no compile-time flag.

Cores need a build with `NDIM=3` and `HYDRO=1`. On GPU, cores run on a single MPI rank, without MHD and with `NENER=0`, and core formation also needs `rho_type_core=4` and `core_saddle_threshold<=0`. The GPU does not write the core clump outputs yet: `output_core_clump` and `output_peak_core`, and with core formation `output_peak_grid`, `_part`, `_star`, `_sink` and `_tree`; support is planned. These combinations stop the code when the namelist is read. Core formation also needs `GRAV=1`: without it, `core_form` has no effect.

## `&CORE_PARAMS`: formation, merging and conversion into stars

| Variable name                          | Fortran type | Default value | Description                                                                                         |
|:---------------------------------------|:-------------|:--------------|:----------------------------------------------------------------------------------------------------|
| `core`                                 | `logical`    | `.false.`     | Turn the core particle module on or off.                                                            |
| `ncoremax`                             | `integer`    | 0             | Maximum number of core particles per MPI process. If 0, set to `ncoretot` divided by the number of processes; if that is also 0, cores are turned off. |
| `ncoretot`                             | `integer(8)` | 0             | Total maximum number of core particles across all MPI processes.                                    |
| `rho_type_core`                        | `integer`    | 4             | Fluid used to build the density field for core formation (1: DM, 2: stars, 3: sinks, 4: gas).       |
| `core_form`                            | `logical`    | `.false.`     | Activate core particle formation. Needs `core_seed_mass > 0`.                                       |
| `core_seed_mass`                       | `real`       | 0             | Seed mass of newly formed cores in solar masses, taken from the host cell.                          |
| `core_relevance_threshold`             | `real`       | 2             | Peak relevance threshold used for core formation sites.                                             |
| `core_density_threshold`               | `real`       | -1            | Formation density threshold in H/cc (proper), like `n_star`. With `core_use_LP_criteria` it is no longer the formation gate, but a tenth of it still sets the clump finder threshold when none of the clump thresholds below is given. |
| `core_clump_density_threshold`         | `real`       | -1            | Clump finder density threshold in H/cc. If <= 0, set to `core_density_threshold`/10.                |
| `core_clump_density_threshold_nH`      | `real`       | -1            | Clump finder density threshold in H/cc at redshift `core_clump_density_threshold_zref`, scaled as (aexp_ref/aexp)^2 in cosmological runs. Takes precedence over `core_clump_density_threshold`. |
| `core_clump_density_threshold_zref`    | `real`       | 3             | Reference redshift of `core_clump_density_threshold_nH`.                                            |
| `core_clump_density_threshold_LP_frac` | `real`       | -1            | Clump finder density threshold as a fraction of the Larson-Penston density at the finest cell size. Takes precedence over the two thresholds above. |
| `core_clump_density_threshold_LP_T`    | `real`       | 100           | Temperature (K) used for that Larson-Penston density.                                               |
| `core_clump_density_threshold_LP_mu`   | `real`       | 1.22          | Mean molecular weight used for that Larson-Penston density.                                         |
| `core_use_LP_criteria`                 | `logical`    | `.false.`     | Form cores in peak cells above the Larson-Penston density instead of `core_density_threshold`.      |
| `core_saddle_threshold`                | `real`       | -1            | Saddle-point density threshold (code units) for merging peaks into halo-patches for cores.          |
| `core_mass_threshold`                  | `real`       | 0             | Minimum clump mass (code units) kept by the clump finder for cores.                                 |
| `core_fraction_threshold`              | `real`       | 2             | Sub-halo mass fraction threshold used to identify main halos as formation sites.                    |
| `core_use_clump_prop_criteria`         | `logical`    | `.false.`     | Apply the mass4, energy and contracting criteria below. Without it (or `core_use_contracting_all_axes`) they have no effect. |
| `core_use_mass4_criteria`              | `logical`    | `.false.`     | Require the mass within 4 finest-cell widths of the peak to exceed `core_mass4_threshold`.          |
| `core_mass4_threshold`                 | `real`       | 0             | Mass threshold (code units) of the mass4 criterion.                                                 |
| `core_use_energy_criteria`             | `logical`    | `.false.`     | Require the clump to be gravitationally bound.                                                      |
| `core_include_radiation_support`       | `logical`    | `.true.`      | Include the radiation term in the support of the energy criterion.                                  |
| `core_include_magnetic_support`        | `logical`    | `.true.`      | Include the magnetic term in the support of the energy criterion.                                   |
| `core_use_contracting_criteria`        | `logical`    | `.false.`     | Require the clump to be contracting (dI/dt < 0).                                                    |
| `core_use_contracting_all_axes`        | `logical`    | `.false.`     | Require the clump to contract along all its principal axes.                                         |
| `core_refine`                          | `logical`    | `.false.`     | Refine the mesh around cores, as `sink_refine` does for sinks, so that they reach the level where `core_to_star` converts them. |
| `output_core_clump`                    | `logical`    | `.false.`     | Write the `core_clump.*` catalogue of the clumps found for core formation. Not on GPU yet.          |
| `core_merge`                           | `logical`    | `.false.`     | Activate core-core merging.                                                                         |
| `core_merger_search_order`             | `integer`    | 4             | Two cores merge if they are closer than this number of cells at their level.                        |
| `core_merging_timescale`               | `real`       | 1000          | Two cores merge only if both are younger than this age in years. If <= 0, cores do not merge.       |
| `core_to_star`                         | `logical`    | `.false.`     | Turn core mass into star particles at the finest level the mesh may reach (`nlevelmax`, or lower in cosmological runs with `aexp_lock_refine`). A core converts only once its cell is at that level: use `core_refine` so cores reach it (without it, a warning is printed at start-up when `levelmax > levelmin`). Needs star particles (`star=.true.` in `&STAR_PARAMS`): without them it does nothing. |
| `core_to_star_mass`                    | `real`       | 100           | Mass of the star particles created by `core_to_star`, in solar masses.                              |

## `&CORE_ACCRETION_PARAMS`: accretion and retirement

| Variable name                           | Fortran type | Default value | Description                                                                                        |
|:----------------------------------------|:-------------|:--------------|:---------------------------------------------------------------------------------------------------|
| `core_accretion_type`                   | `integer`    | 2             | Accretion scheme (0: none, 2: flux accretion, Bleuler et al. 2014).                                |
| `core_accretion_boost`                  | `real`       | 1             | Boost factor of the accretion rate, used as max(abs(value), 1). If < 0 and `n_star` > 0, the boost is (rho/n_star)^2 instead (at least 1). |
| `core_b_spline_order`                   | `integer`    | 2             | Order of the B-spline used to sample the accretion zone (2: CIC, 3: TSC, 4: PCS).                  |
| `core_use_rho_inf`                      | `logical`    | `.true.`      | Extrapolate the density at infinity from the Bondi solution.                                       |
| `core_use_bondi_lambda`                 | `logical`    | `.true.`      | Use the Bondi sonic-point lambda factor.                                                           |
| `core_use_bondi_for_subsonic_accretion` | `logical`    | `.false.`     | Use the Bondi rate instead of the flux rate for subsonic accretion.                                |
| `core_use_fa_fact_B14`                  | `logical`    | `.false.`     | Apply the Bleuler et al. (2014) flux accretion density threshold correction.                       |
| `core_momentum_conserving`              | `logical`    | `.true.`      | Update the core momentum with the accreted gas.                                                    |
| `core_fix_mass`                         | `logical`    | `.false.`     | Keep the core mass fixed during accretion (the gas is still removed).                              |
| `core_retire_if_no_acc`                 | `logical`    | `.false.`     | Turn a core that has stopped accreting into a star particle. Needs star particles.                 |
| `core_zero_acc_intervals`               | `integer`    | 2             | A core retires once more than this number of coarse steps have passed since it formed or last accreted. |

The mass of the cores is included in the total sink mass printed at each coarse step and used in the mass conservation check.
