# Running the thesis code

## 1. Install dependencies

Use Python 3.12 or later and R 4.4 or later. Make `python` and `Rscript` available
in your terminal. Run the commands below from the repository root, which contains
`1.Input data`, `2.data_preprocessing`, and `3.Imputation and evaluation`.

Install the Python packages:

```sh
python -m pip install -r "2.data_preprocessing/2.2 missForest/requirements.txt"
```

Install the R packages in an R console:

```r
install.packages(c("missForest", "data.table", "ranger", "randomForest",
                   "fastshap", "shapviz", "ggplot2", "remotes"))
remotes::install_gitlab(
  "water/ARCHI", host = "code.usgs.gov",
  ref = "d6e952f5d2021303e699a39822b2396886dd1399",
  dependencies = NA, build_vignettes = FALSE, upgrade = "never"
)
```

## 2. Place the input files

Keep this folder structure and these filenames:

```text
1.Input data/
  Location_SGI.xlsx
  EDGE_SE_GWL.csv
  Eobs/
    1980-1994/
    1995-2010/
    2011-2024/
2.data_preprocessing/
  2.1 SGI/
    preprocess_sgi.py
  2.2 missForest/
    preprocess_lerum.py
    preprocess_global.py
    requirements.txt
  2.3 ARCHI/
    preprocess_archi.py
    requirements.txt
3.Imputation and evaluation/
  3.1missForest/
    impute_evaluate_lerum.R
    impute_evaluate_global.R
    np_shap_lerum.R
  3.2ARCHI/
    impute_evaluate_archi.R
```

`Location_SGI.xlsx` must contain the sheet `Location_SGI`. Keep the original
column headers in both input files.

### Download groundwater-level data

Data source: [Groundwater-level data on OSF](https://osf.io/skxg5/files/osfstorage).

Download the Swedish groundwater-level CSV from the linked data collection and
extract it if supplied in an archive. The preprocessing scripts expect the file
`EDGE_SE_GWL.csv` at this location relative to the repository root:

```text
1.Input data/EDGE_SE_GWL.csv
```

Keep the original column names. The scripts read `station_name`, `date`, and
`gwl_(inf)`. Store this large file locally; it does not need to be uploaded to
GitHub. Place it in the indicated folder before running either the missForest
or ARCHI preprocessing scripts.

### Download E-OBS data

Download the large E-OBS files separately and store them locally; they do not
need to be uploaded to GitHub.

Download page: [E-OBS daily gridded meteorological data for Europe — Copernicus Climate Data Store](https://cds.climate.copernicus.eu/datasets/insitu-gridded-observations-europe?tab=overview).

Open the **Download** tab and sign in to a CDS account. Accept the dataset licence
when prompted. Use the following settings to match the input files used here:

| Setting | Selection |
| --- | --- |
| Version | v31.0e |
| Grid resolution | 0.1° × 0.1°, regular latitude–longitude grid |
| Product | Ensemble mean |
| Time step | Daily |
| File format | NetCDF |
| Periods | 1980–1994, 1995–2010, 2011–2024 |

Download each of these five variables for all three periods:

| File prefix | Variable |
| --- | --- |
| `rr` | Precipitation amount |
| `tg` | Mean temperature |
| `hu` | Relative humidity |
| `fg` | Wind speed |
| `qq` | Surface shortwave downwelling radiation |

The preprocessing scripts require these five variables; `tn` and `tx` are not
required. Extract any downloaded archives and place the `.nc` files under
`1.Input data/Eobs/` in the matching period folders:

```text
1.Input data/Eobs/
  1980-1994/
    rr_ens_mean_0.1deg_reg_1980-1994_v31.0e.nc
    tg_ens_mean_0.1deg_reg_1980-1994_v31.0e.nc
    hu_ens_mean_0.1deg_reg_1980-1994_v31.0e.nc
    fg_ens_mean_0.1deg_reg_1980-1994_v31.0e.nc
    qq_ens_mean_0.1deg_reg_1980-1994_v31.0e.nc
  1995-2010/
    (the same five variables for 1995-2010)
  2011-2024/
    (the same five variables for 2011-2024)
```

Preserve the original filenames and NetCDF variable names. Keep only one copy
of each variable and period under `Eobs/`, since the scripts search its
subfolders recursively. The default paths work without editing the code.

## 3. Run preprocessing

Run SGI preprocessing first:

```sh
python "2.data_preprocessing/2.1 SGI/preprocess_sgi.py"
```

This creates `2.data_preprocessing/2.1 SGI/Location_SGI_final.xlsx`.
The following scripts use this generated file automatically:

```sh
python "2.data_preprocessing/2.2 missForest/preprocess_lerum.py"
python "2.data_preprocessing/2.2 missForest/preprocess_global.py"
python "2.data_preprocessing/2.3 ARCHI/preprocess_archi.py"
```

Keep `preprocess_lerum.py` and `preprocess_global.py` together because the global
script imports functions from the Lerum script.

## 4. Run imputation and evaluation

Run each analysis after its corresponding preprocessing has finished:

```sh
Rscript "3.Imputation and evaluation/3.1missForest/impute_evaluate_lerum.R"
Rscript "3.Imputation and evaluation/3.1missForest/impute_evaluate_global.R"
Rscript "3.Imputation and evaluation/3.2ARCHI/impute_evaluate_archi.R"
```

For Lerum SHAP analysis, run the following after Lerum preprocessing. Replace
`Lerum_1` with the desired prepared well:

```sh
Rscript "3.Imputation and evaluation/3.1missForest/np_shap_lerum.R" --well Lerum_1
```

## 5. File locations

All default paths are resolved from the script locations. Keep the directory
structure when moving the repository to another computer. Leave each generated
preprocessing directory intact so the analysis scripts can locate its data and
supporting files.

| Script | Default output location, relative to the repository root |
| --- | --- |
| `preprocess_sgi.py` | `2.data_preprocessing/2.1 SGI/Location_SGI_final.xlsx` |
| `preprocess_lerum.py` | `2.data_preprocessing/2.2 missForest/Lerum/` |
| `preprocess_global.py` | `2.data_preprocessing/2.2 missForest/Global/` |
| `preprocess_archi.py` | `2.data_preprocessing/2.3 ARCHI/` |
| `impute_evaluate_lerum.R` | `3.Imputation and evaluation/3.1missForest/Lerum/` |
| `impute_evaluate_global.R` | `3.Imputation and evaluation/3.1missForest/Global/` |
| `np_shap_lerum.R` | `3.Imputation and evaluation/3.1missForest/Lerum_SHAP/<well>/` |
| `impute_evaluate_archi.R` | `3.Imputation and evaluation/3.2ARCHI/Scenario_A/nref_10/` |

To use different paths or select individual wells, append `--help` to any script
command to see its options. Quote paths that contain spaces.
