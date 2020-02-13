# routine_qiime2_analyses

## Description

Running the microbiome diversity analyses has become easy thanks to tools such as [qiime2](https://docs.qiime2.org/2019.10/install/).
Yet, using the many command line to perform basic steps such as importing files in the right format, 
subsetting to remove rare species, and notably run this on a High-Performance Computer (HPC) can be cumbersome.
This wrapper does the same things as what is done in [qiita](https://qiita.ucsd.edu)-analysis, and is what could be done more 
efficiently using [snakemake](https://snakemake.readthedocs.io/en/stable/). Yet, here's a version tailored for user having access to a HPC using either the
[Torque](http://docs.adaptivecomputing.com/torque/4-0-2/help.htm) or [Slurm](https://slurm.schedmd.com/documentation.html) scheduler, and where a conda and
a qiime2 conda environment are installed. 

Here, it allows the user to pass a folder with .tsv or .biom files corresponding to microbiome-to-sample features tables,
for which another folder containing metadata for these samples is present, and run a series of analyses automatically:
- **features filtering**,
- fetching of the Web of Life **phylogeny** sub tree (only for features generated using the gOTU pipeline), 
- **alpha** diversity analysis (for now, using "observed_otus", "pielou_e", "shannon" metrics),
- **beta** diversity analysis (for now, using "jaccard", "braycurtis", "aitchison" metrics),
- **deicode**
- **permanova**
- **adonis**
- ... _(more to come)_


## Installation
```
pip install --upgrade git+https://github.com/FranckLejzerowicz/routine_qiime2_analyses.git
```
or 
```
git clone https://github.com/FranckLejzerowicz/routine_qiime2_analyses.git
cd routine_qiime2_analyses
python3 setup.py build_ext --inplace --force install
```

*_Note that python and pip should be python3_

### Depencency

- [Xpbs](https://github.com/FranckLejzerowicz/Xpbs): allows automatic preparation of HPC scripts from the basic qiime2 bash scripts
written here. For Xpbs to work, it is necessary that the user provide edit the config.txt file of this tool (simply adding 
the email address for job completion [as explained here](https://github.com/FranckLejzerowicz/Xpbs#requisite)).  

## Input

The input is strict (so that it runs well):
- path to a main folder containing at least two sub-folder named `data` and `metadata` (option `-i`):
  - sub-folder `data` must contain one or more feature table(s):
    - starting with `tab_` and ending either with:
      - `.tsv` (for tables that are tab-separated plain text with samples as columns and features as rows)
      - samples are columns; features are rows; first column is for the "`#OTU ID`"
      - `.biom` (for biom tables possibly generated from these .tsv tables or fetched using e.g. [redbiom](https://github.com/biocore/redbiom))
  - sub-folder `metadata` must contain as many metadata table(s):
    - starting with `meta_` and ending with `.tsv`
    - samples are row; variables are columns; first column is for the "`sample_name`"
- name of the dataset(s) that appear _internally_ in these folder's file names (options `-d`).
There must be a perfect matching of this _internal_ name in the features/metadata file pairs, e.g.
    ```
    datasets_folder
    ├── metadata
    │   └── meta_dataset_number_1.tsv
    │   └── meta_dataset_number_2.tsv
    ├── data
    │   └── tab_dataset_number_1.tsv
    │   └── tab_dataset_number_2.tsv
    ```
    In this case, the matching _internal_ names are `dataset_number_1` and `dataset_number_2`. Note that a found `data` 
    .tsv/.biom files that do not have a matching `metadata` .tsv path will be ignored.

    **The analysis is performed as follows:**
    - If both datasets are to be processed:
    ```
    routine_qiime2_analyses -i datasets_folder -d dataset_number_1 -d dataset_number_2 -n jobs_name -e qiime2-2019.10
    ```
    - If only the `dataset_number_2` dataset is to be processed:
    ```
    routine_qiime2_analyses -i datasets_folder -d dataset_number_2 -n jobs_name -e qiime2-2019.10
    ```
In fact, the tool simply generates scripts files that need to be started manually, and which
output should be scrutinized manually (**highly recommended**). This just a way to help you
obtain the standard qiime2 command lines pre-written for Torque/Slurm and ready to be run on a HPC!

## Outputs

After running this command (you can try):
```
routine_qiime2_analyses -i ./routine_qiime2_analyses/test/files -d dataset_number_1 -n jobs_name -e qiime2-2019.10
```
You would obtain _files_ in the `jobs` folders (scripts to check and run),
and _folders_ in the `qiime` folder (locations for qiime2 outputs).
```
.
├── data
│   ├── tab_dataset_number_1.tsv
│   └── tab_dataset_number_2.tsv
├── jobs
│   ├── alpha
│   ├── alpha_correlations
│   ├── beta
│   ├── emperor
│   ├── import_tables
│   └── pcoa
├── metadata
│   ├── meta_dataset_number_1.tsv
│   └── meta_dataset_number_2.tsv
└── qiime
    ├── alpha
    │   └── dataset_number_1
    ├── alpha_correlations
    │   └── dataset_number_1
    ├── beta
    │   └── dataset_number_1
    ├── emperor
    │   └── dataset_number_1
    └── pcoa
        └── dataset_number_1
```

The jobs to run are in printed in the stdout, i.e. the commands to copy-paste on the
HPC terminal to actually run the jobs are those after the `[TO RUN]` indicators that print
in the terminal as you run `routine_qiime2_analyses`, here, for the above example:  
```
# Fetching data and metadata (in dataset_number_1)
# Import tables to qiime2
[TO RUN] qsub /Data/Programs/routine_qiime2_analyses/routine_qiime2_analyses/test/files/jobs/import_tables/0_run_import.pbs
# Calculate alpha diversity indices
[TO RUN] sh /Data/Programs/routine_qiime2_analyses/routine_qiime2_analyses/test/files/jobs/alpha/1_run_alpha.sh
# Merge alpha diversity indices to metadata
[TO RUN] sh /Data/Programs/routine_qiime2_analyses/routine_qiime2_analyses/test/files/jobs/alpha/2_run_merge_alphas.sh
# Export alpha diversity indices to metadata
[TO RUN] qsub /Data/Programs/routine_qiime2_analyses/routine_qiime2_analyses/test/files/jobs/alpha/3_run_merge_alpha_export.pbs
# Correlate numeric metadata variables with alpha diversity indices
[TO RUN] sh /Data/Programs/routine_qiime2_analyses/routine_qiime2_analyses/test/files/jobs/alpha_correlations/4_run_alpha_correlation.sh
# Calculate beta diversity indices
[TO RUN] sh /Data/Programs/routine_qiime2_analyses/routine_qiime2_analyses/test/files/jobs/beta/2_run_beta.sh
# Export beta diversity matrices
[TO RUN] qsub /Data/Programs/routine_qiime2_analyses/routine_qiime2_analyses/test/files/jobs/beta/2x_run_beta_export.pbs
# Calculate principal coordinates
[TO RUN] sh /Data/Programs/routine_qiime2_analyses/routine_qiime2_analyses/test/files/jobs/pcoa/3_run_pcoa.sh
# Make EMPeror plots

Warning: Make sure you first run alpha -> alpha merge -> alpha export
        (if you want alpha diversity as a variable in the PCoA)!
[TO RUN] sh /Data/Programs/routine_qiime2_analyses/routine_qiime2_analyses/test/files/jobs/emperor/4_run_emperor.sh
```
The job are labeled for you to get a suggestion of the order in which to run them
(this is essentially what `snakemake` would do but it does but the whole point
 here is to have some human verification.)

## **WARNING**: Check the jobs scripts first!

The files to copy-paste in the HPC terminal all have a name that start by a number. 
However, different things happen if the file name ends with: 
- `.pbs`: submit one HPC job (_that may consist in of several qiime2 commands, e.g. import multiple artefacts_)
- `.sh`: script that submits multiple HPC jobs (_each may also consist of several qiime2 commands_)
  - these multiple jobs are the `.pbs` files in the subfolders named `chunks`

**PLEASE ALWAYS DO CHECK OF THE ACTUAL QIIME2 COMMANDS, WRITTEN IN THE `.pbs` FILES**  

Example:
```
├── jobs
│   ├── alpha
│   │   ├── 1_run_alpha.sh
│   │   ├── 2_run_merge_alphas.sh
│   │   ├── 3_run_merge_alpha_export.pbs
│   │   ├── 3_run_merge_alpha_export.sh
│   │   └── chunks
│   │       ├── run_alpha_dataset_number_1.pbs
│   │       ├── run_alpha_dataset_number_1.sh
│   │       ├── run_merge_alpha_dataset_number_1.pbs
│   │       └── run_merge_alpha_dataset_number_1.sh
│   ├── alpha_correlations
│   │   ├── 4_run_alpha_correlation.sh
│   │   └── chunks
│   │       ├── run_alpha_correlation_dataset_number_1.pbs
│   │       └── run_alpha_correlation_dataset_number_1.sh
│   ├── beta
│   │   ├── 2_run_beta.sh
│   │   ├── 2x_run_beta_export.pbs
│   │   ├── 2x_run_beta_export.sh
│   │   └── chunks
│   │       ├── run_beta_dataset_number_1.pbs
│   │       └── run_beta_dataset_number_1.sh
│   ├── emperor
│   │   ├── 4_run_emperor.sh
│   │   └── chunks
│   │       ├── run_emperor_dataset_number_1_tab_dataset_number_1_aitchison.pbs
│   │       ├── run_emperor_dataset_number_1_tab_dataset_number_1_aitchison.sh
│   │       ├── run_emperor_dataset_number_1_tab_dataset_number_1_braycurtis.pbs
│   │       ├── run_emperor_dataset_number_1_tab_dataset_number_1_braycurtis.sh
│   │       ├── run_emperor_dataset_number_1_tab_dataset_number_1_jaccard.pbs
│   │       └── run_emperor_dataset_number_1_tab_dataset_number_1_jaccard.sh
│   ├── import_tables
│   │   ├── 0_run_import.pbs
│   │   └── 0_run_import.sh
│   └── pcoa
│       ├── 3_run_pcoa.sh
│       └── chunks
│           ├── run_PCoA_dataset_number_1.pbs
│           └── run_PCoA_dataset_number_1.sh

```

## PERMANOVA

It is possible to run PERMANOVA for a series of user-defined subsets of the data and to test difference between 
different groups of each subset automatically.

- **permanova tests** (`-t`):

    This use of `-t` will result in one test for each factor to the column `sex`, as well as one subset for each
    factor to the column `age_cat`. As in this example, note that `-t` can be used multiple time, once per group. 
    
- **group subsets** (`-g`): a config file must be provided in the following .yml format:  
    ```
    sex:
    - - Male
    - - Female
    timepoint_months:
    - - '9'
      - '24'
    - - '24'
      - '36'  
    income:
    - - '<15000'
    - - '>15000'
    ```
    which is interpreted as a dictionary which for each metadata variable, lists one or more factor(s) 
    defining a subset:
    ```
    {'sex': [['Male'], ['Female']],
     'timepoint_months': [['9', '24'], ['24', '36']],
     'income': [['<15000'], ['>15000']]}    
    ```
    In this example, there will be one subset for:
     - samples having `Male` in column `sex`
     - samples having `Female` in column `sex`
     - samples having `9` or `24` in column `timepoint_months`
     - samples having `24` or `36` in column `timepoint_months`
     - samples having a value inferior to `15000` in column `income`
     - samples having a value superior to `15000` in column `income`
    For example:
    ```
    routine_qiime2_analyses \
        -t sex \
        -t age_cat \
        -g ./routine_qiime2_analyses/examples/example_PERMANOVA_subsets.yml \
        -i ./routine_qiime2_analyses/test/files \
        -d dataset_number_1 \
        -d test2 \
        -n test \
        -e qiime2-2019.10
    ```
        
- The output is self contained, e.g.: `tab_dataset_number_2_braycurtis_sex_Female__timepoint_months_permanova.qzv`
is for the `Female` subset of metadata variable `sex` (it also does the result for `Male` etc), and using PERMANOVA to perform comparison between 
the groups in columns `timepoint_months`. 


## ADONIS

It is possible to run R's Adonis in Qiime2 for a series of user-defined formulas to test difference as in PERMANOVA
for multivariate data but with continuous and multiple metadata variables as regressors
(see [http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/adonis.html](http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/adonis.html)).
The passed models will perform on each of the subsets defined in the file passed to option `-g`, as above.

- **adonis formula** (`-a`): a config file must be provided in the following .yml format:
    ```
    sexPLUSincomeINTERtime: "sex+income*timepoint_months"
    incomePLUStime: "income+timepoint_months"
    ```
    which is interpreted as a dictionary which for each metadata variable, lists one or more factor(s) 
    defining a subset:
    ```
    {'sexPLUSincomeINTERtime': 'sex+income*timepoint_months',
     'incomePLUStime': 'income+timepoint_months'}
    ```
    In this example, there will be one model for each formula (and for each distance matrix),
    which in R, would correspond to these commands:
     - `adonis(<bray_curtis_distance_matrix-file> ~ sex + income * timepoint_months, <metadata-file>)`
     - `adonis(<bray_curtis_distance_matrix-file> ~ income * timepoint_months, <metadata-file>)`
     
    For example:
    ```
    routine_qiime2_analyses \
        -a ./routine_qiime2_analyses/examples/example_ADONIS_formulas.yml \
        -g ./routine_qiime2_analyses/examples/example_PERMANOVA_subsets.yml \
        -i ./routine_qiime2_analyses/test/files \
        -d dataset_number_1 \
        -d dataset_number_2 \
        -n jobs_name \
        -e qiime2-2019.10
    ```
    This use of `-a` will result in one test for each formula placed as rows in the .yml file. 
    
- **group subsets** (`-g`): a config file must be provided in the following .yml format.
This is the exact same file (and thus format) as for the PERMANOVA above.   
    
- The output is self contained, e.g.: `tab_dataset_number_1_braycurtis_sex_Female__sexPLUSincomeINTERtime_adonis.qzv` is for
the `Female` subset of metadata variable `sex` (it also does the result for `Male` etc), and using ADONIS to perform
testing between the groups in columns `timepoint_months`. 

 
 
## Usage

```
routine_qiime2_analyses -i <input_folder_path> -d <dataset_name> -n <project_name> -e <qiime2_env> [OPTIONS]
```

### Optional arguments

``` 
  -i, --i-datasets-folder TEXT  Path to the folder containing the sub-folders
                                'data' and 'metadata'.  [required]
  -d, --i-datasets TEXT         Identifier(s) of the dataset(s) (e.g. '-d
                                dataset_number_1 -d dataset_number_2' for
                                inputs'data/tab_dataset_number_1.tsv +
                                metadata/meta_dataset_number_1.tsv' as well as
                                'data/tab_dataset_number_2.tsv +
                                metadata/meta_dataset_number_2.tsv')
                                [required]
  -t, --i-wol-tree TEXT         path to the tree containing the genome IDs
                                alse present in the features names (On
                                barnacle, it is there: /projects/wol/profiling
                                /dbs/wol/phylogeny/tree.nwk).
  -n, --p-project-name TEXT     Nick name for your project.  [required]
  -e, --p-qiime2-env TEXT       name of your qiime2 conda environment (e.g.
                                qiime2-2019.10)   [required]
  -t, --p-perm-tests TEXT       Groups to tests between in each PERMANOVA
                                subset (multiple values are possible, e.g. '-t
                                sex -t age_cat').  [default: False]
  -g, --p-perm-groups TEXT      Subsets for PERMANOVA. Must be a yaml file,
                                e.g.
                                (see example in 'examples/example_PERMANOVA_subsets.yml'
                                and the README).  [default: False]
  -a, --p-adonis-formulas TEXT  Formula for Adonis tests for each PERMANOVA
                                subset. Must be a yaml file, e.g.
                                (see example
                                in 'examples/example_ADONIS_formulas.yml' and the README).
                                [default: False]
  -l, --p-longi-column TEXT     If data is longitudinal; provide the time
                                metadata columnfor volatility analysis.
                                [default: False]
  -f, --p-reads-filter INTEGER  Minimum number of reads per sample to be kept.
                                [default: 0]
  --force / --no-force          Force the re-writing of scripts for all
                                commands(default is to not re-run if output
                                file exists).  [default: False]
  --gid / --no-gid              If feature names have the genome ID (to use
                                the Web of Life tree).  [default: False]
  --biom / --no-biom            Use biom files in the input folder(automatic
                                if there's no .tsv but only .biom file(s)).
                                [default: False]
  --version                     Show the version and exit.
  --help                        Show this message and exit.
```

## Example

For the command:
```
routine_qiime2_analyses  \
    -i ./routine_qiime2_analyses/test/files  \
    -d dataset_number_1  \
    -d dataset_number_2  \
    -t ./routine_qiime2_analyses/resources/web_of_life_tree.nwk  \
    -n test_name  \
    -e qiime2-2019.10  \
    -t sex \
    -t age_cat  \
    -g ./routine_qiime2_analyses/examples/example_PERMANOVA_subsets.yml  \
    -a ./routine_qiime2_analyses/examples/example_ADONIS_formulas.yml  \
    --gid  \
    -l timepoint_months
    -f 10000
```
The standard output shows you the scripts that have been written with qiime2 commands and that need to be run:
```
# Fetching data and metadata (in dataset_number_1, dataset_number_2)
# Import tables to qiime2
[TO RUN] qsub /Data/Programs/routine_qiime2_analyses/routine_qiime2_analyses/test/files/jobs/import_tables/0_run_import.pbs
# Filter samples for a min number of 10000 reads
[TO RUN] qsub /Data/Programs/routine_qiime2_analyses/routine_qiime2_analyses/test/files/jobs/import_filtered/1_run_import_filtered.pbs
# Shear Web of Life tree to features' genome IDs
[TO RUN] qsub /Data/Programs/routine_qiime2_analyses/routine_qiime2_analyses/test/files/jobs/import_tree_dataset_number_1/0_import_tree.pbs
[TO RUN] qsub /Data/Programs/routine_qiime2_analyses/routine_qiime2_analyses/test/files/jobs/import_tree_dataset_number_2/0_import_tree.pbs
[TO RUN] qsub /Data/Programs/routine_qiime2_analyses/routine_qiime2_analyses/test/files/jobs/import_tree_dataset_number_1_min10000_339s/0_import_tree.pbs
[TO RUN] qsub /Data/Programs/routine_qiime2_analyses/routine_qiime2_analyses/test/files/jobs/import_tree_dataset_number_2_min10000_339s/0_import_tree.pbs
# Calculate alpha diversity indices
[TO RUN] sh /Data/Programs/routine_qiime2_analyses/routine_qiime2_analyses/test/files/jobs/alpha/1_run_alpha.sh
# Merge alpha diversity indices to metadata
[TO RUN] sh /Data/Programs/routine_qiime2_analyses/routine_qiime2_analyses/test/files/jobs/alpha/2_run_merge_alphas.sh
# Export alpha diversity indices to metadata
[TO RUN] qsub /Data/Programs/routine_qiime2_analyses/routine_qiime2_analyses/test/files/jobs/alpha/3_run_merge_alpha_export.pbs
# Correlate numeric metadata variables with alpha diversity indices
[TO RUN] sh /Data/Programs/routine_qiime2_analyses/routine_qiime2_analyses/test/files/jobs/alpha_correlations/4_run_alpha_correlation.sh
# Longitudinal change in alpha diversity indices

Warning: First make sure you run alpha -> alpha merge -> alpha export before running volatility
        (if you need the alpha as a response variable)!
# Calculate beta diversity indices
[TO RUN] sh /Data/Programs/routine_qiime2_analyses/routine_qiime2_analyses/test/files/jobs/beta/2_run_beta.sh
# Export beta diversity matrices
[TO RUN] qsub /Data/Programs/routine_qiime2_analyses/routine_qiime2_analyses/test/files/jobs/beta/2x_run_beta_export.pbs
# Calculate principal coordinates
[TO RUN] sh /Data/Programs/routine_qiime2_analyses/routine_qiime2_analyses/test/files/jobs/pcoa/3_run_pcoa.sh
# Make EMPeror plots

Warning: Make sure you first run alpha -> alpha merge -> alpha export
        (if you want alpha diversity as a variable in the PCoA)!
[TO RUN] sh /Data/Programs/routine_qiime2_analyses/routine_qiime2_analyses/test/files/jobs/emperor/4_run_emperor.sh
# DEICODE (groups config in ./routine_qiime2_analyses/examples/example_PERMANOVA_subsets.yml)
sh /Data/Programs/routine_qiime2_analyses/routine_qiime2_analyses/test/files/jobs/deicode/3_run_beta_deicode.sh
# Kruskal-Wallis (groups config in ./routine_qiime2_analyses/examples/example_PERMANOVA_subsets.yml)
[TO RUN] sh /Data/Programs/routine_qiime2_analyses/routine_qiime2_analyses/test/files/jobs/alpha_group_significance/6_run_alpha_group_significance.sh
# PERMANOVA (groups config in ./routine_qiime2_analyses/examples/example_PERMANOVA_subsets.yml)
Beta diversity, distances matrices must be generated already to automatise PERMANOVA
        (re-run this after steps "2_run_beta.sh" and "2x_run_beta_export.pbs" are done)
[TO RUN] sh /Data/Programs/routine_qiime2_analyses/routine_qiime2_analyses/test/files/jobs/permanova/3_run_beta_group_significance.sh
# Run Adonis (groups config in ./routine_qiime2_analyses/examples/example_PERMANOVA_subsets.yml)
Beta diversity, distances matrices must be generated already to automatise adonis
        (re-run this after steps "2_run_beta.sh" and "2x_run_beta_export.pbs" are done)
[TO RUN] sh /Data/Programs/routine_qiime2_analyses/routine_qiime2_analyses/test/files/jobs/adonis/3_run_adonis.sh
```


### Bug Reports

contact `flejzerowicz@health.ucsd.edu`