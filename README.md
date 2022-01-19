# microbiome_analyzer

## Description

Running the microbiome diversity analyses has become easy thanks to tools such as [qiime2](https://docs.qiime2.org/2019.10/install/).
Yet, using many command lines to perform basic steps such as importing files in the right format, 
subsetting to remove rare species, subset samples, collapse features based on taxonomix rank or names, 
and for all that to run this on a High-Performance Computer (HPC), can be cumbersome.
This wrapper does the same things as what is available (and not yet available) in [qiita](https://qiita.ucsd.edu)-analysis, and is what could be done more 
efficiently using [snakemake](https://snakemake.readthedocs.io/en/stable/). Yet, here's a version tailored for user having access to a HPC using either the
[Torque](http://docs.adaptivecomputing.com/torque/4-0-2/help.htm) or [Slurm](https://slurm.schedmd.com/documentation.html) scheduler, and where a qiime2 conda environment are installed. 

Here, it allows the user to pass a folder with .tsv or .biom files corresponding to microbiome-to-sample features tables,
for which another folder containing metadata for these samples is present, and run a series of analyses automatically:
- **features filtering**,
- fetching of the Web of Life **phylogeny** sub tree (features generated using the 
- [OGUs](https://www.biorxiv.org/content/10.1101/2021.04.04.438427v1.abstract) pipeline), 
- **alpha** diversity analysis
- **beta** diversity analysis
  - **pcoa**
  - **tsne**
  - **umap**
  - **robust-pca**
- **permanova**
- **procrustes**
- **mantel**
- R's **adonis**
- **nestedness** (incl. NODF calulcation acrss
- **dissimilarity-overlap curves**
- Differential analysis using **(songbird)[https://github.com/biocore/songbird]**
- Co-occurrence estimation between datasets using **(mmvec)[https://github.com/biocore/mmvec]**
- **sourcetracking**
- ... _(more to come)_


## Installation
```
pip install --upgrade git+https://github.com/FranckLejzerowicz/microbiome_analyzer.git
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
    microbiome_analyzer -i datasets_folder -d dataset_number_1 -d dataset_number_2 -n jobs_name -e qiime2-2019.10
    ```
    - If only the `dataset_number_2` dataset is to be processed:
    ```
    microbiome_analyzer -i datasets_folder -d dataset_number_2 -n jobs_name -e qiime2-2019.10
    ```
In fact, the tool simply generates scripts files that need to be started manually, and which
output should be scrutinized manually (**highly recommended**). This just a way to help you
obtain the standard qiime2 command lines pre-written for Torque/Slurm and ready to be run on a HPC!

## Outputs

After running this command (you can try):
```
microbiome_analyzer -i ./microbiome_analyzer/tests/files -d dataset_number_1 -n jobs_name -e qiime2-2021.11
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

Here's the stdout for the simple command above:
```
# import
sh /abs/path/to/microbiome_analyzer/tests/files/jobs/import/jbs_nm.sh
# taxonomy
sh /abs/path/to/microbiome_analyzer/tests/files/jobs/taxonomy/jbs_nm.sh
```

These prints are jobs to run, i.e. these `sh` commands only need to be copy-pasted on the
HPC terminal to actually run the `.pbs` (for Torque) or `.slm` (for Slurm, use option `--slurm`)
scripts within. For example, the first `.sh` file contains:

```
$ cat /abs/path/to/microbiome_analyzer/tests/files/jobs/import/jbs_nm.sh
qsub /abs/path/to/microbiome_analyzer/tests/files/jobs/import/chunks/jbs_nm_dataset_number_1.pbs
```

* Note: If you where to prepare scripts for 5 datasets:
  ```
  microbiome_analyzer -i ./microbiome_analyzer/tests/files \
      -d dataset_number_1 \
      -d dataset_number_2 \
      -d dataset_number_3 \
      -d dataset_number_4 \
      -d dataset_number_5 \
      -n jobs_name -e qiime2-2021.11
  ```
  Then this first `.sh` file would contain:
  ```
  $ cat /abs/path/to/microbiome_analyzer/tests/files/jobs/import/jbs_nm.sh
  qsub /abs/path/to/microbiome_analyzer/tests/files/jobs/import/chunks/jbs_nm_dataset_number_1.pbs
  qsub /abs/path/to/microbiome_analyzer/tests/files/jobs/import/chunks/jbs_nm_dataset_number_2.pbs
  qsub /abs/path/to/microbiome_analyzer/tests/files/jobs/import/chunks/jbs_nm_dataset_number_3.pbs
  qsub /abs/path/to/microbiome_analyzer/tests/files/jobs/import/chunks/jbs_nm_dataset_number_4.pbs
  qsub /abs/path/to/microbiome_analyzer/tests/files/jobs/import/chunks/jbs_nm_dataset_number_5.pbs
  ```
  * *Trick here*: using the option `-chunkit <int>` to group the commands into less jobs. This can be 
  useful if you have say 50 datasets across which distance matrices you plan on doing mantel tests,
  this would send probably too many jobs to the scheduler. Let's see what it does to make For our 5 datasets:
  ```
  microbiome_analyzer -i ./microbiome_analyzer/tests/files \
      -d dataset_number_1 \
      -d dataset_number_2 \
      -d dataset_number_3 \
      -d dataset_number_4 \
      -d dataset_number_5 \
      -n jobs_name -e qiime2-2021.11 -chunkit 3
  ```
  The `.sh` file would now contain only two jobs:
  ```
  $ cat /abs/path/to/microbiome_analyzer/tests/files/jobs/import/jbs_nm.sh
  qsub /abs/path/to/microbiome_analyzer/tests/files/jobs/import/chunks/jbs_nm_jbs_nm_chunk0.pbs
  qsub /abs/path/to/microbiome_analyzer/tests/files/jobs/import/chunks/jbs_nm_jbs_nm_chunk1.pbs
  qsub /abs/path/to/microbiome_analyzer/tests/files/jobs/import/chunks/jbs_nm_jbs_nm_chunk2.pbs
  ```

Only `import` and `taxonomy` are showing up because the former one needs to be run first before most other
usual routine qiime2 analyses, whereas the latter can be run by generating fasta file from the features before
running the taxonomic assignment. 

Note that this feature-to-fasta encoding into taxonomy only would work for features encoded as sequences or that 
correspond to the Web of Life genomes ("G#########" format). *If the features are taxa or OTU names, the taxonomy
will be each taxon or OTU name*: in this case, please edit manually the taxonomy file (column "Taxon").

After running these two jobs, if you re-run the same, simple command above, you get:

```
# alpha
sh /Users/franck/programs/microbiome_analyzer/microbiome_analyzer/tests/files/jobs/alpha/jbs_nm.sh
# tabulate
sh /Users/franck/programs/microbiome_analyzer/microbiome_analyzer/tests/files/jobs/tabulate/jbs_nm.sh
# barplot
sh /Users/franck/programs/microbiome_analyzer/microbiome_analyzer/tests/files/jobs/barplot/jbs_nm.sh
# beta
sh /Users/franck/programs/microbiome_analyzer/microbiome_analyzer/tests/files/jobs/beta/jbs_nm.sh
```
That's more work ready to start, because now the data table was imported to qiime2:

__WARNING: It is strongly recommended to check the jobs scripts first!__

## High-level configurations

Whatever the analyses you will run, `microbiome_analyzer` can run them on 
different version of each dataset, that are all defined using yaml files
passed in command line, including:

### Filtering (option `-f`):

The yaml file tells, for each dataset, which samples to remove (`names`), 
and which thresholds to use to remove rare features (`features`) and poorly
sequenced samples (`samples`), e.g.:

```
dataset_number_1:  # dataset name
  names:
  - "samplename1"
  - "samplename2"
  samples: 100
  features: 0.0001
```
which is interpreted as a dictionary:
```
{"dataset_number_1": {"names": ["samplename1", "samplename2"],
                      features: 0.0001, samples: 100}    
```
The filteting proceeds in this order:
* first remove the samples in the `names` list 
* second remove the `samples` that do not have enough reads 
* third remove the `features` (from each sample) that do not have enough reads
    
For the thresholds, if the number is above 1, then the filtering is based on absolute reads values,
but it it is between 0 and 1, then the filtering is based on percentages, and in effect:
* if `samples: 100`: all samples must have 100 reads or more  
* if `samples: 0.25`: all samples must have 25% of the average amount of reads per sample  
* of `features: 1000`: all features with less than 1000 reads in a sample will be removed from that sample.
* of `features: 0.001`: all features with less than 0.1% of the reads of a sample will be removed from that sample.
    

### Rarefaction depths: (option `-r` and `--raref`):

Note that `--raref` **MUST** be set for the rarefaction to happen. If no yaml file is given along
using option `-r`, then the rarefaction will be done using as rarefaction depth the second percentile
of the distribution of number of reads per sample. 

If the following yaml file is given to option `-r`:
```
dataset_number_1:
  - "100"
  - "200"
dataset_number_2:
  - min
```
There will be two rarefaction depths for `dataset_number_1`, while `dataset_number_2`
will be rarefied to the depth of the minimum of the distribution of number of reads per sample.

### Feature subsets (option `-k`):

The yaml file gives, for each dataset, the list of names **or regex** to find the features 
for which to make different subsets. 
```
dataset_number_1:
  OGUs_selection:
    - "G000146185"
    - "G000153885"
    - "G000153905"
  OGUs_ragex:
    - "G000[12]0000[89]"
dataset_number_5:
  Milk_dataset:
    - "Milk"
```
This will create:
- two subsets for `dataset_number_1`:
  * subset `OGUs_selection` will contain three features given exactly
  * subset `OGUs_ragex` will contain max the four features matching the regex (i.e., `G000100008`, `G000100009`, `G000200009` or `G000200009`) 
- one subset for `dataset_number_5`:
  * subset `Milk_dataset` will contain all features contaiing the work "Milk".

### Taxonomic collapsing (option `--coll`):

The yaml file gives, for each dataset, a given name of taxonomic level as key,
and as value, either:
- the identifier in the taxonomy for this level (e.g.
`"p_Firmicutes"` will be a phylum and hence `"p"`).
- the index of the column to collapse to in the taxonomic table 
(`1` would be for the very first, coarsest taxonomic level).

```
dataset_number_1:
  phylum: "p"
  family: "f"
  genus: "g"
dataset_number_5:
  level1: 1
  level2: 2
  level4: 4
```

### Sample subsets (option `-g`):

The yaml config file must be have  the following format:  
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
which is interpreted as a dictionary which for each metadata variable, 
lists one or more factor(s) defining a subset:
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

Note that each subset will be made on each combination of the above, i.e., for 
each filtering version, each feature subset and each taxonomic level. 

## Analyses of the datasets (and their versions)

### Automatic

A range of analyses are prepared systemacially, including:
* [**Phylogenetic sequence placement**](https://docs.qiime2.org/2021.11/plugins/available/fragment-insertion/)
using **SEPP** (for 16S data, which is detected)
* [**Alpha diversity**](https://docs.qiime2.org/2021.11/plugins/available/diversity/alpha/): by default, the metrics defined in the resource files located in
    `./microbiome_analyzer/resources/alpha_metrics.txt` are used, i.e.:
    `pielou_e`, `shannon`, `faith_pd` and `observed_otus` (or `observed_features`
    in latest qiime2 versions).
* [**Alpha correlations**](https://docs.qiime2.org/2021.11/plugins/available/diversity/alpha-correlation/)
* [**Merging**](https://docs.qiime2.org/2021.11/plugins/available/metadata/tabulate/) of all alpha metrics (i.e., all dataset versions) to the metadata
* [**Alpha rarefaction**](https://docs.qiime2.org/2021.11/plugins/available/diversity/alpha-rarefaction/)
* [**Barplot**](https://docs.qiime2.org/2021.11/plugins/available/taxa/barplot/):
* Beta diversity:
  * [**DEICODE**](https://github.com/biocore/DEICODE) is a robust PCA method.   
  * [**Distance matrices generation**](https://docs.qiime2.org/2021.11/plugins/available/diversity/beta/): by default, the metrics defined in the resource files located in
`./microbiome_analyzer/resources/beta_metrics.txt` are used, i.e.:
`jaccard`, `braycurtis`, `aitchison`, `unweighted_unifrac` and `weighted_unifrac`.
  * Dimensionality reduction:
    * [**PCoA**](https://docs.qiime2.org/2021.11/plugins/available/diversity/pcoa/)
    * [**PCoA (biplot)**](https://docs.qiime2.org/2021.11/plugins/available/diversity/pcoa-biplot/)
    * [**t-SNE**](https://docs.qiime2.org/2021.11/plugins/available/diversity/tsne/)
    * [**UMAP**](https://docs.qiime2.org/2021.11/plugins/available/diversity/umap/)
  * Vizualization:
    * [**Emperor**](https://docs.qiime2.org/2021.11/plugins/available/emperor/plot/)
    * [**Emperor (biplot)**](https://docs.qiime2.org/2021.11/plugins/available/emperor/biplot/)
    * [**Empress**](https://library.qiime2.org/plugins/empress/32/) for all 
    * [**Empress (biplot)**](https://library.qiime2.org/plugins/empress/32/)

### Volatility (option `-l`)

* [**Volatility**](https://docs.qiime2.org/2021.11/plugins/available/longitudinal/volatility/)
will output in `./qiime/volatility` and must be triggered by using the option `-l`, 
which should be a continuous (or numeric ordinal) metadata columns (usually the time points).

### Procrustes (option `-prc`) and Mantel test (option `-mtl`)

* [**Procrustes**](https://docs.qiime2.org/2021.11/plugins/available/diversity/procrustes-analysis/) 
will output in `./qiime/procrustes` under a folder defined for each datasets pair  
* [**Mantel test**](https://docs.qiime2.org/2021.11/plugins/available/diversity/mantel/)
will output in `./qiime/mantel` under a folder defined for each datasets pair

These two-datasets comparison methods are using the same subsetting mechanism,
to run on the sample in common between datasets pairs defined in a yaml file, e.g.:
```
pairs:
  dataset_number_1_5:  # invented name  
    - "dataset_number_1"
    - "dataset_number_5"
  my_datastes_pair:  # invented name  
    - "dataset_number_3"
    - "dataset_number_4"
```

Note that as above, the tests are also computed for sample subsets applied to these pairs.

### PERMANOVA (option(s) `-t`)

* [**PERMANOVA and PERMDISP**](https://docs.qiime2.org/2021.11/plugins/available/diversity/beta-group-significance/):
It is possible to run PERMANOVA for a series of user-defined subsets of the data and to test difference between 
different groups of each subset automatically.

#### **permanova tests** 

This use of `-t` will result in one test for each factor to the column `sex`, as well as one subset for each
factor to the column `age_cat`. As in this example, note that `-t` can be used multiple time, once per group. 

### PHATE (yaml file to option `-phate`)

* [**PHATE**](https://github.com/KrishnaswamyLab/PHATE) is dimensionality reduction method that captures both local and global structures,
that will output in `./qiime/phate`

### Adonis (yaml file to option `-a`)

* [**Adonis**](https://docs.qiime2.org/2021.11/plugins/available/diversity/adonis/):
will output in `./qiime/adonis`

It is possible to run R's Adonis in Qiime2 for a series of user-defined formulas 
to test difference as in PERMANOVA for multivariate data but with continuous and
multiple metadata variables as regressors (see [here](http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/adonis.html)).
The passed models will perform on each of the subsets defined in the file passed to option `-g`, as above.

#### **adonis formula** 

A config file must be provided in the following .yml format:
```
models:
  dataset_number_1:
    timepoint_months_PLUS_income: "timepoint_months+income"
    timepoint_months_INTER_income: "timepoint_months*income"
    timepoint_months_PLUS_income_INTER_sex: "timepoint_months+income*sex"
    sex: "sex"
    timepoint_months: "timepoint_months"
  dataset_number_5:
    timepoint_months: "timepoint_months"
    timepoint_months_PLUS_income: "timepoint_months+income"
strata:
  global: "sex"
  dataset_number_1:
    timepoint_months_PLUS_income:
    - "sex"
    - "group"
    - "notacolumn"
    timepoint_months_INTER_income:
    - "sex"
    timepoint_months_PLUS_income_INTER_sex:
    - "sex"
  vioscreen_foods_consumed_grams_per_day_1800s_noLiquids:
    age_PLUS_height_cm:
    - "sex"
```
In this example, there will be one model for each formula (and for each distance matrix),
which in R, would correspond to these commands:
 - `adonis(<DM> ~ timepoint_months + income, <metadata-file>)`
 - `adonis(<DM> ~ timepoint_months * income, <metadata-file>)`
 - `adonis(<DM> ~ timepoint_months + income * sex, <metadata-file>)`

## Sourcetracking

## 

## Differential abundance (using [songbird](https://github.com/biocore/songbird)) and co-occurrence (using [mmvec](https://github.com/biocore/mmvec)) 

These two methods can "communicate" using several yaml files with different options, in order to:
- [option `-mm`] correlate the songbird differentials with the PC of the mmvec embeddings.
- [option `-hlg`] highlights specific groups of featured in interactive co-occurrence heatmaps using [Xmmvec](https://github.com/FranckLejzerowicz/Xmmvec).

**Note** that both tools can use a custom training set that can be defined based
on a metadata column, using the option `-tt`. This option takes a yaml file (again!)
which is fairly simple. It tells which column name to create for each dataset, and based 
on randomly selecting which proportion of each factor in each other column, e.g.:
```
datasets:
  dataset_number_1:   # a dataset name
    traintest_sex:    # a column to create in metadata of this dataset
      - "sex"         # which columns to use for (balanced) random sampling 
  dataset_number_5:   # another dataset name
    traintest_sex:
      - "sex"
train: 0.7            # proportiton to randomly sample for training
```
These columns (here `traintest_sex`) can be now used in the parameters of 
songbird and mmvec!


### Songbird (option `-s`)

This tool runs [songbird](https://github.com/biocore/songbird), a QIIME2 
plugin for differential abundance measure. It can take as input a number of 
parameters and also, is often used to run several models for several 
datasets. This tool eases the process by reading a single configuration file 
defining all datasets to use, as well as all the filtering and samples 
subsets to make per dataset, and the parameters.

This config file must be provided in the following `.yml` format:
```
models:
  dataset_number_1:
    timeINTERsexPLUSincome: "sex+income*timepoint_months"
    sexPLUSincome: "sex+income"
  dataset_number_2:
    sexPLUSincome: "sex+income"
baselines:
  dataset_number_1:
    sexPLUSincome:
      sex: "sex"
      income: "income"
subsets:
  sex:
  - - Female
  - - Male
filtering:
  dataset_number_1:
    0.1-0.0001:
    - '0.1'
    - '0.0001'
params:
  batches:
    - 20
    - 40
  learns:
    - 1e-4
  epochs:
    - 1000
    - 2000
  thresh_feats:
    - 10
  thresh_samples:
    - 1000
  diff_priors:
    - 0.5
  train:
    - '0.7'
    - 'traintest_sex'
```

The allowed "sections" are `models`, `baselines`, `filtering`, `subsets` and 
`params`:

- `models`: for each dataset, one (or more) model name(s) and associated 
  model formula (which can accommodate categorical variables in formulation, 
  see [here](https://github.com/biocore/songbird#3-specifying-a-formula-)). 
    - In the above example, both `dataset_number_1` and `dataset_number_2` will test 
      for the model named `sexPLUSincome` by the user, which will actually 
      use the formula `sex+income`.


- `baselines`: for each dataset **and for each model name** defined in 
  `models`, one (or more) model name(s) and associated model formula (as in 
  'models'), but here to be run as baseline, i.e., for comparison. For 
  example, with a config having the following: 
  ```
  baselines:
    datasetName2:
      sexPLUSincome:
        sex: "sex"
        income: "income"
  ```
  the model `sexPLUSincome` (which formula was `sex+income`) will be 
  compared for `dataset_number_2` with both the results of model `sex` (which 
  formula is simply `sex`) and `income` (which formula is simply `income`). 
  Note that by default, the baseline is `1` and thus the "section" 
  `baselines` can be missing. It is important to know that Pseudo Q2 values 
  are only reliable to assess models compared against the same baseline 
  (this tool will reuse the same baseline model result when assessed against 
  multiple times, hence saving lots of computation time!) 


- `filtering`: for each dataset, one (or more) filtering name(s), and two 
  filtering values:
    * first: the sample prevalence threshold
    * second: the sample abundance threshold

  In both cases, the value can be between 0 and 1, which will be 
      interpreted as a fraction, e.g. `0.4` would mean _min 40%_ (of samples,
      or the reads per sample), while a value of 1 or more will be 
      interpreted as an absolute number, e.g. `10` would mean _min 10_ 
      (sample occurrences, or reads per sample).

  In the above example, only `dataset_number_1` will be filtered (`dataset_number_2` 
  will be used raw, or filtered using songbird params, see below), to keep 
  only features that have at least 0.01% of the reads of each sample 
  (`0.0001`), for 10% of the samples (`0.1`).   
  ```
  dataset_Name1:
  0.1-0.0001:
  - '0.1'
  - '0.0001'
  ```
  (For the name, I recommend using the filtering values linked by an 
  underscore.)


- `subsets`: the subsets are applied to all datasets (i.e., no sub-header per 
  dataset), e.g.:
  ```
  subsets:
    sex:
    - - Female
    - - Male
    age_cat:
    - - 'baby'
      - 'teen'
    - - '30s'
      - '40s'  
      - '50s'  
    income:
    - - '<15000'
    - - '>15000'
  ```
  which is interpreted as 6 different subsets:
  * Females only (value of `sex` is in `['Female']`
  * Males only (value of `sex` is in `['Male']`
  * Young people only (value of `age_cat` is in `['baby', 'teen']`
  * Older people only (value of `age_cat` is in `['30s', '40s', 50s']`
  * Poor people only (value of `income` is below 15,000
  * Rich people only (value of `income` is above 15,000


  The outputs will have one folder per subset, which will be named `sex_Female` 
for the first subset, etc...

- `params`: just like for "section" `subsets`, the parameters are applied to 
  all datasets (i.e., no sub-header per dataset), and the accepted 
  parameters are:
  - `train`: can be an existing metadata variable to pass to 
    `--p-training-column` (containing only `Train` and `Test` factors), or a 
    number between 0 and 1 to specify the fraction of samples to randomly 
    pick for training (default is 0.7, or 70%).
  - `batches`: `--p-batch-size`
  - `learns`: `--p-learning-rate`
  - `epochs`: `--p-epochs`
  - `diff_priors`: `--differential-prior`
  - `thresh_feats`: `--min-feature-count`
  - `thresh_samples`: `--min-sample-count`
  - `summary_interval`: `--p-summary-interval` 
  
  This "section" is where most of the combinatorial ability of this tool can 
  be leveraged, ass you can pass more than one value per parameters: each
  combination of all parameters will be run, e.g.:
  ```
  params:
    batches:
      - 20
      - 40
    train:
      - 0.6
      - 0.7
  ```
  will run **4** combinations of parameters.

**Note**: it you re-run `microbiome_analyzer`, with any config file, it will parse 
all the outputs from all configs and summarize all model performances (i.e., 
the Pseudo Q2 values) into one main table located in the `qiime/songbird` 
output folder. This table is called `songbird_q2.tsv` and it contains the 
following columns:
  - `pair`: currently only contain `unpaired` (paired datasets in dev...)
  - `dataset`: dataset name
  - `filter`: filtering name (not the filtering values, so be explicit!)
  - `subset`: samples subset (e.g. `sex_Female` fror the explanation above)
  - `model`: model name (not the model formula, so be explicit!)
  - `songbird_filter`: filtering in songbird (`f#_s#` for feature and sample)
  - `parameters`: concatenation of `batchsize`_`learnrate`_`epochs`_`diffprior`_`traintest`_`summary_interval`
  - `baseline`: baseline name of the model (not the model formula, so be explicit!)
  - `differentials`: file path to the feature differentials (main output)
  - `Pseudo_Q_squared`: performance value after cross-validation

**Note2**: The output folders will contain a `readme.txt` file explaining 
the folder name, which can be any of the "sections" setup defined in the 
config file.

### **mmvec** (option `-m`):

It is possible to run Jamie Morton's MMVEC in Qiime2 for a series of user-defined thresholds to get filter multiple
omics datasets to predict co-occurrences for (see [mmvec help page](https://github.com/biocore/mmvec)).
It is also possible to map the previous, Songbird differential ranks onto)

#### **datasets + filtering + parameters** `-m`:

A config file must be provided in the following .yml format:
```
pairs:
  2_3:
    - dataset_number_2
    - dataset_number_3*
  2_4:
   - dataset_number_2
   - dataset_number_4
  3_4:
   - dataset_number_3*
   - dataset_number_4
filtering:
  prevalence:
    - 0
    - 10
  abundance:
    - - 0
      - 0
    - - 1
      - 3
params:
  train_column:
    - 'None'
  n_examples:
    - 10
  batches:
    - 2
  learns:
    - 1e-4
  epochs:
    - 5000
  priors:
    - 0.1
    - 1
  thresh_feats:
    - 0
  latent_dims:
    - 3
```
which is interpreted as a dictionary with the folowing "sections": `pairs`, `filtering` and `params`

- `pairs`: for each named pair of datasets, the list of two datasets.
- `filtering`: for both a `prevalence` and `abundance` filter, the threshold values.
- `params`: parameters to mmvec (see [doc](https://github.com/biocore/mmvec))
 
Use a `*` character after the dataset name to indicate if it is a metabolomics dataset. 

## Usage

```
microbiome_analyzer -i <input_folder_path> -d <dataset_name> -n <project_name> -e <qiime2_env> [OPTIONS]
```

### Optional arguments

``` 
  -i, --i-datasets-folder TEXT    Path to the folder containing the sub-
                                  folders 'data' and 'metadata'.  [required]

  -d, --i-datasets TEXT           Dataset(s) identifier(s). Multiple is
                                  possible: e.g. -d dataset_number_1 and -d
                                  dataset_number_2 for
                                  'tab_dataset_number_1.tsv' and
                                  tab_dataset_number_2.tsv'.  [required]

  -n, --p-project-name TEXT       Nick name for your project.  [required]
  -e, --p-qiime2-env TEXT         name of your qiime2 conda environment (e.g.
                                  qiime2-2021.11).  [required]

  -w, --i-wol-tree TEXT           path to the tree containing the genome IDs
                                  (will check if exist in features names) (On
                                  barnacle, it is there: /projects/wol/profili
                                  ng/dbs/wol/phylogeny/tree.nwk).  [default:
                                  resources/wol_tree.nwk]

  -x, --i-sepp-tree TEXT          Qiime2 SEPP reference database to use for
                                  16S reads placement:
                                  https://docs.qiime2.org/2019.10/data-
                                  resources/#sepp-reference-databases (auto
                                  detection of datasets' tables with sequences
                                  as features).

  -q, --i-qemistree TEXT          Path to a folder containing Qemistree's
                                  feature data (named 'feature-
                                  data_<dataset_identifier>.qza'), and tree
                                  for each metabolomics dataset (named
                                  'qemistree_<dataset_identifier>.qza').

  -z, --i-classifier TEXT         Qiime2 reference taxonomic classifier
                                  database to use for 16Sreads assignment:
                                  https://docs.qiime2.org/2020.2/data-
                                  resources/#taxonomy-classifiers-for-use-
                                  with-q2-feature-classifier.

  -u, --p-run-params TEXT         server run parameters.
  -k, --p-feature-subsets TEXT    Regex to use for subsetting features (yml
                                  file).

  -t, --p-test TEXT               Groups to tests between in each PERMANOVA
                                  subset (multiple values possible, e.g. '-d
                                  sex -d age_cat').  [default: False]

  -g, --p-sample-subsets TEXT     Subsets for DMs, PCoAs, PERMANOVAs, etc (yml
                                  file).  [default: False]

  -nstd, --p-nestedness TEXT      Nestedness analysis config  (yml file).
                                  [default: False]

  -bt, --p-beta-type [permanova|anosim|permdisp]
                                  Type of beta group significance
  -prc, --p-procrustes TEXT       Pairs and subsets for procrustes/protests
                                  (yaml file).  [default: False]

  -mtl, --p-mantel TEXT           Pairs and subsets for mantel test (yaml
                                  file).  [default: False]

  -ddecay, --p-dm-decay TEXT      Parameters for (not geographic) distance
                                  decay analysis (yaml file).  [default:
                                  False]

  -gdecay, --p-geo-decay TEXT     Parameters for geographic distance decay
                                  analysis (yml file).  [default: False]

  -c, --p-collapse TEXT           Nominative or rank-based taxonmic collapse
                                  per dataset (yaml file).  [default: False]

  -tt, --p-train-test TEXT        Train test split per dataset (yaml file).
                                  [default: False]

  -a, --p-adonis TEXT             Formula for Adonis tests for each PERMANOVA
                                  subset (yaml file).  [default: False]

  -phate, --p-phate TEXT          Filters, subsets, parameters and
                                  stratifications for the PHATE latent space
                                  analysis (yaml file).  [default: False]

  -st, --p-sourcetracking TEXT    Filters, subsets, parameters and
                                  sink/sources for sourcetracking (yaml file).
                                  [default: False]

  -doc, --p-doc TEXT              Filters and subsets for the dissimilarity
                                  overlap curves analyses  (yaml file).
                                  [default: False]

  -s, --p-diff-models TEXT        Formulas for multinomial regression-based
                                  differential abundance ranking (songbird)
                                  (yaml file).  [default: False]

  -m, --p-mmvec-pairs TEXT        Pairs of datasets for which to compute co-
                                  occurrences probabilities (mmvec) (yaml
                                  file).  [default: False]

  -hlg, --p-mmvec-highlights TEXT
                                  Features to highlights on mmvec biplot (per
                                  dataset) (yaml file).  [default: False]

  -mm, --p-xmmvec TEXT            Config for Xmmvec (yaml file).  [default:
                                  False]

  -l, --p-longi-column TEXT       If data is longitudinal; provide the time
                                  metadata columnfor volatility analysis.
                                  [default: False]

  -f, --p-filter TEXT             Samples to remove, min. read abundance and
                                  feature prevalence (>1 = based on absolute
                                  reads, 0-1 = based on relative reads). (yaml
                                  file).  [default: False]

  -r, --p-rarefs TEXT             rarefaction depth per dataset (yaml file).
  -chmod, --p-chmod TEXT          Change output files permission (default =
                                  664 [= -rw-rw-r--]).  [default: 664]

  -skip, --p-skip [alpha|merge_alpha|export_alpha|alpha_correlations|alpha_group_significance|wol|taxonomy|barplot|volatility|beta|export_beta|pcoa|biplot|emperor|emperor_biplot|empress|empress_biplot|phate|doc|deicode|sepp|do_pies|alpha_kw|permanova|procrustes|mantel|decay|nestedness|adonis|songbird|mmvec|mmbird]
                                  Steps to skip (e.g. if already done or not
                                  necessary).

  -As, --p-alphas TEXT            Alpha diversity indices to use.
  -Bs, --p-betas TEXT             Beta diversity metrics to use.
  --biplot / --no-biplot          Whether to do the PCoA biplots or not
                                  [default: False]

  --force / --no-force            Force the re-writing of scripts for all
                                  commands(default is to not re-run if output
                                  file exists).  [default: False]

  --gpu / --no-gpu                Use GPUs instead of CPUs for MMVEC.
                                  [default: False]

  --raref / --no-raref            Whether to rarefy and only perform the
                                  routine analyses on the rarefied dataset(s).
                                  [default: False]

  --loc / --no-loc                whether to do compute on scratch.  [default:
                                  True]

  --filt-only / --no-filt-only    Only process the filtered version (and not
                                  also the raw) version of each dataset.
                                  [default: False]

  -filt3d, --p-filt3d TEXT        Levels for the exploration of filtering.
                                  Must be a yaml file.

  --jobs / --no-jobs              Whether to prepare Torque jobs from scripts.
                                  [default: True]

  --slurm / --no-slurm            Whether to prepare Slurm and not Torque
                                  jobs.  [default: False]

  -chunkit, --p-chunkit INTEGER   Maximum number of jobs at which extra jobs
                                  will be added in chunks

  --version                       Show the version and exit.
  --help                          Show this message and exit.
```

## Example

Once all data files imported and using the `--force` option, you can
run the following (using the files in the `tests/files` folder): 
```
microbiome_analyzer \
  -i ./microbiome_analyzer/tests/files \
  -e qiime2-2021.11 \
  -n test \
  -d dataset_number_1 \
  -d dataset_number_5 \
  -d vioscreen_foods_consumed_grams_per_day_1800s_noLiquids \
  -f ./microbiome_analyzer/tests/files/filtering.yml \
  -filt3d ./microbiome_analyzer/tests/files/filtering_3d.yml \
  -r ./microbiome_analyzer/tests/files/rarefactions.yml \
  -k ./microbiome_analyzer/tests/files/clades.yml \
  -c ./microbiome_analyzer/tests/files/collapse.yml \
  -g ./microbiome_analyzer/tests/files/subsets.yml \
  -a ./microbiome_analyzer/tests/files/adonis.yml \
  -mtl ./microbiome_analyzer/tests/files/procrustes.yml \
  -prc ./microbiome_analyzer/tests/files/procrustes.yml \
  -m ./microbiome_analyzer/tests/files/mmvec.yml \
  -s ./microbiome_analyzer/tests/files/songbird.yml \
  -ddecay ./microbiome_analyzer/tests/files/dm_decay.yml \
  -phate ./microbiome_analyzer/tests/files/phate.yml \
  -doc ./microbiome_analyzer/tests/files/doc.yml \
  -t sex -t group \
  -l timepoint_months \
  --no-jobs \
  --raref \
  --biplot \
  --force
```
The standard output shows you the scripts that have been written, i.e.,  
pretty much all the analyses available here so far: 
```
# import
sh /Users/franck/programs/microbiome_analyzer/microbiome_analyzer/tests/files/jobs/import/tst_flt_rrf.sh
# rarefy
sh /Users/franck/programs/microbiome_analyzer/microbiome_analyzer/tests/files/jobs/rarefy/tst_flt_rrf.sh
# taxonomy
sh /Users/franck/programs/microbiome_analyzer/microbiome_analyzer/tests/files/jobs/taxonomy/tst_flt_rrf.sh
# import_trees
sh /Users/franck/programs/microbiome_analyzer/microbiome_analyzer/tests/files/jobs/import_trees/tst_flt_rrf.sh
# subsets
sh /Users/franck/programs/microbiome_analyzer/microbiome_analyzer/tests/files/jobs/subsets/tst_flt_rrf.sh
# alpha
sh /Users/franck/programs/microbiome_analyzer/microbiome_analyzer/tests/files/jobs/alpha/tst_flt_rrf.sh
# alpha_correlations
sh /Users/franck/programs/microbiome_analyzer/microbiome_analyzer/tests/files/jobs/alpha_correlations/tst_flt_rrf.sh
# tabulate
sh /Users/franck/programs/microbiome_analyzer/microbiome_analyzer/tests/files/jobs/tabulate/tst_flt_rrf.sh
# alpha_rarefaction
sh /Users/franck/programs/microbiome_analyzer/microbiome_analyzer/tests/files/jobs/alpha_rarefaction/tst_flt_rrf.sh
# volatility
sh /Users/franck/programs/microbiome_analyzer/microbiome_analyzer/tests/files/jobs/volatility/tst_flt_rrf.sh
# mmvec_single_imports
sh /Users/franck/programs/microbiome_analyzer/microbiome_analyzer/tests/files/jobs/mmvec_single_imports/tst_flt_rrf.sh
# mmvec_paired_imports
sh /Users/franck/programs/microbiome_analyzer/microbiome_analyzer/tests/files/jobs/mmvec_paired_imports/tst_flt_rrf.sh
# phate
sh /Users/franck/programs/microbiome_analyzer/microbiome_analyzer/tests/files/jobs/phate/tst_flt_rrf.sh
# barplot
sh /Users/franck/programs/microbiome_analyzer/microbiome_analyzer/tests/files/jobs/barplot/tst_flt_rrf.sh
# beta
sh /Users/franck/programs/microbiome_analyzer/microbiome_analyzer/tests/files/jobs/beta/tst_flt_rrf.sh
# deicode
sh /Users/franck/programs/microbiome_analyzer/microbiome_analyzer/tests/files/jobs/deicode/tst_flt_rrf.sh
# pcoa
sh /Users/franck/programs/microbiome_analyzer/microbiome_analyzer/tests/files/jobs/pcoa/tst_flt_rrf.sh
# tsne
sh /Users/franck/programs/microbiome_analyzer/microbiome_analyzer/tests/files/jobs/tsne/tst_flt_rrf.sh
# umap
sh /Users/franck/programs/microbiome_analyzer/microbiome_analyzer/tests/files/jobs/umap/tst_flt_rrf.sh
# permanova
sh /Users/franck/programs/microbiome_analyzer/microbiome_analyzer/tests/files/jobs/permanova/tst_flt_rrf.sh
# adonis
sh /Users/franck/programs/microbiome_analyzer/microbiome_analyzer/tests/files/jobs/adonis/tst_flt_rrf.sh
# dm_decay
sh /Users/franck/programs/microbiome_analyzer/microbiome_analyzer/tests/files/jobs/dm_decay/tst_flt_rrf.sh
# dm_decay_plot
sh /Users/franck/programs/microbiome_analyzer/microbiome_analyzer/tests/files/jobs/dm_decay_plot/tst_flt_rrf.sh
# doc_R
sh /Users/franck/programs/microbiome_analyzer/microbiome_analyzer/tests/files/jobs/doc_R/tst_flt_rrf.sh
# songbird_imports
sh /Users/franck/programs/microbiome_analyzer/microbiome_analyzer/tests/files/jobs/songbird_imports/tst_flt_rrf.sh
# songbird_filter
sh /Users/franck/programs/microbiome_analyzer/microbiome_analyzer/tests/files/jobs/songbird_filter/tst_flt_rrf.sh
# songbird_baselines
sh /Users/franck/programs/microbiome_analyzer/microbiome_analyzer/tests/files/jobs/songbird_baselines/tst_flt_rrf.sh
# songbird
sh /Users/franck/programs/microbiome_analyzer/microbiome_analyzer/tests/files/jobs/songbird/tst_flt_rrf.sh
# qurro
sh /Users/franck/programs/microbiome_analyzer/microbiome_analyzer/tests/files/jobs/qurro/tst_flt_rrf.sh
```


### Bug Reports

contact `franck.lejzerowicz@gmail.com`