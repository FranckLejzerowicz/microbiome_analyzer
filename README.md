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
- ... _(more to come)_


## Installation

```
git clone https://github.com/FranckLejzerowicz/routine_qiime2_analyses.git
cd routine_qiime2_analyses
python setup.py build_ext --inplace --force install
```
or
```
pip install git+https://github.com/FranckLejzerowicz/routine_qiime2_analyses.git
```
*_Note that python should be python3_

## Input

The input is strict (so that it runs well):
- path to a folder named `XX_datasets` containing the features tables (option `-i`), which must contain files:
  - starting with `tab_` and ending either with:
    - `.tsv` (for tables that are tab-separated plain text with samples as columns and features as rows)
    - `.biom` (for biom tables possibly generated from these .tsv tables or fetched using e.g. [redbiom](https://github.com/biocore/redbiom))
- next to this folder must exist another folder named `00_metadata`,  which must contain metadata table files:
  - starting with `meta_` and ending with `.tsv`
  - samples must be as row and metadata variables as columns (first column is for the `sample_name`)     
- there must be a perfect matching of the internal names of the features tables _vs._ metadata tables, e.g.
    ```
    projectname_folder
    ├── 00_metadata
    │   └── meta_mGtax_gOTU_382s.tsv
    ├── XX_datasets
    │   └── tab_mGtax_gOTU_382s.tsv
    ```
 
    In this case, the analysis is done as follows:
    ```
    routine_qiime2_analyses -i projectname_folder -n projectname -e qiime2-2019.10
    ```

In fact, the tool simply generates the scripts to by started manually and which output should be scrutinized manually,
as a way to help user have default qiime command lines written in Torque's / Slurm's scripts ready to be run on a HPC.

For example, after running the above command, one would obtain the following files.

## Outputs

  
## Usage

```
routine_qiime2_analyses -i <input_folder_path> -o <output_path> -e <conda_env> [OPTIONS]
```

### Optional arguments

``` 
  -i, --i-datasets-folder TEXT  Path to the folder containing the .tsv
                                datasets  [required]
  -t, --i-wol-tree TEXT         default on barnacle /projects/wol/profiling/db
                                s/wol/phylogeny/tree.nwk  [default: /projects/
                                wol/profiling/dbs/wol/phylogeny/tree.nwk]
  -n, --p-project-name TEXT     Nick name for your project.  [required]
  -e, --p-qiime2-env TEXT       name of your qiime2 conda environment (e.g.
                                qiime2-2019.10).  [required]
  -p, --p-perm-subsets TEXT     Subsets for PERMANOVA.  [default: False]
  -g, --p-perm-groups TEXT      Groups to test between in each PERMANOVA
                                subset. Must be a yaml file, e.g.
                                (see example
                                in 'example_PERMANOVA.yml' and the README)
                                [default: False]
  -l, --p-longi-column TEXT     If data is longitudinal; provide the time
                                metadata columnfor volatility analysis.
                                [default: False]
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


### Bug Reports

contact `flejzerowicz@health.ucsd.edu`