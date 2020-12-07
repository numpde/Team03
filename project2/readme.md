## idiva: identification of disease variants

Team03 solution manual.


### Installation

Create/update/activate the conda environment `Team03`:
```{shell script}
cd solution
conda env create -f environment.yml
conda env update -f environment.yml
conda activate Team03
```


### Usage

#### Small example

The folder [solution/example_data/chrom17/input](solution/example_data/chrom17/input)
contains a small excerpt from 
the VCF files provided
(see [here](solution/example_data/chrom17/input/readme.md)).

To run the code on those files:

```
# It is important to start from this folder to resolve dependencies
cd solution

# Input path (for convenience)
I="example_data/chrom17/input"

# Output folder
O="example_data/chrom17/output"

PYTHONPATH=. python3 bin/main.py "$I/case_processed_v2.vcf" "$I/control_v2.vcf" "$O"
```

The code first runs
the whole program
on the first few datalines of the input,
and the results are written to the subfolder 
[head](solution/example_data/chrom17/output/head).
The results include a compressed VCF file
and a few figures (with the corresponding data).

Then the full files are processed
into the subfolder 
[full](solution/example_data/chrom17/output/full).



#### Your example

Mutatis mutandis.


### License

See [here](solution/license.txt).


## Team members

[RA](https://github.com/numpde/), LB, HK
