# Team03 solution manual

### Humpty Dumpty

Our aligner is called Humpy Dumpty (humdum).

    Humpty Dumpty sat on a wall,
    Humpty Dumpty had a great fall;
    All the King's horses
    And all the King's men,
    Couldn't put Humpty together again.


### Usage

#### Creating an index and aligning reads

To put the reads 
on the provided **small dataset**
together again, 
run the aligner as follows:

```
cd solution

I="input/data_small"
O="output/data_small"

mkdir -p "${O}"

PYTHONPATH=. \
    python3 bin/humdum_aligner.py \
    "${I}"/genome*.fa \
    "${I}"/*30xCov1.fq "${I}"/*30xCov2.fq \
    > "${O}/alignment.sam"
```

This will first create an index file
in the same folder as the reference genome,
then output the SAM file.

The aligner reports exactly one alignment per read.
This alignment may be poor.


Run on the **large dataset** as follows.
It will take some 15min to create the index
before outputting the alignments.

```{shell script}
cd solution

I="../input/data"
O="../output/data"
Cov=30

mkdir -p "${O}"

PYTHONPATH=. \
    python3 bin/humdum_aligner.py \
    "${I}"/genome*.fa.gz \
    "${I}"/*${Cov}xCov1.fq.gz "${I}"/*${Cov}xCov2.fq.gz \
    > "${O}/alignment.sam"
```



#### Metrics

To analyze the resulting SAM file, run:

```{shell script}
cd solution

O="../output/data_small"

PYTHONPATH=. \
    python3 bin/humdum_qc.py \
    "${O}"/*.sam "${O}"
```

This will create diagnostic files in 
the output folder.

Mutatis mutandis for the large dataset. 

## Reference aligner

To start developing and compare the quality metrics, 
we used the [bwa tool](http://bio-bwa.sourceforge.net/)
on the given data.
This is in the folder [output_ref](./output_ref). 


## Team members

[RA](https://github.com/numpde/), LB, HK
