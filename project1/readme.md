## Team03 solution manual

This is the readme file on how to run our solution.


### Humpty Dumpty

Our aligner is called Humpy Dumpty (humdum).

    Humpty Dumpty sat on a wall,
    Humpty Dumpty had a great fall;
    All the King's horses
    And all the King's men,
    Couldn't put Humpty together again.


### Usage

To put the reads 
on the provided `data_small` dataset
together again, 
run the aligner as follows:

```
cd solution

I="../input/data_small"
O="../output/data_small"

mkdir -p "${O}"

PYTHONPATH=. \
    python bin/humdum_aligner.py \
    "${I}/genome.chr22.5K.fa" \
    "${I}"/*30xCov1.fq "${I}"/*30xCov2.fq \
    > "${O}/alignment.sam"
```



## Reference aligner

To start developing and compare the quality metrics, 
we used the [bwa tool](http://bio-bwa.sourceforge.net/)
on the given data.
This is in the folder [output_ref](./output_ref). 



## Team members

RA, LB, HK
