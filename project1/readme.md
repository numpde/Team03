## Team03 solution manual

TODO: This is the readme file on how to run our solution

Run aligner on `data_small`:

```
I="./input/data_small" \
O="./output/data_small" \
    python solution/bin/aligner03.py \
    ${I}/genome.chr22.5K.fa \
    ${I}/*30xCov1.fq ${I}/*30xCov2.fq \
    > ${O}/alignment.sam
```


## Reference aligner

To start developing and compare the quality metrics, 
we used the [bwa tool](http://bio-bwa.sourceforge.net/)
on the given data.
This is in the folder [output_ref](./output_ref). 



## Team03 members

TODO
