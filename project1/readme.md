## Team03 solution manual

TODO: This is the readme file on how to run our solution

Run aligner on `data_small`:

```
I="./input/data_small" \
O="./output/data_small" \
    python solution/bin/humdum.py \
    ${I}/genome.chr22.5K.fa \
    ${I}/*30xCov1.fq ${I}/*30xCov2.fq \
    > ${O}/alignment.sam
```


## Humpty Dumpty

    Humpty Dumpty sat on a wall,
    Humpty Dumpty had a great fall;
    All the King's horses
    And all the King's men,
    Couldn't put Humpty together again.



## Reference aligner

To start developing and compare the quality metrics, 
we used the [bwa tool](http://bio-bwa.sourceforge.net/)
on the given data.
This is in the folder [output_ref](./output_ref). 



## Team03 members

TODO
