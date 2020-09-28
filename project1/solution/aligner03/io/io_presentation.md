# IO Presentation

## Input

### Read Fastq
The fastq file consists of many examples which are four lines each.
- The first line is the identifier, there does not exist an optional description, which is mentioned in the [spec](https://en.wikipedia.org/wiki/FASTQ_format).
- Raw sequence of bases, just write through.
- Only a +.
- Raw sequence of probabilities/ quality values in [phred quality](https://en.wikipedia.org/wiki/Phred_quality_score).
There are two functions, *read_fastq_sample(file)* returns only the next instance of (identifier,bases,probabilities).

*read_fastq_file(path)* uses the previous function in it's body, but reads the whole file on *path* and returns an array/list of (identifier,bases,probabilities) tuples.

Reading the whole file is obviously more resource intensive compared to only the samples.

The proposed way of handling the samples is:
- Read one sample.
- Align the sample.
- Measure how good it performs.
- Write it to file (\*).

*read_fastq_file(path)* gives a good implementation how to structure the loop so that the whole file gets read.

### Read FASTA
The [FASTA](https://en.wikipedia.org/wiki/FASTA_format) file has an identifier and then many bases which are separated by new lines for readability.

FASTA allows to store many different sequences in the same file, but in the data given by the exercise, each file only contains one sequence each.

Also there is the possibility to have many combiniations of possibilities e.g. the letter "S" stands for "C or G", thankfully only the bases "ACGTU" and the wildcard "N" occur in our data.

Each reported sequence has size of the parameter *k*, which may be small or bigger than lines in the file.
To accomodate that there is a bit more code for the first line and also an error message if the FASTA data is smaller than *k*.

As a balance of I/O reads and memory reads, a buffer is introduced, that has the size of one line. Everytime a sequence gets used the buffer pops out the next base and the current string loses the first base,
such that every returned instance has size *k*. The returned object is a tuple of (Sequence, starting_position).

Since the longer data has many lines of just "N" (Wildcard), which hold no valuable information for alignment and act as buffer, there is a process to get rid of those strings without reading the whole sequence.
As a result the returned sequence is never all "N"s.