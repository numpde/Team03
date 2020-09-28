import resource
print("FASTQ_mem: %s" (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))