bins = 100  # Number of nucleotides per "transcripts"

# Read genome
with open("hiv_genome_1M7.shape", "r") as f:
    _ = f.readline()
    line = f.readline()
    shape = line.strip().split()

# Partition it
tcnt = 0
with open("hiv_1M7_partitioned.shape", "w") as f:
    f.write("> {:d}\n".format(tcnt))

    cnt = 0
    for i in shape:

        if cnt < bins:
            f.write(i + " ")
            cnt += 1
        else:
            cnt = 0
            tcnt += 1
            f.write("\n")
            f.write("> {:d}\n".format(tcnt))
            f.write(i + " ")
