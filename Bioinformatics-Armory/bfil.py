import os, subprocess
import shlex
# you can install trimmomatic using conda:
# conda install -c bioconda trimmomatic
with open("rosalind_bfil.txt", "r") as inp, open("input.fastq", "w") as out:
    q = int(inp.readline().strip())
    for line in inp:
        out.write(line)
command = f"trimmomatic SE input.fastq output.fastq LEADING:{q} TRAILING:{q}"
command = shlex.split(command)
subprocess.run(command)
os.remove("input.fastq")