import gzip
import statistics
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

illumina_lines = 0
illumina_files = ['~/data/C/IL_1.fastq.gz', '~/data/C/IL_2.fastq.gz']

for file_path in illumina_files:
    full_path = os.path.expanduser(file_path)
    with gzip.open(full_path, 'rt') as f:
        for line in f:
            illumina_lines += 1

total_illumina_reads = illumina_lines / 4
print(f"Całkowita liczba odczytów Illumina: {int(total_illumina_reads)}")

nanopore_lengths = []
nanopore_gc = 0
nanopore_bases = 0

nanopore_path = os.path.expanduser('~/data/C/NP.fastq.gz')

with gzip.open(nanopore_path, 'rt') as f:
    i = 0
    for line in f:
        if i % 4 == 1:
            seq = line.strip()
            length = len(seq)
            nanopore_lengths.append(length)
            nanopore_bases += length
            nanopore_gc += seq.count('G') + seq.count('C')
        i += 1

print(f"Całkowita liczba nukleotydów NanoPore: {nanopore_bases}")

if len(nanopore_lengths) > 0:
    mean_val = statistics.mean(nanopore_lengths)
    median_val = statistics.median(nanopore_lengths)
    gc_percent = (nanopore_gc / nanopore_bases) * 100

    print(f"Średnia długość odczytu NanoPore: {mean_val:.2f}")
    print(f"Mediana długości odczytu NanoPore: {median_val}")
    print(f"Średnia zawartość %GC: {gc_percent:.2f}%")

    plt.hist(nanopore_lengths, bins=50, color='skyblue', edgecolor='black')
    plt.title("Rozkład długości odczytów NanoPore")
    plt.xlabel("Długość odczytu")
    plt.ylabel("Liczba odczytów")
    plt.savefig("histogram.png")
