# Benchmarks:

Benchmarks were run 5x each with 1 warmup run - using [hyperfine](https://github.com/sharkdp/hyperfine).

They were run on my own PC running Ubuntu 23.04:

- i7 12700k - 12 cores
- 32GB memory
- M.2 SSD drive

Hyperfine command:
```bash
hyperfine --warmup 1 -r 5 --export-json --export-markdown
```
## Unpaired

Single tests were based on the [Zymo mock community](https://github.com/LomanLab/mockcommunity) Zymo-GridION-EVEN-BB-SN dataset.

### Benchmark 1:
ONT long read sequencing, input `fastq` format, outputting a non-compressed `fastq` file. Extracting all reads classified as *Bacillus spizizenii*.

*Inputs:*

| File type | Platform | Total reads | Reads to extract | Kraken output size | Output   |
|-----------|----------|-------------|------------------|--------------------|----------|
| `.fastq`  | ONT      | 3,491,078   | 490,984          | 885MB              | `.fastq` |

*Commands run:*

| Tool | Command   |
|------|---------|
| `KrakenTools` | `extract_kraken_reads.py -s Zymo-GridION-EVEN-BB-SN.fq -k out.kraken  -o krakentools.fq -t 96241 --fastq-output` |
| `Kractor` | `kractor -i Zymo-GridION-EVEN-BB-SN.fq -k out.kraken -o kractor.fq -t 96241` |


*Results:*
| Tool | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `krakenTools` | 254.881 ± 9.482 | 242.553  | 263.158 | 1.00 |
| `Kractor` | 58.253 ± 5.651 | 48.358 | 62.222  | 4.38 ± 0.45 |

### Benchmark 2:
ONT long read sequencing, input `fastq.gz` format, outputting a non-compressed `fastq` file. Extracting all reads classified as *Bacillus spizizenii*.


*Inputs:*

| File type   | Platform | Total reads | Reads to extract | Kraken output size | Output   |
|-------------|----------|-------------|------------------|--------------------|----------|
| `.fastq.gz` | ONT      | 3,491,078   | 490,984          | 885MB              | `.fastq` |

*Commands run:*

| Tool        | Command                                                                                                             |
|-------------|---------------------------------------------------------------------------------------------------------------------|
| KrakenTools | `extract_kraken_reads.py -s Zymo-GridION-EVEN-BB-SN.fq.gz -k out.kraken  -o krakentools.fq -t 96241 --fastq-output` |
| Kractor     | `kractor -i Zymo-GridION-EVEN-BB-SN.fq.gz -k out.kraken -o kractor.fq -t 96241`                                     |


*Results:*
| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `krakenTools` | 376.592 ± 3.905 | 373.343 | 383.315 | 1.00 |
| `Kractor` | 100.044 ± 3.599 | 98.044 | 106.449 | 3.76 ± 0.14 |

## Paired

Paired end tests were based on SRA accession: SRR19995508

### Benchmark 3:
Illumina paired end sequencing, input `fastq` format, outputting a non-compressed `fastq` file.

*Inputs:*

| File type | Platform        | Total reads | Reads to extract | Kraken output size | Output   |
|-----------|-----------------|-------------|------------------|--------------------|----------|
| `.fastq`  | Illumina paired | 53,526,611  | 1,646,117        | 9.3GB              | `.fastq` |

*Commands run:*

| Tool        | Command                                                                                                                                       |
|-------------|-----------------------------------------------------------------------------------------------------------------------------------------------|
| KrakenTools | `extract_kraken_reads.py -s SRR19995508_R1.fastq -s2 SRR19995508_R2.fastq -o R1_tools.fq -o2 R2_tools.fq -k out.kraken -t 590 --fastq-output` |
| Kractor     | `kractor -i SRR19995508_R1.fastq -i SRR19995508_R2.fastq -k out.kraken -t 590 -o R1_kractor.fq -o R2_kractor.fq`                              |


*Results:*

| Command       | Mean [s]         | Min [s] | Max [s] | Relative    |
|---------------|------------------|---------|---------|-------------|
| `krakenTools` | 898.306 ± 14.203 | 884.229 | 920.653 | 1.00        |
| `Kractor`     | 94.198 ± 2.317   | 90.852  | 96.474  | 9.54 ± 0.28 |

### Benchmark 4:
Illumina paired end sequencing, input `fastq.gz` format, outputting a non-compressed `fastq` file.

*Inputs:*

| File type   | Platform        | Total reads | Reads to extract | Kraken output size | Output   |
|-------------|-----------------|-------------|------------------|--------------------|----------|
| `.fastq.gz` | Illumina paired | 53,526,611  | 1,646,117        | 9.3GB              | `.fastq` |

*Commands run:*

| Tool        | Command                                                                                                                                           |
|-------------|---------------------------------------------------------------------------------------------------------------------------------------------------|
| KrakenTools | `extract_kraken_reads.py -s SRR19995508_R1.fastq.gz -s2 SRR19995508_R2.fastq.gz -o R1_tools.fq -o2 R2_tools.fq -k out.kraken -t 2 --fastq-output` |
| Kractor     | `kractor -i SRR19995508_R1.fastq.gz -i SRR19995508_R2.fastq.gz -k out.kraken -t 2 -o R1_kractor.fq -o R2_kractor.fq`                              |


*Results:*

| Command       |          Mean [s] |  Min [s] |  Max [s] |     Relative |
|:--------------|------------------:|---------:|---------:|-------------:|
| `krakenTools` | 1033.379 ± 25.238 | 1005.720 | 1068.522 |         1.00 |
| `Kractor`     |    49.071 ± 0.179 |   48.857 |   49.334 | 21.06 ± 0.52 |