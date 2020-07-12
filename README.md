# Agora: **A** **G**raph-**O**riented **R**emote homology detection and **A**linger

> The agora (/ˈæɡərə/; Ancient Greek: ἀγορά agorá) was a central public space in ancient Greek city-states. It is the best representation of city form’s response to accommodate the social and political order of the polis. The literal meaning of the word is "gathering place" or "assembly". The agora was the center of the athletic, artistic, spiritual and political life in the city. The Ancient Agora of Athens is the best-known example.

https://en.wikipedia.org/wiki/Agora

## Prerequisites

### Python >= 3.6

Agora requires the Python >= 3.6, which is available at: https://www.python.org/downloads/.

For MacOS users, the Homebrew provied a formula:

```bash
$ brew install python
```

### DELTA-BLAST

Agora uses the DELTA-BLAST internally. It's available at [NCBI](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download).

For MacOS users, the Homebrew providesa formula:

```
$ brew install blast
```

### BLAST DB

Agora requires some BLAST databases:

* `cdd_delta`, Conserved Domain Database: For DELTA-BLAST
* `pdbaa`, Protein Databank: For template database
* `uniref50`, UniRef by UniProt: For intermediate databaase

`cdd_delta` and `pdbaa` are provied by NCBI and can be retrieved by simple commands:

```bash
$ cd /data/BLASTDB   ### The directory is just an example.
$ update_blastdb.pl --decompress cdd_delta
$ update_blastdb.pl --decompress pdbaa
```

`uniref50` should be downloaded from UniProt firstly:

```bash
$ cd /data/BLASTDB
$ wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz       
$ gunzip uniref50.fasta.gz
```

The decompressed FASTA file must be converted to BLAST database format:

```bash
$ makeblastdb -dbtype prot -in uniref50.fasta -hash_index -parse_seqids -out uniref50 -title uniref50 
```

## Installation

```bash
$ git clone https://github.com/shuichiro-makigaki/agora
$ cd agora
$ pip3 install [--upgrade] [--user] -r requirement.txt
```

Set `BLASTDB` environemnt variable

```bash
$ export BLASTDB=<e.g. /data/BLASTDB>
```

## How to use

`--in-fasta` and `--out-dir` are required options. The input FASTA file can contain multiple sequences identified by `>` header. All sequences in the file will be treated as query sequences.

```bash
$ ls target1
target1/target1.fasta
$ python3 agora_cli.py iss --in-fasta target1/target.fasta --out-dir target1
```
