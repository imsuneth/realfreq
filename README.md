# <em>realfreq</em>
Real-time base modification frequency monitoring tool along with a real-time modification calling pipeline.

# Installation
## Building from source
```bash
git clone https://github.com/imsuneth/realfreq
cd realfreq
./scripts/install-hts.sh
make
```

Append parent directory where binary <em>realfreq</em> is in to PATH variable by appending following line to the <em>~/.bashrc</em> file.
```bash
export PATH=$PATH:/parent/dir/of/realfreq
```
Make sure to run ```source ~/.bashrc``` in order to use <em>realfreq</em> on already opened terminals.

# Usages
## Running <em>realfreq.sh</em> pipeline on a experiment directory
### 01. Execute the following command on a terminal.
```bash
./scripts/realfreq.sh -m [directory] -g [guppy_bin] -f [reference] -x [reference_index] -e [model] [options ...]
```
Make sure to set <em>[directory]</em> to the absolute path of the experiment directory set in Minknow. (ex:<em>/opt/ont/minknow/data/exp001</em>).

<em>reference_index</em> can be generated from <em>reference.fasta</em> using samtools running the command below
```bash
samtools faidx <reference.fasta>
```

Example command
```bash
export GUPPY_BIN=/tools/ont-dorado-server/bin
export REF=/data/ref/hg38noAlt.fa
export REFIDX=/data/ref/hg38noAlt.idx
export MODEL="dna_r10.4.1_e8.2_400bps_5khz_modbases_5hmc_5mc_cg_hac.cfg"

./scripts/realfreq/realfreq.sh -r -g /$GUPPY_BIN -m /data/minknow/test3 -f $REF -x $REFIDX -e $MODEL
```

## Running <em>realfreq</em> alone
<em>realfreq</em> takes the input BAM file path by reading the stdin. Therefore, BAM file path can be either piped to <em>realfreq</em> or a list of BAM file paths can be given in a file. Command for the two scenarios are as follows.

Example commands
```bash
# using pipe
echo /path/to/reads.bam | ./realfreq -r ref.fa -o freq.tsv
# input a list of bam files
./realfreq -r ref.fa -o freq.tsv < bams_list.txt
```

# <em>realfreq</em> server

<em>realfreq</em> server provides an interface to access the real-time modification frequency information using simple socket connections.

<em>realfreq.sh</em> pipeline and <em>realfreq</em> doesn't start it unless specified by setting -c \<port> or --server \<port> flag. The server supports several query commands.

Available query commands
```bash
help
    show this help message
get_contig:<contig>
    query by contig
get_range:<start_pos>:<end_pos>
    query data between start and end positions (both inclusive)
get_contig_range:<contig>:<start_pos>:<end_pos>
    query by contig and between start and end positions
get_contig_range_mod:<contig>:<start_pos>:<end_pos>:<mod_code>
    query by contig and between start and end positions and by mod code
```

Sample query commands
```bash
nc localhost 8080 <<< help
nc localhost 8080 <<< get_contig:chr1
nc localhost 8080 <<< get_range:1:100
nc localhost 8080 <<< get_contig_range:chr22:18850302:49514860
nc localhost 8080 <<< get_contig_range_mod:chr22:18850302:49514860:m
```

## <em>realfreq</em> output
### freq.tsv
Each field of <em>freq.tsv</em> is listed below with their definition.

| Field    | Definition    |
|----------|-------------|
| 1. chrom | choromosome |
| 2. start | position (0-based) of the base |
| 3. end   | position (0-based) of the base |
| 4. depth | number of reads covering the base |
| 5. n_mod | number of reads with probability (>0.2) for base modification |
| 6. n_called | number of reads called for base modification |
| 7. n_skipped | number of reads skipped the base as having less likelihood for modification |
| 8. freq | n_mod/n_called ratio |
| 9. mod_code | modification code as in [SAMtags: 1.7 Base modifications](https://github.com/samtools/hts-specs/blob/master/SAMtags.pdf) |
| 10. strand | strand (+/-) where the base modification is observed |

#### Sample freq.tsv
```
contig	start	end	strand	n_called	n_mod	freq	mod_code
chr22	20016337	20016337	+	5	0	0.000000	m
chr22	20016594	20016594	+	2	0	0.000000	m
chr22	20017045	20017045	+	1	0	0.000000	m
chr22	19970705	19970705	+	1	0	0.000000	m
chr22	19981716	19981716	+	1	1	1.000000	m
chr22	20020909	20020909	+	3	0	0.000000	m
chr22	19995719	19995719	+	4	2	0.500000	m
chr22	20017060	20017060	+	1	0	0.000000	m
chr22	19971259	19971259	+	1	1	1.000000	m
```


**Sample freq.bedmethyl**
```
chr22	20016337	20016338	m	5	+	20016337	20016337	255,0,0	5	0.000000
chr22	20016594	20016595	m	2	+	20016594	20016594	255,0,0	2	0.000000
chr22	20017045	20017046	m	1	+	20017045	20017045	255,0,0	1	0.000000
chr22	19970705	19970706	m	1	+	19970705	19970705	255,0,0	1	0.000000
chr22	19981716	19981717	m	1	+	19981716	19981716	255,0,0	1	1.000000
chr22	20020909	20020910	m	3	+	20020909	20020909	255,0,0	3	0.000000
chr22	19995719	19995720	m	4	+	19995719	19995719	255,0,0	4	0.500000
chr22	20017060	20017061	m	1	+	20017060	20017060	255,0,0	1	0.000000
chr22	19971259	19971260	m	1	+	19971259	19971259	255,0,0	1	1.000000
chr22	19973437	19973438	m	1	+	19973437	19973437	255,0,0	1	1.000000
```