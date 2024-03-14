# <em>realfreq</em>
Realtime methylation frequency tool.

## Installation
### Building from source
```bash
git clone https://github.com/imsuneth/realfreq
cd realfreq
./scripts/install-hts.sh
make
```

## Usage
Append parent directory where binary <em>realfreq</em> is in to PATH variable by appending following line to the <em>~/.bashrc</em> file.
```bash
export PATH=$PATH:/parent/dir/of/realfreq
```
Make sure to run ```source ~/.bashrc``` in order to use <em>realfreq</em> on already opened terminals.

### Run <em>realfreq</em> on a Minknow experiment directory
#### 01. Execute the following command on a terminal **before** starting sequencing run on Minknow.
```bash
./scripts/realfreq.sh <exp_dir>
```
Make sure to set <em><exp_dir></em> to the absolute path of the experiment directory set in Minknow. (ex:<em>/opt/ont/minknow/data/exp001</em>). Note that <em>realfreq.sh</em> will create the <em><exp_dir></em> if it doesn't exist.

#### 04. Start sequencing run on Minknow
When <em>realfreq</em> finishes processing a new batch of reads, it writes the updated methylation frequencies data to <em>freq.tsv</em> file inside <em><exp_dir></em>.

## <em>realfreq</em> output
### Fields
Each field of <em>freq.tsv</em> is listed below with their definition.

| Field    | Definition    |
|----------|:-------------:|
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

### Sample output
```
chrom	start	end	depth	n_mod	n_called	n_skipped	freq	mod_code	strand
chr22	20026776	20026776	1	1	1	0	1.000000	m	-
chr22	20016594	20016594	2	0	2	0	0.000000	m	+
chr22	20019069	20019069	1	1	1	0	1.000000	m	+
chr22	19970705	19970705	1	0	1	0	0.000000	m	+
chr22	19981716	19981716	1	1	1	0	1.000000	m	+
chr22	20020909	20020909	3	0	3	0	0.000000	m	-
chr22	19988672	19988672	2	2	2	0	1.000000	m	-
chr22	20017060	20017060	1	0	1	0	0.000000	m	+
chr22	20016854	20016854	5	0	2	0	0.000000	m	-
```