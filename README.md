# cas9_nanopore
A software package for analyzing genetic and epigenetic variations in the data of nanopore Cas9-targeted sequencing (nCATS).    
## Introduction
Cas9_nanopore package includes three bioinformatics analysis pipelines and is used to detect the variations of the regions of interest with high coverage depth in the data of nCATS, mainly involving structural variations (SVs), single nucleotide variations and 5-methylcytosine (5mC) modification.
#### Dependency
- [megalodon](https://github.com/nanoporetech/megalodon.git)
- [bcftools](https://github.com/samtools/bcftools.git)
- [NanoSV](https://github.com/mroosmalen/nanosv.git)
- [guppy](https://community.nanoporetech.com/downloads/guppy/release_notes)
- [minimap2](https://github.com/lh3/minimap2)
- [samtools](https://github.com/samtools/samtools)
- [SnpEff](https://pcingola.github.io/SnpEff/)
- [SnpSite](https://pcingola.github.io/SnpEff/ss_introduction/)

## Usage
Only one running command is needed, and various parameters are adjusted according to the data and results in time to output the analysis results of the three variations.
```
usage: cas_nanopore <subprogram> [options]

Commands:
megalodon               megalodon analysis
bcftools                bcftools analysis
nanosv                  nanosv analysis

cas_nanopore

positional arguments:
  {megalodon,bcftools,NanoSV}

optional arguments:
  -h, --help            show this help message and exit
```

### megalodon analysis
```angular2html
usage: cas_nanopore megalodon [options]

Megalodon analysis

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input the FAST5 folder
  -o OUTPUT, --output OUTPUT
                        output folder
  -g GENOME, --genome GENOME
                        reference genome sequences in fasta format
  --guppy_config GUPPY_CONFIG
                        Guppy config.Change the path in the config.ini file.
  --guppy_server_path GUPPY_SERVER_PATH
                        Absolute path to the guppy_basecall_server. The default value is in the config.ini.
  -b, --trim_barcode    Whether to remove the barcode sequence. default = False
  --remora              Whether to use the remora model. default: rerio model
  -t T                  Number of parallel processes.default:16
  --guppy_param GUPPY_PARAM
                        The default value is in the config.ini. This parameter cannot be used in remora mode
```
#### Commands
```
python /path/cas9_nanopore.py megalodon -i example_fast5/ -o rerio_output -g reference.fa -t 16
```
Or use the remora mode:
```
python /path/cas9_nanopore.py megalodon -i example_fast5/ -o remora_output -g reference.fa -t 16 -b --remora
```

### bcftools analysis
```
usage: cas_nanopore bcftools [options]

bcftools analysis

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input reads in fastq format
  -o OUTPUT, --output OUTPUT
                        Name your project output directory.
  -g GENOME, --genome GENOME
                        Reference genome sequence in FASTA format.
  --snpEff SNPEFF       The path of snpEff.jar
  --snpSift SNPSIFT     The path of snpSift.jar
  --snpEff_config SNPEFF_CONFIG
                        The path of snpEff.config
  --snpEff_mode SNPEFF_MODE
                        Verbose mode, default:GRCh38.p14
  --annotate ANNOTATE   Address of the vcf annotation file
  -t T                  Number of parallel processes. default:16
  -p                    Whether to use Porechop to remove the adapter sequence. default： False
  -q Q                  minimum mapping quality of sub-reads. default:None
  -l [MINLENGTH ...], --minLength [MINLENGTH ...]
                        Filter minimum length,Multiple values can be entered ,No filtering by default
  -qual QUAL            The vcf file minimum quality. default:20
  --minDP MINDP         The vcf file minimum depth to call variants. default：10
  --type TYPE           Type of variation analyzed. default:'snp'
```
#### Command
```
python cas9_nanopore.py bcftools -i example.fastq -o bcftools_output -g reference.fa
```

### NanoSV analysis
```
usage: cas_nanopore nanosv [options]

NanoSV analysis

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input reads in fastq format
  -o OUTPUT, --output OUTPUT
                        Name your project output directory.
  -g GENOME, --genome GENOME
                        Reference genome sequence in FASTA format.
  -t T                  Number of parallel processes. default:16
  -p                    Whether to use Porechop to remove the adapter sequence. default： False
  -q Q                  minimum mapping quality of sub-reads. default:None
  -l [MINLENGTH ...], --minLength [MINLENGTH ...]
                        Filter minimum length,Multiple values can be entered ,No filtering by default
```
#### Command
```
python cas9_nanopore.py nanosv -i example.fastq -o nanosv_output -g reference.fa
```


