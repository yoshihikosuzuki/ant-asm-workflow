# ant-asm-workflow

Scripts for genome assembly with HiFi + Omni-C. Currently, our goal is haplotype-merged assembly.

- [Workflow google doc](https://docs.google.com/document/d/12-pf9O7lHs2xxj6XQZjtEVPWICHrmc37GXqrOut-2AI/edit)
- [Stats goole sheet](https://docs.google.com/spreadsheets/d/1d9j-m88aHG6rK9bzatN4KyjZq6pDtYoDzsfgAdjD3Uk/edit)
- [Figures and tables doc](https://docs.google.com/document/d/1RwdPGJw9Yg86xsIoVGdsSVKAvT0edEfBFF1_wNOHt7A/edit)

## Directory structure

For each species, we assume the following directory structure to generate different types of assemblies:

```
.
├── 00-data/
│   ├── hifi/
│   │    └── hifi.fastq
│   │
│   └── omnic/
│        ├── omnic_R1_001.fastq
│        └── omnic_R2_001.fastq
│
├── 01-asm/
│   ├── <assembly-1>/
│   ├── <assembly-2>/
|       ...
│
├── 10-contigs/
│   ├── <contig-set-1>/
│   ├── <contig-set-2>/
|       ...
│
├── 11-scaf/
│   ├── <scaffolding-1>/
│   ├── <scaffolding-2>/
|       ...
│
└── 20-scaffolds/
    ├── <scaffold-set-1>/
    ├── <scaffold-set-2>/
        ...
```

- `00-data/`
  - Input read datasets (Use symlinks), which must be:
    - a single FASTQ file named `hifi.fastq` for HiFi; and
    - two FASTQ files named `omnic_R1_001.fastq` and `omnic_R2_001.fastq` for Omni-C.
  - GenomeScope+Smudgeplot and FASTK+GeneScope are performed for each of them.
- `01-asm/`
  - In each subdirectory, we perform a contig assembly task (e.g. hifiasm, HiCanu) using symlinks to `00-data/`.
- `10-contigs/`
  - In each subdirectory,
    - A single set of contigs is placed (using symlink to `01-asm/`).
    - With the contig set, we generate index files (e.g. `samtools index` and `bwa index`) and evaluate the assembly.
- `11-scaf/`
  - In each subdirectory, we perform a scaffolding task (e.g. SALSA, 3D-DNA) using symlinks to `00-data/` and `10-contigs/`.
- `20-scaffolds/`
  - In each subdirectory,
    - A single set of scaffolds is placed (using symlink to `11-scaf/`).
    - With the scaffold set, we generate index files (e.g. `samtools index` and `bwa index`) and evaluate the assembly.

## Using template directory

To automatically generate the directory structure above and to provide scripts to run several tools, we have a ready-made template directory `template/`. The contents in the template directory look like the following:

```
template/
├── 00-data
│   ├── hifi
│   │   ├── run_fastk.sh
│   │   ├── run_genescope.sh
│   │   └── run_genomescope.sh
│   └── omnic
│       ├── run_fastk.sh
│       ├── run_genescope.sh
│       └── run_genomescope.sh
├── 01-asm
│   ├── hicanu
│   │   └── run_hicanu.sh
│   ├── hifiasm
│   │   └── run_hifiasm.sh
│   ├── ipa
│   │   └── run_ipa.sh
│   ├── peregrine
│   │   └── run_peregrine.sh
│   └── template-purge-dups
│       ├── run_purge_dups_plot.sh
│       └── run_purge_dups.sh
├── 10-contigs
│   └── template
│       ├── 01-busco
│       │   └── run_busco.sh
│       ├── 02-merqury
│       │   └── run_merqury.sh
│       ├── 04-winnowmap
│       │   └── run_winnowmap.sh
│       ├── 05-deepvariant
│       │   └── run_deepvariant.sh
│       ├── 06-mapqv
│       │   └── run_mapqv.sh
│       ├── 07-asset
│       │   └── run_asset.sh
│       ├── 08-mosdepth
│       │   └── run_mosdepth.sh
│       ├── 09-telomere
│       │   └── run_make_telomere_bed.sh
│       └── make_index.sh
├── 11-scaf
│   ├── template-3ddna
│   │   └── run_3ddna.sh
│   └── template-salsa
│       └── run_salsa.sh
└── 20-scaffolds
    └── template
        ├── 01-busco
        │   └── run_busco.sh
        ├── 02-merqury
        │   └── run_merqury.sh
        ├── 03-bwa
        │   └── run_bwa.sh
        ├── 04-winnowmap
        │   └── run_winnowmap.sh
        ├── 05-deepvariant
        │   └── run_deepvariant.sh
        ├── 06-mapqv
        │   └── run_mapqv.sh
        ├── 07-asset
        │   └── run_asset.sh
        ├── 08-mosdepth
        │   └── run_mosdepth.sh
        ├── 09-telomere
        │   └── run_make_telomere_bed.sh
        └── make_index.sh
```

## How to use

The workflow is not (and will not be) fully automatic for many reasons. You need to do the followings step-by-step.

### 0. Copy the template directory and scripts

```bash
# Copy the template directory
cp -r /path/to/ant-asm-workflow/template/ <your_working_dir> &&
    cd <your_working_dir>

# Make symlinks to input read datasets
cd 00-data/hifi/ &&
    ln -sf /path/to/<your-hifi-reads>.fastq ./hifi.fastq &&
    cd ../..
cd 00-data/omnic/ &&
    ln -sf /path/to/<your-omnic-reads>_R1_001.fastq ./omnic_R1_001.fastq &&
    ln -sf /path/to/<your-omnic-reads>_R2_001.fastq ./omnic_R2_001.fastq &&
    cd ../..
```

### 1. Run GenomeScope and GeneScope

```bash
cd 00-data/
cd hifi/ &&
    ./run_fastk.sh &&   # NOTE: You can also use `sbatch` to submit the script. Run these scripts in this order while modifying variables if necessary.
    ./run_genescope.sh &&
    ./run_genomescope.sh &&
    cd ..
cd omnic/ &&
    ./run_fastk.sh &&
    ./run_genescope.sh &&
    ./run_genomescope.sh &&
    cd ..
cd ..
```

### 2. Run assemblers

```bash
cd 01-asm/
cd XXX/ &&   # NOTE: XXX = hifiasm, hicanu, etc.
    ./run_XXX.sh &&
    cd ..
# To perform purge_dups after the assembly with XXX, run the following
cp -r template-purge-dups/ XXX-pd &&
    cd XXX-pd &&
    ln -sf ../XXX/YYY.fasta ./contigs.fasta &&   # NOTE: YYY = prefix of contig FASTA file
    ./run_purge_dups_plot.sh &&
    ./run_purge_dups.sh <L> <M> <U> &&   # NOTE: Replace <L>, <M>, <U> to the values of `-l`, `-m`, `-u` options, based on the plot generated by run_purge_dups_plot.sh
    cd ..
cd ..
```

### 3. Evaluate contigs

```bash
cd 10-contigs/
cp -r templete/ XXX &&   # NOTE: XXX = directory name in 01-asm/ (i.e. hifiasm, hifiasm-pd, hicanu, etc.)
    cd XXX &&
    ln -sf ../../01-asm/XXX/YYY.fasta ./contigs.fasta &&   # NOTE: YYY = prefix of contig FASTA file
    ./make_index.sh &&
    cd 01-busco/ && ./run_busco.sh && cd .. &&
    cd 02-merqury/ && ./run_merqury.sh && cd .. &&
    # same applies to the remainings (i.e. 03-xxx ... 09-xxx)   [TODO: make a script for this?]
    cd ..
cd ..
```

### 4. Run scaffolding tools

```bash
cd 11-scaf/
cp -r template-ZZZ/ XXX-ZZZ &&   # NOTE: ZZZ = {3ddna, salsa}, XXX = directory name in 01-asm/ (i.e. hifiasm, hifiasm-pd, hicanu, etc.)
    cd XXX-ZZZ &&
    ln -sf ../../01-asm/XXX/YYY.fasta ./contigs.fasta &&   # NOTE: YYY = prefix of contig FASTA file
    ./run_ZZZ.sh &&
    cd ..
cd ..
```

### 5. Evaluate scaffolds

```bash
cd 20-scaffolds/
cp -r templete/ WWW &&   # NOTE: WWW = directory name in 11-scaf/ (i.e. hifiasm-salsa, hifiasm-pd-3ddna, etc.)
    cd WWW &&
    ln -sf ../../11-scaf/WWW/YYY.fasta ./scaffolds.fasta &&   # NOTE: YYY = prefix of scaffold FASTA file
    ./make_index.sh &&
    cd 01-busco/ && ./run_busco.sh && cd .. &&
    cd 02-merqury/ && ./run_merqury.sh && cd .. &&
    # same applies to the remainings (i.e. 03-xxx ... 09-xxx)   [TODO: make a script for this?]
    cd ..
cd ..
```

## Yoshi TODO memo

- Write run_genescope.sh
- Create make_centromere_bed
