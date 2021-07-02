# ant-asm-workflow

Scripts and template working directory for genome assembly with HiFi + Omni-C. Currently, our goal is haplotype-merged assembly.

- [Workflow google doc](https://docs.google.com/document/d/12-pf9O7lHs2xxj6XQZjtEVPWICHrmc37GXqrOut-2AI/edit)
- [Stats goole sheet](https://docs.google.com/spreadsheets/d/1d9j-m88aHG6rK9bzatN4KyjZq6pDtYoDzsfgAdjD3Uk/edit)
- [Figures and tables doc](https://docs.google.com/document/d/1RwdPGJw9Yg86xsIoVGdsSVKAvT0edEfBFF1_wNOHt7A/edit)

## Directory structure

For each species, we assume the following directory structure to generate different types of assemblies:

```
<working-dir-for-your-species>
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

Short description about the role of each directory:

- `00-data/`
  - We store input read datasets, which must consist of the followings:
    - a single FASTQ file exactly named `hifi.fastq` for HiFi; and
    - two FASTQ files exactly named `omnic_R1_001.fastq` and `omnic_R2_001.fastq` for Omni-C.
    - (NOTE: Use symlinks as described below.)
  - QC using GenomeScope+Smudgeplot and FASTK+GeneScope are performed for each of them.
- `01-asm/`
  - In each subdirectory, we perform a single contig assembly task (e.g. hifiasm, HiCanu) using symlinks to input reads at `00-data/`.
- `10-contigs/`
  - In each subdirectory,
    - A FASTA file of contigs generated by a single assembly method is placed (using symlink to `01-asm/`).
    - With the contig set, we generate index files (e.g. `samtools index` and `bwa index`) and evaluate the assembly based on the criteria proposed by [VGP](https://www.nature.com/articles/s41586-021-03451-0).
- `11-scaf/`
  - In each subdirectory, we perform a single scaffolding task (e.g. SALSA, 3D-DNA) using symlinks to input reads at `00-data/` and a draft assembly at `10-contigs/`.
- `20-scaffolds/`
  - Almost same as `10-contigs/`.
  - In each subdirectory,
    - A FASTA file of scaffolds generated by a single scaffolding method is placed (using symlink to `11-scaf/`).
    - With the scaffold set, we generate index files and evaluate the assembly.

## Using template directory

To automatically generate the directory structure above and to provide scripts to run several tools, we have a ready-made template directory named `template/`. The contents in the template directory are the following:

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

NOTE: Variables and parameters written in the scripts should be modified as necessary.

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
    ./run_fastk.sh &&   # NOTE: You can also use `sbatch` to submit the script.
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

### 2. Run assemblers (without or with purge_dups)

```bash
cd 01-asm/
cd XXX/ &&   # NOTE: XXX = hifiasm, hicanu, etc.
    ./run_XXX.sh &&
    cd ..
# To run purge_dups, do the followings (WARN: Only one purging task per subdirectory)
cp -r template-purge-dups/ XXX-pd &&   # NOTE: XXX = hifiasm, hicanu, etc.
    cd XXX-pd &&
    ln -sf ../XXX/<contig-fasta-file> ./contigs.fasta &&   # NOTE: <contig-fasta-file> depends on the assembler
    ./run_winnowmap.sh &&
    ./run_purge_dups_plot.sh &&
    # Here, replace the values of `L`, `M`, `U` in run_purge_dups.sh to the values of `-l`, `-m`, `-u` options of purge_dups, based on the plot generated by run_purge_dups_plot.sh
    ./run_purge_dups.sh &&
    cd ..
cd ..
```

### 3. Evaluate contigs

WARN: Only one contig FASTA file to be assessed must be placed for each subdirectory. For example, if you wish to evaluate both `*.p_ctg.fasta` and `p_utg.fasta` generated from the same sample with hifiasm, then you need to make a distinct directory for each of them.

```bash
cd 10-contigs/
cp -r templete/ YYY &&   # NOTE: YYY = name of directory made in step 2.1 (i.e. hifiasm, hifiasm-pd, hicanu, etc.)
    cd YYY &&
    ln -sf ../../01-asm/YYY/<contig-fasta-file> ./contigs.fasta &&   # NOTE: <contig-fasta-file> = "purged.fa" for those after purge_dups
    ./make_index.sh &&
    cd 01-busco/ && ./run_busco.sh && cd .. &&
    cd 02-merqury/ && ./run_merqury.sh && cd .. &&
    # Same applies to the remainings (i.e. 03-xxx, 04-xxx, ...)
    cd ..
cd ..
```

### 4. Run scaffolding tools

WARN: One subdirectory must be for a single scaffolding task. That is, you must NOT make nested directories like `hifiasm/salsa/` and `hifiasm/3ddna`. Instead, you need to make `hifiasm-salsa` and `hifiasm-3ddna` right under `11-scaf/`.

```bash
cd 11-scaf/
cp -r template-ZZZ/ YYY-ZZZ &&   # NOTE: ZZZ = {3ddna, salsa}, YYY = directory name in 10-contigs/ (i.e. hifiasm, hifiasm-pd, hicanu, etc.)
    cd YYY-ZZZ &&
    ln -sf ../../10-contigs/YYY/contigs.fasta . &&
    ./run_ZZZ.sh &&
    cd ..
cd ..
```

### 5. Evaluate scaffolds

WARN: One scaffold FASTA per subdirectory.

```bash
cd 20-scaffolds/
cp -r templete/ WWW &&   # NOTE: WWW = directory name in 11-scaf/ (i.e. hifiasm-salsa, hifiasm-pd-3ddna, etc.)
    cd WWW &&
    ln -sf ../../11-scaf/WWW/<scaf-fasta-file> ./scaffolds.fasta &&   # NOTE: <scaf-fasta-file> depends on the scaffolding tool
    ./make_index.sh &&
    cd 01-busco/ && ./run_busco.sh && cd .. &&
    cd 02-merqury/ && ./run_merqury.sh && cd .. &&
    # Same applies to the remainings (i.e. 03-xxx, 04-xxx, ...)
    cd ..
cd ..
```

## Yoshi TODO memo

- Create make_centromere_bed
