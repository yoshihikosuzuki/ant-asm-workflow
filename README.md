# ant-asm-workflow

Scripts and template working directory that can be used for semi-automated genome assembly with HiFi + Omni-C reads.
Currently, our goal is haplotype-merged assembly.

- [Workflow google doc](https://docs.google.com/document/d/12-pf9O7lHs2xxj6XQZjtEVPWICHrmc37GXqrOut-2AI/edit)
- [Stats goole sheet](https://docs.google.com/spreadsheets/d/1d9j-m88aHG6rK9bzatN4KyjZq6pDtYoDzsfgAdjD3Uk/edit)
- [Figures and tables doc](https://docs.google.com/document/d/1RwdPGJw9Yg86xsIoVGdsSVKAvT0edEfBFF1_wNOHt7A/edit)

Below we start with an abstract description about the overall structure of the working directory for assembly, scaffolding, and evaluation that we propose. Then, we provide practical information on how to work with the template directory containing ready-made scripts (`template/` of this repository).

**IMPORTANT: The structure of the directories and the names of the files MUST be exactly the same as described below except those surrounded by `<` `>`. DO NOT change them unless you know how everything works.**

**WARNING: Some values hard-coded in scripts in the template directory are specific to ants (e.g BUSCO's lineage is `hymenoptera_odb10` and telomeric motif sequence is `TTAGG`, not `TTAGGG`). Change them as necessary if you wish to use it for different creatures.**

## Directory structure

For each species, we assume the following directory structure to generate different types of assemblies:

```text
<your-working-dir-name>
├── 00-data/
│   ├── hifi/
│   │    └── hifi.fastq -> <path-to-your-hifi-reads>.fastq
│   │
│   └── omnic/
│        ├── omnic_R1_001.fastq -> <path-to-your-omnic-reads>_R1_001.fastq
│        └── omnic_R2_001.fastq -> <path-to-your-omnic-reads>_R2_001.fastq
│
├── 01-asm/
│   ├── <assembler-1>/
│   │    ├── hifi.fastq -> ../../00-data/hifi/hifi.fastq
│   │        ...
│   │
│   ├── <assembler-2>/
|       ...
│
├── 10-contigs/
│   ├── <assembly-name-1>/
│   │    ├── contigs.fasta -> ../../01-asm/<assembler-N>/<contig-fasta-file>
│   │        ...
│   │
│   ├── <assembly-name-2>/
|       ...
│
├── 11-scaf/
│   ├── <scaf-tool-1>/
│   │    ├── omnic_R1_001.fastq -> ../../00-data/omnic/omnic_R1_001.fastq
│   │    ├── omnic_R2_001.fastq -> ../../00-data/omnic/omnic_R2_001.fastq
│   │    ├── contigs.fasta* -> ../../10-contigs/<assembly-name-N>/contigs.fasta*
│   │        ...
│   │
│   ├── <scaf-tool-2>/
|       ...
│
└── 20-scaffolds/
    ├── <scaf-name-1>/
    │    ├── scaffolds.fasta -> ../../11-scaf/<scaf-tool-N>/<scaf-fasta-file>
    │        ...
    │
    ├── <scaf-name-2>/
        ...
```

A short description about the role of each directory is as follows:

- `00-data/`
  - We store input read datasets (using symlinks as described below), which must consist of the followings:
    - a single FASTQ file exactly named `hifi.fastq` for HiFi; and
    - two FASTQ files exactly named `omnic_R1_001.fastq` and `omnic_R2_001.fastq` for Omni-C.
  - We perform QC using GenomeScope+Smudgeplot and FASTK+GeneScope.
- `01-asm/`
  - In each subdirectory, we perform one of the followings:
    - a single contig assembly task (e.g. hifiasm, HiCanu) using a symlink to HiFi reads at `00-data/`, or
    - a single purge_dups task using a symlink to an assembled contig FASTA file.
- `10-contigs/`
  - In each subdirectory, we
    - put a single contig FASTA file (i.e. a symlink to `01-asm/`),
    - assess the quality of the contigs basically based on the criteria proposed by [VGP](https://www.nature.com/articles/s41586-021-03451-0), and
    - generate some .bed files useful for manual inspection of the assembly.
- `11-scaf/`
  - In each subdirectory, we perform a single scaffolding task (e.g. SALSA, 3D-DNA) using symlinks to Omni-C reads at `00-data/` and a draft assembly at `10-contigs/`.
- `20-scaffolds/`
  - Almost same as `10-contigs/`.
  - In each subdirectory, we
    - put a single scaffold FASTA file (i.e. a symlink to `01-asm/`),
    - assess the quality of the scaffolds basically based on the criteria proposed by [VGP](https://www.nature.com/articles/s41586-021-03451-0), and
    - generate some .bed files useful for manual curation of the scaffolds.

## Contents in the template directory

To automatically generate the directory structure above and to provide scripts for several tools, we have a ready-made template directory named `template/` in this GitHub repository. The contents in the template directory (except symlinks of data files) are the following:

```text
template/
├── 00-data
│   ├── hifi
│   │   ├── run_all.sh
│   │   ├── run_fastk.sh
│   │   ├── run_genescope.sh
│   │   └── run_genomescope.sh
│   └── omnic
│       ├── run_all.sh
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
│       ├── 00-make_index.sh
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
│       ├── 09-telomere
│       │   └── run_make_telomere_bed.sh
│       └── run_all.sh
├── 11-scaf
│   ├── template-3ddna
│   │   └── run_3ddna.sh
│   └── template-salsa
│       └── run_salsa.sh
└── 20-scaffolds
    └── template
        ├── 00-make_index.sh
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
        ├── 09-telomere
        │   └── run_make_telomere_bed.sh
        └── run_all.sh
```

## How to work with the template

The workflow is **NOT** (and will never be) fully automatic for many reasons.
That is, you need to run the followings step-by-step although most of the steps can be easily executed using ready-made shell scripts.

**IMPORTANT: In the commnds below, you need to replace each name surrounded by `<` `>` (e.g. path to an input data file, name of an assembly) with something up to you. On the other hand, however, you must NOT change the names of the other files/directories.**

**NOTE: Variables and parameters written in the scripts can be (and should be) modified as necessary (Mandatory changes are mentioned in the comments between the commands below). We strongly recommend reading through each script before running to know what it actually does.**

**NOTE: `&&` at the end of a line below means "execute the next command if the previous command finishes successfully." You do not have to add it if you are running each line interactively in a shell.**

### 0. Copy the template directory and make symlinks to input read datasets

```bash
cp -r <path-to>/ant-asm-workflow/template/ <your-working-dir-name> &&
    cd <your-working-dir-name>

cd 00-data/ &&
    cd hifi/ &&
        ln -sf <path-to-your-hifi-reads>.fastq ./hifi.fastq &&
        cd ..
    cd omnic/ &&
        ln -sf <path-to-your-omnic-reads>_R1_001.fastq ./omnic_R1_001.fastq &&
        ln -sf <path-to-your-omnic-reads>_R2_001.fastq ./omnic_R2_001.fastq &&
        cd ..
cd ..
```

### 1. Run GenomeScope and GeneScope for each of HiFi reads and Omni-C reads

```bash
cd 00-data/ &&
    cd hifi/ &&
        ./run_all.sh &&   # Will submit several jobs using `sbatch`
        cd ..
    cd omnic/ &&
        ./run_all.sh &&
        cd ..
cd ..
```

After `run_all.sh` finishes, check the GeneScope plot, GenomeScope plot, and Smudgeplot for QC of HiFi/Omni-C reads.

- Important output files for QC:
  - `00-data/hifi/`
    - `hifi.genescope/`: GeneScope for HiFi reads
      - `summary.txt`: Fitting result
      - `transformed_linear_plot.png`: k-mer (k=40 by default) count histogram with fitting
    - `hifi.genomescope_L*_U*_smudgeplot.png`: Smudgeplot for HiFi reads
  - `00-data/omnic/`
    - `omnic.genescope/`: GeneScope for Omni-C reads
      - `summary.txt`: Fitting result
      - `transformed_linear_plot.png`: k-mer (k=21 by default) count histogram with fitting
    - `omnic.genomescope_L*_U*_smudgeplot.png`: Smudgeplot for Omni-C reads

### 2. Run assemblers (and purge_dups if necessary)

**WARNING: A purge_dups task needs to be performed (not in the directory you ran hifiasm/HiCanu/etc. but) in a separated, new subdirectory right under `01-asm/`.**

```bash
cd 01-asm/ &&
    cd <assembler>/ &&   # NOTE: <assembler> = hifiasm, hicanu, etc. If you wish to try different parameters with the same assembler, then make a new directory.
        ./run_<assembler>.sh &&
        cd ..
    # To run purge_dups, do the followings after the assembly finishes
    cp -r template-purge-dups/ <assembler>-pd &&
        cd <assembler>-pd &&
        ln -sf ../<assembler>/<contig-fasta-file> ./contigs.fasta &&   # NOTE: <contig-fasta-file> depends on the assembler
        ./run_purge_dups_plot.sh &&
        # After run_purge_dups_plot.sh finishes, here write the values of `L`, `M`, `U` in run_purge_dups.sh based on `PB.hist.plot`
        ./run_purge_dups.sh &&
        cd ..
cd ..
```

### 3. Evaluate contigs

**WARNING: Only a single contig FASTA file to be assessed must be placed under each subdirectory. That is, for example, if you wish to evaluate all of the i) primary contigs, ii) purged primary contigs, and iii) primary unitigs generated from the same assembler, then you need to make a distinct directory for each of them.**

```bash
cd 10-contigs/ &&
    cp -r templete/ <assembly-name> &&   # NOTE: <assembly-name> = an arbitrary name representing an assembly generated in 01-asm/ (i.e. hifiasm, hifiasm-p_ctg, hicanu-pd, etc.)
        cd <assembly-name> &&
        ln -sf ../../01-asm/<assembly-dir>/<contig-fasta-file> ./contigs.fasta &&   # NOTE: <contig-fasta-file> = "purged.fa" for those after purge_dups; otherwise, same as that in step 2
        ./run_all.sh &&
        cd ..
cd ..
```

For each contig assembly, after `run_all.sh` finishes, check assembly stats (such as contig N50), BUSCO score, Merqury QV, mapping QV, and reliable block N50, and pick up some contig assemblies that have the highest qualities.

- Important output files for quality metrics:
  - `10-contigs/<assembly-name>/`
    - `01-busco/busco.log`: BUSCO result
    - `02-merqury/`: Merqury result
      - `contigs.hifi.merqury.qv`: K-mer based QV
      - `contigs.hifi.merqury.contigs.spectra-cn.fl.png`: Copy-number spectrum plot
    - `06-mapqv/mapping.qv`: Mapping based QV
    - `07-asset/reliable_blocks.n50`: Reliable block N50 length

### 4. Run scaffolding tools

**WARNING: A single subdirectory must be only for a single scaffolding task. That is, you must NOT make nested directories like `hifiasm/salsa/` nor `hifiasm/3ddna`. Instead, you need to make `hifiasm-salsa` and `hifiasm-3ddna` right under `11-scaf/`.**

```bash
cd 11-scaf/ &&
    cp -r template-<scaf-tool>/ <assembly-name>-<scaf-tool> &&   # NOTE: <scaf-tool> = {3ddna, salsa}, <assembly-name> = a directory name in 10-contigs/. If you wish to try different parameters with the same scaffolding tool, then make a new directory.
        cd <assembly-name>-<scaf-tool> &&
        ln -sf ../../10-contigs/<assembly-name>/contigs.fasta* . &&   # WARNING: Do not forget * !!!
        ./run_<scaf-tool>.sh &&
        cd ..
cd ..
```

- Path to the final scaffold FASTA file:
  1. SALSA: `11-scaf/<scaf-name>/contigs.omnic_salsa/scaffolds_FINAL.fasta`
  2. 3D-DNA: `11-scaf/<scaf-name>/scaffolding/contigs.FINAL.fasta`

### 5. Evaluate scaffolds

**WARNING: Only one scaffold FASTA per subdirectory, just like above.**

```bash
cd 20-scaffolds/ &&
    cp -r templete/ <scaf-name> &&   # NOTE: <scaf-name> = an arbitrary name representing an assembly generated in 11-scaf/ (i.e. hifiasm-salsa, hifiasm-pd-3ddna, etc.)
        cd <scaf-name> &&
        ln -sf ../../11-scaf/<scaf-dir>/<scaf-fasta-file> ./scaffolds.fasta &&   # NOTE: <scaf-fasta-file> depends on the scaffolding tool
        ./run_all.sh &&
        cd ..
cd ..
```

For each scaffold assembly, after `run_all.sh` finishes, check assembly stats (such as scaffold N50), BUSCO score, Merqury QV, mapping QV, and reliable block N50. Along with these quantitative qualty metrics, by looking at the Omni-C contact matrix and .bed files, choose the best scaffold assembly and proceed to manual curation using the .bed files as a guide.

- Important output files for quality metrics:
  - `20-scaffolds/<scaf-name>/`
    - `01-busco/busco.log`: BUSCO result
    - `02-merqury/`: Merqury result
      - `scaffolds.hifi.merqury.qv`: K-mer based QV
      - `scaffolds.hifi.merqury.contigs.spectra-cn.fl.png`: Copy-number spectrum plot
    - `06-mapqv/mapping.qv`: Mapping based QV
    - `07-asset/reliable_blocks.n50`: Reliable block N50 length
- Important output files for manual curation:
  - `11-scaf/<scaf-name>/`
    - If the scaffolding tool is **SALSA**:
      - `contigs.omnic.dedup.sorted.bam`: Omni-C read mappings used for scaffolding
      - `contigs.omnic_salsa/`
        - `scaffolds_FINAL.hic`: Omni-C contact matrix for Juicebox
        - `scaffolds_FINAL.assembly`: Index file for Juicebox
    - If the scaffolding tool is **3D-DNA**:
      - `hic/contigs.hic`: Initial Omni-C contact matrix using the contigs, not scaffolds
      - `references/contigs.assembly`: For Juicebox with the initial contact matrix
      - `scaffolding/`
        - `contigs.final.hic`: Omni-C contact matrix for Juicebox
        - `contigs.final.assembly`: Index file for Juicebox
  - `20-scaffolds/<scaf-name>/`
    - `scaffolds.fasta.fai`: Index file for e.g. IGV
    - `02-merqury/scaffolds_only.bed`: K-mers that appear only in scaffolds (potential misassemblies)
    - `03-bwa/`: Omni-C read mapping
      - `scaffolds.omnic.sorted.bam`: Mappings
      - `scaffolds.omnic.sorted.bam.regions.bedgraph`: Mapping depth (per 1 kb by default)
    - `04-winnowmap/`: HiFi read mapping
      - `scaffolds.hifi.winnowmap.sorted.bam`: Mappings
      - `scaffolds.hifi.winnowmap.sorted.bam.regions.bedgraph`: Mapping depth (per 1 kb by default)
    - `07-asset`: Asset result
      - `scaffolds.gaps.bed`: Gaps (= N bases) between contigs
      - `scaffolds.unreliable.bed`: Unreliable regions (= not supported by at least one of HiFi and Omni-C; unreliable regions that are not gaps are potential misassemblies)
      - `scaffolds.hic.unreliable.bed`: Regions not supported by Omni-C reads
      - `scaffolds.pb.unreliable.bed`: Regions not supported by HiFi reads
    - `09-telomere/scaffolds.telomere.bed`: Long tandem arrays of telomeric motifs

## Visual dependencies among the files and commands in the workflow

### 1. Assembly and scaffolding

![](assets/dependency-asm-light.png)

### 2. Contig assembly evaluation

![](assets/dependency-eval-contig-light.png)

### 3. Scafold assembly evaluation

![](assets/dependency-eval-scaf-light.png)
