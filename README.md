# ant-asm-workflow

A template working directory (`template/` of this repository) containing scripts that can be used for semi-automated genome assembly and quality assessment with HiFi + Omni-C reads.
Currently, our goal is **haplotype-merged** assembly.
(We will offer another template for haplotype-phased assembly in the near future.)

- [Workflow doc](https://docs.google.com/document/d/12-pf9O7lHs2xxj6XQZjtEVPWICHrmc37GXqrOut-2AI/edit)
- [Stats spreadsheet](https://docs.google.com/spreadsheets/d/1d9j-m88aHG6rK9bzatN4KyjZq6pDtYoDzsfgAdjD3Uk/edit)
- [Figures and tables doc](https://docs.google.com/document/d/1RwdPGJw9Yg86xsIoVGdsSVKAvT0edEfBFF1_wNOHt7A/edit)

Below we start with an abstract description about the overall structure of the working directory. Then, we provide practical information and commands on how to work with it.

> :information_source: **IMPORTANT:**
> 
> The structure of the directories and the names of the files MUST be exactly the same as described below **except** those surrounded by `<` `>`. DO NOT change them unless you know how everything works.

> :warning: **WARNING if you wish to apply this workflow to species other than ants:**
> 
> Some values hard-coded in scripts are specific to ants (e.g. `hymenoptera_odb10` for BUSCO's lineage and `TTAGG` for telomeric motif sequence whereas human is `TTAGGG`). You need to change them as necessary for different creatures.

## Prerequisites

This workflow is supposed to be run on [the Deigo HPC cluster](https://groups.oist.jp/scs/documentation) at OIST with [the Bioinfo User Group module set](https://github.com/oist/BioinfoUgrp).
If you have an account of Deigo, you should be able to run it.

> :information_source: **IMPORTANT:**
> 
> To make the Bioinfo User Group module set available from the scripts in this workflow, you MUST add the following lines to your `$HOME/.bashrc` on Deigo:
> ```
> module use /apps/.bioinfo-ugrp-modulefiles81
> module use /apps/unit/BioinfoUgrp/DebianMed/10.7/modulefilesge
> ```

If you wish to run outside Deigo, you need to change the lines for loading modules (e.g. `ml samtools`) in the scripts as necessary.

## Directory structure

For each species, we propose the following directory structure to generate different types of assemblies.

- Those surrounded by `<` `>` are the names up to you.
- `->` means a symbolic link (made with `$ ln -s`).
- Multiple directories such as `<assembler-1>/`, `<assembler-2>/`, ..., under `01-asm/` mean multiple assemblies with different assemblers or different parameters to be compared.

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

#### `00-data/`: **Storage and QC of reads**

- Here we store input read datasets (using symlinks as described below), which must consist of the followings:
  - a single FASTQ file (not gzipped) exactly named `hifi.fastq` for HiFi; and
  - two FASTQ files (not gzipped) exactly named `omnic_R1_001.fastq` and `omnic_R2_001.fastq` for Omni-C.
- We also perform QC using GenomeScope+Smudgeplot and FASTK+GeneScope.

#### `01-asm/`: **Contig assembly**

- In each subdirectory, one of the followings is supposed to be performed:
  - a single contig assembly task (e.g. hifiasm, HiCanu) using a symlink to HiFi reads at `00-data/`; or
  - a single purge_dups task using a symlink to an assembled contig FASTA file.

#### `10-contigs/`: **Quality assessment of contigs**

- In each subdirectory, we
  - put a symlink to a single contig FASTA file at `01-asm/`, and
  - assess the quality of the contigs basically based on the criteria proposed by [VGP](https://www.nature.com/articles/s41586-021-03451-0).

#### `11-scaf/`: **Scaffolding**

- In each subdirectory, a single scaffolding task (e.g. SALSA, 3D-DNA) is performed using symlinks to Omni-C reads at `00-data/` and a smylink to a draft contig assembly at `10-contigs/`.
- Contact matrix files (i.e. .hic & .assembly files for Juicebox and .mcool & .chrom_sizes for HiGlass) are also generated for manual inspection.

#### `20-scaffolds/`: **Quality assessment of scaffolds**

- In each subdirectory, we
  - put a symlink to a single scaffold FASTA file at `11-scaf/`,
  - assess the quality of the scaffolds with the VGP criteria, and
  - generate some .bed files useful for manual curation of the assembly.

## Contents in the template directory

To automatically generate the directory structure above and to provide scripts for several tools, we have a ready-made template directory named `template/` in this GitHub repository. The contents in the template directory (except symlinks of data files) are as follows:

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
That is, you need to run the following commands step-by-step, although most of the steps can be easily executed using ready-made shell scripts.

> :information_source: **IMPORTANT:**
> 
> In the commnds below, you need to replace each name surrounded by `<` `>` (e.g. path to an input data file, name of an assembly) with something up to you. On the other hand, however, you must NOT change the names of the other files/directories.

> :information_source: **IMPORTANT:**
> 
> Variables and parameters written in the scripts can be (and should be) modified as necessary (Mandatory changes are mentioned in the comments between the commands below). We strongly recommend reading through each script before running to know what it actually does.

> :memo: **NOTE:**
> 
> `&&` at the end of a line below means "execute the next command if the previous command finishes successfully." You do not have to add it if you are running each line interactively in a shell.

### 0. Copy the template and make symlinks to input read datasets

```bash
cp -r <path-to>/ant-asm-workflow/template/ <your-working-dir-name> &&   # Assuming <your-working-dir-name> does not exist yet
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

> :memo: **NOTE:**
> 
> You can update the scripts in your working directory while keeping the other files as they are (which is useful when the workflow itself is updated) by:
> ```bash
> rsync -auv <path-to>/ant-asm-workflow/template/ <your-working-dir-name>/   # <your-working-dir-name> needs to exist
> # For example, to update scripts in `10-contigs/<assembly-name>`:
> cd <your-working-dir-name>/10-contigs/
> rsync -auv template/ <assembly-name>/
> ```
> :warning: **WARNING**: Do NOT forget `/` after `template` in the command above.

### 1. Run GenomeScope and GeneScope for each of HiFi reads and Omni-C reads

```bash
cd 00-data/ &&
    cd hifi/ &&
        sbatch run_all.sh &&   # Will submit several child jobs using SLURM
        cd ..
    cd omnic/ &&
        sbatch run_all.sh &&
        cd ..
cd ..
```

After submitting these jobs, you can **immediately** proceed to the next step 2 (i.e. You do not have to wait for them to finish).

After these jobs (and their child jobs) finish, check the GeneScope plot and the Smudgeplot for each of HiFi and Omni-C:

> :information_source: **IMPORTANT output files in this step for QC:**
> 
>  - `00-data/hifi/`
>    - `hifi.genescope/`: GeneScope for HiFi reads
>      - `summary.txt`: Fitting result
>      - `transformed_linear_plot.png`: k-mer (k=40 by default) count histogram with the fitting
>    - `hifi.genomescope_L*_U*_smudgeplot.png`: Smudgeplot for HiFi reads
>  - `00-data/omnic/`
>    - `omnic.genescope/`: GeneScope for Omni-C reads
>      - `summary.txt`: Fitting result
>      - `transformed_linear_plot.png`: k-mer (k=21 by default) count histogram with the fitting
>    - `omnic.genomescope_L*_U*_smudgeplot.png`: Smudgeplot for Omni-C reads

### 2. Run assemblers (and purge_dups if necessary)

> :warning: **WARNING about purge_dups:**
> 
> You need to run a purge_dups task (not in the directory you ran hifiasm/HiCanu/etc. but) in a separated, new subdirectory right under `01-asm/`.
> The commands below will do so.

```bash
cd 01-asm/ &&
    cd <assembler>/ &&   # NOTE: <assembler> = hifiasm, hicanu, etc. If you wish to try different parameters with the same assembler, then make a new directory.
        sbatch run_<assembler>.sh &&
        cd ..
    # To run purge_dups, do the followings after the assembly finishes
    cp -r template-purge-dups/ <assembler>-pd &&
        cd <assembler>-pd &&
        ln -sf ../<assembler>/<contig-fasta-file> ./contigs.fasta &&   # NOTE: <contig-fasta-file> depends on the assembler
        sbatch run_purge_dups_plot.sh &&
        # After run_purge_dups_plot.sh finishes, here write the values of `L`, `M`, `U` in run_purge_dups.sh based on `PB.hist.plot`
        sbatch run_purge_dups.sh &&
        cd ..
cd ..
```

### 3. Evaluate contigs

> :warning: **WARNING:**
> 
> Only a single contig FASTA file to be assessed must be placed under each subdirectory. That is, for example, if you wish to evaluate all of the i) primary contigs, ii) purged primary contigs, and iii) primary unitigs generated by the same assembler, then you need to make a distinct directory for each of them.

```bash
cd 10-contigs/ &&
    cp -r templete/ <assembly-name> &&   # NOTE: <assembly-name> = an arbitrary name representing an assembly generated in 01-asm/ (i.e. hifiasm, hifiasm-p_ctg, hicanu-pd, etc.)
        cd <assembly-name> &&
        ln -sf ../../01-asm/<assembly-dir>/<contig-fasta-file> ./contigs.fasta &&   # NOTE: <contig-fasta-file> = "purged.fa" for those after purge_dups; otherwise, something specific to an assembler
        sbatch run_all.sh &&
        cd ..
cd ..
```

After index files (`contigs.fasta.fai`, `contigs.fasta.bwt`, etc.) are created right under the `<assembly-name>` directory (i.e. after the SLURM job named `make_index` finishes), you can proceed to the next step 4 (i.e. You do not have to wait for the remainings to finish).
However, if you have many assembly candidates, then you should probably just wait for them to finish and choose only a few of them based on the quality metrics.

For each contig assembly, after `run_all.sh` finishes, check i) assembly stats such as contig N50 length, ii) BUSCO score, iii) Merqury QV, iv) mapping QV, and v) reliable block N50:

> :information_source: **IMPORTANT output files in this step for quality metrics (Use `$ cat` to see the results):**
> 
>  - `10-contigs/<assembly-name>/`
>    - `01-busco/busco.log`: BUSCO result
>    - `02-merqury/`: Merqury result
>      - `contigs.hifi.merqury.qv`: K-mer based QV
>      - `contigs.hifi.merqury.contigs.spectra-cn.fl.png`: Copy-number spectrum plot
>    - `06-mapqv/mapping.qv`: Mapping based QV
>    - `07-asset/reliable_blocks.n50`: Reliable block N50 length

### 4. Run scaffolding tools

> :warning: **WARNING:**
> 
> A single subdirectory must be only for a single scaffolding task. That is, you must NOT make nested directories like `hifiasm/salsa/` nor `hifiasm/3ddna`. Instead, you need to make `hifiasm-salsa` and `hifiasm-3ddna` right under `11-scaf/`.
> Likewise, if you wish to try different parameters with the same scaffolding tool, then make a new directory.

```bash
cd 11-scaf/ &&
    cp -r template-<scaf-tool>/ <assembly-name>-<scaf-tool> &&   # NOTE: <scaf-tool> = {3ddna, salsa}, <assembly-name> = a directory name in 10-contigs/.
        cd <assembly-name>-<scaf-tool> &&
        ln -sf ../../10-contigs/<assembly-name>/contigs.fasta* . &&   # WARNING: Do not forget * !!!
        sbatch run_<scaf-tool>.sh &&
        cd ..
cd ..
```

> :information_source: **IMPORTANT output files in this step (SALSA):**
> 
> - `11-scaf/<scaf-name>/contigs.omnic_salsa/`
>   - `scaffolds_FINAL.fasta`: Final scaffolds
>   - `scaffolds_FINAL.hic`: Omni-C contact matrix editable in Juicebox
>   - `scaffolds_FINAL.assembly`: Index file for Juicebox
>   - `scaffolds_FINAL.mcool`: Omni-C contact matrix for HiGlass
>   - `scaffolds_FINAL.chrom_sizes`: Index file for HiGlass

> :information_source: **IMPORTANT output files in this step (3D-DNA):**
> 
> - `11-scaf/<scaf-name>/scaffolding/`
>   - `contigs.FINAL.fasta`: Final scaffolds
>   - `contigs.final.hic`: Omni-C contact matrix editable in Juicebox
>   - `contigs.final.assembly`: Index file for Juicebox
>   - `contigs.final.mcool`: Omni-C contact matrix for HiGlass
>   - `contigs.final.chrom_sizes`: Index file for HiGlass
> - Be careful with the difference between `*.FINAL.*` and `*.final.*`.

### 5. Evaluate scaffolds

> :warning: **WARNING:**
> 
> Only one scaffold FASTA per subdirectory, just like above.

```bash
cd 20-scaffolds/ &&
    cp -r templete/ <scaf-name> &&   # NOTE: <scaf-name> = an arbitrary name representing an assembly generated in 11-scaf/ (i.e. hifiasm-salsa, hifiasm-pd-3ddna, etc.)
        cd <scaf-name> &&
        ln -sf ../../11-scaf/<scaf-dir>/<scaf-fasta-file> ./scaffolds.fasta &&   # NOTE: <scaf-fasta-file> depends on the scaffolding tool
        sbatch run_all.sh &&
        cd ..
cd ..
```

For each scaffold assembly, after `run_all.sh` finishes, check i) assembly stats such as scaffold N50 length, ii) BUSCO score, iii) Merqury QV, iv) mapping QV, and v) reliable block N50:

> :information_source: **IMPORTANT output files in this step for quality metrics (Use `$ cat` to see the results):**
> 
> - `20-scaffolds/<scaf-name>/`
>    - `01-busco/busco.log`: BUSCO result
>    - `02-merqury/`: Merqury result
>      - `scaffolds.hifi.merqury.qv`: K-mer based QV
>      - `scaffolds.hifi.merqury.contigs.spectra-cn.fl.png`: Copy-number spectrum plot
>    - `06-mapqv/mapping.qv`: Mapping based QV
>    - `07-asset/reliable_blocks.n50`: Reliable block N50 length

Along with these qualty metrics, we generate several files used for manual assessment and curation of the assembly (see the next section for how to do manual cuartion):
  
> :information_source: **IMPORTANT output files in this step for manual curation:**
> 
>  - `20-scaffolds/<scaf-name>/`
>    - `03-bwa/`:
>      - `scaffolds.omnic.sorted.bam`: Omni-C read mappings
>      - `scaffolds.omnic.sorted.bam.regions.bedgraph`: Mapping depth (per 1 kb by default)
>    - `04-winnowmap/`:
>      - `scaffolds.hifi.winnowmap.sorted.bam`: HiFi read mappings
>      - `scaffolds.hifi.winnowmap.sorted.bam.regions.bedgraph`: Mapping depth (per 1 kb by default)
>    - `07-asset/`:
>      - `scaffolds.gaps[.bed|.JBAT.bed]`: Locations of sequence gaps (= `N` bases) between contigs (`.JBAT.bed` is for Juicebox)
>      - `scaffolds.unreliable[.bed|.JBAT.bed]`: Locations of unreliable regions (where read coverage is very low, meaning potential misassemblies; `.JBAT.bed` is for Juicebox)
>    - `09-telomere/`:
>      - `scaffolds.telomere.filtered[.bed|.JBAT.bed]`: Locations of long tandem arrays of telomeric motifs (`.JBAT.bed` is for Juicebox)

## Visual dependencies among the files and commands in the workflow

### 1. Assembly and scaffolding

![](assets/dependency-asm-light.png)

### 2. Contig assembly evaluation

![](assets/dependency-eval-contig-light.png)

### 3. Scafold assembly evaluation

![](assets/dependency-eval-scaf-light.png)

# How to do manual curation using Juicebox (JBAT)

**Input files**:

- Hi-C contact matrix
  - `11-scaf/<scaf-name>/contigs.omnic_salsa/scaffolds_FINAL.hic` (SALSA)
  - `11-scaf/<scaf-name>/scaffolding/contigs.final.hic` (3D-DNA)
- Hi-C index file
  - `11-scaf/<scaf-name>/contigs.omnic_salsa/scaffolds_FINAL.assembly` (SALSA)
  - `11-scaf/<scaf-name>/scaffolding/contigs.final.assembly` (3D-DNA)
- Sequence gaps
  - `20-scaffolds/<scaf-name>/scaffolds.gaps.JBAT.bed`
- Unreliable regions
  - `20-scaffolds/<scaf-name>/scaffolds.unreliable.JBAT.bed`
- Long telomere arrays
  - `20-scaffolds/<scaf-name>/scaffolds.telomere.filtered.JBAT.bed`

**Steps:**

1. Load the `.hic` and `.assembly` files in JBAT.
2. Load the `.bed` files of sequence gaps, unreliable regions, and tandem arrays as 1D tracks (via `View` → `Show Annotation Panel` → `Add Local`).
3. Turn on `Enable straight edge` via Right mouse click.
4. Fix the assembly if there is a breakpoint of contacts. If the breakpoint is a sequence gap and/or unreliable region, the cutting location should be at that position in most cases.
5. After modifying the assembly, save the reviewed `.assembly` file and make a curated `.fasta` file with the following command:

```bash
ml samtools bwa Other/3d-dna
3d-dna-post-review -r <reviewed-asm>.assembly scaffolds.fasta merged.nodups.txt
```

where

- `scaffolds.fasta` = `20-scaffolds/<scaf-name>/scaffolds.fasta`
- `merged.nodups.txt` =
  - `11-scaf/<scaf-name>/contigs.omnic_salsa/aligned/merged_nodups.txt` (SALSA)
  - `11-scaf/<scaf-name>/aligned/merged_nodups.txt` (3D-DNA)

**Example:**

<img src="assets/jbat_1.png" width="500">

1. There exist two contact breakpoints at sequence gaps.

<img src="assets/jbat_2.png" width="500">

2. Cut the scaffold at the first breakpoint so that it is exactly at the sequence gap.

<img src="assets/jbat_3.png" width="500">

3. Cut the other location in the same manner.

<img src="assets/jbat_4.png" width="500">

4. Translocate the third block to the second.

<img src="assets/jbat_5.png" width="500">

5. Concatenate the blocks. Telomere motifs are now not located at the chromosome end. Telomere motifs can be found in the middle of a chromosome as well, so this might be no problem, but there also might be a further improvement.

# TODO

- [ ] Is it better to keep only HiFi reads mapped to purged contigs (rather than discarded haplotigs), before proceeding to assembly evaluation?
