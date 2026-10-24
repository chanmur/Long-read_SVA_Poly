# Long-read structural variant and SVA/TE annotation pipeline

A Nextflow (DSL2) pipeline that calls structural variants (SVs) from ONT long-read DNAseq BAM files with Sniffles2, merges them across all samples, and annotates inserted/deleted sequences for transposable elements (RepeatMasker with Dfam), followed by SVA classification with xTea (long-read version).

Germline and somatic (mosaic) calls are made separately for every sample. A phased mode (`--phase`) is available through the `phased` parameter.

---

## Directory layout

Overall, main.nf contains the workflows used in the pipeline. All the processes used are defined in the files under modules/. Input parameters are added to the params.json file, while the execution and environment settings are included in the nextflow.config file.
```
.
├── main.nf
├── nextflow.config
├── params.json
└── modules/
    ├── germline_call.nf
    ├── somatic_call.nf
    ├── population_call.nf
    ├── vcf_processing.nf
    └── repeat_mask.nf
```

---
## Input parameters

All parameters are set in `params.json` and passed with `-params-file params.json`. Any of them can be overridden on the command line (for example `--phased true`).

| Parameter | Type | Description |
|---|---|---|
| `out_dir` | path | Output directory. All results are published here. |
| `reference_tag` | string | Short label for the reference genome (for example `CommonRefhg38`). It is part of every output file name, and the pipeline throws an error if files from different reference tags are mixed. |
| `reference` | path | Reference genome FASTA (for example hg38). The BAM files must be aligned to this reference. |
| `phased` | boolean | `true` runs Sniffles with `--phase`; `false` runs the standard pipeline. Must be a real JSON boolean (`false`, not `"false"`), otherwise the pipeline stops at startup. |
| `metadata` | path | Sample sheet in CSV format (see below). |
| `population_tag` | string | Name of the merged cohort (for example `AllSamples`). Used in the population file names. |
| `dfam_db_dir` | path | Directory with Dfam's consensus parition files (see notes below). Mounted read-only into the RepeatMasker container at `/usr/local/share/famdb-3.0.0/Libraries/famdb`. |
| `te_library_dir` | path | Directory that contains the TE library file (see notes below). Mounted read-only into the container at `/repeat_library`. |
| `te_library` | string | File name of the TE library inside `te_library_dir` in FASTA format (see notes below). |
| `te_library_tag` | string | Short label for the TE library (for example `HsDfam40`). It is added to the RepeatMasker output names. |
| `mast2_fa` | path | MAST2 gene FASTA used by xTea (see notes below). |
| `hexamer_fa` | path | Hexamer file used by xTea (see notes below). |
| `queue_size` | integer | Maximum number of jobs Nextflow submits to SLURM at the same time (see notes below). |

### Notes:
1. The dfam consensus files can be downloaded from the dfam database. And, the corresponding TE library, as .fasta, can be created using famdb.
2. `dfam_db_dir`, or `te_library_dir`, are only required for the RepeatMasker singularity container, downloaded from Galaxy repository. By default, it lacks the famdb library files. Hence, `dfam_db_dir` and `te_library_dir` are mounted read-only into the container at `/usr/local/share/famdb-3.0.0/Libraries/famdb` and `/repeat_library`, respectively. In case the RepeatMasker container is not used, the settings and parameters need to be changed accordingly in the nextflow.config.
3. Regardless of point 2 above, `te_library_dir` and `te_library` are required for the process, masking the SVs using RepeatMasker.
4. `mast2_fa` and `hexamer_fa`, can be downloaded from the xTEA long-read GitHub page.
5. `queue_size`is included for a SLURM setup. In case a different job scheduler is used, or the absence of one, this parameter can be removed, after adjusting the corresponding settings on the nextflow.config file.

### Example `params.json`

```json
{
    "out_dir": "/path/to/results",
    "reference_tag": "CommonRefhg38",
    "reference": "/path/to/hg38.fa",
    "phased": false,
    "metadata": "/path/to/metadata.csv",
    "population_tag": "AllSamples",
    "dfam_db_dir": "/path/to/dfam40_lib",
    "te_library_dir": "/path/to/data",
    "te_library": "Homo_sapiens_Dfam.fa",
    "te_library_tag": "HsDfam40",
    "mast2_fa": "/path/to/MAST2_gene.fa",
    "hexamer_fa": "/path/to/heximer4.txt",
    "queue_size": 10
}
```

### Metadata file

One row per BAM file, with a header line. The pipeline stops if a required field is empty. The BAM and index paths are checked for existence before anything runs.

Note: Ensure that all bam files included in a given metadata file, have been
aligned to the same reference genome. This otherwise causes an error, as the
population call requires it.

| Column | Required | Description |
|---|---|---|
| `sample_name` | yes | Used in the per-sample output file names, so it must be unique per sample. |
| `sample_id` | yes | Short sample identifier, shown in task names. |
| `donor_id` | yes | Donor identifier. |
| `project` | yes | Dataset or project name. |
| `condition` | yes | Used as a sub-folder in the output (for example `control`, `case`). |
| `brain_region` | yes | Brain region of the sample. |
| `age` | yes | Age of the donor. |
| `path` | yes | Directory that contains the BAM and index files. |
| `bam_file` | yes | BAM file name inside `path`. |
| `index_file` | yes | BAM index file name (`.bai`) inside `path`. |
| `sex` | no | Sex information of the sample. Use the same format in all rows. |
| `pmi` | no | Post-mortem interval, stored in the sample metadata. Use the same unit throughout. |

Example:

```csv
sample_name,sample_id,donor_id,project,condition,brain_region,age,sex,pmi,path,bam_file,index_file
S001,ID001,D01,ProjectA,control,cortex,72,F,18,/data/bams,S001.bam,S001.bam.bai
S002,ID002,D02,ProjectA,case,cortex,68,M,22,/data/bams,S002.bam,S002.bam.bai
```

For phased mode (`phased: true`) the BAM files must be haplotagged (contain HP tags), since Sniffles `--phase` relies on them.

---

## Software and Environment

- Nextflow (DSL2), sniffles2, htslib, tabix++, bcftools, RepeatMasker, xTea (long-read version)
- Default settings in the nextflow.config file are for a SLURM system, and employing singularity containers, or environment modules where applicable. These need to be re-configured based on the execution environment and available resources.
- Singularity containers, set in `nextflow.config`: Sniffles 2.8.0, htslib 1.24 (bgzip), tabix++ 1.1.2, bcftools 1.24 and RepeatMasker 4.2.4.
- xTea (long-read version) through environment modules: `GCC/11.3.0`, `OpenMPI/4.1.4` and `xTea/v0.1.0-long`.
- A SLURM account, set globally in `nextflow.config` with `-A`.

---

## Pipeline overview

### Processes
Following are the various processes included in the pipeline. Existing workflows utilize them. The processes can also be adapted into new workflows.

| Step | Process | Tool | What it does |
|---|---|---|---|
| 1 | `GERMLINE_SV` / `GERMLINE_SV_PHASED` | Sniffles2 | Per-sample germline SV calling (`--minsvlen 10`). Writes `.vcf` and `.snf`. |
| 1 | `SOMATIC_SV` / `SOMATIC_SV_PHASED` | Sniffles2 | Per-sample mosaic SV calling (`--mosaic --minsvlen 10`). Writes `.vcf` and `.snf`. |
| 2 | `POPULATION_SV` | Sniffles2 | Merges all `.snf` files of one call type into one multi-sample VCF (one run for germline, one for somatic). |
| 3 | `COMPRESS_VCF`, `INDEX_VCF` | bgzip, tabix | Compresses and indexes the unfiltered population VCF. |
| 4 | `FILTER_VCF` | bcftools | Keeps INS/DEL with a known length of at least 50 bp and a non-symbolic ALT allele. |
| 5 | `COMPRESS_VCF_FILTERED`, `INDEX_VCF_FILTERED` | bgzip, tabix | Compresses and indexes the filtered VCF. |
| 6 | `VCF_TO_FASTA` | bcftools | Writes the inserted (ALT) or deleted (REF) sequence of each SV to a FASTA file. |
| 7 | `REPEATMASK` | RepeatMasker | Annotates the SV sequences against the Dfam library. |
| 8 | `REFINE_SVA_ANNOTATIONS` | xTea (long) | Classifies SVA insertions from the RepeatMasker output. |

The population steps run once per call type, so a standard run produces `germline` and `somatic` population results. In phased mode the call types are `germline_phased` and `somatic_phased`.

---

### Workflows
`main.nf` is built from small workflows that are chained together. The unnamed workflow acts as the entry point, checks `phased` and starts one of the two pipeline workflows.

| Workflow | Type | Inputs | What it does | Emits |
|---|---|---|---|---|
| Unnamed | Entry point | Pipeline parameters | Checks that `phased` is a boolean, loads the reference FASTA, and starts `sv_sva_pipeline` (`phased: false`) or `sv_sva_phased_pipeline` (`phased: true`). | Nothing |
| `sv_sva_pipeline` | Full pipeline | `metadata`, `ref_tag`, `reference`, `mast2_fa`, `hexamer_fa`, `population_tag`, `te_lib`, `te_lib_tag` | Runs the standard analysis from the sample sheet to the xTea SVA classification by chaining the workflows below. | Nothing |
| `sv_sva_phased_pipeline` | Full pipeline | Same as `sv_sva_pipeline` | Same chain, but uses `structural_variant_call_phased`, so Sniffles runs with `--phase`. | Nothing |
| `read_metadata` | Sub-workflow | `metadata` | Reads the sample sheet, checks the required columns, and checks that the BAM and index files exist. | `bam_files`: `(meta, bam, bai)` |
| `structural_variant_call` | Sub-workflow | `ref_tag`, `reference`, `bam_files` | Adds the reference to each sample and runs per-sample germline and somatic Sniffles calling. | `germline_vcf`, `germline_snf`, `somatic_vcf`, `somatic_snf` |
| `structural_variant_call_phased` | Sub-workflow | `ref_tag`, `reference`, `bam_files` | Same as `structural_variant_call`, with the phased processes (call types `germline_phased` and `somatic_phased`). | Same four channels |
| `population_variant_call` | Sub-workflow | `population_tag`, `ref_tag`, `reference`, `snf_files` (SNFs grouped by call type) | Stops if an SNF was made with another reference tag, then merges the SNFs of each call type into one population VCF. | `population_sv`: `(population_tag, ref_tag, sv_type, vcf)` |
| `index_vcf` | Sub-workflow | `vcf_ch`: `(..., vcf)` | Compresses and indexes the VCF (last tuple element) and keeps the other tuple elements. Used for the unfiltered population VCF. | `indexed_vcf_ch`: `(..., vcf.gz, tbi)` |
| `index_vcf_filtered` | Sub-workflow | `vcf_ch`: `(..., vcf)` | The same as `index_vcf`, but with the aliased processes `COMPRESS_VCF_FILTERED` and `INDEX_VCF_FILTERED`, so the filtered VCF can be indexed in the same run. | `indexed_vcf_ch`: `(..., vcf.gz, tbi)` |
| `filter_vcf` | Sub-workflow | `vcf_ch`: `(..., vcf)` | Filters the VCF to INS/DEL of at least 50 bp with a non-symbolic ALT and keeps the other tuple elements. | `filtered_vcf_ch`: `(..., filtered_vcf)` |
| `te_annot_sv` | Sub-workflow | `te_lib`, `te_lib_tag`, `vcf_ch`: `(..., vcf.gz, tbi)` | Writes the SV sequences to a FASTA file and annotates them with RepeatMasker. | `te_annot_sv_ch`: `(..., vcf.gz, tbi, sv_fasta, rmsk_out)` |
| `refine_sva` | Sub-workflow | `ref_tag`, `reference`, `mast2_fa`, `hexamer_fa`, `te_annot_sv_ch` | Stops on a reference tag mismatch, then classifies SVA insertions with xTea. | `refined_sva_ch`: `(..., sv_fasta, rmsk_out, xtea_sva.txt)` |

Call order inside both pipeline workflows:

```
read_metadata
  └─ structural_variant_call (or _phased)
       └─ population_variant_call
            ├─ index_vcf                      (unfiltered VCF)
            └─ filter_vcf ─ index_vcf_filtered
                 └─ te_annot_sv
                      └─ refine_sva
```

#### Notes:
1. The workflows `index_vcf`, `index_vcf_filtered`, `filter_vcf` and `te_annot_sv` are written to be reusable: they take any channel whose tuple ends in a VCF file, and the first tuple element must be unique per tuple.

2. Every file remains tagged to their corresponding metadata elements throughout the workflow, until the population call. So, custom grouping can be done on files prior to the population call. For example, files can be grouped by condition, or brain_region, or sex, prior to the population call.

3. Above workflows can be-reused or customized as required.
---


## Running the pipeline

1. Edit `params.json` with your paths.
2. Check that the config resolves as expected:
   ```bash
   nextflow config -params-file params.json
   ```
3. Start the pipeline. Run it inside `tmux` or `screen`, or as a SLURM job, because the main Nextflow process stays alive for the whole run:
   ```bash
   nextflow run main.nf -params-file params.json
   ```
4. If a run is interrupted, or you changed something downstream, continue from the cached results:
   ```bash
   nextflow run main.nf -params-file params.json -resume
   ```
---

## Output structure

```
<out_dir>/
├── 01_structural_variants/
│   ├── germline/<condition>/
│   │   ├── vcf/<sample_name>_<reference_tag>_germline.vcf
│   │   └── snf/<sample_name>_<reference_tag>_germline.snf
│   ├── somatic/<condition>/
│   │   ├── vcf/<sample_name>_<reference_tag>_somatic.vcf
│   │   └── snf/<sample_name>_<reference_tag>_somatic.snf
│   └── population/
│       ├── <population_tag>_<reference_tag>_<call_type>.vcf          merged, unfiltered
│       ├── <population_tag>_<reference_tag>_<call_type>.vcf.gz       compressed
│       ├── <population_tag>_<reference_tag>_<call_type>.vcf.gz.tbi   index
│       ├── <population_tag>_<reference_tag>_<call_type>_filtered.vcf
│       ├── <population_tag>_<reference_tag>_<call_type>_filtered.vcf.gz
│       ├── <population_tag>_<reference_tag>_<call_type>_filtered.vcf.gz.tbi
│       └── <population_tag>_<reference_tag>_<call_type>_filtered.fa  SV sequences
└── 02_te_annotation/
    ├── <...>_filtered_<te_library_tag>/       RepeatMasker output directory (.out, .tbl, .masked)
    └── <...>_filtered_xtea_sva.txt            xTea SVA classification
```

`<call_type>` is `germline` and `somatic`, or `germline_phased` and `somatic_phased` in phased mode. The `.tbl` and `.masked` RepeatMasker files may be missing when no repeats are found; the `.out` file is required.

---

## Notes

- **Failed tasks:** the config sets `errorStrategy = 'ignore'`, so a failed task does not stop the pipeline. Samples can then silently drop out of the population call. Check the `.err` files and the Nextflow report after a run.
- **Reference consistency:** the pipeline stops with an error if the reference tag of the SNF files differs from `reference_tag`.
- **Changing the reference:** set a new `reference_tag` and `reference`. The BAM files must be aligned to that reference.
- **Name uniqueness:** `sample_name` appears in file names and must not repeat in the metadata file.