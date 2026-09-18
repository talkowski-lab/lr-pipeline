# Conventions
This document defines the comment, WDL and Python style conventions to follow when writing or editing code in this repository.


## Comments
These rules apply to WDL, Python, R and Bash alike, including code embedded in a WDL `command` block.
- Comments should mark key sections of code and explain at a high level what that section does. They should never restate what a line of code already says, nor annotate individual lines.
- A comment should be a single line. Use a second line only for a genuine gotcha that cannot be compressed without losing the reason it exists - for example an upstream tool quirk or a correctness constraint that is not visible in the code.
- Comments should be imperative and verb-first - e.g. `# Concatenate the TSVs, keeping a single header when requested`, not `# Concatenation` or `# This concatenates the TSVs`.
- Comments should start with a capital letter and should not end with a period. Multi-sentence comments should be rewritten as a single statement.
- Avoid vague verbs such as `Handle`, `Process` and `Run script to`, matching the workflow naming rules below. Name the actual operation.
- Comments should not exceed 120 characters, should be ASCII only, and should never use decorative characters such as `->`, `→` or ` -- `. Prefer `and`, `so`, `because` or a semicolon.
- Comments should be placed on their own line directly above the code they describe, never trailing at the end of a line of code. The exception is a machine-read directive such as `# noqa: E302`.
- Commented-out code should be deleted rather than left in place.
- Comments inside a WDL `command` block must never contain `~{`, `${` or backticks. Cromwell and Bash expand these even within a comment.
- License, copyright and upstream-provenance headers are exempt from all of the above and must be preserved verbatim. Provenance headers use the form `# Derived from <repo> <path or URL>`.
- Terminology should be used consistently: `gVCF` (not `GVCF`), `GLNexus`, `contig` (not `chr` or `chrom`), `locus`/`loci`, `FORMAT`/`INFO` field names in uppercase, and `chrX`/`chrY` rather than `chrX/Y`.


## WDL
- Workflows should be structured in the following order, with each of the below separated by a blank line:
	1. Imports.
	2. Inputs.
	3. Definition of variables dynamically generated in the workflow itself.
	4. Calls to tasks.
	5. Outputs.
- Tasks should be structured in the following order, with each of the below separated by a blank line:
	1. Inputs.
	2. Definition of variables dynamically generated in the task itself.
	3. Command.
	4. Outputs.
	5. Runtime settings - default parameters, followed by a select first with the runtime override, then the actual runtime block.
- Inputs should be structured in the following order, with each of the below separated by a blank line:
	1. Core input files that will be run through the workflow - e.g. VCFs being annotated, BAMs being analyzed etc (as well as their indexes if applicable). Also the contigs to be run on as well as the prefix.
	2. Parameters that govern how the file will be processed - e.g. prefixes, modes, input arguments to tools being called, PEDs, metadata files etc.
	3. Reference files - e.g. reference fasta, their indexes, catalogs used for annotations, etc.
	4. Runtime-related information that are not of type RuntimeAttr - e.g. docker paths, cores if applicable, sharding information if applicable.
	5. All RuntimeAttr? inputs - there should be one per task called, with its name reflective of the task's function.
- Workflows should take in an input `prefix` that is passed to every task that creates output files, which should be used in conjunction with a descriptive suffix when creating outputs.
- Workflow imports should not be renamed using the `as` operator.
- Workflows should never contain any blank comments - e.g. `#########################`.
- Workflows should never contain any consecutive blank lines - i.e. they should have a maximum of one blank line at a time.
- Inputs passed to a task should not have blank lines between inputs.
- The order of inputs passed to a task should reflect their order in the inputs on the workflow level.
- Inputs passed to a task should have a space on either side of the `=` character.
- The inputs section of a task should not have blank lines between inputs.
- Tasks should always have input fields `docker` and `runtime_attr_override` defined, though what is passed to each one of these when calling the task should be explicitly named - e.g. `docker = utils_docker` and `runtime_attr_override = runtime_attr_annotate_svan` respectively.
- Tasks should also have a prefix input defined, which is passed and set at the workflow level - the outputs from the task should simply use the prefix along with the file type.
- Every command block within a task should begin with `set -euo pipefail` followed by a blank line.
- The default `disk_gb` for a task should be calculated dynamically based on the largest sized input file - or multiple if there are several large inputs, like multiple reference fastas or input catalogs. It should be defined in-line in the default runtime attributes section, unless it is a complicated function in which it can have a dedicated variable `disk_gb`.
- The default `mem_gb`, `boot_disk_gb` and `cpu_cores` for a task should be explicitly defined rather than based on an input file - it should be set based on the intensity of compute needed by that task.
- The default `preemptible_tries` for a task should always be 1.
- The default `max_retries` for a task should always be 0.
- The names of workflows and tasks should never include a `_` character within them - rather, they should always be in Pascal case.
- Top-level workflow names should describe their primary operation using consistent action verbs:
	- `Annotate` adds fields or tags to existing records.
	- `Create` derives a new artifact, such as a matrix, metadata table, interval file, plot set or summary.
	- `Convert` changes the representation or file format of existing data.
	- `Extract` emits selected records separately, while `Subset` retains the same representation with fewer records or columns.
	- `Concatenate` joins ordered, non-overlapping shards or contigs.
	- `Merge` reconciles files, records or sample sets of the same kind.
	- `Combine` applies domain-specific logic across caller-specific or type-specific inputs.
	- `Integrate` combines different variant classes into a unified callset.
	- `Normalize`, `Resolve`, `Filter`, `Summarize` and `Evaluate` should be used when they accurately describe the primary operation.
- Do not use `Create` merely because a workflow produces an output; choose the verb that describes its primary operation.
- Thin tool wrappers should use the underlying tool name. Add a purpose suffix when multiple wrappers or orchestration variants exist, such as `PALMERAssembly` and `PALMERDiploid`.
- Avoid vague verbs such as `Generate`, `Process` and `Plot` when a more specific operation or output can be named.
- The names of inputs, variables and outputs should include a `_` to separate words, and be entirely lowercase unless they refer to a noun that is capitalized (e.g. PALMER or L1MEAID) - i.e. they should always be in snake case.
- There should never be any additional indentation in order to better align parts of the code to the length or horizontal/vertical spacing of other components in its section - indentation should only be applied at the start of a line.
- All mentions of `fasta` should instead use `fa` - e.g. `ref_fa` instead of `ref_fasta`.
- All mentions of `fasta_index`, `fasta_fai` or `fa_fai` should instead use `fai` - e.g. `ref_fai` instead of `ref_fasta_index`, `ref_fasta_fai` or `ref_fa_fai`.
- All mentions of `vcf_index` or `vcf_tbi` should instead use `vcf_idx`.
- All VCFs should have suffix `_vcf`, and be coupled with a VCF index file that has a suffix `_vcf_idx`.
- Reusable tasks should live in `Helpers.wdl` and be imported rather than duplicated across workflows.
- Before adding a reusable task, check existing helper modules for equivalent behavior. Generalize an existing task only when current callers' interfaces, commands, runtime behavior and outputs can be preserved or explicitly migrated.
- Do not extract inline WDL command logic into a repository script without explicit approval. Such extraction changes Docker image dependencies and deployment requirements and must be disclosed before implementation.
- Workflow file names must always match the workflow defined within them.
- Annotation workflows should always output a TSV file rather than a VCF, unless its annotations are done for every single variant in the input VCF or if the underlying workflow is designed to annotate variants in a VCF.


## Python
- All code should be compliant with `flake8`, using the repository configuration: a 130-character maximum line length and `E203` and `W503` ignored.
- Executable scripts should use `#!/usr/bin/env python3`.
- Command-line scripts should use `argparse` to define and parse their arguments.
- Scripts with an application entry point should define a `main()` function and call it from an `if __name__ == "__main__":` guard.
- Functions and variables should use snake_case; module-level constants should use uppercase names.
- File I/O should use context managers.
- Use f-strings for string interpolation.
- `pyproject.toml` sets Black's line length to 88. Use Black when formatting a script, but do not assume existing scripts are already Black-formatted.


## Markdown
- Do not leave blank lines between a heading and its content.
- Use two blank lines before every `##` heading and before the first `###` heading in a `##` section.
- Use one blank line between subsequent `###` sections.
- Do not use thematic section dividers such as `---`.


## Workspace
- All reference files - i.e. those not specific to an input callset - should be passed in via workspace data.
- All dockers should be passed in via workspace data.
