# Conventions

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
- Inputs should be structured in the following order, with each of the below separated by a blank lines:
	1. Core input files that will be run through the workflow - e.g. VCFs being annotated, BAMs being analyzed etc (as well as their indexes if applicable). Also the contigs to be run on as well as the prefix.
	2. Parameters that govern how the file will be processed - e.g. prefixes, modes, input arguments to tools being called, PEDs, metadata files etc.
	3. Reference files - e.g. reference fasta, their indexes, catalogs used for annotations, etc.
	4. Runtime-related information that are not of type RuntimeAttr - e.g. docker paths, cores if applicable, sharding information if applicable.
	5. All RuntimeAttr? inputs - there should be one per task called, with its name reflective of the task's function.
- Workflows should take in an input `prefix` that is passed to every task that creates output files, which should be used in conjunction with a descriptive suffix when creating outputs.
- Workflow imports should not be renamed using the `as` operator.
- Workflows should never contain any blank comments - e.g. `#########################`.
- Workflows should never contain be any consecutive blank lines - i.e. they should have a maximum of one blank line at a time.
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
- All code should be compliant with `flake8`.


## Workspace
- All reference files - i.e. those not specific to an input callset - should be passed in via workspace data.
- All dockers should be passed in via workspace data.
