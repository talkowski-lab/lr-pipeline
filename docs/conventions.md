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
- Terminology should be used consistently: `gVCF` (not `GVCF`), `contig` (not `chr` or `chrom`), `locus`/`loci` and `FORMAT`/`INFO` field names in uppercase.


## WDL
- Every file should be indented with 4 spaces per level of nesting, from the `workflow` or `task` declaration down through every nested block. Tabs should never be used, and a file should never mix indentation widths.
- Every file should begin with `version 1.0` followed by a blank line, preceded only by an exempt provenance or license header.
- Every file should end with exactly one newline, and no line should carry trailing whitespace - including blank lines inside a `command` block.
- Each task should contain exactly one `input`, `command`, `output` and `runtime` block.
- Workflows should be structured in the following order, with each of the below separated by a blank line:
  1. Imports.
  2. Documentation - a `meta` block, then a `parameter_meta` block, as described in [Workflow documentation](#workflow-documentation).
  3. Inputs.
  4. Definition of variables dynamically generated in the workflow itself.
  5. Calls to tasks.
  6. Outputs.
- Tasks should be structured in the following order, with each of the below separated by a blank line:
  1. Inputs.
  2. Definition of variables dynamically generated in the task itself.
  3. Command.
  4. Outputs.
  5. Runtime settings - default parameters, followed by a select first with the runtime override, then the actual runtime block. The `RuntimeAttr runtime_attr = select_first(...)` line should be followed immediately by `runtime {`, with no blank line between them.
- Inputs should be structured in the following order, with each of the below separated by a blank line:
  1. Core input files that will be run through the workflow - e.g. VCFs being annotated, BAMs being analyzed etc (as well as their indexes if applicable). Also the contigs to be run on as well as the prefix.
  2. Parameters that govern how the file will be processed - e.g. prefixes, modes, input arguments to tools being called, PEDs, metadata files etc.
  3. Reference files - e.g. reference fasta, their indexes, catalogs used for annotations, etc.
  4. Runtime-related information that are not of type RuntimeAttr - e.g. docker paths, cores if applicable, sharding information if applicable.
  5. All RuntimeAttr? inputs - there should be one per task called, with its name reflective of the task's function.
- Workflows should take in an input `prefix` that is passed to every task that creates output files, which should be used in conjunction with a descriptive suffix when creating outputs.
- Workflow imports should not be renamed using the `as` operator. A `call Namespace.Workflow` therefore produces a call whose name matches the namespace; this is expected for sub-workflow calls and should not be worked around.
- Call aliases introduced with `as` should be in Pascal case, matching the task and workflow naming rules - e.g. `as RunTruvari09`, not `as RunTruvari_09` or `as run_truvari_09`. Do not alias a call to the name it already has.
- A task should not declare an input it never references, except for index and companion files - e.g. `vcf_idxs` alongside `vcfs`, or `ref_fai` alongside `ref_fa`. Those inputs exist so Cromwell localizes the index next to the file it belongs to, which tools such as `bcftools concat --allow-overlaps`, `bcftools merge` and `tabix` require. They are load-bearing and must not be removed as unused.
- A workflow should not declare an input or variable it never references. Every `RuntimeAttr?` input should be passed as the `runtime_attr_override` of a call.
- A call whose outputs nothing references is acceptable only when the call exists for its failure behavior - e.g. `CheckSampleConsistency` and `ValidateContigOrder`, which gate the workflow by exiting non-zero.
- Workflows should never contain any blank comments - e.g. `#########################`.
- Workflows should never contain any consecutive blank lines outside a `command` block - i.e. they should have a maximum of one blank line at a time. Embedded Python inside a `command` block follows PEP 8 instead, so two blank lines between definitions are expected there.
- Inputs passed to a task should not have blank lines between inputs.
- The order of inputs passed to a task should reflect their order in the inputs on the workflow level.
- Inputs passed to a task should have a space on either side of the `=` character.
- The inputs section of a task should not have blank lines between inputs.
- Tasks should always have input fields `docker` and `runtime_attr_override` defined, though what is passed to each one of these when calling the task should be explicitly named - e.g. `docker = utils_docker` and `runtime_attr_override = runtime_attr_annotate_svan` respectively.
- Tasks should also have a prefix input defined, which is passed and set at the workflow level - the outputs from the task should simply use the prefix along with the file type.
- Every command block within a task should begin with `set -euo pipefail` followed by a blank line, using exactly those flags in that order. Tracing variants such as `set -euxo pipefail` and bare `set -x` should not be used.
- Every command block should use the `command <<< ... >>>` heredoc form, with the `command <<<` and `>>>` delimiters at the task's indentation, a blank line before `command <<<` and a blank line after `>>>`.
- The Bash body of a command block should be indented 4 spaces past `command <<<`. The payload of a nested heredoc such as `python3 <<'CODE'` should start at column 0, and is exempt from the indentation rules since it is Python rather than WDL.
- The default `disk_gb` for a task should be calculated dynamically based on the largest sized input file - or multiple if there are several large inputs, like multiple reference fastas or input catalogs. It should be defined in-line in the default runtime attributes section, unless it is a complicated function in which it can have a dedicated variable `disk_gb`.
- The default `mem_gb`, `boot_disk_gb` and `cpu_cores` for a task should be explicitly defined rather than based on an input file - it should be set based on the intensity of compute needed by that task.
- The default `preemptible_tries` for a task should always be 1.
- The default `max_retries` for a task should always be 0.
- The names of workflows and tasks should never include a `_` character within them - rather, they should always be in Pascal case.
- A task should not share its name with the workflow in the same file. Where a thin wrapper workflow and its task would collide, prefix the task with `Run` - e.g. workflow `TRGTLPS` with task `RunTRGTLPS`.
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
- In a task library - a file that defines tasks but no workflow, such as `Helpers.wdl` - tasks should be declared in alphabetical order, so a task can be located without searching. A provenance or license comment directly above a task belongs to that task and moves with it. Tasks in a workflow file should instead follow the order in which the workflow calls them.
- Before adding a reusable task, check existing helper modules for equivalent behavior. Generalize an existing task only when current callers' interfaces, commands, runtime behavior and outputs can be preserved or explicitly migrated.
- Do not extract inline WDL command logic into a repository script without explicit approval. Such extraction changes Docker image dependencies and deployment requirements and must be disclosed before implementation.
- Workflow file names must always match the workflow defined within them.
- Annotation workflows should always output a TSV file rather than a VCF, unless its annotations are done for every single variant in the input VCF or if the underlying workflow is designed to annotate variants in a VCF.
- The mechanically checkable rules above are enforced by `.github/scripts/check_wdl_style.py`, which should be run after editing any file under `wdl/`. The rules it cannot check - input grouping, naming verbs, runtime sizing, comment wording - still apply and are left to review.
- `miniwdl check --strict wdl/<file>.wdl` is a useful deeper audit for unused declarations and name collisions, but it is not part of CI: it also flags the index localization inputs and the sub-workflow namespace collisions described above, both of which are intentional here.


## Workflow documentation
[`docs/workflows.md`](workflows.md) is generated from the `meta` and `parameter_meta` blocks of every workflow by [`generate_workflows_doc.py`](../.github/scripts/generate_workflows_doc.py). It should never be edited by hand: write the documentation in the WDL and the document follows. The rules below are enforced by [`check_wdl_style.py`](../.github/scripts/check_wdl_style.py), so a workflow that passes the style check is guaranteed to generate a complete section.
- Every workflow should open with a `meta` block, then a `parameter_meta` block, then its `input` block, before any other statement.
- `meta` should declare exactly one key, `description`, whose value is an array of strings holding one string per paragraph - an array even when there is only one paragraph. The first paragraph should say what the workflow does; later paragraphs should cover method detail or caveats.
- **Every input and every output should have a `parameter_meta` entry**, written as `name: "Description."`, and the entries should appear in declaration order: all documented inputs in `input` order, then all outputs in `output` order.
- An entry should give only the description. The parameter's type, whether it is optional, and its default are all read from the declaration and added when the document is generated, so restating them in the description duplicates what the reader already sees.
- `RuntimeAttr?` inputs, Docker image inputs (`docker` and `*_docker`) and `String prefix` should never appear in `parameter_meta`. They are boilerplate, so the generator emits a standard bullet for each instead.
- Descriptions should be plain prose. They should never contain Markdown links, emphasis, headings, list markers or block quotes, since the generated document supplies its own structure. Refer to another workflow, task, input or file by bare name in backticks - e.g. `LongReadCNVs`.
- The one URL a description may cite is a bare parenthesized URL immediately after the first mention of an external tool - e.g. `L1ME-AID (https://github.com/Markloftus/L1ME-AID)`. Never write it as a Markdown link.
- Descriptions should use single quotes rather than escaped double quotes, which keeps the WDL string readable.
- An index input should be described as `Index for <name>.`, and a shared reference file from [references.md](references.md) as `From references.`, which the generated document expands into a pointer to that document.
- After changing any workflow's `meta`, `parameter_meta`, `input` or `output` block, run `python .github/scripts/check_wdl_style.py`. Regenerating the document is not necessary: CI regenerates and commits it on push to `main`.


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
- Indent a nested list item to its parent's content column - two spaces under a `-` parent, three under a `1.` parent - never with tabs. This is what CommonMark requires for the item to nest, and four or more spaces past the content column would turn it into a code block instead.
- Every file should end with exactly one newline, and no line should carry trailing whitespace.
- The rules above are enforced by `.github/scripts/check_markdown_style.py`, which should be run after editing any Markdown file. Fenced code blocks are exempt from the heading and list rules, since they hold code rather than Markdown.


## Workspace
- All reference files - i.e. those not specific to an input callset - should be passed in via workspace data.
- All dockers should be passed in via workspace data.
