# Single-cell Expression Atlas control workflow for single-cell expression data analysis

This is a Nextflow workflow which:

 * Derives configuration from an SDRF file
 * Determines protocol and associated settings
 * Runs quantifications appropriate to the protocol 
 * Aggregates resulting matrices
 * If specified, Runs downstream analysis using [Scanpy](https://scanpy.readthedocs.io/en/latest/), leveraging the [scanpy-scripts](https://github.com/ebi-gene-expression-group/scanpy-scripts) package to run individual steps of the Scanpy workflow.

It's actually a workflow-of-workflows comprising:

 * [scxa-smartseq-quantification-workflow](https://github.com/ebi-gene-expression-group/scxa-smartseq-quantification-workflow)
 * [scxa-aggregation-workflow](https://github.com/ebi-gene-expression-group/scxa-aggregation-workflow)
 * Clustering using Scanpy, via Galaxy workflows encoded [here](https://github.com/ebi-gene-expression-group/scxa-workflows)
 * [scxa-bundle-workflow](https://github.com/ebi-gene-expression-group/scxa-bundle-workflow)

## Setup

### Conda/ Bioconda

Workflow dependencies are managed via Conda and Bioconda, so you'll need to set that up, see instructions [here](https://bioconda.github.io/#install-conda). 

### Nextflow

Obviously you'll need Nextflow itself. If you don't have it already you can install via Conda:

```
conda install nextflow
```

You may well want to do this within a Conda environment you create for the purpose.

### Workflows

Clone this repo into a directory called 'workflow' under the directory from where you will execute it, such that you have (RUN DIR)/workflow/scxa-control-workflow.

You will also need to clone [scxa-workflows](https://github.com/ebi-gene-expression-group/scxa-workflows), which contains child workflows as submodules, and acts as a central place to store parameters. Clone it under the same directory, so that you have (RUN DIR)/workflow/scxa-workflows. 

## Run the workflow

Routine analysis is triggered (from the above directories) like:

```
./workflow/scxa-control-workflow/bin/submitControlWorkflow.sh -t scanpy-galaxy
```

This will look for SDRF files in the directory specified by the environment variable SCXA_SDRF_DIR, triggering analyses for any new experiments found there, running quantifications via Nextflow child workflows, and tertiary analysis via the API to a Galaxy setup. 

### Outputs

Outputs will be placed in the directory defined as SCXA_RESULTS under 'env' in nextflow.config ('results' by default).
