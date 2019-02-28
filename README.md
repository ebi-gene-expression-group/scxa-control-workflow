# Single-cell Expression Atlas SMART-seq workflow

This is a Nextflow workflow which:

 * Derives configuration from an SDRF file
 * Generates quantifications using Kallisto
 * Runs downstream analysis using [Scanpy](https://scanpy.readthedocs.io/en/latest/), leveraging the [scanpy-scripts](https://github.com/ebi-gene-expression-group/scanpy-scripts) package to run individual steps of the Scanpy workflow.

It's actually a workflow-of-workflows comprising:

 * [scxa-smartseq-quantification-workflow](https://github.com/ebi-gene-expression-group/scxa-smartseq-quantification-workflow)
 * [scxa-aggregation-workflow](https://github.com/ebi-gene-expression-group/scxa-aggregation-workflow)
 * [scanpy-workflow](https://github.com/ebi-gene-expression-group/scanpy-workflow)
 * [scxa-bundle-workflow](https://github.com/ebi-gene-expression-group/scxa-bundle-workflow) - IN PROGRESS


## Setup

### Conda/ Bioconda

Workflow dependencies are managed via Conda and Bioconda, so you'll need to set that up, see instructions [here](https://bioconda.github.io/#install-conda). 

### Nextflow

Obviously you'll need Nexflow itself. If you don't have it already you can install via Conda:

```
conda install nextflow
```

You may well want to do this within a Conda environment you create for the purpose.

## Run the workflow

### Inputs

Expected inputs are:

 * A path to a directory where SDRF and IDF files are stored
 * An expriment ID for files in that directory. e.g. specifying 'foo' will imply the existence of foo.sdrf.txt and foo.idf.txt in the above directory 


### Parameters

Default base configuration is in [nextflow.config](nextflow.config). Further configuration is derived from the SDRF files and passed to child workflows.

### Execution

The workflow can be run directly from this repository:

```
nextflow run -resume scxa-smartseq-workflow --exp_name <exp name> --sdrf_dir <sdrf dir> 
```

See Nextflow documentation for further command line options.

This will download the workflow, create any necessary environments, and run the workflow with the specified inputs. Future executions will use a cached copy of the pipeline, should you wish to update the code in future, you can do so like:

```
nextflow pull ebi-gene-expression-group/scanpy-workflow
```

### Outputs

Outputs will be placed in the directory defined as SCXA_RESULTS under 'env' in nextflow.config ('results' by default). Outputs include:
