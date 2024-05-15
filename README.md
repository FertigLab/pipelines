# CO-SPACE

Run Cogaps followed by Spacemarkers to find spatially interacting genes in Visium data.

Locally with docker installed (see samplesheet.csv for reference):

```
nextflow run main.nf --input ./[samplesheet] -w [workdir] -profile docker
```

Other running options - check out profiles in nextflow.config.