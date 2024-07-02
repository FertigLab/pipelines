# CO-SPACE

Run CoGAPS followed by SpaceMarkers to find spatially interacting genes in Visium data.

Locally with docker installed (see [samplesheet.csv](samplesheet.csv) for reference):

```
nextflow run main.nf --input [samplesheet] -w [workdir] -profile docker
```

Other running options - check out profiles in [nextflow.config](nextflow.config).

>[!TIP]
>Use `--max_memory --max_cpus --max_time` to control the maximum resources to grant on the host.
