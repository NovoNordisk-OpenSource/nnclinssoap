# nnclinssoap

**Novo Nordisk Clinical Single-cell Spatial Omics Analytical Pipeline**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.10.0-brightgreen.svg)](https://www.nextflow.io/)
[![Singularity](https://img.shields.io/badge/singularity-%E2%89%A53.8-blue.svg)](https://sylabs.io/singularity/)
[![R](https://img.shields.io/badge/R-4.5.2-blue.svg)](https://www.r-project.org/)

---

## Overview

A comprehensive analysis framework for spatial transcriptomics data, with primary focus on 10x Genomics Xenium platform. The repository provides both Nextflow-based and R-based containerized workflows for processing and analyzing spatial single-cell data.

### Repository Structure

```
nnclinssoap/
├── docker/                    # Container definitions (Dockerfiles & Apptainer .def)
│   ├── preprocessing/         # QC, doublet detection, ambient RNA removal
│   ├── analysis/              # Cell annotation, differential expression
│   └── spatialxenium/         # Xenium-specific spatial analysis
│   └── qc-report/             # QC reporting for SpatialXenium
├── containers/                # Built Singularity images (.sif files)
├── SpatialXenium/             # Nextflow pipeline (nf-core style)
├── scrnaseq/                  # R-based pipeline using whirl
├── tests/                     # Testing scripts and validation
└── docs/                      # Documentation and guides
```

---

## Quick Start (Test case executed locally)

### 1. Download Test Data
```bash
./setup_test_data.sh
```

### 2. Build Containers
```bash
cd docker
./build_all_images.sh -v docker

# Builds 4 containers: preprocessing, analysis, spatialxenium, qc-report
# Total time: ~3-6 hours (parallelizable)
```

### 3. Run sc/snRNA-seq Pipeline
```bash
cd scrnaseq
./run_preprocessing_containerized.sh
./run_analysis_containerized.sh
```

### 4. Run SpatialXenium Pipeline
```bash
cd SpatialXenium
nextflow run main.nf -profile test,docker --outdir test_results
```

---

## Pipelines

### scrnaseq (R-based)

**Purpose:** Flexible R workflow for single-cell analysis  
**Framework:** whirl package  
**Workflows:**
- Preprocessing: QC, doublet detection, batch correction
- Analysis: Cell annotation, differential expression

**Documentation:** [`scrnaseq/README.md`](scrnaseq/README.md)

### SpatialXenium (Nextflow)

**Purpose:** Automated spatial transcriptomics analysis  
**Platform:** 10x Genomics Xenium  
**Features:**
- Quality control and filtering
- Spatial feature analysis
- Optional reference-based cell type annotation
- Comprehensive QC reporting

**Documentation:** [`SpatialXenium/README.md`](SpatialXenium/README.md)

---

## Container Images

All R-based containers use **R 4.5.2** with **Bioconductor 3.22**

| Container | Purpose | Key Packages | Size |
|-----------|---------|--------------|------|
| **preprocessing** | QC & batch correction | Seurat, SoupX, scDblFinder, harmony | ~8.09 GB |
| **analysis** | Annotation & DE | Azimuth, SingleR, edgeR, Signac | ~10.97 GB |
| **spatialxenium** | Spatial analysis | presto, Seurat, argparser, arrow, spatstat | ~7.95 GB |
| **qc-report** | QC PDF generation | Python 3.11, fpdf2 | ~0.3 GB |

**Build documentation:** [`docker/BUILD_GUIDE.md`](docker/BUILD_GUIDE.md)

---

## Documentation

| Document | Description |
|----------|-------------|
| [`docs/QUICKSTART.md`](docs/QUICKSTART.md) | One-page quick reference |
| [`docs/TESTING_GUIDE.md`](docs/TESTING_GUIDE.md) | Complete testing procedures |
| [`docker/BUILD_GUIDE.md`](docker/BUILD_GUIDE.md) | Container build guide |

---

## Requirements

### System
- Linux/Windows/Mac (tested on MacOS Tahoe 26.1 and SLURM HPC)
- 18+ GB RAM (32+ GB recommended)
- 40+ GB disk space for containers and test data

### Software
- Docker or Apptainer/Singularity ≥3.8
- Nextflow ≥23.10.0
- Standard Linux tools (wget/curl, unzip, tar)

---

## Data Sources

### Test Data (Auto-downloaded)
- **scRNA-seq:** 4 single-cell h5 matrixes (10x Genomics)
- **Non-diseased kidney:** Xenium V1 human kidney (10x Genomics)
- **Cancer kidney:** Xenium V1 kidney cancer (10x Genomics)
- **Total size:** ~6 GB

### Optional Reference Data
- Kidney reference: [Azimuth HubMAP](https://azimuth.hubmapconsortium.org/references/#Human%20-%20Kidney)
- Additional datasets: [10x Genomics Dataset Explorer](https://www.10xgenomics.com/products/xenium-in-situ/dataset-explorer)

---

## Development

### Contributing

This is a technical research repository. For improvements:
1. Test changes thoroughly
2. Update relevant documentation
3. Maintain portability (no external paths)
4. Follow existing code style

---

## Troubleshooting

### Common Issues

**Container not found**
```bash
ls containers/*.sif
cd docker && bash build_all_images.sh
```

**Test data incomplete**
```bash
./setup_test_data.sh
find test_data/ -name "*.tar.gz"  # Should be empty
```

**Pipeline fails**
```bash
tail SpatialXenium/.nextflow.log
nextflow run main.nf -profile test,singularity -resume
```

**Full troubleshooting:** See [`TESTING_GUIDE.md`](docker/TESTING_GUIDE.md#common-issues-and-solutions)

---

## Citation

If you use this pipeline in your research, please cite:

```
[Citation information to be added upon publication]
```
Pre-print available at [bioRxiv](https://www.biorxiv.org/content/10.64898/2026.01.23.701261v1)

---

## License

- Code - Apache License 2.0
- Documentation - CC BY 4.0

---

## Contact

**Maintainer:** Martin E. Garcia Sola, Amaya Zaratiegui, Timothy Burfield  
**Email:**   martin.garciasola@zs.com, azeu@novonordisk.com, tmbf@novonordisk.com
**Institution:** ZS Associates & Novo Nordisk A/S

---

## Acknowledgments

- Built using [nf-core](https://nf-co.re/) framework
- Test data from [10x Genomics](https://www.10xgenomics.com/)
- Container infrastructure based on [Rocker Project](https://rocker-project.org/)

---

**Version:** 1.0  
**Last Updated:** February 4, 2026
