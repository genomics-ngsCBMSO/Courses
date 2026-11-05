

# CONDA ENVIRONMENT INSTALLATION

**Servicio de Análisis Biocomputacional (SABio); CBM** \
**Edition**: October, 2026 \
**Last update**: 05/11/2026

This course uses Conda environments which are defined in YAML files.

## Prerequisites

Make sure that Conda is installed and available in your PATH:

```bash
conda --version
```

## Create the environment

To create the environment from the provided YAML file:

```bash
conda env create -f <environment.yml>
```

Replace _environment.yml_ with the path to your environment file.

## Verify the installation

Check that the environment was created successfully:

```bash
conda env list
```
or
```bash
conda info --envs
```

## Activate the environment

After installation, activate the environment:

```bash
conda activate <environment_name>
```

# Troubleshooting

#### Error: "_packages excluded by strict channel priority_"

Example:

```bash
LibMambaUnsatisfiableError: Encountered problems while solving:
  - package mafft-7.508-hec16e2b_0 is excluded by strict repo priority
  - package sortmerna-2.0-he860b03_4 is excluded by strict repo priority
```

This is commonly caused by Conda being configured with:

```bash
channel_priority: strict
```
which may prevent dependency resolution for environments that require packages from multiple channels.

**SOLUTION**

Set channel priority to flexible:

```bash
conda config --set channel_priority flexible
```

Alternatively, disable channel priority completely:

```bash
conda config --set channel_priority disabled
```

Then recreate the environment:

```bash
conda env create -f <environment.yml>
```
