# TRIAD: Multimodal T-Cell Epitope Prioritization & Immunology AI Agent Platform

[![Docker Pulls](https://img.shields.io/docker/pulls/lsbnb/triad.svg)](https://hub.docker.io/r/lsbnb/triad)
[![Docker Image](https://img.shields.io/badge/docker-lsbnb%2Ftriad-blue.svg)](https://hub.docker.io/r/lsbnb/triad)
[![Platform Version](https://img.shields.io/badge/version-v1.2-green.svg)](https://hub.docker.io/r/lsbnb/triad)
[![Port Requirement](https://img.shields.io/badge/port%20mapping-5001%3A5001-important.svg)](https://hub.docker.io/r/lsbnb/triad)

> ⚡ **CRITICAL DOCKER PORT MAPPING REQUIREMENT (`5001:5001`):**  
> To access TRIAD's Web Dashboard, you **MUST** configure port mapping **`5001:5001`**:  
> • **Command Line (CLI)**: Run `docker run -d -p 5001:5001 lsbnb/triad:latest`  
> • **Docker Desktop GUI**: Expand **Optional settings** → **Ports** and set **Host port: `5001`** (mapping to `:5001/tcp`).  
> *Without setting Host Port 5001, `http://localhost:5001` will fail to load!*

---

## 🧬 Overview

**TRIAD** (*Multimodal T-Cell Epitope Prioritization and Autonomous Immunology AI Agent Platform*) is an integrated, high-performance web platform designed for rapid CD8+ T-cell epitope discovery, antigen presentation prediction, physicochemical T-cell immunogenicity scoring, and global population coverage analysis.

Developed by the **Laboratory of Systems Biology and Bioinformatics (LSBNB)**, Institute of Information Science, Academia Sinica, Taiwan.

---

## ⚡ Key Features

- **High-Throughput Epitope Discovery**: Supports raw peptide lists or full-length protein FASTA sequences with configurable k-mer sliding windows (8–14 aa).
- **MHC-I Antigen Presentation ANN**: Integrates NetMHCpan-4.2 trained on over 1,000,000 mass spectrometry eluted ligands (MS-EL) and quantitative binding affinity (BA) datasets.
- **Physicochemical T-Cell Immunogenicity**: Implements Calis et al. (2013) log-odds scoring model based on position-weighted amino acid properties at TCR-contact positions (positions 4–8).
- **Physicochemical TCR Contact Feature Evaluation**: Evaluates T-cell receptor contact amino acid properties to predict immune activation potential and candidate concordance.
- **Population-Aware Coverage**: Computes non-redundant population coverage across **22 selectable populations** — Taiwan island-wide plus each of its Han and Indigenous groups, East Asia, Europe and World — using a locus-aware Hardy–Weinberg model.
- **Multimodal Candidate Tiering**: Automatically categorizes candidates into **Tier 1** (High priority), **Tier 2** (Secondary), and **Tier 3** based on a harmonized 4-dimensional Composite Priority Index.
- **In-Memory RAM Disk Acceleration Engine**: Parallel multi-core batch processing using `/dev/shm` for zero-disk-latency temporary sequence ingestion.

---

## 🆕 What's New in v1.2

### 🇹🇼 Taiwan Demographic Genetics & Locus-Aware Population Coverage

v1.2 upgrades the population genetics framework:

| Addition | What it does |
| :--- | :--- |
| **🇹🇼 Taiwan Populations** | Taiwan is now **15 selectable populations** instead of a silent fallback to East Asia: island-wide, Han, Minnan, Hakka, Indigenous, the 11 recognised peoples individually, and Pingpu. |
| **Locus-Aware Model** | Alleles at the same locus are summed before the Hardy–Weinberg square instead of multiplied as independent events, preventing under-estimation within loci. |
| **Lower Bound Reporting** | Alleles with no frequency record contribute 0 and are explicitly reported rather than receiving fabricated frequencies. |
| **AFND Data Quality Fixes** | Corrected thousand-separator parsing in AFND sample sizes, recovering 99.0% of sampled individuals (including DKMS and NMDP registries). |

### Frequency-Table Provenance

All 22 population tables live in `data/allele_frequencies.json`, regenerated and re-verified from one public input by `tools/build_afnd_tables.py`. Source: [slowkow/allelefrequencies](https://github.com/slowkow/allelefrequencies) mirror of AFND.
