<p align="center">
  <img src="https://raw.githubusercontent.com/cylin2022/triad/main/static/triad_logo.png" alt="TRIAD Logo" width="560" />
</p>

# TRIAD: Multimodal T-Cell Epitope Prioritization & Immunology AI Agent Platform

[![Docker Pulls](https://img.shields.io/docker/pulls/lsbnb/triad.svg)](https://hub.docker.io/r/lsbnb/triad)
[![Docker Image](https://img.shields.io/badge/docker-lsbnb%2Ftriad-blue.svg)](https://hub.docker.io/r/lsbnb/triad)
[![Platform Version](https://img.shields.io/badge/version-v1.1-green.svg)](https://hub.docker.io/r/lsbnb/triad)
[![Port Requirement](<https://img.shields.io/badge/port%20mapping-5001%3A5001-important.svg>)](https://hub.docker.io/r/lsbnb/triad)

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
- **Population-Aware Coverage**: Computes non-redundant population coverage across global and regional cohorts (Taiwan, East Asia, Europe, World).
- **Multimodal Candidate Tiering**: Automatically categorizes candidates into **Tier 1** (High priority), **Tier 2** (Secondary), and **Tier 3** based on a harmonized 4-dimensional Composite Priority Index.
- **In-Memory RAM Disk Acceleration Engine**: Parallel multi-core batch processing using `/dev/shm` for zero-disk-latency temporary sequence ingestion.

---

## 🚀 Quick Start Guide

### 1. Pull Image from Docker Hub

```bash
docker pull lsbnb/triad:latest
# or
docker pull lsbnb/mhcpan_shell:latest
```

---

### 2. Launch Container with Port Mapping `5001:5001`

#### Option A: Command Line Interface (CLI)

Run the container mapping host port `5001` to container port `5001`:

```bash
docker run -d \
  --name triad_app \
  -p 5001:5001 \
  --restart unless-stopped \
  lsbnb/triad:latest
```

Open your browser and navigate to:
👉 **`http://localhost:5001`** (or `http://127.0.0.1:5001`)

---

#### Option B: Docker Desktop Graphical Interface (Windows / macOS GUI Users)

If you use **Docker Desktop** GUI:

1. Click **Run** on the `lsbnb/triad:latest` image.
2. Expand **Optional settings** &rarr; **Ports**.
3. Fill in **`5001`** in the **Host port** field (which maps to `:5001/tcp`).
4. Click **Run**.

> ⚠️ **IMPORTANT:** Leaving **Host port** empty causes Docker Desktop to assign a random host port, making `http://localhost:5001` inaccessible.

![Docker Desktop Host Port 5001 Setup](https://raw.githubusercontent.com/cylin2022/triad/main/static/docker_desktop_port_setting.png)

---

## 🔑 NetMHCpan-4.2 Academic Licensing & Setup Guide

`netMHCpan-4.2` is academic-licensed software copyrighted by **DTU Health Tech (Technical University of Denmark)**. Its academic license explicitly prohibits third-party redistribution of binary executables and model parameters (`"not give the program to third parties"`).

For legal compliance, the public Docker image (`lsbnb/triad`) ships as a clean application shell **without** pre-bundled DTU files. Each user must obtain their own copy directly from DTU.

### Step 1: Request Linux Package from DTU

1. Visit the [Official DTU Download Portal](https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/).
2. Request the Linux tarball package (filename format: `netMHCpan-4.2.Linux.tar.gz`).

### Step 2: Configure Container (Choose Method 1 or Method 2)

#### Method 1: Host Directory Volume Mount (Server / CLI)

Extract the package on your host machine and mount it to `/opt/netMHCpan-4.2`:

```bash
# 1. Unpack DTU package on host
tar -xzvf netMHCpan-4.2.Linux.tar.gz

# 2. Start container with volume mount and port mapping 5001:5001
docker run -d \
  --name triad_app \
  -p 5001:5001 \
  -v /dev/shm:/dev/shm \
  -v "$(pwd)/netMHCpan-4.2:/opt/netMHCpan-4.2" \
  lsbnb/triad:latest
```

#### Method 2: Web Upload Wizard (No Terminal Access Needed)

1. Start the container with port mapping `5001:5001`:

   ```bash
   docker run -d --name triad_app -p 5001:5001 lsbnb/triad:latest
   ```

2. Open **`http://localhost:5001/setup`** in your web browser.
3. Upload your `netMHCpan-4.2.Linux.tar.gz` file.
4. Confirm license compliance and click **Upload & Install**. The platform automatically extracts, verifies, tests, and activates the engine.

---

## 🐳 Docker Compose Deployment

Create a `docker-compose.yml` file:

```yaml
version: '3.8'

services:
  triad:
    image: lsbnb/triad:latest
    container_name: triad_app
    ports:
      - "5001:5001"
    volumes:
      - /dev/shm:/dev/shm
      - ./netMHCpan-4.2:/opt/netMHCpan-4.2
    restart: unless-stopped
```

Run with:

```bash
docker compose up -d
```

---

## 📊 Output Data Format & Column Ordering

Output tables and exported CSV / Excel reports strictly follow standard NetMHCpan ordering with **Protein ID / Identifier** positioned in the first column:

| Col # | Field                    | Label                             | Description                                                                                           |
| :---- | :----------------------- | :-------------------------------- | :---------------------------------------------------------------------------------------------------- |
| 1     | `Identity`             | **Protein ID / Identifier** | FASTA header sequence ID or source protein identifier (e.g.`sp\|P0DTC2\|SPIKE_SARS2` or `PEPLIST`). |
| 2     | `Pos`                  | Position                          | Amino acid starting position in the source protein.                                                   |
| 3     | `MHC`                  | HLA Allele                        | Targeted HLA allele (e.g.,`HLA-A*02:01`).                                                           |
| 4     | `Peptide`              | Peptide Sequence                  | Predicted k-mer peptide amino acid sequence.                                                          |
| 5     | `Core`                 | Core Motif                        | Binding core motif predicted by NetMHCpan.                                                            |
| 6     | `Score_EL`             | Presentation Score                | Raw eluted ligand presentation probability (`0.0000 ~ 1.0000`).                                     |
| 7     | `Rank_EL`              | %Rank EL                          | Presentation percentile rank ($\le 0.5\%$: Strong Binder, $\le 2.0\%$: Weak Binder).              |
| 8     | `Affinity_nM`          | Binding IC50 (nM)                 | Quantitative IC50 binding affinity in nanomolar ($<50\text{ nM}$: High affinity).                   |
| 9     | `Immunogenicity_Score` | T-Cell Immunogenicity             | Calis et al. physicochemical TCR contact activation score ($>0.0$: Active).                         |
| 10    | `Composite_Score`      | Composite Priority Index          | 4-Dimensional harmonized priority index (`0.0000 ~ 1.0000`).                                        |
| 11    | `Tier`                 | Decision Priority Tier            | **Tier 1** (High priority), **Tier 2** (Secondary), **Tier 3** (Low priority).      |
| 12    | `BindLevel`            | Binder Category                   | `SB` (Strong Binder), `WB` (Weak Binder), or empty.                                               |

---

## 🌐 REST API Usage Examples

### 1. Check System Specs & RAM Disk Status

```bash
curl -s http://localhost:5001/api/system_resources | jq .
```

### 2. Submit Prediction Payload

```bash
curl -s -X POST http://localhost:5001/api/predict \
  -H "Content-Type: application/json" \
  -d '{
    "mode": "peptide",
    "input": "AAAWYLWEV\nAEFGPWQTV\nYLLPAIVHI\nGILGFVFTL",
    "alleles": ["HLA-A*02:01", "HLA-B*07:02"],
    "include_ba": true
  }'
```

### 3. Download Filtered CSV Report

```bash
curl -O http://localhost:5001/api/download/<JOB_ID>/csv
```

---

## 🔬 Literature & References

1. **NetMHCpan Presentation Model**:Reynisson B, et al. *NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation.* **Nucleic Acids Res.** 2020;48(W1):W449-W454. [doi:10.1093/nar/gkaa379](https://doi.org/10.1093/nar/gkaa379)
2. **T-Cell Immunogenicity Model**:Calis JJ, et al. *Properties of MHC Class I Presented Peptides That Inspire Immunogenicity.* **PLoS Comput Biol.** 2013;9(10):e1003266. [doi:10.1371/journal.pcbi.1003266](https://doi.org/10.1371/journal.pcbi.1003266)
3. **Population Coverage Database**:
   Gonzalez-Galarza FF, et al. *Allele frequency net database (AFND) 2020 update.* **Nucleic Acids Res.** 2020;48(D1):D783-D788. [doi:10.1093/nar/gkz1029](https://doi.org/10.1093/nar/gkz1029)

---

## 🏛️ Maintained By

**Laboratory of Systems Biology and Bioinformatics (LSBNB)**
Institute of Information Science, Academia Sinica, Taipei, TAIWAN.
Web: [https://hub.docker.io/r/lsbnb/triad](https://hub.docker.io/r/lsbnb/triad)
