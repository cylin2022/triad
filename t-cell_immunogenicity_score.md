# TRIAD: T-Cell Immunogenicity Score Reference

Comprehensive documentation of the **T-Cell Immunogenicity Scoring Model** implemented in the TRIAD platform, covering theoretical foundations, mathematical formulations, parameter tables, length adaptation, and clinical candidate ranking.

---

## 1. Scientific Origin & Core Concept

* **Reference**:  
  Calis JJ, Maybeno M, Greenbaum JA, Weiskopf D, De Silva AD, Sette A, Keşmir C, Peters B.  
  *Properties of MHC Class I Presented Peptides That Inspire Immunogenicity.*  
  **PLoS Computational Biology** (2013), 9(10): e1003266. [doi:10.1371/journal.pcbi.1003266](https://doi.org/10.1371/journal.pcbi.1003266)

* **Problem Addressed**:  
  Conventional epitope prediction pipelines (such as NetMHCpan) predict **MHC-I binding affinity and presentation** ($\text{IC}_{50}$ or $\% \text{Rank\_EL}$). However, **over 50% of peptides that bind strongly to MHC fail to trigger a CD8+ T-cell response in vivo**.

* **Underlying Mechanism**:  
  MHC-presented peptides must interact productively with the **T-cell receptor (TCR)**. TCR activation is primarily dictated by the **physicochemical properties of amino acids exposed at TCR-facing contact positions** (the solvent-exposed central bulge of the peptide), rather than the anchor residues buried inside the MHC cleft.

---

## 2. Mathematical Formulation

The raw immunogenicity score is computed as a position-weighted sum of amino acid log-odds enrichment values, excluding anchor positions buried within MHC pockets:

$$\text{Immunogenicity\_Score} = \sum_{i \notin \text{Mask}} w_i \times s(\text{aa}_i)$$

Where:
* $i$: 0-indexed position along the peptide sequence ($0 \dots L-1$).
* $\text{Mask}$: Set of anchor positions buried inside the MHC binding cleft (interacting with MHC, not TCR).
* $w_i$: Positional contribution weight reflecting how prominently position $i$ protrudes toward the TCR CDR3 loops.
* $s(\text{aa}_i)$: Physicochemical log-odds enrichment score for amino acid $\text{aa}_i$.

---

## 3. Amino Acid Propensity Scale ($s(\text{aa})$)

Calis et al. derived these values by comparing thousands of experimentally verified immunogenic vs. non-immunogenic pMHC complexes from the Immune Epitope Database (IEDB):

| Category | Amino Acid | 1-Letter Code | Log-Odds Weight ($s$) | Structural & Biophysical Rationale |
| :--- | :--- | :---: | :---: | :--- |
| **Highly Immunogenic** | **Tryptophan**<br>**Isoleucine**<br>**Phenylalanine**<br>**Glutamate** | **W**<br>**I**<br>**F**<br>**E** | **+0.719**<br>**+0.432**<br>**+0.380**<br>**+0.325** | **Bulky, aromatic, and rigid aliphatic side-chains** protrude deeply into the TCR CDR3 binding pocket, establishing strong van der Waals contacts and hydrophobic interactions. Negatively charged glutamate (E) readily forms salt bridges. |
| **Moderately Positive** | **Arginine**<br>**Valine**<br>**Alanine**<br>**Threonine**<br>**Glycine**<br>**Histidine**<br>**Aspartate** | **R**<br>**V**<br>**A**<br>**T**<br>**G**<br>**H**<br>**D** | +0.168<br>+0.134<br>+0.127<br>+0.126<br>+0.110<br>+0.105<br>+0.072 | Favorable to neutral contributions without significant steric repulsion. |
| **Weak / Neutral** | **Tyrosine**<br>**Asparagine**<br>**Leucine**<br>**Proline** | **Y**<br>**N**<br>**L**<br>**P** | -0.012<br>-0.021<br>-0.036<br>-0.036 | Minimal net contribution to T-cell receptor engagement. |
| **Unfavorable / Disfavored** | **Cysteine**<br>**Glutamine**<br>**Serine**<br>**Methionine**<br>**Lysine** | **C**<br>**Q**<br>**S**<br>**M**<br>**K** | -0.175<br>-0.376<br>-0.537<br>-0.570<br>**-0.700** | **Small, overly flexible, or highly hydrated hydrophilic residues**. They fail to provide sufficient desolvation free energy or steric engagement needed to trigger TCR conformational activation. |

---

## 4. Positional Weighting ($w_i$) & Length Adaptation

### Standard 9-mers
For standard 9-mer peptides, positional weights follow the canonical structural protrusion profile:

```python
# Defined in netmhcpan_platform/immunogenicity_scorer.py
POSITION_WEIGHTS = [0.00, 0.00, 0.10, 0.31, 0.30, 0.29, 0.26, 0.18, 0.00]
#                     P1    P2    P3    P4    P5    P6    P7    P8    P9
```

* **P1 & P2 ($w = 0.00$)**: Buried in the N-terminal A/B pockets of the MHC-I groove.
* **P4, P5, P6, P7 ($w = 0.26 \sim 0.31$)**: The apex of the peptide bulge facing directly into the CDR3$\alpha$ and CDR3$\beta$ loops.
* **P9 ($w = 0.00$)**: C-terminal anchor buried in the F pocket.

### Longer Peptides (> 9-mers: 10-mer, 11-mer, etc.)
Because the MHC Class I binding groove is closed at both ends, extra residues are forced into a **central zigzag arch**. TRIAD dynamically adapts by inserting a weight of `0.30` for each additional residue into the central bulge:

```python
if peplen > 9:
    weights = POSITION_WEIGHTS[:5] + [0.30] * (peplen - 9) + POSITION_WEIGHTS[5:]
```

---

## 5. Allele-Specific Anchor Masking (`ALLELE_POSITION_MASK`)

Different HLA alleles utilize distinct internal binding pockets (Pockets A through F). TRIAD dynamically masks non-accessible anchor positions depending on the presenting allele:

| Presenting Allele | Masked Positions (1-indexed) | Anchor Pocket Mechanism |
| :--- | :---: | :--- |
| **Default / Canonical** (e.g., HLA-A\*02:01, HLA-A\*03:01) | `1, 2, C-term` | Standard N-terminal P2 and C-terminal P9 anchors. |
| **HLA-B\*08:01** | `2, 5, 9` | Pocket D anchors position P5 internally; masked from TCR contact. |
| **HLA-A\*24:02 / HLA-A\*23:01** | `2, 7, 9` | Auxiliary anchor at P7. |
| **HLA-B\*27:05 / HLA-B\*44:02 / HLA-B\*44:03** | `2, 3, 9` | Pocket B/D anchors P3 as a major secondary anchor. |
| **HLA-A\*01:01** | `2, 3, 9` | P3 aspartate/glutamate anchor preference. |

---

## 6. Score Interpretation & Platform Integration

### 1. Raw `Immunogenicity_Score`
* Typical range: **$-0.4000$ to $+0.4000$**.
* **$> 0.00$**: **Positive immunogenic propensity** (enriched in bulky/aromatic residues at contact sites; CD8+ T-cell activation potential).
* **$< 0.00$**: **Low immunogenicity likelihood** (high risk of being a non-activating MHC binder).

### 2. Normalized $S_{\text{imm}}$ ($[0.0, 1.0]$)
For multi-criteria decision making, TRIAD rescales the raw score into a normalized interval $[0, 1]$:

$$S_{\text{imm}} = \text{clip}\left(\frac{\text{Immunogenicity\_Score} + 0.4}{0.8}, 0.0, 1.0\right)$$

### 3. Integration into the TRIAD `Composite_Score`
TRIAD synthesizes four complementary biological dimensions into a single golden index to rank epitope candidates:

$$\text{Composite\_Score} = \underbrace{0.40 \times S_{\text{pres}}}_{\text{MHC-I Presentation (%Rank)}} + \underbrace{0.25 \times S_{\text{aff}}}_{\text{Binding Affinity (IC}_{50}\text{)}} + \underbrace{\mathbf{0.25 \times S_{\text{imm}}}}_{\mathbf{TCR\ Immunogenicity}} + \underbrace{0.10 \times S_{\text{breadth}}}_{\text{HLA Promiscuity}}$$

---

## 7. Implementation Files in TRIAD
* Core Scorer: [`netmhcpan_platform/immunogenicity_scorer.py`](file:///home/cylin/NETMHC/netmhcpan_platform/immunogenicity_scorer.py)
* Batch Execution & Composite Ranking: [`netmhcpan_platform/batch_runner.py`](file:///home/cylin/NETMHC/netmhcpan_platform/batch_runner.py)
* Interactive Web UI Dashboard: [`netmhcpan_platform/templates/index.html`](file:///home/cylin/NETMHC/netmhcpan_platform/templates/index.html)
