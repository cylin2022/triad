# Multimodal T-cell Epitope Prioritization Platform Plan

## 1. Project Goal

建立一套整合式 T-cell epitope prioritization 平台，從 protein/peptide sequence 出發，依序整合：

1. MHC-I binding / presentation
2. Antigen processing
3. Immunogenicity
4. HLA breadth
5. Population coverage
6. 3D pMHC structural analysis
7. Optional TCR–pMHC analysis
8. Integrated candidate ranking and visualization

目標不是只預測「哪一條 peptide 會 binding MHC」，而是回答：

> 哪些 peptide 最可能被有效處理、呈現、被 T cell 辨識，並且在目標族群中具有較高的生物學與實驗驗證價值？

---

## 2. Biological Rationale

T-cell epitope 與 MHC 密切相關，因為 T cell receptor（TCR）實際辨識的是：

**peptide–MHC complex（pMHC）**

而不是單獨的 peptide。

因此一條 peptide 要成為有效的 T-cell epitope，通常需經過下列流程：

```text
Protein antigen
    ↓
Peptide generation
    ↓
MHC binding
    ↓
Antigen processing / presentation
    ↓
pMHC displayed on cell surface
    ↓
TCR recognition
    ↓
T-cell response
```

重要觀念：

```text
MHC binding
    ≠
MHC presentation
    ≠
T-cell immunogenicity
```

因此不應只用 MHC binding score 判斷最終候選。

---

## 3. Proposed Analysis Pipeline

```text
Protein / Peptide Input
        ↓
[Stage 1]
NetMHCpan
MHC-I Binding / Eluted Ligand Prediction
        ↓
[Stage 2]
Antigen Processing
Proteasome + TAP + MHC Presentation
        ↓
[Stage 3]
IEDB Class I Immunogenicity
        ↓
[Stage 4]
HLA Breadth Analysis
        ↓
[Stage 5]
Population Coverage
        ↓
[Stage 6]
Integrated Candidate Ranking
        ↓
Top 20–50 candidates
        ↓
[Stage 7]
3D pMHC Structural Analysis
        ↓
Top 5–20 candidates
        ↓
[Stage 8, Optional]
TCR–pMHC Structural / Specificity Analysis
        ↓
Experimental Validation
```

---

## 4. Stage 1 — MHC-I Binding and Presentation

### Recommended tool
- NetMHCpan

### Main outputs
- `%Rank_EL`
- `%Rank_BA`
- predicted affinity
- binding level

### Preferred metric
優先使用 `%Rank_EL`，因為 EL 模型較接近實際被 MHC 呈現的 ligand，而不只是純 binding affinity。

### Suggested initial thresholds

```text
Strong binder:   %Rank_EL < 0.5%
Weak binder:     0.5–2%
Non-binder:      >2%
```

建議：

- candidates 很多時：先保留 `%Rank_EL < 0.5%`
- candidates 較少時：可放寬到 `%Rank_EL < 2%`

---

## 5. Stage 2 — Antigen Processing

### Purpose

評估 peptide 是否可能經歷：

1. proteasomal cleavage
2. TAP transport
3. MHC loading / presentation

### Options

#### Option A
直接使用 NetMHCpan EL 作為 practical approximation。

#### Option B
額外加入 IEDB MHC-I Processing prediction。

### Suggested outputs

- proteasome score
- TAP score
- MHC score
- processing score
- total score

這一層主要回答：

> 此 peptide 是否不只「能 bind」，而且有機會實際被產生並呈現在細胞表面？

---

## 6. Stage 3 — Immunogenicity Prediction

### Recommended tool
IEDB T Cell Prediction – Class I 中的：

**Immunogenicity Prediction**

### Interpretation

此步驟是在已形成 pMHC 的前提下，評估其被 T cell 辨識並引發 response 的可能性。

簡化理解：

```text
NetMHCpan
→ 能不能被呈現？

IEDB Immunogenicity
→ 呈現之後，是否較可能被 T cell 辨識？
```

### Important note

immunogenicity score 不應單獨做 hard cutoff。

較合理的使用方式：

- 先以 presentation / binding 篩選
- 再以 immunogenicity 作 secondary ranking

---

## 7. Stage 4 — HLA Breadth

每條 peptide 不應只看單一 HLA allele。

建議計算：

- 可結合的 HLA-A 數目
- 可結合的 HLA-B 數目
- 可結合的 HLA-C 數目
- strong-binding allele count
- total supported HLA alleles

### Example

```text
Peptide X

HLA-A*02:01   %Rank_EL = 0.14
HLA-A*24:02   %Rank_EL = 0.41
HLA-A*11:01   %Rank_EL = 0.62
```

此類 peptide 通常比只對單一 rare HLA 有極強 binding 的 peptide 更值得優先驗證。

---

## 8. Stage 5 — Population Coverage

### Purpose

評估候選 peptide 對目標族群的適用性。

可計算：

- global population coverage
- Asian / East Asian coverage
- Taiwan-relevant HLA coverage
- disease-specific population if available

### Suggested interpretation

高優先 candidate 應同時具有：

```text
Strong presentation
+
Positive immunogenicity
+
Broad HLA support
+
High population coverage
```

---

## 9. Stage 6 — Integrated Candidate Ranking

不建議直接把所有 raw score 相加，因為不同工具的 score scale 不一致。

### Recommended approach

先將各項指標 normalize 到 0–1，再建立 composite score。

例如：

```text
Candidate Score =
0.40 × Presentation
+
0.30 × Immunogenicity
+
0.15 × HLA Breadth
+
0.15 × Population Coverage
```

其中：

- Presentation 權重最高
- Immunogenicity 次之
- HLA breadth 與 population coverage 用於增加 translational relevance

### Alternative

初期可以先不建立單一 score，而使用 rule-based tier system。

---

## 10. Candidate Tier System

### Tier 1 — High priority

條件建議：

```text
%Rank_EL < 0.5%
+
Positive / high immunogenicity
+
Common HLA support
+
Good population coverage
```

用途：

- 優先進入 3D analysis
- 優先 wet-lab validation

### Tier 2 — Medium priority

```text
%Rank_EL = 0.5–2%
+
Positive immunogenicity
```

用途：

- secondary candidates

### Tier 3 — Backup

可能包含：

- strong binding，但 immunogenicity 較低
- HLA coverage 較窄
- population frequency 較低

---

## 11. Stage 7 — 3D pMHC Structural Analysis

3D structure 建議放在後期，而不是一開始就對所有 candidates 執行。

### Input

- peptide
- HLA allele

### Goal

分析 peptide 在 MHC groove 中的結構合理性與 TCR-facing surface。

### Suggested structural features

- peptide–MHC interface quality
- anchor residue placement
- hydrogen bonds
- hydrophobic contacts
- peptide protrusion
- solvent exposure
- TCR-facing residue exposure
- structural confidence
- predicted complex stability

### Why this matters

T-cell epitope 雖然通常由 linear peptide 定義，但 TCR 實際看到的是 peptide 在 MHC groove 中形成的 3D surface。

因此：

```text
Sequence-compatible
≠
Structurally equivalent
```

3D analysis 可以作為 Top candidates 的 reranking layer。

---

## 12. Stage 8 — Optional TCR–pMHC Analysis

如果有已知 TCR sequence 或 TCR repertoire，可進一步分析：

```text
TCR
+
Peptide
+
MHC
```

此層才真正接近：

**T-cell specificity prediction**

### Limitation

如果只有 peptide + HLA，而沒有 TCR sequence，就不能進行真正個人化的 TCR recognition prediction。

因此建議將此功能設為 optional advanced module。

---

## 13. Recommended Final Candidate Profile

每個 candidate 建議產生一個多維度 profile：

```text
Candidate P17

MHC presentation       0.91
Immunogenicity         0.82
HLA breadth            0.74
Population coverage    0.81
Structural confidence  0.89

Overall priority       Tier 1
```

這比只提供單一 binding score 更容易解讀。

---

## 14. Recommended Result Table

| Peptide | HLA | %Rank_EL | %Rank_BA | Processing | Immunogenicity | HLA Breadth | Population Coverage | Structural Score | Tier |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| Pep-01 | HLA-A*02:01 | 0.12 | 0.35 | 1.80 | 0.42 | High | 0.81 | 0.89 | Tier 1 |
| Pep-02 | HLA-A*11:01 | 0.31 | 0.70 | 1.40 | 0.58 | Medium | 0.73 | 0.82 | Tier 1 |
| Pep-03 | HLA-B*07:02 | 0.95 | 1.60 | 1.10 | 0.71 | Low | 0.42 | 0.77 | Tier 2 |

---

## 15. Recommended Visualizations

### Figure 1 — Integrated Pipeline

```text
Sequence
   ↓
Presentation
   ↓
Processing
   ↓
Immunogenicity
   ↓
Population
   ↓
Structure
   ↓
Candidate
```

### Figure 2 — Candidate Priority Map

X-axis:

```text
−log10(%Rank_EL)
```

Y-axis:

```text
IEDB Immunogenicity Score
```

Additional encodings:

- bubble size = population coverage
- label = peptide ID
- shape = HLA-A / HLA-B / HLA-C
- optional color = candidate tier

右上角代表：

```text
Strong presentation
+
High immunogenicity
```

通常為最高優先 candidate。

---

## 16. Platform Positioning

目前 IEDB 已可涵蓋：

- MHC-I binding
- processing
- immunogenicity
- population coverage

NetMHCpan 則提供很成熟的 MHC-I binding / ligand presentation prediction。

然而目前仍缺乏一個完整平台，將以下資訊自動整合：

```text
NetMHCpan
+
IEDB Processing
+
IEDB Immunogenicity
+
HLA Breadth
+
Population Coverage
+
3D pMHC Structure
+
Automated Ranking
+
Interactive Visualization
```

因此平台 novelty 應避免定位成單純的 tool wrapper，而應強調：

**Multimodal T-cell Epitope Prioritization**

---

## 17. Proposed Platform Architecture

```text
                    ┌─────────────────────┐
                    │ Protein / Peptide   │
                    └─────────┬───────────┘
                              ↓
                    ┌─────────────────────┐
                    │ NetMHCpan           │
                    │ MHC Presentation    │
                    └─────────┬───────────┘
                              ↓
                    ┌─────────────────────┐
                    │ Processing Module   │
                    │ Proteasome / TAP    │
                    └─────────┬───────────┘
                              ↓
                    ┌─────────────────────┐
                    │ Immunogenicity      │
                    │ IEDB Class I        │
                    └─────────┬───────────┘
                              ↓
              ┌───────────────┴───────────────┐
              ↓                               ↓
     ┌──────────────────┐           ┌──────────────────┐
     │ HLA Breadth      │           │ Population       │
     │ Analysis         │           │ Coverage         │
     └─────────┬────────┘           └─────────┬────────┘
               └──────────────┬───────────────┘
                              ↓
                    ┌─────────────────────┐
                    │ Integrated Ranking  │
                    │ Tier 1 / 2 / 3      │
                    └─────────┬───────────┘
                              ↓
                    ┌─────────────────────┐
                    │ 3D pMHC Analysis    │
                    └─────────┬───────────┘
                              ↓
                    ┌─────────────────────┐
                    │ Final Candidates    │
                    └─────────┬───────────┘
                              ↓
                    ┌─────────────────────┐
                    │ Experimental        │
                    │ Validation          │
                    └─────────────────────┘
```

---

## 18. Experimental Validation

最終 Top 5–20 candidates 可進入實驗驗證：

```text
Top candidates
      ↓
pMHC binding / stability
      ↓
PBMC stimulation
      ↓
IFN-γ ELISpot
or
Intracellular Cytokine Staining
      ↓
Tetramer validation
      ↓
Validated T-cell epitope
```

最終 workflow：

```text
Prediction
    ↓
Prioritization
    ↓
Structural Reranking
    ↓
Experimental Validation
```

---

## 19. Suggested Development Phases

### Phase I — Core Sequence-based Platform

整合：

- NetMHCpan
- IEDB processing
- IEDB immunogenicity
- HLA breadth
- population coverage

完成：

- candidate table
- candidate ranking
- Tier 1/2/3 classification
- basic visualization

### Phase II — Structural Module

加入：

- pMHC structure prediction
- interface analysis
- TCR-facing residue analysis
- structural reranking

### Phase III — TCR-aware Module

若有 TCR sequence：

- TCR–pMHC modeling
- TCR specificity analysis
- repertoire-aware candidate prioritization

---

## 20. Proposed Research Positioning

### Conventional approach

```text
Epitope Prediction
≈
MHC Binding Prediction
```

### Proposed approach

```text
Multimodal T-cell Epitope Prioritization
=
Presentation
+
Processing
+
Immunogenicity
+
HLA Breadth
+
Population Coverage
+
3D Structure
+
Optional TCR Recognition
```

核心概念：

> 從「哪條 peptide 可以 bind MHC」提升為「哪條 peptide 最值得進入 T-cell experimental validation」。

---

## 21. Potential Project Title

**Multimodal T-cell Epitope Prioritization through Integrated MHC Presentation, Immunogenicity, Population Coverage, and Structural Analysis**

較簡短版本：

**An Integrated Platform for Multimodal T-cell Epitope Prioritization**

或：

**Structure-Aware Multimodal Prioritization of T-cell Epitopes**

---

## 22. Key Novelty

平台的核心創新不應只是整合既有 predictors，而是：

1. 將多個免疫學層次整合為單一 candidate profile
2. 建立可解釋的 prioritization framework
3. 同時考慮 HLA breadth 與 population coverage
4. 將 3D pMHC structure 作為後期 reranking layer
5. 可進一步加入 TCR-aware analysis
6. 最終直接輸出適合 experimental validation 的候選清單

最終核心定位：

**Prediction → Integration → Prioritization → Structural Refinement → Validation**
