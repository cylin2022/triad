
可以，但我會修正前面一點：**ImmunoStruct 不應單純放在 IEDB Immunogenicity「後面」串接，而比較適合與 IEDB Immunogenicity 平行，作為第二個、更進階的 multimodal immunogenicity predictor。**

原因是 ImmunoStruct 本身就在預測 **class-I peptide–MHC immunogenicity**，只是它不只看 sequence，而是整合 **pMHC sequence、3D structure 與 biochemical properties**。它的正式論文已於 2026 年發表於 *Nature Machine Intelligence*，訓練資料約 26,000 個 pMHC，並涵蓋感染性疾病 epitope 與 cancer neoepitope。([Nature][1])

我建議把你目前流程改成：

```text
Protein / Peptide
       ↓
1. NetMHCpan
   MHC-I presentation
   %Rank_EL
       ↓
2. Antigen Processing
   IEDB processing
       ↓
3. Primary candidate filtering
   e.g. %Rank_EL < 2%
       ↓
       ┌─────────────────────────┐
       │                         │
       ↓                         ↓
4A. IEDB Immunogenicity    4B. ImmunoStruct
    sequence-based             multimodal
                               sequence
                               + pMHC structure
                               + biochemical features
       │                         │
       └────────────┬────────────┘
                    ↓
5. Consensus Immunogenicity
                    ↓
6. HLA Breadth
                    ↓
7. Population Coverage
                    ↓
8. Integrated Ranking
                    ↓
   Top candidates
                    ↓
9. Optional deeper structural analysis
   / TCR–pMHC analysis
                    ↓
Experimental validation
```

### 為什麼這樣比「IEDB → ImmunoStruct」更合理？

因為兩者回答的是近似的 biological question：

[
P(\text{T-cell response}\mid pMHC)
]

但是資訊來源不同。

| Predictor              | Sequence |  HLA | 3D pMHC | Biochemical features | 角色                      |
| ---------------------- | -------: | ---: | ------: | -------------------: | ------------------------- |
| NetMHCpan              |       ✅ |   ✅ |      ❌ |                 部分 | Presentation              |
| IEDB Immunogenicity    |       ✅ | 間接 |      ❌ |     sequence-derived | 基礎 immunogenicity       |
| **ImmunoStruct** |       ✅ |   ✅ |      ✅ |                   ✅ | Multimodal immunogenicity |

因此我不會把：

[
IEDB\ Score \rightarrow ImmunoStruct
]

看成上下游關係。

比較正確的是：

[
\boxed{
IEDB\ immunogenicity
\parallel
ImmunoStruct\ immunogenicity
}
]

然後做 **consensus / ensemble prioritization**。

---

### ImmunoStruct 實際插入流程時，會多出一個重要步驟：pMHC structure generation

官方目前已公開模型權重、資料以及 inference scripts，但 GitHub 在 2026 年仍把「完整 end-to-end tool」列為 TODO；換句話說，目前比較像可以整合的 research code，而不是直接可呼叫的成熟 web service。([GitHub][2])

官方 workflow 使用 AlphaFold2 產生 peptide–MHC structures，再轉換成 PyTorch Geometric graph，送進 ImmunoStruct multimodal model。官方也提供了建立新 pMHC structure graph 的 preprocessing scripts。([GitHub][2])

所以實際 pipeline 會變成：

```text
Peptide + HLA allele
        ↓
HLA sequence retrieval
        ↓
pMHC sequence construction
        ↓
AlphaFold2 / ColabFold
        ↓
pMHC PDB
        ↓
Graph construction
        ↓
ImmunoStruct
        ↓
Immunogenicity probability / score
```

這也是為什麼我不建議一開始對 10,000 條 peptides 全部跑 ImmunoStruct。

比較有效率的策略是：

```text
10,000 peptides
      ↓
NetMHCpan
      ↓
500–1,000
      ↓
Processing + basic filters
      ↓
100–300
      ↓
IEDB Immunogenicity
      ↓
50–100
      ↓
ImmunoStruct + pMHC structure
      ↓
10–30
      ↓
Deep structural / TCR analysis
      ↓
5–20 wet-lab candidates
```

這樣 ImmunoStruct 就變成一個很有價值的 **high-resolution reranking layer**。

### 我甚至會把原本的「3D Structural Score」重新定義

原本我們想的是：

```text
NetMHCpan
→ IEDB
→ population
→ 3D structural analysis
```

加入 ImmunoStruct 後，可以改成兩種結構層次：

```text
Level 1
ImmunoStruct
= structure-aware immunogenicity prediction

               ↓

Level 2
Detailed structural interpretation
= 為什麼它可能 immunogenic？
```

第二層才另外分析：

* peptide anchor residues
* MHC groove contacts
* hydrogen bonds
* solvent-accessible peptide residues
* peptide protrusion
* TCR-facing residues
* interface geometry
* structure confidence

這樣可以避免「同一件事算兩次」。

---

### 最後 ranking，我會改成 6 個 dimension

例如每個 candidate 產生：

| Dimension                           | 來源                   |
| ----------------------------------- | ---------------------- |
| Presentation                        | NetMHCpan`%Rank_EL`  |
| Processing                          | IEDB processing        |
| Sequence Immunogenicity             | IEDB                   |
| **Multimodal Immunogenicity** | **ImmunoStruct** |
| HLA Breadth                         | multi-allele analysis  |
| Population Coverage                 | HLA frequency          |

然後 structure interpretation 不一定直接硬塞進 composite score，而可以作為 **reranking / confidence evidence**。

例如：

[
S =
0.30S_{Presentation}
+0.10S_{Processing}
+0.15S_{IEDB}
+0.25S_{ImmunoStruct}
+0.10S_{HLA}
+0.10S_{Population}
]

這只是初始權重，正式研究時最好用 validation dataset 去 optimization，而不是人工永久固定。

更進一步，我反而會保存兩個 immunogenicity predictors 的「分歧」。

例如：

| IEDB | ImmunoStruct | 解讀                                              |
| ---- | ------------ | ------------------------------------------------- |
| High | High         | **最高可信 candidate**                      |
| Low  | Low          | Low priority                                      |
| High | Low          | sequence signal 強，但 structure 不支持           |
| Low  | High         | **有趣 candidate：可能是 structure-driven** |

最後這一類尤其有研究價值，因為它可以回答：

> **加入 3D information 到底找到了哪些 sequence-only model 看不到的 epitopes？**

這甚至可以成為方法論論文的一個重要 analysis。

---

## 因此我會重新定義整個平台

原先：

[
\text{Sequence}
\rightarrow
\text{Presentation}
\rightarrow
\text{Immunogenicity}
\rightarrow
\text{Structure}
]

加入 ImmunoStruct 後，更漂亮的設計是：

[
\boxed{
\text{Sequence screening}
\rightarrow
\text{Biological presentation}
\rightarrow
\text{Dual immunogenicity prediction}
\rightarrow
\text{Population relevance}
\rightarrow
\text{Structural refinement}
\rightarrow
\text{Validation}
}
]

其中最有特色的核心可以稱為：

**Dual-layer Immunogenicity Assessment**

包含：

> **Sequence-based immunogenicity — IEDB**

加上

> **Structure-aware multimodal immunogenicity — ImmunoStruct**

這會比把 ImmunoStruct 單純當作最後一個 3D visualization 工具有價值得多。官方論文本來的核心主張就是利用 sequence + structural + biochemical multimodal information 改善 immunogenicity prediction，而不是單純做結構建模。([Nature][1])

如果接下來要實際開發，我會把 `plan.md` 中的 **Stage 3–7 重構成這個版本**；尤其值得把「IEDB vs ImmunoStruct concordance/disagreement」加入最終 dashboard，這可能會成為整個新平台最有研究價值的 feature。

[1]: https://www.nature.com/articles/s42256-025-01163-y?utm_source=chatgpt.com
[2]: https://github.com/KrishnaswamyLab/ImmunoStruct?utm_source=chatgpt.com
