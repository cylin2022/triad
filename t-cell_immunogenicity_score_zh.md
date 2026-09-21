# TRIAD: T 細胞免疫原性評分（T-Cell Immunogenicity Score）技術說明指南

本文件為 TRIAD 平台中採用的 **T 細胞免疫原性評分模型（T-Cell Immunogenicity Scoring Model）** 完整技術說明，涵蓋科學背景、數學公式、參數權重表、胜肽長度自適應演算法、HLA 錨定遮蔽以及在臨床候選標靶排序中的整合應用。

---

## 1. 科學背景與核心機制

* **經典文獻出處**：  
  Calis JJ, Maybeno M, Greenbaum JA, Weiskopf D, De Silva AD, Sette A, Keşmir C, Peters B.  
  *Properties of MHC Class I Presented Peptides That Inspire Immunogenicity.*  
  **PLoS Computational Biology** (2013), 9(10): e1003266. [doi:10.1371/journal.pcbi.1003266](https://doi.org/10.1371/journal.pcbi.1003266)

* **解決的臨床痛點**：  
  傳統表位預測工具（如 NetMHCpan）主要評估 **MHC-I 結合力與呈遞概率**（$\text{IC}_{50}$ 或 $\% \text{Rank\_EL}$）。然而實驗研究顯示：**超過 50% 能夠緊密結合並呈遞於 MHC 分子上的胜肽，在體內（in vivo）根本無法誘發 CD8+ T 細胞免疫反應（假陽性）**。

* **分子生物學機轉**：  
  胜肽呈遞至細胞表面後，必須與 **T 細胞受體（TCR）** 產生具備足夠親和力與構型誘發力的實質接觸。T 細胞的活化與否，高度取決於**暴露於溶劑中、朝向 TCR 的接觸殘基（主要為胜肽中央隆起處）之物理化學特性（Physicochemical Properties）**，而非深埋於 MHC 凹槽內部的錨定殘基。

---

## 2. 數學公式與演算法定義

原始免疫原性評分（`Immunogenicity_Score`）是透過對非錨定接觸位點之胺基酸對數優勢比（log-odds enrichment）進行位置加權求和：

$$\text{Immunogenicity\_Score} = \sum_{i \notin \text{Mask}} w_i \times s(\text{aa}_i)$$

* $i$：胜肽序列中的 0-indexed 位置（$0 \dots L-1$）。
* $\text{Mask}$：深埋於 MHC 凹槽口袋內部的錨定位點集合（這幾處與 MHC 相互作用，TCR 無法碰觸）。
* $w_i$：位置加權係數（Positional Contribution Weight），反映該位置凸出朝向 TCR CDR3 結合環的空間暴露度。
* $s(\text{aa}_i)$：胺基酸 $\text{aa}_i$ 的物理化學對數優勢免疫原性數值。

---

## 3. 胺基酸理化偏好量表（$s(\text{aa})$）

Calis 等人比對了免疫表位資料庫（IEDB）中數千筆經實驗驗證「能誘發免疫反應」與「無免疫活性」的 pMHC 數據，統計得出 20 種天然胺基酸的免疫活化偏好權重：

| 活化分級 | 胺基酸名稱 | 單字母代碼 | 對數權重 ($s$) | 結構生物物理機制說明 |
| :--- | :--- | :---: | :---: | :--- |
| **高度免疫原性<br>(Highly Immunogenic)** | **色胺酸 (Tryptophan)**<br>**異白胺酸 (Isoleucine)**<br>**苯丙胺酸 (Phenylalanine)**<br>**麩胺酸 (Glutamate)** | **W**<br>**I**<br>**F**<br>**E** | **+0.719**<br>**+0.432**<br>**+0.380**<br>**+0.325** | **大體積、剛性芳香環與強疏水性側鏈**能深植入 TCR 的 CDR3 凹槽，提供極強的凡得瓦爾力（van der Waals）與疏水作用力；帶負電的麩胺酸（E）則易與 TCR 形成穩定鹽橋，促成強烈的受體活化訊號。 |
| **中度正向貢獻<br>(Moderately Positive)** | **精胺酸 (Arginine)**<br>**纈胺酸 (Valine)**<br>**丙胺酸 (Alanine)**<br>**酥胺酸 (Threonine)**<br>**甘胺酸 (Glycine)**<br>**組胺酸 (Histidine)**<br>**天門冬胺酸 (Aspartate)** | **R**<br>**V**<br>**A**<br>**T**<br>**G**<br>**H**<br>**D** | +0.168<br>+0.134<br>+0.127<br>+0.126<br>+0.110<br>+0.105<br>+0.072 | 提供中性至溫和的正向結合貢獻，不會產生嚴重的空間位阻碰撞。 |
| **微弱 / 中性貢獻<br>(Weak / Neutral)** | **酪胺酸 (Tyrosine)**<br>**天門冬醯胺 (Asparagine)**<br>**白胺酸 (Leucine)**<br>**脯胺酸 (Proline)** | **Y**<br>**N**<br>**L**<br>**P** | -0.012<br>-0.021<br>-0.036<br>-0.036 | 對 TCR 辨識活化的淨貢獻極小或呈中性。 |
| **不利 / 抑制活化<br>(Unfavorable / Disfavored)** | **半胱胺酸 (Cysteine)**<br>**麩醯胺酸 (Glutamine)**<br>**絲胺酸 (Serine)**<br>**甲硫胺酸 (Methionine)**<br>**離胺酸 (Lysine)** | **C**<br>**Q**<br>**S**<br>**M**<br>**K** | -0.175<br>-0.376<br>-0.537<br>-0.570<br>**-0.700** | **分子過小、彈性過高、或高度親水/水合化殘基**。此類側鏈周圍被水分子緊密包覆，無法提供足夠的去溶劑化自由能或立體抓握力來觸發 TCR 構型變化。 |

---

## 4. 空間位置加權值（$w_i$）與胜肽長度自適應

### 標準 9-mer 胜肽
9-mer 胜肽的位置加權精確模擬了胜肽在 MHC 溝槽內的幾何隆起輪廓：

```python
# 定義於 netmhcpan_platform/immunogenicity_scorer.py
POSITION_WEIGHTS = [0.00, 0.00, 0.10, 0.31, 0.30, 0.29, 0.26, 0.18, 0.00]
#                     P1    P2    P3    P4    P5    P6    P7    P8    P9
```

* **P1 與 P2 ($w = 0.00$)**：深埋於 MHC-I 分子 N 端 A/B 口袋內，TCR 無法接觸。
* **P4、P5、P6、P7 ($w = 0.26 \sim 0.31$)**：胜肽中央拱形隆起的最高頂點，直接與 TCR CDR3$\alpha$ / CDR3$\beta$ 結合環交會。
* **P9 ($w = 0.00$)**：C 端主要錨定點，深插在 F 口袋內。

### 較長胜肽（> 9-mer，如 10-mer、11-mer 等）
MHC Class I 分子溝槽兩端皆為封閉結構，多出的 1～2 個胺基酸無法向兩端延展，只能在中央形成更顯著的**「拱橋效應（Central Zigzag Arch）」**。TRIAD 在演算法中設計了動態擴展機制，在中央插入加權係數 `0.30`：

```python
if peplen > 9:
    weights = POSITION_WEIGHTS[:5] + [0.30] * (peplen - 9) + POSITION_WEIGHTS[5:]
```

---

## 5. HLA 等位基因專屬錨定遮蔽（`ALLELE_POSITION_MASK`）

不同的 HLA 等位基因具有不同的口袋結構（Pockets A 至 F）。當特定位置被 HLA 吸納作為錨定點時，該位置即喪失與 TCR 接觸的機會。TRIAD 內建動態遮蔽機制：

| 呈現之 HLA 型別 | 遮蔽位點 (1-indexed) | 結構生物學機制 |
| :--- | :---: | :--- |
| **標準預設** (如 HLA-A\*02:01, HLA-A\*03:01) | `1, 2, C-term` | 標準 N 端 P2 與 C 端 P9 錨定模式。 |
| **HLA-B\*08:01** | `2, 5, 9` | Pocket D 將 P5 深深吸入作為次要錨定點；因此 P5 不納入 TCR 計分。 |
| **HLA-A\*24:02 / HLA-A\*23:01** | `2, 7, 9` | 具備 P7 輔助錨定特徵。 |
| **HLA-B\*27:05 / HLA-B\*44:02 / HLA-B\*44:03** | `2, 3, 9` | Pocket B/D 強烈吸附 P3 作為主要錨定點。 |
| **HLA-A\*01:01** | `2, 3, 9` | 偏好以 P3 的天門冬胺酸/麩胺酸作為錨定。 |

---

## 6. 分數解讀與 TRIAD 綜合指標整合

### 1. 原始分數 `Immunogenicity_Score`
* 典型數值範圍：**$-0.4000 \sim +0.4000$**。
* **$> 0.00$**：**具備正向免疫原性**（外露位點富含大體積疏水/芳香族殘基，具備誘發 CD8+ T 細胞毒殺之高潛力）。
* **$< 0.00$**：**低免疫原性風險**（外露位點多為小分子或親水殘基，臨床上常屬「結合呈遞但無法活化 T 細胞」之非有效候選者）。

### 2. 標準化指標 $S_{\text{imm}}$ ($[0.0, 1.0]$)
為了便於多維度指標融合，TRIAD 將原始分數線性縮放至 $[0, 1]$ 區間：

$$S_{\text{imm}} = \text{clip}\left(\frac{\text{Immunogenicity\_Score} + 0.4}{0.8}, 0.0, 1.0\right)$$

### 3. 整合至 TRIAD 黃金綜合指數（`Composite_Score`）
TRIAD 將四大關鍵免疫維度合成為單一排序指標，用於精準疫苗與免疫療法候選篩選：

$$\text{Composite\_Score} = \underbrace{0.40 \times S_{\text{pres}}}_{\text{MHC-I 呈遞概率 (%Rank)}} + \underbrace{0.25 \times S_{\text{aff}}}_{\text{結合親和力 (IC}_{50}\text{)}} + \underbrace{\mathbf{0.25 \times S_{\text{imm}}}}_{\mathbf{TCR\ 免疫原性}} + \underbrace{0.10 \times S_{\text{breadth}}}_{\text{HLA 廣效覆蓋度}}$$

---

## 7. 專案核心程式碼對照
* 核心評分器：[`netmhcpan_platform/immunogenicity_scorer.py`](file:///home/cylin/NETMHC/netmhcpan_platform/immunogenicity_scorer.py)
* 批次運算與綜合評分引擎：[`netmhcpan_platform/batch_runner.py`](file:///home/cylin/NETMHC/netmhcpan_platform/batch_runner.py)
* 前端互動視覺化看板：[`netmhcpan_platform/templates/index.html`](file:///home/cylin/NETMHC/netmhcpan_platform/templates/index.html)
