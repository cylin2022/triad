
如果目標是「從一批 peptide 中有效篩出最可能成為 CD8+ T-cell epitope 的 candidates」，我建議不要只呈現單一 score，而是做成一個多層篩選表 + 最終 ranking score。這樣最容易解讀，也最適合後續實驗驗證。

可以把結果分成四層：

MHC binding / presentation
NetMHCpan %Rank_EL
NetMHCpan %Rank_BA
affinity（nM，可選）
建議先以 %Rank_EL 為主，因為比較接近實際 presentation。
Antigen processing
如果有跑 IEDB processing，就加入：
proteasome score
TAP score
processing score
若沒有，NetMHCpan EL 本身已經可以作為較實用的第一層 approximation。
Immunogenicity
IEDB Class I immunogenicity score
分數越高，代表 pMHC 被 T cell 辨識的潛力越高，但不要單獨拿來做 cutoff。
Population relevance
peptide 可結合多少常見 HLA allele
population coverage
是否跨多個 HLA-A / HLA-B types

最實用的結果表可以長這樣：

Peptide HLA %Rank_EL %Rank_BA Immunogenicity Processing Population freq. Final rank
PEPTIDE1 HLA-A*02:01 0.12 0.35 0.42 1.8 25% 1
PEPTIDE2 HLA-A*11:01 0.31 0.70 0.58 1.4 18% 2
PEPTIDE3 HLA-B*07:02 0.95 1.6 0.71 1.1 12% 3

篩選邏輯我會建議採用「先 hard filter、再 rank」，而不是一開始就把所有 score 混在一起。

第一層先去掉明顯不好的：

NetMHCpan %Rank_EL
    ↓
Strong binder: <0.5%
Weak binder:   0.5–2%
Non-binder:    >2%

通常可以先保留：

%Rank_EL < 2%

如果 candidates 很多，可以更嚴格：

%Rank_EL < 0.5%

接著才用 immunogenicity 做排序。例如：

Strong presentation
        +
Positive immunogenicity score
        +
Common HLA allele
        ↓
High-priority candidate

因此不是：

highest immunogenicity score=best candidate

而應該是：

good presentation+good immunogenicity+population relevance

我會建議最後分成三個 tier，這在論文或報告中很直觀：

Tier 建議條件 解讀
Tier 1 %Rank_EL <0.5% + positive immunogenicity + common HLA 最優先驗證
Tier 2 %Rank_EL 0.5–2% + positive immunogenicity 次優先
Tier 3 binding 強但 immunogenicity 較低，或 HLA 較少見 備選

如果一條 peptide 可以同時 bind 多個 HLA，應該額外加分，例如：

Peptide A
HLA-A*02:01   EL rank 0.14
HLA-A*24:02   EL rank 0.41
HLA-A*11:01   EL rank 0.62

通常會比只對單一 rare HLA 有極強 binding 的 peptide 更值得優先考慮，尤其你要做的是 vaccine 或廣泛族群適用的 epitope。

因此最終可以建立一個簡單 composite score，例如：

Candidate Score=w
1
 
(Presentation)+w
2
 
(Immunogenicity)+w
3
 
(HLA breadth)+w
4
 
(Population coverage)

不過我不建議直接拿 raw scores 相加，因為它們的尺度不同。比較好的做法是先各自 normalize 成 0–1，再加權。

例如：

S=0.4S
EL
 
+0.3S
Immunogenicity
 
+0.15S
HLAbreadth
 
+0.15S
Population
 
其中 EL 建議權重最高，因為 peptide 如果根本不會被呈現，後面的 immunogenicity 再高也沒有太大意義。

圖形呈現方面，我會非常推薦再加一張 2D candidate map：

                    High immunogenicity
                           ↑
                           │
        Medium priority   │   ★ Best candidates
                           │
                           │
───────────────────────────┼────────────→
                           │       Stronger
                           │       presentation
                           │
        Low priority      │   Presentation only

X 軸用：

−log
10
 
(%Rank
EL
 
)

所以越右邊代表 presentation 越強。

Y 軸：

IEDB Immunogenicity Score

再讓：

bubble size = population coverage
label = peptide ID
shape = HLA-A / HLA-B / HLA-C

這樣一眼就能找到右上角的 candidates。

如果是我替你設計最後的分析結果，我會做成：

Table 1：完整 NetMHCpan + IEDB 結果

→ Table 2：Tier 1 / Tier 2 shortlist

→ Figure：Presentation vs Immunogenicity bubble plot

→ Top 5–20 peptides 進入實驗驗證

而實驗驗證則再接：

Top candidates
      ↓
pMHC binding / stability
      ↓
PBMC stimulation
      ↓
IFN-γ ELISpot / intracellular cytokine staining
      ↓
Tetramer validation
      ↓
Validated T-cell epitope

這樣就形成一條很完整的：

Prediction→Prioritization→Experimental validation
 
pipeline。

如果你的目標不是「找 epitope」，而是相反地要從 therapeutic peptide 中排除可能造成 human immunogenicity 的 peptide，同一套結果也能直接反過來使用：右上角就變成你最想排除的高風險 candidates。
