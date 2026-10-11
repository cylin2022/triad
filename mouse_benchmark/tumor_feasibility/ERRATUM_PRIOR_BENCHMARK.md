# 前次廣泛基準的資料限制更正（2026-10-10）

在腫瘤情境清理時確認：`tcell_search.linear_sequence` 指的是 epitope，未必是
該筆 assay 真正測試的 antigen。應同時核對 `tcell_export.assay_antigen__name`。

例如 PMID 40586324 的 assay 23044448 是 mutant LTVRYGRFLNL 的陽性反應；
assay 23044468 的實際 antigen 則是 WT MTVRYGRFLNL，結果陰性。摘要卻把兩者
都列在 mutant epitope 下。前次按 epitope–allele 合併時會把這類紀錄誤判為衝突。
這個研究共確認 7 組此類 WT 對照；新流程保留真正 mutant 的標籤，排除 WT assay。

另如 PMID 33879241，同一 mutant 在未接種／僅帶瘤時沒有反應，接種後有反應，
不應在未區分情境時直接合併。本輪已限定 vaccination endpoint。

因此前次 `mouse_benchmark/results/metrics.csv` 與報告仍只是較粗略摘要整理的
探索結果，**尚未經 assay-antigen 層級全面重算**；不能用來宣稱可靠的 mouse
immunogenicity 準確度或作為新模型訓練標籤。原始檔案保留以便追溯，沒有偷偷
改寫舊結果。本次腫瘤模型從原始 assay 及詳細 export 重新整理，不使用前次合併標籤。
