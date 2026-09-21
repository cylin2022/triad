#!/usr/bin/env python3
"""
Convert markdown documentation to professional publication-quality PDF
using Pandoc MathML and headless Chrome.
"""
import os
import sys
import base64
import subprocess

def convert_md_to_pdf(md_path, pdf_path, doc_title="TRIAD Technical Documentation"):
    if not os.path.exists(md_path):
        print(f"Error: {md_path} does not exist.")
        return False

    # Logo base64
    logo_path = "/home/cylin/NETMHC/static/triad_logo.png"
    logo_b64 = ""
    if os.path.exists(logo_path):
        with open(logo_path, "rb") as f:
            logo_b64 = base64.b64encode(f.read()).decode("utf-8")

    # Run pandoc to convert Markdown to HTML fragment with native MathML
    cmd = ["pandoc", md_path, "--mathml"]
    res = subprocess.run(cmd, capture_output=True, text=True, check=True)
    body_html = res.stdout

    full_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>{doc_title}</title>
<style>
  @page {{
    size: A4;
    margin: 16mm 15mm 16mm 15mm;
  }}
  *, *::before, *::after {{
    box-sizing: border-box;
  }}
  body {{
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Noto Sans TC", "PingFang SC", "Microsoft JhengHei", Arial, sans-serif;
    color: #1e293b;
    background: #ffffff;
    line-height: 1.55;
    font-size: 10pt;
    margin: 0;
    padding: 0;
  }}
  .header-banner {{
    display: flex;
    align-items: center;
    justify-content: space-between;
    border-bottom: 2px solid #e2e8f0;
    padding-bottom: 10px;
    margin-bottom: 18px;
  }}
  .header-logo {{
    height: 44px;
    width: auto;
  }}
  .header-meta {{
    font-size: 8.5pt;
    color: #64748b;
    text-align: right;
    line-height: 1.4;
  }}
  h1 {{
    font-size: 17pt;
    font-weight: 700;
    color: #0f172a;
    border-bottom: 2.5px solid #2563eb;
    padding-bottom: 6px;
    margin-top: 0;
    margin-bottom: 14px;
    page-break-after: avoid;
  }}
  h2 {{
    font-size: 12.5pt;
    font-weight: 600;
    color: #1e3a8a;
    border-bottom: 1px solid #e2e8f0;
    padding-bottom: 4px;
    margin-top: 20px;
    margin-bottom: 10px;
    page-break-after: avoid;
  }}
  h3 {{
    font-size: 11pt;
    font-weight: 600;
    color: #2563eb;
    margin-top: 14px;
    margin-bottom: 6px;
    page-break-after: avoid;
  }}
  p, li {{
    color: #334155;
  }}
  ul, ol {{
    padding-left: 20px;
    margin-top: 4px;
    margin-bottom: 10px;
  }}
  li {{
    margin-bottom: 3px;
  }}
  strong {{
    color: #0f172a;
  }}
  table {{
    width: 100%;
    border-collapse: collapse;
    margin: 12px 0;
    font-size: 9pt;
    page-break-inside: avoid;
  }}
  th, td {{
    border: 1px solid #cbd5e1;
    padding: 6px 9px;
    text-align: left;
    vertical-align: top;
  }}
  th {{
    background-color: #f1f5f9;
    color: #0f172a;
    font-weight: 600;
  }}
  tr:nth-child(even) {{
    background-color: #f8fafc;
  }}
  pre, code {{
    font-family: ui-monospace, SFMono-Regular, "JetBrains Mono", Menlo, Consolas, monospace;
    font-size: 8.5pt;
  }}
  code {{
    background-color: #f1f5f9;
    padding: 1px 4px;
    border-radius: 3px;
    color: #0f172a;
    border: 1px solid #e2e8f0;
  }}
  pre {{
    background-color: #f8fafc;
    border: 1px solid #cbd5e1;
    border-radius: 6px;
    padding: 8px 12px;
    overflow-x: auto;
    page-break-inside: avoid;
  }}
  pre code {{
    background: none;
    border: none;
    padding: 0;
  }}
  hr {{
    border: 0;
    height: 1px;
    background: #e2e8f0;
    margin: 16px 0;
  }}
  math {{
    font-size: 1.12em;
  }}
  p > math[display="block"] {{
    margin: 12px 0;
    text-align: center;
  }}
</style>
</head>
<body>
<div class="header-banner">
  <img class="header-logo" src="data:image/png;base64,{logo_b64}" alt="TRIAD Logo">
  <div class="header-meta">
    <strong>TRIAD Technical Documentation</strong><br>
    Laboratory of Systems Biology and Bioinformatics (LSBNB)
  </div>
</div>
{body_html}
</body>
</html>
"""
    tmp_html = md_path.replace(".md", "_export_tmp.html")
    with open(tmp_html, "w", encoding="utf-8") as f:
        f.write(full_html)

    chrome_cmd = [
        "google-chrome",
        "--headless=new",
        "--disable-gpu",
        "--no-sandbox",
        "--no-pdf-header-footer",
        f"--print-to-pdf={pdf_path}",
        tmp_html
    ]
    print(f"Converting {md_path} -> {pdf_path} via Chrome...")
    p = subprocess.run(chrome_cmd, capture_output=True, text=True)
    if os.path.exists(tmp_html):
        os.remove(tmp_html)

    if p.returncode == 0 and os.path.exists(pdf_path):
        size = os.path.getsize(pdf_path)
        print(f"Successfully generated {pdf_path} ({size:,} bytes)")
        return True
    else:
        print(f"Error converting to PDF: {p.stderr}")
        return False

if __name__ == "__main__":
    convert_md_to_pdf(
        "/home/cylin/NETMHC/t-cell_immunogenicity_score.md",
        "/home/cylin/NETMHC/t-cell_immunogenicity_score.pdf",
        "TRIAD: T-Cell Immunogenicity Score Reference"
    )
    convert_md_to_pdf(
        "/home/cylin/NETMHC/t-cell_immunogenicity_score_zh.md",
        "/home/cylin/NETMHC/t-cell_immunogenicity_score_zh.pdf",
        "TRIAD: T細胞免疫原性評分技術說明指南"
    )
