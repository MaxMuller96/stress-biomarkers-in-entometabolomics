#!/usr/bin/env python3
"""
Robustness analysis for Stressbiomarkers_insects.py

Runs the extraction script repeatedly with different tokenization/fuzzy parameters,
compares predictions against a gold-standard file, and writes ROC points + summary.

Expected gold CSV columns:
- pdf
- metabolite

Matching is done case-insensitively and uses InChIKey when available in predictions,
otherwise metabolite name.
"""

import argparse
import itertools
import json
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

import pandas as pd
from Stressbiomarkers_insects import (
    PDF_DIR as DEFAULT_PDF_DIR,
    HMDB_PATH as DEFAULT_HMDB_PATH,
    TAXON_PATH as DEFAULT_TAXON_PATH,
)


@dataclass(frozen=True)
class RunConfig:
    max_n: int
    token_min_len: int
    fuzzy_threshold: int
    fuzzy_len_diff_ratio: float
    fuzzy_min_len: int
    fuzzy_max_len: int
    fuzzy_k: int
    fuzzy_m: int


def parse_int_list(s: str) -> List[int]:
    return [int(x.strip()) for x in s.split(",") if x.strip()]


def parse_float_list(s: str) -> List[float]:
    return [float(x.strip()) for x in s.split(",") if x.strip()]


def load_gold(gold_csv: Path) -> Dict[str, Set[str]]:
    df = pd.read_csv(gold_csv, sep=';', encoding='utf-8-sig')
    required = {"pdf", "metabolite"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Gold CSV is missing required columns: {sorted(missing)}")

    gold: Dict[str, Set[str]] = {}
    for _, row in df.iterrows():
        pdf = str(row["pdf"]).strip()
        met = str(row["metabolite"]).strip().lower()
        if not pdf or not met or met == "nan":
            continue
        gold.setdefault(pdf, set()).add(met)
    return gold


def normalize_prediction_value(v: object) -> str:
    if v is None:
        return ""
    s = str(v).strip()
    if not s or s.lower() == "nan":
        return ""
    return s.lower()


def load_predictions(output_xlsx: Path) -> Dict[str, Set[str]]:
    df = pd.read_excel(output_xlsx)

    required_any = [("InChIKey", "Metabolite name")]
    if not all(any(c in df.columns for c in group) for group in required_any):
        raise ValueError(
            "Prediction output must contain at least one of InChIKey or Metabolite name columns."
        )
    if "DOI references" not in df.columns:
        raise ValueError(
            "Prediction output is missing 'DOI references'. Please update Stressbiomarkers_insects.py to include it in results."
        )

    preds: Dict[str, Set[str]] = {}
    for _, row in df.iterrows():
        pdf = str(row.get("DOI references", "")).strip()
        if not pdf:
            continue

        inchikey = normalize_prediction_value(row.get("InChIKey", ""))
        name = normalize_prediction_value(row.get("Metabolite name", ""))
        key = inchikey or name
        if not key:
            continue

        preds.setdefault(pdf, set()).add(key)
    return preds


def compute_metrics(
    gold: Dict[str, Set[str]],
    preds: Dict[str, Set[str]],
    universe: Set[str],
) -> Dict[str, float]:
    tp = fp = fn = tn = 0

    all_pdfs = set(gold.keys()) | set(preds.keys())
    for pdf in all_pdfs:
        g = gold.get(pdf, set())
        p = preds.get(pdf, set())

        g_u = g & universe
        p_u = p & universe

        tp_set = p_u & g_u
        fp_set = p_u - g_u
        fn_set = g_u - p_u
        tn_set = universe - (tp_set | fp_set | fn_set)

        tp += len(tp_set)
        fp += len(fp_set)
        fn += len(fn_set)
        tn += len(tn_set)

    tpr = tp / (tp + fn) if (tp + fn) else 0.0
    fpr = fp / (fp + tn) if (fp + tn) else 0.0
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tpr

    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "tn": tn,
        "tpr": tpr,
        "fpr": fpr,
        "precision": precision,
        "recall": recall,
    }


def run_once(
    script_path: Path,
    pdf_dir: Path,
    hmdb_path: Path,
    taxon_path: Path,
    output_xlsx: Path,
    cfg: RunConfig,
    timeout_sec: int,
) -> Tuple[int, str, str]:
    cmd = [
        sys.executable,
        str(script_path),
        "--pdf-dir",
        str(pdf_dir),
        "--hmdb-path",
        str(hmdb_path),
        "--taxon-path",
        str(taxon_path),
        "--output-xlsx",
        str(output_xlsx),
        "--max-n",
        str(cfg.max_n),
        "--token-min-len",
        str(cfg.token_min_len),
        "--fuzzy-threshold",
        str(cfg.fuzzy_threshold),
        "--fuzzy-len-diff-ratio",
        str(cfg.fuzzy_len_diff_ratio),
        "--fuzzy-min-len",
        str(cfg.fuzzy_min_len),
        "--fuzzy-max-len",
        str(cfg.fuzzy_max_len),
        "--fuzzy-k",
        str(cfg.fuzzy_k),
        "--fuzzy-m",
        str(cfg.fuzzy_m),
    ]

    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=timeout_sec,
    )
    return proc.returncode, proc.stdout, proc.stderr


def main() -> None:
    ap = argparse.ArgumentParser(description="Run tokenization robustness grid and compute ROC points")
    ap.add_argument("--script-path", default="Stressbiomarkers_insects.py")
    ap.add_argument("--pdf-dir", default=DEFAULT_PDF_DIR)
    ap.add_argument("--hmdb-path", default=DEFAULT_HMDB_PATH)
    ap.add_argument("--taxon-path", default=DEFAULT_TAXON_PATH)
    ap.add_argument(
        "--gold-csv",
        default=str(Path(__file__).resolve().parent / "gold.csv"),
        help="CSV with columns: pdf, metabolite",
    )
    ap.add_argument("--work-dir", default="robustness_runs")
    ap.add_argument("--timeout-sec", type=int, default=14000)

    # Sweep ranges
    ap.add_argument("--max-n-values", default="2")
    ap.add_argument("--token-min-len-values", default="3")
    ap.add_argument("--fuzzy-threshold-values", default="90")
    ap.add_argument("--fuzzy-len-diff-ratio-values", default="0.3")
    ap.add_argument("--fuzzy-min-len-values", default="7")
    ap.add_argument("--fuzzy-max-len-values", default="80")
    ap.add_argument("--fuzzy-k-values", default="4000")
    ap.add_argument("--fuzzy-m-values", default="50")

    args = ap.parse_args()

    script_path = Path(args.script_path).resolve()
    pdf_dir = Path(args.pdf_dir).resolve()
    hmdb_path = Path(args.hmdb_path).resolve()
    taxon_path = Path(args.taxon_path).resolve()
    gold_csv = Path(args.gold_csv).resolve()
    work_dir = Path(args.work_dir).resolve()
    work_dir.mkdir(parents=True, exist_ok=True)

    if not script_path.exists():
        raise FileNotFoundError(f"Script not found: {script_path}")

    gold = load_gold(gold_csv)
    universe = set()
    for mets in gold.values():
        universe |= mets

    max_n_vals = parse_int_list(args.max_n_values)
    token_min_len_vals = parse_int_list(args.token_min_len_values)
    fuzzy_threshold_vals = parse_int_list(args.fuzzy_threshold_values)
    fuzzy_len_diff_ratio_vals = parse_float_list(args.fuzzy_len_diff_ratio_values)
    fuzzy_min_len_vals = parse_int_list(args.fuzzy_min_len_values)
    fuzzy_max_len_vals = parse_int_list(args.fuzzy_max_len_values)
    fuzzy_k_vals = parse_int_list(args.fuzzy_k_values)
    fuzzy_m_vals = parse_int_list(args.fuzzy_m_values)

    grid = [
        RunConfig(*vals)
        for vals in itertools.product(
            max_n_vals,
            token_min_len_vals,
            fuzzy_threshold_vals,
            fuzzy_len_diff_ratio_vals,
            fuzzy_min_len_vals,
            fuzzy_max_len_vals,
            fuzzy_k_vals,
            fuzzy_m_vals,
        )
    ]

    summary_rows = []

    for i, cfg in enumerate(grid, start=1):
        run_name = (
            f"run_{i:04d}"
            f"_n{cfg.max_n}"
            f"_tml{cfg.token_min_len}"
            f"_ft{cfg.fuzzy_threshold}"
            f"_ldr{cfg.fuzzy_len_diff_ratio}"
            f"_fmin{cfg.fuzzy_min_len}"
            f"_fmax{cfg.fuzzy_max_len}"
            f"_k{cfg.fuzzy_k}"
            f"_m{cfg.fuzzy_m}"
        )
        run_dir = work_dir / run_name
        run_dir.mkdir(parents=True, exist_ok=True)

        output_xlsx = run_dir / "metabolite_results.xlsx"
        rc, out, err = run_once(
            script_path=script_path,
            pdf_dir=pdf_dir,
            hmdb_path=hmdb_path,
            taxon_path=taxon_path,
            output_xlsx=output_xlsx,
            cfg=cfg,
            timeout_sec=args.timeout_sec,
        )

        (run_dir / "stdout.txt").write_text(out or "", encoding="utf-8")
        (run_dir / "stderr.txt").write_text(err or "", encoding="utf-8")

        row = {
            "run": run_name,
            "return_code": rc,
            **cfg.__dict__,
        }

        if rc != 0 or not output_xlsx.exists():
            row.update(
                {
                    "tp": None,
                    "fp": None,
                    "fn": None,
                    "tn": None,
                    "tpr": None,
                    "fpr": None,
                    "precision": None,
                    "recall": None,
                    "error": "script_failed_or_no_output",
                }
            )
            summary_rows.append(row)
            continue

        try:
            preds = load_predictions(output_xlsx)
            metrics = compute_metrics(gold, preds, universe)
            row.update(metrics)
            row["error"] = ""
        except Exception as ex:
            row.update(
                {
                    "tp": None,
                    "fp": None,
                    "fn": None,
                    "tn": None,
                    "tpr": None,
                    "fpr": None,
                    "precision": None,
                    "recall": None,
                    "error": f"metrics_error: {ex}",
                }
            )

        summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows)
    summary_csv = work_dir / "robustness_summary.csv"
    summary_df.to_csv(summary_csv, index=False)

    roc_df = (
        summary_df.dropna(subset=["tpr", "fpr"]) 
        .sort_values(["fpr", "tpr"]) 
        [["run", "fpr", "tpr", "precision", "recall", "max_n", "token_min_len", "fuzzy_threshold", "fuzzy_len_diff_ratio", "fuzzy_min_len", "fuzzy_max_len", "fuzzy_k", "fuzzy_m"]]
    )
    roc_csv = work_dir / "roc_points.csv"
    roc_df.to_csv(roc_csv, index=False)

    best = summary_df.dropna(subset=["tpr", "fpr"]).sort_values(["tpr", "precision"], ascending=[False, False]).head(20)
    best_csv = work_dir / "best_runs.csv"
    best.to_csv(best_csv, index=False)

    print(f"Done. Wrote:\n- {summary_csv}\n- {roc_csv}\n- {best_csv}")


if __name__ == "__main__":
    main()
