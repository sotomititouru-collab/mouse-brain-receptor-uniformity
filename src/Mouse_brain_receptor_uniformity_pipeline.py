#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
WMB-10Xv3 receptor uniformity analysis
FINAL PAPER-READY PIPELINE (+ hierarchical bootstrap CIs)

Primary outputs (point estimates):
- region_detection_fractions.csv
- region_cell_counts_primary.csv
- region_compositions.csv
- Figure1_heatmap.png
- uniformity_U75_U95.csv
- aitchison_distance_matrix.csv
- Figure2A_aitchison_heatmap.png
- Figure2B_aitchison_hist.png

Bootstrap outputs (95% CI; hierarchical donor->cell; 1000 reps by default):
- bootstrap_summary.csv              (Aitchison median/p95, U75/U95 medians over genes)
- bootstrap_distributions.csv        (optional: replicate-wise values used for CIs)

Notes:
- Detected expression rule: X > 0.
  NOTE: Your input shards are '*-log2.h5ad'. If X stores log2(UMI+1), then X>0 is exactly UMI>0,
  matching the manuscript definition (UMI > 0).
- Hierarchical bootstrap:
  donors are sampled with replacement;
  within each donor×region stratum, cells are resampled with replacement.
  For Bernoulli detection indicators, cell-bootstrap is equivalent to a Binomial(n, p_hat).
"""

from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
import json
from datetime import datetime

# =============================================================================
# CONFIG
# =============================================================================

H5AD_DIR = Path(r"C:\Users\abeke\abc_cache\expression_matrices\WMB-10Xv3\20230630")
LOG2_PATTERN = "*-log2.h5ad"

CANDIDATE_RECEPTOR_GENES = [
    "Gabra1", "Gabra2", "Gabra3", "Gabrb1", "Gabrb2", "Gabrg2",
    "Htr1a", "Htr2a", "Htr2c", "Htr7",
    "Drd1", "Drd2", "Drd3", "Drd4",
    "Oprk1", "Oprm1", "Oprd1",
    "Grin1", "Grin2a", "Grin2b", "Gria1", "Gria2",
    "Adra2a", "Adra2c", "Chrna4", "Chrnb2",
]

REGION_MAP = {
    "ISO": "Isocortex",
    "ISOCORTEX": "Isocortex",
    "CTXSP": "Cortical subplate",
    "HPF": "Hippocampal formation",
    "HY": "Hypothalamus",
    "STR": "Striatum",
    "TH": "Thalamus",
    "MB": "Midbrain",
    "P": "Pons",
    "MY": "Medulla",
    "OLF": "Olfactory areas",
    "PAL": "Pallidum",
    "CB": None,  # excluded from primary
}

PRIMARY_REGIONS = [
    "Isocortex",
    "Hippocampal formation",
    "Hypothalamus",
    "Striatum",
    "Thalamus",
    "Midbrain",
    "Pons",
    "Medulla",
    "Olfactory areas",
    "Cortical subplate",
    "Pallidum",
]

REGION_KEY_CANDIDATES = [
    "region_label",
    "anatomical_region_label",
    "anatomical_division_label",
    "struct_acronym",
    "structure",
    "structure_acronym",
    "structure_label",
]

# donor/specimen candidates (dataset-dependent; we autodetect)
DONOR_KEY_CANDIDATES = [
    "donor_id",
    "donor",
    "specimen_id",
    "specimen",
    "individual_id",
    "individual",
    "mouse_id",
    "sample_id",
]

OUT_DIR = Path.cwd() / "WMB_receptor_uniformity_results"
OUT_DIR.mkdir(exist_ok=True)

LOG_FILE = OUT_DIR / "analysis_log.txt"
if LOG_FILE.exists():
    LOG_FILE.unlink()

# CLR pseudo-count
DELTA = 1e-6

# Bootstrap config
N_BOOTSTRAP = 1000
RANDOM_SEED = 42
BOOTSTRAP_WRITE_DISTRIBUTIONS = True  # replicate-wise CSV; set False if you want smaller outputs

# =============================================================================
# UTIL
# =============================================================================

def log(msg: str) -> None:
    print(msg)
    with open(LOG_FILE, "a", encoding="utf-8") as f:
        f.write(msg + "\n")

def find_region_key(adata) -> str:
    for k in REGION_KEY_CANDIDATES:
        if k in adata.obs.columns:
            return k
    raise KeyError("No region label column found in .obs")

def find_donor_key_or_fallback(adata):
    """
    Try to find a donor/specimen/individual column.
    If none is found, return None (single-donor fallback).
    """
    for k in DONOR_KEY_CANDIDATES:
        if k in adata.obs.columns:
            return k
    return None

def get_gene_symbol_indexer(adata) -> dict:
    if "gene_symbol" in adata.var.columns:
        symbols = list(map(str, adata.var["gene_symbol"]))
    else:
        symbols = list(map(str, adata.var_names))
    idx = {}
    for i, g in enumerate(symbols):
        if g not in idx:
            idx[g] = i
    return idx

def closure(mat: np.ndarray) -> np.ndarray:
    s = mat.sum(axis=1, keepdims=True)
    s[s == 0] = 1.0
    return mat / s

def clr_transform(mat: np.ndarray) -> np.ndarray:
    mat = np.where(mat <= 0, DELTA, mat)
    logm = np.log(mat)
    return logm - logm.mean(axis=1, keepdims=True)

def compute_uniformity_metrics(comp_df: pd.DataFrame) -> pd.DataFrame:
    """
    For each gene column:
      - compute pairwise absolute differences across regions in percentage points
      - summarize with U75 and U95 (upper quantiles)
    """
    rows = []
    for g in comp_df.columns:
        v = comp_df[g].values
        diffs = pdist(v[:, None], metric=lambda a, b: abs(a - b)) * 100.0  # pp
        rows.append({
            "gene": g,
            "U75_pp": float(np.quantile(diffs, 0.75)),
            "U95_pp": float(np.quantile(diffs, 0.95)),
        })
    return pd.DataFrame(rows).set_index("gene")

def summarize_aitchison_from_distance_matrix(D: np.ndarray) -> dict:
    """
    D is square matrix across PRIMARY_REGIONS (11x11).
    Summaries are computed across the 55 unique unordered pairs.
    """
    vals = squareform(D)
    return {
        "aitchison_median": float(np.median(vals)),
        "aitchison_p95": float(np.quantile(vals, 0.95)),
        "aitchison_max": float(np.max(vals)),
    }

# =============================================================================
# MAIN
# =============================================================================

def main():
    log("=== WMB receptor uniformity analysis (FINAL PIPELINE + BOOTSTRAP) ===")
    log(f"Timestamp: {datetime.now().isoformat()}")
    log(f"Bootstrap: N={N_BOOTSTRAP}, seed={RANDOM_SEED}, delta={DELTA}")

    files = sorted(H5AD_DIR.glob(LOG2_PATTERN))
    log(f"Found {len(files)} h5ad files")

    if len(files) == 0:
        raise FileNotFoundError(f"No files matched pattern {LOG2_PATTERN} in {H5AD_DIR}")

    # -------------------------------------------------------------------------
    # PASS 1: gene set union
    # -------------------------------------------------------------------------
    union = set()
    for f in files:
        ad_b = ad.read_h5ad(f, backed="r")
        union.update(get_gene_symbol_indexer(ad_b).keys())
        ad_b.file.close()

    genes = [g for g in CANDIDATE_RECEPTOR_GENES if g in union]
    if not genes:
        raise RuntimeError("No candidate receptor genes found in the dataset union")

    log(f"Genes used (n={len(genes)}): {genes}")

    # -------------------------------------------------------------------------
    # PASS 2: donor×region aggregation for detection indicators
    #   We store:
    #     donor_region_ncells[donor][region] = n_cells
    #     donor_region_detcounts[donor][region] = array(len(genes)) of detected cell counts
    #
    # This supports:
    #   - point estimates (summing counts across donors)
    #   - hierarchical bootstrap donor->cell using Binomial(n, p_hat) equivalence for Bernoulli
    # -------------------------------------------------------------------------
    donor_region_ncells = {}      # dict[str][str] -> int
    donor_region_detcounts = {}   # dict[str][str] -> np.ndarray(len(genes))

    region_key_used = None
    donor_key_used = None

    for f in files:
        log(f"Processing {f.name}")
        ad_b = ad.read_h5ad(f, backed="r")

        region_key = find_region_key(ad_b)
        donor_key = find_donor_key_or_fallback(ad_b)

        if donor_key is None:
            log(f"WARNING: No donor column found in {f.name}; using single-donor fallback.")
            donor_raw = pd.Series(["_single_donor_"] * ad_b.n_obs, index=ad_b.obs.index)
        else:
            donor_raw = ad_b.obs[donor_key].astype(str).str.strip()
        region_key_used = region_key_used or region_key
        donor_key_used = donor_key_used or donor_key or "_single_donor_fallback_"

        # standardize region codes
        region_raw = ad_b.obs[region_key].astype(str).str.strip().str.upper()

        # map to macro-region
        macro = region_raw.map(REGION_MAP)
        valid = macro.notna() & macro.isin(PRIMARY_REGIONS)

        if valid.sum() == 0:
            ad_b.file.close()
            continue

        macro_valid = macro[valid].values
        donor_valid = donor_raw[valid].values

        # gene index mapping for this file
        sym2idx = get_gene_symbol_indexer(ad_b)
        present = [g for g in genes if g in sym2idx]

        if not present:
            log("  No target genes present in this shard; skipping detections.")
            ad_b.file.close()
            continue

        # load only present genes to memory for detection (boolean)
        sub = ad_b[valid.values, [sym2idx[g] for g in present]].to_memory()

        # NOTE: original pipeline uses log2(expression) > 0 for detection.
        # That implies the stored matrix is log2-normalized and >0 indicates detected.
        # We keep it consistent with your existing code & output expectations.
        detected = (sub.X > 0)

        # ensure we can slice & sum; support sparse matrices
        # (sum over rows returns 1 x G matrix)
        try:
            import scipy.sparse as sp
            is_sparse = sp.issparse(detected)
        except Exception:
            is_sparse = False

        # group indices by (macro, donor)
        df_key = pd.DataFrame({"macro": macro_valid, "donor": donor_valid})
        groups = df_key.groupby(["macro", "donor"]).indices

        # update donor-region aggregates
        for (reg, donor), idx in groups.items():
            n = int(len(idx))
            if n == 0:
                continue

            # sum detections for present genes within this donor×region×shard
            if is_sparse:
                det_sum_present = np.asarray(detected[idx, :].sum(axis=0)).reshape(-1)
            else:
                det_sum_present = np.asarray(detected[idx, :].sum(axis=0)).reshape(-1)

            if donor not in donor_region_ncells:
                donor_region_ncells[donor] = {}
                donor_region_detcounts[donor] = {}

            if reg not in donor_region_ncells[donor]:
                donor_region_ncells[donor][reg] = 0
                donor_region_detcounts[donor][reg] = np.zeros(len(genes), dtype=float)

            donor_region_ncells[donor][reg] += n

            # map present-gene sums into full gene vector
            for j, g in enumerate(present):
                full_idx = genes.index(g)
                donor_region_detcounts[donor][reg][full_idx] += float(det_sum_present[j])

        ad_b.file.close()

    if not donor_region_ncells:
        raise RuntimeError("No valid donor×region data accumulated. Check REGION_MAP and donor/region keys.")

    donors = sorted(donor_region_ncells.keys())
    log(f"Detected donor count: {len(donors)}")
    log(f"Region key used: {region_key_used}")
    log(f"Donor key used: {donor_key_used}")

    # -------------------------------------------------------------------------
    # POINT ESTIMATES (primary)
    # -------------------------------------------------------------------------
    region_cell_counts = {r: 0 for r in PRIMARY_REGIONS}
    detected_counts = {r: np.zeros(len(genes), dtype=float) for r in PRIMARY_REGIONS}

    for donor in donors:
        for reg in donor_region_ncells[donor].keys():
            n = donor_region_ncells[donor][reg]
            region_cell_counts[reg] += int(n)
            detected_counts[reg] += donor_region_detcounts[donor][reg]

    frac = np.vstack([
        detected_counts[r] / region_cell_counts[r] if region_cell_counts[r] > 0 else np.zeros(len(genes))
        for r in PRIMARY_REGIONS
    ])

    frac_df = pd.DataFrame(frac, index=PRIMARY_REGIONS, columns=genes)
    frac_df.to_csv(OUT_DIR / "region_detection_fractions.csv")

    pd.DataFrame(
        {"region": PRIMARY_REGIONS, "n_cells": [region_cell_counts[r] for r in PRIMARY_REGIONS]}
    ).set_index("region").to_csv(OUT_DIR / "region_cell_counts_primary.csv")

    comp = closure(frac)
    comp_df = pd.DataFrame(comp, index=PRIMARY_REGIONS, columns=genes)
    comp_df.to_csv(OUT_DIR / "region_compositions.csv")

    # Figure 1
    plt.figure(figsize=(12, 6))
    plt.imshow(comp_df.values, aspect="auto")
    plt.yticks(range(len(PRIMARY_REGIONS)), PRIMARY_REGIONS)
    plt.xticks(range(len(genes)), genes, rotation=90, fontsize=7)
    plt.colorbar(label="Within-region composition")
    plt.title("Receptor subtype compositions across brain macro-regions")
    plt.tight_layout()
    plt.savefig(OUT_DIR / "Figure1_heatmap.png", dpi=300)
    plt.close()

    # Uniformity metrics (point)
    U = compute_uniformity_metrics(comp_df)
    U.to_csv(OUT_DIR / "uniformity_U75_U95.csv")

    # Aitchison distance matrix (point)
    D = squareform(pdist(clr_transform(comp), metric="euclidean"))
    D_df = pd.DataFrame(D, index=PRIMARY_REGIONS, columns=PRIMARY_REGIONS)
    D_df.to_csv(OUT_DIR / "aitchison_distance_matrix.csv")

    # Figure 2a: Aitchison heatmap
    plt.figure(figsize=(6, 5))
    im = plt.imshow(D_df.values, cmap="viridis")
    plt.xticks(range(len(PRIMARY_REGIONS)), PRIMARY_REGIONS, rotation=90)
    plt.yticks(range(len(PRIMARY_REGIONS)), PRIMARY_REGIONS)
    cbar = plt.colorbar(im)
    cbar.set_label("Aitchison distance (CLR space)")
    plt.title("Cross-regional Aitchison distances")
    plt.tight_layout()
    plt.savefig(OUT_DIR / "Figure2A_aitchison_heatmap.png", dpi=300)
    plt.close()

    # Figure 2b: Distance histogram (point)
    vals_point = squareform(D)
    plt.figure(figsize=(5, 4))
    plt.hist(vals_point, bins=20, color="gray", edgecolor="black")
    plt.xlabel("Aitchison distance")
    plt.ylabel("Number of region pairs")
    plt.title("Distribution of cross-regional Aitchison distances")
    plt.tight_layout()
    plt.savefig(OUT_DIR / "Figure2B_aitchison_hist.png", dpi=300)
    plt.close()

    # -------------------------------------------------------------------------
    # HIERARCHICAL BOOTSTRAP (donor -> cell)
    # For Bernoulli detection indicators, resampling cells within donor×region is
    # equivalent to Binomial(n, p_hat) with p_hat = detected_count/n.
    # -------------------------------------------------------------------------
    rng = np.random.default_rng(RANDOM_SEED)

    boot_ait_median = np.zeros(N_BOOTSTRAP, dtype=float)
    boot_ait_p95 = np.zeros(N_BOOTSTRAP, dtype=float)
    boot_U75_median = np.zeros(N_BOOTSTRAP, dtype=float)
    boot_U95_median = np.zeros(N_BOOTSTRAP, dtype=float)

    # Pre-pack donor data for speed
    donor_payload = []
    for donor in donors:
        regs = list(donor_region_ncells[donor].keys())
        donor_payload.append((donor, regs))

    log("Starting hierarchical bootstrap...")
    for b in range(N_BOOTSTRAP):
        sampled_donors = rng.choice(donors, size=len(donors), replace=True)

        # aggregate bootstrap counts per region
        reg_n = {r: 0 for r in PRIMARY_REGIONS}
        reg_det = {r: np.zeros(len(genes), dtype=float) for r in PRIMARY_REGIONS}

        for donor in sampled_donors:
            for reg in donor_region_ncells[donor].keys():
                n = int(donor_region_ncells[donor][reg])
                if n <= 0:
                    continue

                det = donor_region_detcounts[donor][reg].astype(float)
                # empirical p_hat per gene
                p = det / n
                p = np.clip(p, 0.0, 1.0)

                # cell-bootstrap within donor×region (Bernoulli -> Binomial equivalence)
                det_b = rng.binomial(n=n, p=p)

                reg_n[reg] += n
                reg_det[reg] += det_b

        # compute bootstrap fractions/compositions
        frac_b = np.vstack([
            reg_det[r] / reg_n[r] if reg_n[r] > 0 else np.zeros(len(genes))
            for r in PRIMARY_REGIONS
        ])
        comp_b = closure(frac_b)

        # bootstrap U75/U95 (gene-wise), then median across genes (as in PDF text)
        comp_b_df = pd.DataFrame(comp_b, index=PRIMARY_REGIONS, columns=genes)
        U_b = compute_uniformity_metrics(comp_b_df)
        boot_U75_median[b] = float(U_b["U75_pp"].median())
        boot_U95_median[b] = float(U_b["U95_pp"].median())

        # bootstrap Aitchison summaries
        D_b = squareform(pdist(clr_transform(comp_b), metric="euclidean"))
        vals_b = squareform(D_b)
        boot_ait_median[b] = float(np.median(vals_b))
        boot_ait_p95[b] = float(np.quantile(vals_b, 0.95))

        if (b + 1) % 50 == 0:
            log(f"  bootstrap {b + 1}/{N_BOOTSTRAP}")

    # -------------------------------------------------------------------------
    # BOOTSTRAP SUMMARY (use point estimate as "estimate", CI from bootstrap distribution)
    # -------------------------------------------------------------------------
    point_ait = summarize_aitchison_from_distance_matrix(D)
    point_U75_med = float(U["U75_pp"].median())
    point_U95_med = float(U["U95_pp"].median())

    def ci95(x: np.ndarray):
        return float(np.quantile(x, 0.025)), float(np.quantile(x, 0.975))

    ait_med_lo, ait_med_hi = ci95(boot_ait_median)
    ait_p95_lo, ait_p95_hi = ci95(boot_ait_p95)
    U75_lo, U75_hi = ci95(boot_U75_median)
    U95_lo, U95_hi = ci95(boot_U95_median)

    summary = pd.DataFrame([
        {
            "metric": "Aitchison_distance_median",
            "estimate": point_ait["aitchison_median"],
            "ci_lower": ait_med_lo,
            "ci_upper": ait_med_hi,
        },
        {
            "metric": "Aitchison_distance_95th_percentile",
            "estimate": point_ait["aitchison_p95"],
            "ci_lower": ait_p95_lo,
            "ci_upper": ait_p95_hi,
        },
        {
            "metric": "Aitchison_distance_max",
            "estimate": point_ait["aitchison_max"],
            "ci_lower": "",
            "ci_upper": "",
        },
        {
            "metric": "Gene-wise_U75_median_pp",
            "estimate": point_U75_med,
            "ci_lower": U75_lo,
            "ci_upper": U75_hi,
        },
        {
            "metric": "Gene-wise_U95_median_pp",
            "estimate": point_U95_med,
            "ci_lower": U95_lo,
            "ci_upper": U95_hi,
        },
    ])

    summary.to_csv(OUT_DIR / "bootstrap_summary.csv", index=False)

    if BOOTSTRAP_WRITE_DISTRIBUTIONS:
        dist = pd.DataFrame({
            "bootstrap_rep": np.arange(1, N_BOOTSTRAP + 1),
            "aitchison_median": boot_ait_median,
            "aitchison_p95": boot_ait_p95,
            "U75_median_pp": boot_U75_median,
            "U95_median_pp": boot_U95_median,
        })
        dist.to_csv(OUT_DIR / "bootstrap_distributions.csv", index=False)

    # -------------------------------------------------------------------------
    # Metadata
    # -------------------------------------------------------------------------
    meta = {
        "dataset": "Allen Cell Census v2023 (WMB-10Xv3)",
        "snapshot": "2023-06-30",
        "files_processed": [f.name for f in files],
        "genes_used": genes,
        "regions": PRIMARY_REGIONS,
        "detection_rule": "log2(expression) > 0",
        "num_cells_total_primary": int(sum(region_cell_counts.values())),
        "region_key_used": region_key_used,
        "donor_key_used": donor_key_used,
        "clr_pseudocount_delta": DELTA,
        "bootstrap": {
            "enabled": True,
            "n_replicates": N_BOOTSTRAP,
            "random_seed": RANDOM_SEED,
            "ci_method": "percentile",
            "hierarchy": "donor->cell (within donor×region binomial equivalence for Bernoulli detection)",
        },
        "outputs": [
            "region_detection_fractions.csv",
            "region_cell_counts_primary.csv",
            "region_compositions.csv",
            "Figure1_heatmap.png",
            "uniformity_U75_U95.csv",
            "aitchison_distance_matrix.csv",
            "Figure2A_aitchison_heatmap.png",
            "Figure2B_aitchison_hist.png",
            "bootstrap_summary.csv",
            "bootstrap_distributions.csv" if BOOTSTRAP_WRITE_DISTRIBUTIONS else None,
        ],
    }
    meta["outputs"] = [x for x in meta["outputs"] if x is not None]

    with open(OUT_DIR / "run_metadata.json", "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2)

    log("=== Analysis completed successfully ===")
    log(f"Outputs written to: {OUT_DIR}")

if __name__ == "__main__":
    main()
