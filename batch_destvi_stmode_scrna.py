
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import gc
import re
import sys
import json
import time
import glob
import argparse
import traceback
import logging
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp


# -----------------------------
# Parse args early
# -----------------------------
def parse_args():
    example = r"""
Examples:
 nohup python batch_destvi_bin50_cli_optimized.py \
    --st-dir /your/st/bin50_dir \
    --obs-csv /your/scanvi_obs.csv \
    --out-dir /your/output_dir \
    --scrna-10x-dir /your/scrna_10x_dir \
    --obs-index-is-barcode \
    --label-col coarse_label_scanvi \
    --batch-col sample \
    --st-symbol-col real_gene_name \
    --retrain-reference \
    --num-threads 1 \
    --num-interop-threads 1 \
    > /your/output_dir/master.stdout.log 2>&1 &

 nohup python batch_destvi_bin50_cli_optimized.py \
    --st-dir /your/st/bin50_dir \
    --obs-csv /your/scanvi_obs.csv \
    --out-dir /your/output_dir_run2 \
    --scrna-10x-dir /your/scrna_10x_dir \
    --obs-index-is-barcode \
    --label-col coarse_label_scanvi \
    --batch-col sample \
    --st-symbol-col real_gene_name \
    --reference-h5ad /your/output_dir/reference/scrna_condscvi_reference.h5ad \
    --reference-model-dir /your/output_dir/reference/condscvi_model \
    --num-threads 1 \
    --num-interop-threads 1 \
    > /your/output_dir_run2/master.stdout.log 2>&1 &

 Cellbin example using an existing reference:
 nohup python batch_destvi_bin50_cellbin_cli.py \
    --st-mode cellbin \
    --st-dir /your/st/cellbin_dir \
    --obs-csv /your/scanvi_obs.csv \
    --out-dir /your/output_cellbin \
    --obs-index-is-barcode \
    --label-col coarse_label_scanvi \
    --batch-col sample \
    --st-gene-name-source var_names \
    --reference-h5ad /your/output_dir/reference/scrna_condscvi_reference.h5ad \
    --num-threads 1 \
    --num-interop-threads 1 \
    > /your/output_cellbin/master.stdout.log 2>&1 &
"""

    parser = argparse.ArgumentParser(
        description=(
            "Batch DestVI annotation for Stereo-seq bin50/cellbin h5ad files. "
            "This optimized version isolates each sample in a fresh subprocess "
            "to avoid GPU memory residue across looped training runs."
        ),
        epilog=example,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # regular args (validated later so worker mode can reuse the same parser)
    parser.add_argument("--st-dir", default=None, help="Directory containing Stereo-seq ST h5ad files")
    parser.add_argument("--obs-csv", default=None, help="scANVI-corrected obs.csv for scRNA")
    parser.add_argument("--out-dir", default=None, help="Output directory")

    # scRNA input
    parser.add_argument("--scrna-input", default="10x", choices=["10x", "h5ad"],
                        help="Input format for retraining the scRNA reference: 10x mtx directory or prepared scRNA h5ad")
    parser.add_argument("--scrna-10x-dir", default=None, help="Directory of scRNA 10x mtx; can also supplement h5ad counts when requested")
    parser.add_argument("--scrna-h5ad", default=None, help="Prepared scRNA AnnData h5ad used when --scrna-input h5ad")
    parser.add_argument("--scrna-var-names-mode", default="gene_symbols", choices=["gene_symbols", "gene_ids"])
    parser.add_argument("--scrna-h5ad-gene-name-source", default="var_names", choices=["var_names", "symbol_col"],
                        help="Where to get scRNA h5ad gene symbols: var_names or a column in adata.var")
    parser.add_argument("--scrna-h5ad-symbol-col", default="gene_symbols",
                        help="Gene-symbol column in scRNA h5ad var when --scrna-h5ad-gene-name-source symbol_col")
    parser.add_argument("--supplement-10x-counts", action="store_true", default=False,
                        help="When --scrna-input h5ad, replace the h5ad count layer with matching raw counts from --scrna-10x-dir for overlapping barcodes/genes")

    # obs.csv alignment
    parser.add_argument("--obs-index-is-barcode", action="store_true", default=False,
                        help="obs.csv index is barcode")
    parser.add_argument("--barcode-col", default=None, help="barcode column in obs.csv if index is not barcode")

    # columns
    parser.add_argument("--label-col", default="coarse_label_scanvi", help="label column in obs.csv")
    parser.add_argument("--batch-col", default="sample", help="batch column in obs.csv")
    parser.add_argument("--st-symbol-col", default="real_gene_name", help="gene symbol column in ST var, used when --st-gene-name-source symbol_col")
    parser.add_argument("--st-gene-name-source", default="symbol_col", choices=["symbol_col", "var_names"],
                        help="Where to get ST gene symbols: symbol_col for Ensembl IDs plus symbol column; var_names if adata.var_names are already gene symbols")

    # label handling for CondSCVI / DestVI output proportions
    parser.add_argument("--cond-label-source", default="scanvi_pred", choices=["scanvi_pred", "original"],
                        help=(
                            "Which labels to train CondSCVI/DestVI on. "
                            "scanvi_pred: train SCANVI using --label-col, then use predicted labels. "
                            "original: use --label-col directly, e.g. leiden clusters, without replacing them by SCANVI predictions."
                        ))
    parser.add_argument("--cond-label-col", default=None,
                        help=(
                            "Column name to store/use as the CondSCVI label key. "
                            "Default: <label-col>_scanvi_pred when --cond-label-source scanvi_pred; "
                            "otherwise the same as --label-col."
                        ))

    # retrain / reuse reference
    parser.add_argument("--retrain-reference", action="store_true",
                        help="Retrain SCVI/SCANVI/CondSCVI reference")
    parser.add_argument("--reference-h5ad", default=None,
                        help="Existing CondSCVI reference h5ad, used when not retraining")
    parser.add_argument("--reference-model-dir", default=None,
                        help="Existing CondSCVI model dir, kept for compatibility")

    # thread settings
    parser.add_argument("--num-threads", type=int, default=1)
    parser.add_argument("--num-interop-threads", type=int, default=1)
    parser.add_argument("--auto-benchmark-threads", action="store_true", default=False)

    # gene selection
    parser.add_argument("--use-hvg", action="store_true", default=False,
                        help="Use HVG + custom genes for reference")
    parser.add_argument("--n-top-genes", type=int, default=3000)
    parser.add_argument("--custom-genes", nargs="*", default=[],
                        help="Custom genes to force keep in reference")
    parser.add_argument("--hep-bg-file", default=None,
                        help="Optional txt file, one background gene per line")

    # training params
    parser.add_argument("--scvi-max-epochs", type=int, default=200)
    parser.add_argument("--scanvi-max-epochs", type=int, default=50)
    parser.add_argument("--condscvi-max-epochs", type=int, default=200)
    parser.add_argument("--destvi-max-epochs", type=int, default=1200)

    parser.add_argument("--scvi-batch-size", type=int, default=2048)
    parser.add_argument("--scanvi-batch-size", type=int, default=2048)
    parser.add_argument("--condscvi-batch-size", type=int, default=2048)
    parser.add_argument("--destvi-batch-size", type=int, default=256)

    parser.add_argument("--n-latent", type=int, default=30)
    parser.add_argument("--n-layers", type=int, default=2)
    parser.add_argument("--gene-likelihood", default="nb")

    # qc
    parser.add_argument("--min-total-counts", type=int, default=30)
    parser.add_argument("--min-genes-by-counts", type=int, default=30)

    # misc
    parser.add_argument("--counts-layer", default="counts")
    parser.add_argument("--unlabeled", default="Unknown")
    parser.add_argument("--st-mode", default="bin50", choices=["bin50", "cellbin", "custom"],
                        help="ST input type. bin50 uses *bin50*.h5ad; cellbin uses *cellbin*.h5ad; custom requires/uses --st-glob")
    parser.add_argument("--st-glob", default=None,
                        help="Optional custom glob for ST files. Overrides --st-mode default glob when provided")
    parser.add_argument("--accelerator", default="auto")
    parser.add_argument("--devices", default="auto")
    parser.add_argument("--save-destvi-model", action="store_true", default=False,
                        help="Save DestVI model dir for each sample (default: False, saves memory/disk pressure)")

    # hidden worker args
    parser.add_argument("--worker-mode", action="store_true", default=False, help=argparse.SUPPRESS)
    parser.add_argument("--worker-st-path", default=None, help=argparse.SUPPRESS)
    parser.add_argument("--worker-reference-h5ad", default=None, help=argparse.SUPPRESS)
    parser.add_argument("--worker-result-json", default=None, help=argparse.SUPPRESS)

    return parser.parse_args()


ARGS = parse_args()


# -----------------------------
# Set env before torch import
# -----------------------------
def set_thread_env(n_threads=1):
    os.environ["OMP_NUM_THREADS"] = str(n_threads)
    os.environ["MKL_NUM_THREADS"] = str(n_threads)
    os.environ["OPENBLAS_NUM_THREADS"] = str(n_threads)
    os.environ["NUMEXPR_NUM_THREADS"] = str(n_threads)
    os.environ["VECLIB_MAXIMUM_THREADS"] = str(n_threads)
    os.environ["BLIS_NUM_THREADS"] = str(n_threads)


set_thread_env(ARGS.num_threads)

import torch
import scvi


# -----------------------------
# Logging
# -----------------------------
def setup_logger(log_path, mode="w"):
    logger_name = f"destvi::{os.path.abspath(log_path)}::{mode}"
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)
    logger.propagate = False
    logger.handlers = []
    fmt = logging.Formatter("[%(asctime)s] %(levelname)s - %(message)s")

    fh = logging.FileHandler(log_path, mode=mode, encoding="utf-8")
    fh.setLevel(logging.INFO)
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(logging.INFO)
    sh.setFormatter(fmt)
    logger.addHandler(sh)
    return logger


# -----------------------------
# Validation / thread control
# -----------------------------
def validate_args():
    required = ["st_dir", "out_dir"]
    if ARGS.worker_mode:
        required += ["worker_st_path", "worker_reference_h5ad", "worker_result_json"]
    else:
        if ARGS.retrain_reference:
            if ARGS.scrna_input == "10x":
                required += ["scrna_10x_dir", "obs_csv"]
            elif ARGS.scrna_input == "h5ad":
                required += ["scrna_h5ad"]
                # obs_csv is optional for h5ad only if label/batch columns already exist in adata.obs.
                if ARGS.obs_csv not in [None, ""] and (not ARGS.obs_index_is_barcode) and (ARGS.barcode_col is None):
                    raise ValueError("Need --barcode-col when --obs-index-is-barcode is not set and --obs-csv is provided")
            if ARGS.supplement_10x_counts and ARGS.scrna_10x_dir in [None, ""]:
                raise ValueError("--supplement-10x-counts requires --scrna-10x-dir")
        else:
            required += ["reference_h5ad"]

    missing = [k for k in required if getattr(ARGS, k) in [None, ""]]
    if missing:
        raise ValueError(f"Missing required arguments: {missing}")

    if ARGS.scrna_input == "10x" or ARGS.obs_csv not in [None, ""]:
        if (not ARGS.obs_index_is_barcode) and (ARGS.barcode_col is None):
            raise ValueError("Need --barcode-col when --obs-index-is-barcode is not set")



def apply_torch_threads(n_threads=1, n_interop=1):
    set_thread_env(n_threads)
    torch.set_num_threads(n_threads)
    try:
        torch.set_num_interop_threads(n_interop)
    except RuntimeError:
        pass


# -----------------------------
# Memory utilities
# -----------------------------
def cuda_cleanup(logger=None, note=None):
    gc.collect()
    if torch.cuda.is_available():
        try:
            torch.cuda.synchronize()
        except Exception:
            pass
        try:
            torch.cuda.empty_cache()
        except Exception:
            pass
        try:
            torch.cuda.ipc_collect()
        except Exception:
            pass

        if logger is not None:
            try:
                allocated = torch.cuda.memory_allocated() / 1024**2
                reserved = torch.cuda.memory_reserved() / 1024**2
                prefix = f"[{note}] " if note else ""
                logger.info(f"{prefix}CUDA memory | allocated={allocated:.2f} MB, reserved={reserved:.2f} MB")
            except Exception:
                pass



def release_model(model):
    if model is None:
        return
    try:
        module = getattr(model, "module", None)
        if module is not None:
            module.cpu()
    except Exception:
        pass
    try:
        trainer = getattr(model, "trainer", None)
        if trainer is not None:
            try:
                trainer.strategy.teardown()
            except Exception:
                pass
    except Exception:
        pass


# -----------------------------
# Utilities
# -----------------------------
def ensure_sparse_float32(x):
    if sp.issparse(x):
        return x.astype(np.float32).tocsr()
    return sp.csr_matrix(np.asarray(x, dtype=np.float32))



def clean_categories(series, unlabeled):
    # Important for Leiden: pandas may read clusters as int/float.
    # scvi-tools needs labels to be categorical/string-like, not numeric.
    s = series.astype(str).copy()
    s = s.fillna(unlabeled)
    s.loc[s.isna()] = unlabeled
    s.loc[s.str.lower().isin(["nan", "none", ""])] = unlabeled
    # Turn 1.0, 2.0 into 1, 2 when a CSV numeric column was read as float.
    s = s.str.replace(r"\.0$", "", regex=True)
    return pd.Categorical(s)



def get_cond_label_col():
    if ARGS.cond_label_col not in [None, ""]:
        return ARGS.cond_label_col
    if ARGS.cond_label_source == "original":
        return ARGS.label_col
    return f"{ARGS.label_col}_scanvi_pred"



def get_reference_label_col_or_raise(adata, logger=None):
    desired = get_cond_label_col()
    if desired in adata.obs.columns:
        return desired

    # Helpful failure mode for old references made by the previous script.
    if "coarse_label_scanvi" in adata.obs.columns and desired != "coarse_label_scanvi":
        raise ValueError(
            f"Reference h5ad contains old label column 'coarse_label_scanvi' but this run expects '{desired}'. "
            f"Changing --label-col only affects a newly trained reference. Re-run with --retrain-reference "
            f"and --label-col {ARGS.label_col}, or set --cond-label-col coarse_label_scanvi if you intentionally want to reuse the old reference."
        )

    candidate_cols = [c for c in adata.obs.columns if c.endswith("_scanvi_pred") or c in [ARGS.label_col, "coarse_label_scanvi"]]
    raise ValueError(
        f"Could not find CondSCVI label column '{desired}' in reference h5ad. "
        f"Available likely label columns: {candidate_cols}. "
        f"If you changed --label-col, rebuild the reference with --retrain-reference."
    )



def sample_name_from_path(path):
    name = Path(path).name
    name = re.sub(r"\.h5ad$", "", name)
    return name


def resolve_st_glob():
    """Resolve the ST file pattern from --st-mode, unless --st-glob is provided."""
    if ARGS.st_glob not in [None, ""]:
        return ARGS.st_glob
    if ARGS.st_mode == "bin50":
        return "*bin50*.h5ad"
    if ARGS.st_mode == "cellbin":
        return "*cellbin*.h5ad"
    raise ValueError("--st-mode custom requires --st-glob, e.g. --st-glob '*cellbin*.h5ad'")


def read_gene_list(path):
    if path is None:
        return []
    with open(path, "r", encoding="utf-8") as f:
        genes = [x.strip() for x in f if x.strip()]
    return genes


def _make_h5ad_serializable(value):
    """
    Convert metadata values that AnnData/HDF5 cannot write directly.

    This mainly protects against bytes objects carried over from input h5ad
    files, often inside .uns. AnnData may try to write a bytes object as if it
    were a group member and raise:
      IORegistryError: No method registered for writing <class 'bytes'> into <class 'h5py._hl.group.Group'>
    """
    if isinstance(value, bytes):
        try:
            return value.decode("utf-8")
        except UnicodeDecodeError:
            return value.hex()

    if isinstance(value, np.bytes_):
        try:
            return value.astype(str)
        except Exception:
            return bytes(value).hex()

    if isinstance(value, dict):
        return {str(_make_h5ad_serializable(k)): _make_h5ad_serializable(v) for k, v in value.items()}

    if isinstance(value, tuple):
        return tuple(_make_h5ad_serializable(v) for v in value)

    if isinstance(value, list):
        return [_make_h5ad_serializable(v) for v in value]

    if isinstance(value, np.ndarray) and value.dtype == object:
        return np.asarray([_make_h5ad_serializable(v) for v in value], dtype=object)

    return value



def sanitize_anndata_for_h5ad(adata, logger=None, name="adata"):
    """
    In-place cleanup before AnnData.write_h5ad().

    AnnData/HDF5 cannot always serialize Python bytes objects, especially when
    they are nested in .uns or object-typed obs/var columns inherited from an
    input h5ad. Decode bytes to strings where possible.
    """
    # Recursive metadata slots
    for attr in ["uns"]:
        try:
            slot = getattr(adata, attr)
            cleaned = _make_h5ad_serializable(slot)
            slot.clear()
            slot.update(cleaned)
        except Exception as e:
            if logger is not None:
                logger.warning(f"Could not fully sanitize {name}.{attr}: {repr(e)}")

    # DataFrames: obs/var can contain object columns with bytes.
    for attr in ["obs", "var"]:
        try:
            df = getattr(adata, attr)
            for col in df.columns:
                if df[col].dtype == object:
                    has_bytes = df[col].map(lambda x: isinstance(x, (bytes, np.bytes_))).any()
                    if has_bytes:
                        df[col] = df[col].map(_make_h5ad_serializable)
                        if logger is not None:
                            logger.info(f"Sanitized bytes in {name}.{attr}['{col}']")
        except Exception as e:
            if logger is not None:
                logger.warning(f"Could not fully sanitize {name}.{attr}: {repr(e)}")

    return adata


# -----------------------------
# Matrix selection helpers
# -----------------------------
def get_counts_from_st(adata, counts_layer="counts"):
    if counts_layer in adata.layers:
        return ensure_sparse_float32(adata.layers[counts_layer])
    if adata.raw is not None and adata.raw.X is not None:
        return ensure_sparse_float32(adata.raw.X)
    return ensure_sparse_float32(adata.X)



def get_counts_from_scrna(adata, counts_layer="counts"):
    if counts_layer in adata.layers:
        return ensure_sparse_float32(adata.layers[counts_layer])
    if adata.raw is not None and adata.raw.X is not None:
        return ensure_sparse_float32(adata.raw.X)
    return ensure_sparse_float32(adata.X)


# -----------------------------
# Load scRNA + obs.csv
# -----------------------------
def _apply_scrna_h5ad_gene_names(adata, logger):
    if ARGS.scrna_h5ad_gene_name_source == "symbol_col":
        col = ARGS.scrna_h5ad_symbol_col
        if col not in adata.var.columns:
            raise ValueError(f"--scrna-h5ad-symbol-col {col!r} not found in scRNA h5ad var.columns")
        symbols = adata.var[col].astype(str).str.strip()
        valid = symbols.notna() & (symbols != "") & (symbols.str.lower() != "nan") & (symbols.str.lower() != "none")
        if valid.sum() == 0:
            raise ValueError(f"No valid gene symbols found in scRNA h5ad var[{col!r}]")
        if valid.sum() < adata.n_vars:
            logger.info(f"scRNA h5ad gene-symbol filter: kept {valid.sum()}/{adata.n_vars} genes")
            adata = adata[:, valid.values].copy()
            symbols = symbols.loc[valid]
        adata.var["original_var_names"] = adata.var_names.astype(str)
        adata.var_names = symbols.astype(str).values
        adata.var_names_make_unique()
    else:
        adata.var_names = adata.var_names.astype(str)
        adata.var_names_make_unique()
    return adata


def _load_10x_counts_as_adata(logger):
    logger.info("Loading scRNA 10x matrix ...")
    adata = sc.read_10x_mtx(
        ARGS.scrna_10x_dir,
        var_names=ARGS.scrna_var_names_mode,
        make_unique=True,
    )
    adata.var_names_make_unique()
    adata.X = ensure_sparse_float32(adata.X)
    adata.layers[ARGS.counts_layer] = adata.X.copy()
    logger.info(f"10x scRNA loaded: cells={adata.n_obs}, genes={adata.n_vars}")
    return adata


def _merge_obs_csv_into_scrna(adata, logger):
    if ARGS.obs_csv in [None, ""]:
        missing = [c for c in [ARGS.label_col, ARGS.batch_col] if c not in adata.obs.columns]
        if missing:
            raise ValueError(
                f"No --obs-csv provided, and scRNA h5ad obs is missing required columns: {missing}. "
                f"Either provide --obs-csv or make sure the h5ad has --label-col and --batch-col in .obs."
            )
        logger.info("No obs.csv provided; using label/batch columns already present in scRNA h5ad obs.")
        return adata

    logger.info("Loading obs.csv ...")
    obs = pd.read_csv(ARGS.obs_csv, index_col=0 if ARGS.obs_index_is_barcode else None)

    if not ARGS.obs_index_is_barcode:
        if ARGS.barcode_col not in obs.columns:
            raise ValueError(f"{ARGS.barcode_col} not found in obs.csv")
        obs = obs.set_index(ARGS.barcode_col)

    obs.index = obs.index.astype(str)
    adata.obs_names = adata.obs_names.astype(str)

    common = adata.obs_names.intersection(obs.index)
    if len(common) == 0:
        raise ValueError("No overlapping barcodes between scRNA AnnData and obs.csv")

    adata = adata[common].copy()
    obs = obs.loc[common].copy()

    for col in [ARGS.label_col, ARGS.batch_col]:
        if col not in obs.columns:
            raise ValueError(f"{col} not found in obs.csv")
        adata.obs[col] = obs[col].values

    logger.info(f"After obs.csv barcode align: cells={adata.n_obs}, genes={adata.n_vars}")
    return adata


def _supplement_h5ad_counts_from_10x(adata, logger):
    if not ARGS.supplement_10x_counts:
        return adata

    tenx = _load_10x_counts_as_adata(logger)
    adata.obs_names = adata.obs_names.astype(str)
    tenx.obs_names = tenx.obs_names.astype(str)
    adata.var_names = adata.var_names.astype(str)
    tenx.var_names = tenx.var_names.astype(str)

    common_cells = adata.obs_names.intersection(tenx.obs_names)
    common_genes = adata.var_names.intersection(tenx.var_names)
    if len(common_cells) == 0 or len(common_genes) == 0:
        raise ValueError(
            f"Cannot supplement h5ad counts from 10x: overlapping cells={len(common_cells)}, overlapping genes={len(common_genes)}"
        )

    logger.info(
        f"Supplementing scRNA h5ad counts from 10x using overlaps: cells={len(common_cells)}, genes={len(common_genes)}"
    )
    adata = adata[common_cells, common_genes].copy()
    tenx = tenx[common_cells, common_genes].copy()
    tenx = tenx[:, adata.var_names].copy()

    adata.X = ensure_sparse_float32(tenx.X)
    adata.layers[ARGS.counts_layer] = adata.X.copy()
    return adata


# -----------------------------
# Load scRNA + obs.csv
# -----------------------------
def load_scrna_reference(logger):
    if ARGS.scrna_input == "h5ad":
        logger.info(f"Loading scRNA h5ad reference input: {ARGS.scrna_h5ad}")
        adata = sc.read_h5ad(ARGS.scrna_h5ad)
        adata = _apply_scrna_h5ad_gene_names(adata, logger)
        adata = _merge_obs_csv_into_scrna(adata, logger)
        adata = _supplement_h5ad_counts_from_10x(adata, logger)
        adata.X = get_counts_from_scrna(adata, ARGS.counts_layer)
        adata.layers[ARGS.counts_layer] = get_counts_from_scrna(adata, ARGS.counts_layer)
        logger.info(f"scRNA h5ad loaded for reference training: cells={adata.n_obs}, genes={adata.n_vars}")
    else:
        adata = _load_10x_counts_as_adata(logger)
        adata = _merge_obs_csv_into_scrna(adata, logger)

    adata.obs[ARGS.label_col] = clean_categories(adata.obs[ARGS.label_col], ARGS.unlabeled)
    adata.obs[ARGS.batch_col] = pd.Categorical(adata.obs[ARGS.batch_col].astype(str))
    logger.info(f"Final scRNA reference input: cells={adata.n_obs}, genes={adata.n_vars}")
    return adata


# -----------------------------
# Gene selection: HVG + custom
# -----------------------------
def select_reference_genes(adata, logger):
    custom_genes = pd.Index([str(g) for g in ARGS.custom_genes])
    custom_set_upper = {g.upper() for g in custom_genes}
    hep_bg = read_gene_list(ARGS.hep_bg_file)
    hep_bg_upper = {str(x).upper() for x in hep_bg}

    if ARGS.use_hvg:
        logger.info("Selecting HVG on scRNA counts ...")
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=ARGS.n_top_genes,
            batch_key=ARGS.batch_col,
            layer=ARGS.counts_layer,
            flavor="seurat",
        )

        exclude = (
            adata.var_names.str.upper().str.match(r"^MT-") |
            adata.var_names.str.upper().str.match(r"^RP[SL]") |
            adata.var_names.str.upper().str.match(r"^HB(?!P)") |
            adata.var_names.str.upper().isin(hep_bg_upper)
        ) & (~adata.var_names.str.upper().isin(custom_set_upper))

        adata.var["highly_variable"] = adata.var["highly_variable"] & (~exclude)
        adata.var.loc[adata.var_names.isin(custom_genes), "highly_variable"] = True
        final_genes = adata.var_names[adata.var["highly_variable"]].copy()
        logger.info(f"Reference genes (HVG + custom): {len(final_genes)}")
    else:
        final_genes = adata.var_names.copy()
        logger.info(f"Reference genes (all genes): {len(final_genes)}")

    return final_genes


# -----------------------------
# Train reference once and persist to disk
# -----------------------------
def train_reference_models(scrna_adata, logger):
    logger.info("Preparing scRNA reference ...")
    final_genes = select_reference_genes(scrna_adata, logger)
    adata_ref = scrna_adata[:, final_genes].copy()

    adata_ref.uns = {}
    adata_ref.obsm.clear()
    adata_ref.varm.clear()
    adata_ref.obsp.clear()
    adata_ref.varp.clear()

    adata_ref.X = get_counts_from_scrna(adata_ref, ARGS.counts_layer)
    adata_ref.layers[ARGS.counts_layer] = get_counts_from_scrna(adata_ref, ARGS.counts_layer)
    adata_ref.var_names_make_unique()

    scvi_model = None
    scanvi_model = None
    condscvi_model = None

    try:
        logger.info("Training SCVI ...")
        scvi.model.SCVI.setup_anndata(
            adata_ref,
            layer=ARGS.counts_layer,
            batch_key=ARGS.batch_col,
        )
        scvi_model = scvi.model.SCVI(
            adata_ref,
            n_layers=ARGS.n_layers,
            n_latent=ARGS.n_latent,
            gene_likelihood=ARGS.gene_likelihood,
        )
        scvi_model.train(
            max_epochs=ARGS.scvi_max_epochs,
            early_stopping=True,
            accelerator=ARGS.accelerator,
            devices=ARGS.devices,
            batch_size=ARGS.scvi_batch_size,
        )
        cuda_cleanup(logger, "after SCVI")

        cond_label_col = get_cond_label_col()

        if ARGS.cond_label_source == "scanvi_pred":
            logger.info(f"Training SCANVI using labels_key={ARGS.label_col!r} ...")
            scanvi_model = scvi.model.SCANVI.from_scvi_model(
                scvi_model,
                adata=adata_ref,
                labels_key=ARGS.label_col,
                unlabeled_category=ARGS.unlabeled,
            )
            scanvi_model.train(
                max_epochs=ARGS.scanvi_max_epochs,
                early_stopping=True,
                accelerator=ARGS.accelerator,
                devices=ARGS.devices,
                batch_size=ARGS.scanvi_batch_size,
            )
            adata_ref.obs[cond_label_col] = clean_categories(
                pd.Series(scanvi_model.predict(adata_ref), index=adata_ref.obs_names),
                ARGS.unlabeled,
            )
            cuda_cleanup(logger, "after SCANVI")
        else:
            logger.info(f"Skipping SCANVI-predicted label replacement; using original labels from {ARGS.label_col!r} for CondSCVI.")
            adata_ref.obs[cond_label_col] = clean_categories(adata_ref.obs[ARGS.label_col], ARGS.unlabeled)

        logger.info(f"Training CondSCVI with labels_key={cond_label_col!r}")
        logger.info(f"CondSCVI label categories: {list(adata_ref.obs[cond_label_col].cat.categories)}")
        adata_cond = sc.AnnData(
            X=get_counts_from_scrna(adata_ref, ARGS.counts_layer),
            obs=adata_ref.obs[[cond_label_col]].copy(),
            var=adata_ref.var.copy(),
        )
        adata_cond.layers[ARGS.counts_layer] = get_counts_from_scrna(adata_ref, ARGS.counts_layer)
        adata_cond.var_names_make_unique()

        scvi.model.CondSCVI.setup_anndata(
            adata_cond,
            layer=ARGS.counts_layer,
            labels_key=cond_label_col,
        )

        condscvi_model = scvi.model.CondSCVI(
            adata_cond,
            weight_obs=False,
        )
        condscvi_model.train(
            max_epochs=ARGS.condscvi_max_epochs,
            accelerator=ARGS.accelerator,
            devices=ARGS.devices,
            batch_size=ARGS.condscvi_batch_size,
        )

        ref_dir = os.path.join(ARGS.out_dir, "reference")
        os.makedirs(ref_dir, exist_ok=True)
        ref_scanvi_h5ad = os.path.join(ref_dir, "scrna_scanvi_reference.h5ad")
        ref_cond_h5ad = os.path.join(ref_dir, "scrna_condscvi_reference.h5ad")
        ref_model_dir = os.path.join(ref_dir, "condscvi_model")

        sanitize_anndata_for_h5ad(adata_ref, logger, "adata_ref")
        sanitize_anndata_for_h5ad(adata_cond, logger, "adata_cond")
        adata_ref.write_h5ad(ref_scanvi_h5ad, compression="gzip")
        adata_cond.write_h5ad(ref_cond_h5ad, compression="gzip")
        condscvi_model.save(ref_model_dir, overwrite=True, save_anndata=True)

        logger.info("Reference training done.")
        return ref_cond_h5ad, ref_model_dir
    finally:
        release_model(condscvi_model)
        release_model(scanvi_model)
        release_model(scvi_model)
        del condscvi_model, scanvi_model, scvi_model
        try:
            del adata_ref
        except Exception:
            pass
        cuda_cleanup(logger, "after reference cleanup")


# -----------------------------
# ST preprocessing
# -----------------------------
def keep_top_gene_per_symbol(st_adata, logger):
    """Convert ST genes to symbols and keep one feature per symbol.

    For bin50-like files, use --st-gene-name-source symbol_col and --st-symbol-col real_gene_name.
    For cellbin files whose var_names are already symbols, use --st-gene-name-source var_names.
    When duplicate symbols exist, keep the feature with the largest total count.
    """
    if ARGS.st_gene_name_source == "symbol_col":
        symbol_col = ARGS.st_symbol_col
        if symbol_col not in st_adata.var.columns:
            raise ValueError(
                f"'{symbol_col}' not found in st_adata.var.columns. "
                "If this file already has gene symbols in var_names, rerun with "
                "--st-gene-name-source var_names"
            )
        symbols = pd.Series(st_adata.var[symbol_col].astype(str).str.strip().values)
        source_note = f"var['{symbol_col}']"
    else:
        symbols = pd.Series(pd.Index(st_adata.var_names).astype(str).str.strip().values)
        source_note = "var_names"

    valid = (
        symbols.notna()
        & (symbols != "")
        & (symbols.str.lower() != "nan")
        & (symbols.str.lower() != "none")
    ).to_numpy()
    st = st_adata[:, valid].copy()
    valid_symbols = symbols.loc[valid].reset_index(drop=True)

    X = get_counts_from_st(st, ARGS.counts_layer)
    colsum = np.asarray(X.sum(axis=0)).ravel()

    tmp = pd.DataFrame({
        "original_gene_id": st.var_names.astype(str),
        "symbol": valid_symbols.values,
        "colsum": colsum,
        "idx": np.arange(st.n_vars),
    }).sort_values(["symbol", "colsum"], ascending=[True, False])

    keep = tmp.drop_duplicates(subset="symbol", keep="first")["idx"].to_numpy()
    keep = np.sort(keep)

    st = st[:, keep].copy()
    kept_symbols = valid_symbols.iloc[keep].astype(str).values
    st.var["original_gene_id"] = st.var_names.astype(str)
    st.var_names = kept_symbols
    st.var_names_make_unique()

    st.layers[ARGS.counts_layer] = get_counts_from_st(st, ARGS.counts_layer)
    st.X = st.layers[ARGS.counts_layer]

    logger.info(f"ST genes after symbol cleanup from {source_note}: {st.n_vars}")
    return st

def prepare_st_sample(st_path, logger):
    logger.info(f"Loading ST sample: {st_path}")
    st = sc.read_h5ad(st_path)

    st.layers[ARGS.counts_layer] = get_counts_from_st(st, ARGS.counts_layer)
    st.X = st.layers[ARGS.counts_layer]

    keep1 = np.ones(st.n_obs, dtype=bool)
    keep2 = np.ones(st.n_obs, dtype=bool)

    if "total_counts" in st.obs.columns:
        keep1 = st.obs["total_counts"].fillna(0).values >= ARGS.min_total_counts
    if "n_genes_by_counts" in st.obs.columns:
        keep2 = st.obs["n_genes_by_counts"].fillna(0).values >= ARGS.min_genes_by_counts

    keep = keep1 & keep2
    if keep.sum() == 0:
        raise ValueError("All ST bins filtered out by QC thresholds.")
    if keep.sum() < st.n_obs:
        logger.info(f"ST QC filter: kept {keep.sum()}/{st.n_obs} bins")
        st = st[keep].copy()

    st = keep_top_gene_per_symbol(st, logger)
    return st


# -----------------------------
# DestVI per sample (worker)
# -----------------------------
def run_destvi_on_sample(st_path, ref_h5ad_path, logger):
    sample_name = sample_name_from_path(st_path)
    sample_outdir = os.path.join(ARGS.out_dir, "samples", sample_name)
    os.makedirs(sample_outdir, exist_ok=True)

    adata_cond_ref = None
    st = None
    st_sub = None
    adata_cond_sub = None
    condscvi_model_sub = None
    destvi_model = None
    props = None

    try:
        logger.info(f"Loading reference h5ad: {ref_h5ad_path}")
        adata_cond_ref = sc.read_h5ad(ref_h5ad_path)
        logger.info(f"Reference loaded: cells={adata_cond_ref.n_obs}, genes={adata_cond_ref.n_vars}")

        st = prepare_st_sample(st_path, logger)

        common_genes = adata_cond_ref.var_names.intersection(st.var_names)
        if len(common_genes) < 500:
            raise ValueError(f"Too few common genes after alignment: {len(common_genes)}")

        adata_cond_sub = adata_cond_ref[:, common_genes].copy()
        st_sub = st[:, common_genes].copy()
        st_sub = st_sub[:, adata_cond_sub.var_names].copy()

        st_sub.layers[ARGS.counts_layer] = get_counts_from_st(st_sub, ARGS.counts_layer)
        st_sub.X = st_sub.layers[ARGS.counts_layer]

        logger.info(f"Common genes for DestVI: {len(common_genes)}")
        logger.info(f"ST bins after prep: {st_sub.n_obs}, genes: {st_sub.n_vars}")
        cuda_cleanup(logger, "before subset CondSCVI")

        cond_label_col = get_reference_label_col_or_raise(adata_cond_sub, logger)
        logger.info(f"Subset CondSCVI labels_key={cond_label_col!r}")
        scvi.model.CondSCVI.setup_anndata(
            adata_cond_sub,
            layer=ARGS.counts_layer,
            labels_key=cond_label_col,
        )
        condscvi_model_sub = scvi.model.CondSCVI(
            adata_cond_sub,
            weight_obs=False,
        )
        condscvi_model_sub.train(
            max_epochs=max(50, ARGS.condscvi_max_epochs // 2),
            accelerator=ARGS.accelerator,
            devices=ARGS.devices,
            batch_size=ARGS.condscvi_batch_size,
        )
        cuda_cleanup(logger, "after subset CondSCVI")

        scvi.model.DestVI.setup_anndata(
            st_sub,
            layer=ARGS.counts_layer,
        )
        destvi_model = scvi.model.DestVI.from_rna_model(
            st_sub,
            condscvi_model_sub,
        )
        destvi_model.train(
            max_epochs=ARGS.destvi_max_epochs,
            accelerator=ARGS.accelerator,
            devices=ARGS.devices,
            batch_size=ARGS.destvi_batch_size,
        )
        cuda_cleanup(logger, "after DestVI train")

        props = destvi_model.get_proportions()
        for col in props.columns:
            st_sub.obs[f"prop_{col}"] = props[col].values

        out_h5ad = os.path.join(sample_outdir, f"{sample_name}.destvi_annotated.h5ad")
        out_csv = os.path.join(sample_outdir, f"{sample_name}.destvi_proportions.csv")
        out_model = os.path.join(sample_outdir, f"{sample_name}.destvi_model")

        sanitize_anndata_for_h5ad(st_sub, logger, "st_sub")
        st_sub.write_h5ad(out_h5ad, compression="gzip")
        props.to_csv(out_csv)
        if ARGS.save_destvi_model:
            destvi_model.save(out_model, overwrite=True, save_anndata=True)
        else:
            out_model = ""

        return {
            "sample_name": sample_name,
            "status": "success",
            "n_bins": int(st_sub.n_obs),
            "n_genes": int(st_sub.n_vars),
            "n_common_genes": int(len(common_genes)),
            "out_h5ad": out_h5ad,
            "props_csv": out_csv,
            "model_dir": out_model,
            "error": "",
        }
    finally:
        release_model(destvi_model)
        release_model(condscvi_model_sub)
        del destvi_model, condscvi_model_sub, props, st_sub, adata_cond_sub, st, adata_cond_ref
        cuda_cleanup(logger, "worker final cleanup")


# -----------------------------
# Subprocess launcher (master)
# -----------------------------
def as_cli_value(v):
    return str(v)



def build_worker_cmd(st_path, reference_h5ad, result_json):
    cmd = [
        sys.executable,
        os.path.abspath(__file__),
        "--worker-mode",
        "--worker-st-path", st_path,
        "--worker-reference-h5ad", reference_h5ad,
        "--worker-result-json", result_json,
        "--st-dir", ARGS.st_dir,
        "--obs-csv", ARGS.obs_csv,
        "--out-dir", ARGS.out_dir,
        "--scrna-input", ARGS.scrna_input,
        "--scrna-var-names-mode", ARGS.scrna_var_names_mode,
        "--scrna-h5ad-gene-name-source", ARGS.scrna_h5ad_gene_name_source,
        "--scrna-h5ad-symbol-col", ARGS.scrna_h5ad_symbol_col,
        "--label-col", ARGS.label_col,
        "--batch-col", ARGS.batch_col,
        "--st-symbol-col", ARGS.st_symbol_col,
        "--st-gene-name-source", ARGS.st_gene_name_source,
        "--st-mode", ARGS.st_mode,
        "--cond-label-source", ARGS.cond_label_source,
        "--num-threads", as_cli_value(ARGS.num_threads),
        "--num-interop-threads", as_cli_value(ARGS.num_interop_threads),
        "--n-top-genes", as_cli_value(ARGS.n_top_genes),
        "--scvi-max-epochs", as_cli_value(ARGS.scvi_max_epochs),
        "--scanvi-max-epochs", as_cli_value(ARGS.scanvi_max_epochs),
        "--condscvi-max-epochs", as_cli_value(ARGS.condscvi_max_epochs),
        "--destvi-max-epochs", as_cli_value(ARGS.destvi_max_epochs),
        "--scvi-batch-size", as_cli_value(ARGS.scvi_batch_size),
        "--scanvi-batch-size", as_cli_value(ARGS.scanvi_batch_size),
        "--condscvi-batch-size", as_cli_value(ARGS.condscvi_batch_size),
        "--destvi-batch-size", as_cli_value(ARGS.destvi_batch_size),
        "--n-latent", as_cli_value(ARGS.n_latent),
        "--n-layers", as_cli_value(ARGS.n_layers),
        "--gene-likelihood", ARGS.gene_likelihood,
        "--min-total-counts", as_cli_value(ARGS.min_total_counts),
        "--min-genes-by-counts", as_cli_value(ARGS.min_genes_by_counts),
        "--counts-layer", ARGS.counts_layer,
        "--unlabeled", ARGS.unlabeled,
        "--accelerator", ARGS.accelerator,
        "--devices", ARGS.devices,
    ]

    if ARGS.st_glob not in [None, ""]:
        cmd.extend(["--st-glob", ARGS.st_glob])
    if ARGS.obs_index_is_barcode:
        cmd.append("--obs-index-is-barcode")
    if ARGS.scrna_10x_dir is not None:
        cmd.extend(["--scrna-10x-dir", ARGS.scrna_10x_dir])
    if ARGS.scrna_h5ad is not None:
        cmd.extend(["--scrna-h5ad", ARGS.scrna_h5ad])
    if ARGS.supplement_10x_counts:
        cmd.append("--supplement-10x-counts")
    if ARGS.barcode_col is not None:
        cmd.extend(["--barcode-col", ARGS.barcode_col])
    if ARGS.cond_label_col is not None:
        cmd.extend(["--cond-label-col", ARGS.cond_label_col])
    if ARGS.use_hvg:
        cmd.append("--use-hvg")
    if ARGS.hep_bg_file is not None:
        cmd.extend(["--hep-bg-file", ARGS.hep_bg_file])
    if ARGS.custom_genes:
        cmd.extend(["--custom-genes", *ARGS.custom_genes])
    if ARGS.save_destvi_model:
        cmd.append("--save-destvi-model")

    return cmd



def write_json(path, payload):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)



def read_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


# -----------------------------
# Worker / master entrypoints
# -----------------------------
def worker_main():
    sample_name = sample_name_from_path(ARGS.worker_st_path)
    log_dir = os.path.join(ARGS.out_dir, "logs")
    os.makedirs(log_dir, exist_ok=True)
    logger = setup_logger(os.path.join(log_dir, f"{sample_name}.log"), mode="w")

    t0 = time.time()
    logger.info(f"=== Worker started for sample: {sample_name} ===")
    logger.info(f"worker_st_path = {ARGS.worker_st_path}")
    logger.info(f"worker_reference_h5ad = {ARGS.worker_reference_h5ad}")

    apply_torch_threads(ARGS.num_threads, ARGS.num_interop_threads)
    cuda_cleanup(logger, "worker startup")

    try:
        res = run_destvi_on_sample(ARGS.worker_st_path, ARGS.worker_reference_h5ad, logger)
        res["seconds"] = round(time.time() - t0, 2)
        write_json(ARGS.worker_result_json, res)
        logger.info(f"Worker finished sample {sample_name} in {res['seconds']} sec")
        return 0
    except Exception as e:
        logger.error(traceback.format_exc())
        res = {
            "sample_name": sample_name,
            "status": "failed",
            "n_bins": "",
            "n_genes": "",
            "n_common_genes": "",
            "out_h5ad": "",
            "props_csv": "",
            "model_dir": "",
            "seconds": round(time.time() - t0, 2),
            "error": repr(e),
        }
        write_json(ARGS.worker_result_json, res)
        return 1
    finally:
        cuda_cleanup(logger, "worker exit")



def prepare_reference(master_logger):
    if ARGS.retrain_reference:
        master_logger.info("Mode: retrain reference")
        scrna = load_scrna_reference(master_logger)
        try:
            ref_h5ad, ref_model_dir = train_reference_models(scrna, master_logger)
        finally:
            del scrna
            cuda_cleanup(master_logger, "after scrna cleanup")
        return ref_h5ad, ref_model_dir

    master_logger.info("Mode: reuse existing reference")
    if not os.path.exists(ARGS.reference_h5ad):
        raise FileNotFoundError(f"reference_h5ad not found: {ARGS.reference_h5ad}")
    if ARGS.reference_model_dir and (not os.path.exists(ARGS.reference_model_dir)):
        master_logger.warning(f"reference_model_dir not found (ignored in optimized loop): {ARGS.reference_model_dir}")
    return ARGS.reference_h5ad, ARGS.reference_model_dir



def master_main():
    os.makedirs(ARGS.out_dir, exist_ok=True)
    log_dir = os.path.join(ARGS.out_dir, "logs")
    os.makedirs(log_dir, exist_ok=True)

    master_logger = setup_logger(os.path.join(log_dir, "master.log"), mode="w")
    master_logger.info("=== Batch DestVI started (optimized subprocess mode) ===")
    master_logger.info("Args:")
    master_logger.info(vars(ARGS))

    apply_torch_threads(ARGS.num_threads, ARGS.num_interop_threads)
    cuda_cleanup(master_logger, "master startup")

    st_glob = resolve_st_glob()
    st_files = sorted(glob.glob(os.path.join(ARGS.st_dir, st_glob)))
    if len(st_files) == 0:
        raise FileNotFoundError(f"No files matched: {os.path.join(ARGS.st_dir, st_glob)}")
    master_logger.info(f"ST mode: {ARGS.st_mode}; file glob: {st_glob}")
    master_logger.info(f"Found {len(st_files)} ST files.")

    reference_h5ad, reference_model_dir = prepare_reference(master_logger)
    master_logger.info(f"Reference ready: {reference_h5ad}")
    if reference_model_dir:
        master_logger.info(f"Reference model dir: {reference_model_dir}")

    results = []
    for st_path in st_files:
        sample_name = sample_name_from_path(st_path)
        result_json = os.path.join(ARGS.out_dir, "samples", sample_name, f"{sample_name}.worker_result.json")
        os.makedirs(os.path.dirname(result_json), exist_ok=True)

        master_logger.info(f"=== Launch worker for sample: {sample_name} ===")
        master_logger.info(f"File: {st_path}")
        t0 = time.time()

        cmd = build_worker_cmd(st_path, reference_h5ad, result_json)
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(ARGS.num_threads)
        env["MKL_NUM_THREADS"] = str(ARGS.num_threads)
        env["OPENBLAS_NUM_THREADS"] = str(ARGS.num_threads)
        env["NUMEXPR_NUM_THREADS"] = str(ARGS.num_threads)
        env["VECLIB_MAXIMUM_THREADS"] = str(ARGS.num_threads)
        env["BLIS_NUM_THREADS"] = str(ARGS.num_threads)

        try:
            completed = subprocess.run(cmd, env=env, check=False)
            if os.path.exists(result_json):
                res = read_json(result_json)
            else:
                res = {
                    "sample_name": sample_name,
                    "status": "failed",
                    "n_bins": "",
                    "n_genes": "",
                    "n_common_genes": "",
                    "out_h5ad": "",
                    "props_csv": "",
                    "model_dir": "",
                    "seconds": round(time.time() - t0, 2),
                    "error": f"Worker exited with code {completed.returncode} and no result json.",
                }

            if completed.returncode != 0 and res.get("status") != "failed":
                res["status"] = "failed"
                res["error"] = f"Worker exited with code {completed.returncode}"

            res["seconds"] = round(time.time() - t0, 2)
            results.append(res)

            if res["status"] == "success":
                master_logger.info(f"Finished sample {sample_name} in {res['seconds']} sec")
            else:
                master_logger.error(f"Sample failed and skipped: {sample_name} | {res.get('error', '')}")
        except Exception as e:
            master_logger.error(traceback.format_exc())
            results.append({
                "sample_name": sample_name,
                "status": "failed",
                "n_bins": "",
                "n_genes": "",
                "n_common_genes": "",
                "out_h5ad": "",
                "props_csv": "",
                "model_dir": "",
                "seconds": round(time.time() - t0, 2),
                "error": repr(e),
            })
        finally:
            cuda_cleanup(master_logger, f"after sample {sample_name}")

    summary = pd.DataFrame(results)
    summary_path = os.path.join(ARGS.out_dir, "destvi_batch_summary.csv")
    summary.to_csv(summary_path, index=False)

    master_logger.info(f"Summary saved: {summary_path}")
    master_logger.info(f"Success: {(summary['status'] == 'success').sum()}")
    master_logger.info(f"Failed: {(summary['status'] == 'failed').sum()}")



def main():
    validate_args()
    if ARGS.worker_mode:
        raise SystemExit(worker_main())
    master_main()


if __name__ == "__main__":
    main()
