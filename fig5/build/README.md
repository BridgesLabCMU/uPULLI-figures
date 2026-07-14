# build/: regenerate the embedding-UMAP coordinate tables

Two generators produce the 5B-left and 5C-left coordinate tables in `../data/` from the DINOv2 CLS
caches (KiltHub). The other Figure-5 source tables (5A barcode + projection, 5B/5C confusion matrices)
are produced by their datasets' analysis pipelines and provided as source data (see `../README.md`).

## Inputs (KiltHub; override via env)

- `kleb_cls.npy` + `kleb_embIndex.csv` and the kleb `master_frame_features.csv` (growth filter) —
  `FIG5_KLEB_CLS`, `FIG5_KLEB_EMBIDX`, `FIG5_KLEB_FRAME`.
- `multispecies_10X_cls.npy` + `multispecies_10X_embIndex.csv` — `FIG5_MULTI_CLS`, `FIG5_MULTI_EMBIDX`.

## Run

```bash
python build/B_left_klebUmap.py            # -> data/kleb_embeddingUmap_coords.csv (cosine + euclid)
python build/C_left_multispeciesUmap.py    # -> data/multispecies_100pctLB_10X_umap_coords.csv
```

UMAP output is library-version sensitive, so the figures are pinned to the tables in `../data/`; use the
pinned `environment.yml`.
