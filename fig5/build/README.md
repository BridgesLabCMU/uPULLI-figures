# build/: regenerate the embedding-UMAP coordinate tables

Two generators produce the 5B-left and 5C-left coordinate tables in `../data/` from the DINOv2 CLS
caches (KiltHub). The other Figure-5 source tables (5A barcode + projection, 5B/5C confusion matrices)
are produced by their datasets' analysis pipelines and provided as source data (see `../README.md`).

## Inputs

```bash
python fetch_data.py fig5        # from the repository root
```

Logical names resolved via `config.input(...)` — see [`../../INPUTS.md`](../../INPUTS.md):

- `kleb/embeddings/cls.npy`, `kleb/embeddings/index.csv`, `kleb/master_frame_features.csv` (growth filter).
- `multispecies/embeddings/cls_10X.npy`, `multispecies/embeddings/index_10X.csv`.

Already have the deposit unpacked elsewhere? `export UPULLI_DATA_ROOT=/path/to/deposit` instead.

## Run

```bash
python build/B_left_klebUmap.py            # -> data/kleb_embeddingUmap_coords.csv (cosine + euclid)
python build/C_left_multispeciesUmap.py    # -> data/multispecies_100pctLB_10X_umap_coords.csv
```

UMAP output is library-version sensitive, so the figures are pinned to the tables in `../data/`; use the
pinned `environment.yml`.
