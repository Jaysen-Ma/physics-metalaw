# physics-metalaw

The HTML ledger in `data/catalog/` and the numbers on its contributions page.
Operator counts and the exponential fit live in `physics_metalaw/`.

## Serve the ledger

```bash
.venv/bin/python -m physics_metalaw.catalog --serve
```

Then open `http://127.0.0.1:8765/`.

## Rebuild

Feynman, Wikipedia, and Encyclopaedia Inflationaris are seeded from the
released CSVs cloned at `.scratch/physics-formulae-patterns` (commit `2aa199b`).
That clone is not part of this repository.

```bash
.venv/bin/python -m physics_metalaw.catalog --seed
.venv/bin/python scripts/build_leemis_catalog.py
.venv/bin/python scripts/build_fr_catalog.py
.venv/bin/python scripts/make_plots.py
```

`scripts/make_plots.py` redraws every figure the contributions page embeds.
The sample of 10 Encyclopaedia Inflationaris expressions uses seed `20260930`
(`physics_metalaw.sampling.SAMPLE_SEED`). The other fits are deterministic.
