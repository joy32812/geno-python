# geno-python

## Project layout

```
geno-python/
├─ src/
│  └─ geno.py
└─ tests/
   └─ test_geno.py
```

## Requirements

- macOS/Linux
- Python 3.9+
- `pytest` (for tests)

## Setup (optional but recommended)

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install pytest
```

## Run the program

```bash
python src/geno.py
```

> If `geno.py` expects input files (e.g., under `data/`), make sure they exist before running.

## Run the tests

```bash
PYTHONPATH=. pytest -q
```

---

