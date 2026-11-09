#!/bin/bash
# Python environment on the Mac (models, figures, paper builders).
# Usage:  bash ambiente/setup_mac.sh      then   source .venv/bin/activate
set -eu
cd "$(dirname "$0")/.."
python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -r ambiente/requirements.txt
.venv/bin/python -c "import numpy, scipy, matplotlib, docx, lxml, latex2mathml; print('ambiente ok')"
