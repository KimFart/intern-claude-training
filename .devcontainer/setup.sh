#!/bin/bash
set -e

# Register Jupyter kernel (sbml env already exists in the image)
/opt/conda/envs/sbml/bin/python -m ipykernel install --user \
  --name sbml --display-name "Python (sbml)"

echo "==> Setup complete."
