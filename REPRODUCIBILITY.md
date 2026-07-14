# Reproducibility

## Requirements

- Python 3.10 or newer
- `mpmath`
- `pdflatex` to rebuild the paper

Run the full numerical audit:

```bash
python scripts/verify_scale_optimized_polygonal_pi.py --precision 1300 --output-dir output
```

The script verifies:

- unity and moment conditions for the scale weights;
- the corrected stability constants for `b = 2, 3, 6, 9, 27`;
- slopes through order `r = 24`;
- the explicit analytic remainder bound;
- the exact dyadic chord recurrence;
- the `O(N^-3)` linear tripling seed asymptotic.

Rebuild the paper:

```bash
pdflatex -interaction=nonstopmode -halt-on-error -output-directory paper paper/scale_optimized_polygonal_pi_acceleration.tex
pdflatex -interaction=nonstopmode -halt-on-error -output-directory paper paper/scale_optimized_polygonal_pi_acceleration.tex
```
