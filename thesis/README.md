# Masterarbeit thesis seed

This is a starter LaTeX thesis structure for the Laplacian-eigenmode static-source and gluonic flux-tube project.

## Build

```bash
cd thesis
latexmk -pdf -interaction=nonstopmode main.tex
```

or:

```bash
make
```

## Layout

- `main.tex` is the root file.
- `tex/preamble.tex` contains packages and macros.
- `frontmatter/` contains title page, abstract, acknowledgements, declaration, and notation.
- `chapters/` contains the thesis body.
- `appendices/` contains reproducibility notes and extra plots.
- `figures/` contains the BUW logo and preliminary project plots.
- `bib/references.bib` contains the starting bibliography.

## Notes

The current text is intentionally a working draft. It contains validated static-potential benchmark material and placeholders for the final local gluonic probe / flux-tube profile results.
