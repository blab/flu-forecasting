# README

Compile the blog post to HTML with pandoc.

```bash
pandoc \
    --filter pandoc-citeproc \
    --bibliography=references.bib \
    -s predicting-seasonal-influenza-evolution.md \
    -o predicting-seasonal-influenza-evolution.html
```
