# README

Install [pandoc](https://pandoc.org/) and [pandoc-crossref](http://lierdakil.github.io/pandoc-crossref/).
Compile the blog post to HTML with pandoc.

```bash
pandoc \
    --filter pandoc-crossref \
    --filter pandoc-citeproc \
    --bibliography=references.bib \
    -s predicting-seasonal-influenza-evolution.md \
    -o predicting-seasonal-influenza-evolution.html
```
