TITLE=manuscript

rm -f *.log
rm -f *.aux
rm -f *.out
rm -f *.bbl
rm -f *.bcf
rm -f *.blg
rm -f .DS_Store
rm -f ${TITLE}.pdf

pdflatex ${TITLE}.tex
bibtex ${TITLE}
pdflatex ${TITLE}.tex
pdflatex ${TITLE}.tex
pdflatex ${TITLE}.tex
