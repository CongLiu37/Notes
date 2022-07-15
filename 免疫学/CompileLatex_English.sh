pdflatex $1.tex
pdflatex $1.tex
bibtex $1
pdflatex $1.tex
pdflatex $1.tex

rm $1.aux
rm $1.log
rm $1.out
rm $1.bbl
rm $1.blg
rm $1.synctex.gz
rm $1.fdb_latexmk
rm $1.fls
rm $1.nav
rm $1.snm
rm $1.toc
