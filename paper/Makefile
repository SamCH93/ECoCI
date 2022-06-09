all: pdf clean

# names of .Rnw .bib. and .bst files
FILE=ecoci
BIB=bibliography
BIBSTYLE=apalikedoiurl

## generate .tex file from .Rnw file with knitr
tex: $(FILE).Rnw $(BIB).bib $(BIBSTYLE).bst
	Rscript -e "knitr::knit('$(FILE).Rnw')" --vanilla

## generate .tex file from .Rnw file with knitr using journal template
texjournal: $(FILE)Journal.Rnw
	Rscript -e "knitr::knit('$(FILE)Journal.Rnw')" --vanilla

## generate .pdf file from .tex file
pdf: tex
	pdflatex $(FILE)
	bibtex $(FILE)
	pdflatex $(FILE)
	pdflatex $(FILE)

## generate .pdf file from .tex file using journal template
pdfjournal: texjournal
	pdflatex $(FILE)Journal
	bibtex $(FILE)Journal
	pdflatex $(FILE)Journal
	pdflatex $(FILE)Journal

## generate .pdf file from .tex (without first running the tex recipe)
pdf2:
	pdflatex $(FILE)
	bibtex $(FILE)
	pdflatex $(FILE)
	pdflatex $(FILE)

## generate docx file with pandoc
## FIXME: references do not work
docx: tex
## create copy of tex file as input for pandoc
	cp ./$(FILE).tex ./$(FILE)Pandoc.tex		
## use sed to add .pdf extension to file-names in tex file (ugly hack)
	figures=`ls figure/*.pdf | sed -e 's/\.pdf$///'`; \
	for f in $$figures; do sed -i "s|$$f|$$f\.pdf|g" ./$(FILE)Pandoc.tex; done	
## use pandoc to generate .docx from .tex file
	pandoc $(FILE)Pandoc.tex --bibliography=$(BIB).bib --natbib -o $(FILE).docx

## clean
clean:  
	rm $(FILE).aux  $(FILE).blg  $(FILE).log  $(FILE).tex  $(FILE).bbl \
	$(FILE).out  $(FILE).brf  $(FILE).synctex.gz \
	$(FILE)Journal.aux  $(FILE)Journal.blg  $(FILE)Journal.log  $(FILE)Journal.tex \
	$(FILE)Journal.bbl  $(FILE)Journal.out  $(FILE)Journal.brf \
	$(FILE)Journal.synctex.gz; \
	rm -r ./figure/
