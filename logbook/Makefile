logbookentries := $(shell grep 'entries' logbook.tex | grep -v '^%' | grep -o 'entries/[A-Za-z0-9_]*\.tex')
img_inputs := $(shell grep '\\includegraphics' $(logbookentries) | grep -v '%.*\\includegraphics' | sed 's/^.*\\includegraphics\[.*\]{\([\./a-zA-Z0-9_]*\)}.*/\1/')
tex_inputs := $(shell grep --no-filename '\\input' $(logbookentries) | grep -v '^%' | grep -o 'input{[A-Za-z0-9_\/\.]*}'  | sed 's/input{//' | sed 's/}$$//')


all: logbook.pdf

clean:
	rm *.log *.out *.toc *.aux

logbook.bbl: ../bib/notes.bib | logbook.tex
	pdflatex -draftmode logbook.tex
	rm logbook.log logbook.out logbook.toc
	bibtex logbook.aux
	rm logbook.blg

logbook.pdf: logbook.tex $(logbookentries) $(img_inputs) $(tex_inputs) logbook.bbl
	pdflatex -draftmode -shell-escape $<
	pdflatex $<
	rm $(addprefix $(basename $<),.log .out .aux .toc)

input/:
	mkdir $@
