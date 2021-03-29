# ---------------------------------------------------------
# type the "make" command to create report.pdf 
# ---------------------------------------------------------

#doc name
FILE=report
EPS  := $(wildcard ../mat/*.eps)

#list of LibreOffice Draw files
LODFIGS = $(patsubst %.odg,%.pdf,$(wildcard *.odg)) $(patsubst %.odp,%.pdf,$(wildcard *.odp))

all: $(FILE).pdf
	evince $< &

$(FILE).pdf: $(wildcard *.tex) $(LODFIGS)
ifeq ($(EPS),)
	pdflatex  $(FILE).tex
	pdflatex  $(FILE).tex
	pdflatex  $(FILE).tex
else
	cp -f $(EPS) .
	pdflatex  $(FILE).tex
	pdflatex  $(FILE).tex
	pdflatex  $(FILE).tex
endif

#convert libreoffice draw figures to pdf
%.pdf: %.odg
	soffice --headless --convert-to pdf $<

%.pdf: %.odp
	soffice --headless --convert-to pdf $<

clean:
	@rm -rf *.aux *.bbl *.blg *.glg *.glo *.gls *.ilg *.ist *.lof
	@rm -rf *.log *.lot *.nlo *.nls *.out *.toc *~ *.*% ./*.pdf ./*.eps

.PHONY: all clean


