NET := $(wildcard *.net)
LOG := $(wildcard *.log)
OBJ := $(addsuffix .tex, $(basename $(NET)))

all: ngspice.log
	$(foreach var, $(wildcard *.net), make $(basename $(var)).executavel;)


%.executavel:
	$(eval NGSFIGS=$(shell grep _FIG $(basename $@).log | sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGPS=$(addsuffix .ps, $(NGSFIGS)))
	$(eval NGSTABS=$(shell grep _TAB $(basename $@).log | sed 's/_TAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(NGSTABS), sed -n '/^$i_TAB/,/^$i_END/{p;/^$i_END/q}' $(basename $@).log | grep -v $i_TAB | grep -v $i_END | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab.tex;)
	$(foreach i, $(NGSFIGPS), ps2pdf $i;)
	touch $(basename $@).tex
	cp $(basename $@).log $(basename $@).tex
	$(eval OBJTABS=$(shell grep _TAB $(basename $@).tex | sed 's/_TAB/\_TAB/g' | sed ':a;N;$!ba;s/\n/\par/g'))


ngspice.log:
	@echo olaaaaaaaaaaaa
	$(foreach var, $(NET), ngspice -b $(var) -o $(basename $(var)).log;)

clean:
	@rm -f $(LOG) *.tex *.pdf *.ps *~

.PHONY: all clean