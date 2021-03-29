
NET := $(wildcard *.net)
LOG := $(wildcard *.log)
OBJ := $(addsuffix .tex, $(basename $(NET)))

all: ngspice.log
	$(eval NGSFIGS=$(shell grep _FIG $< | sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGPS=$(addsuffix .ps, $(NGSFIGS)))
	$(eval NGSTABS=$(shell grep _TAB $< | sed 's/_TAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(NGSTABS), sed -n '/^$i_TAB/,/^$i_END/{p;/^$i_END/q}' $< | grep -v $i_TAB | grep -v $i_END | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab.tex;)
	$(foreach i, $(NGSFIGPS), ps2pdf $i;)
	touch $(OBJ)
	cp $(NET) $(OBJ)
	$(eval OBJTABS=$(shell grep _TAB $(OBJ) | sed 's/_TAB/\_TAB/g' | sed ':a;N;$!ba;s/\n/\par/g'))

ngspice.log:
	$(foreach var, $(NET), ngspice -b $(var) -o ngspice.log)

clean:
	@rm -f $(LOG) *.tex *.pdf *.ps *~

.PHONY: all clean
