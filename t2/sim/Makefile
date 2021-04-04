
NET := $(wildcard *.net)
PRG := $(patsubst %.net,%,$(NET))
LOG := $(wildcard *.log)
OBJ := $(addsuffix .tex, $(basename $(NET)))

all: ze1 ze2 ze3 ze4 ze5 ze6 ze7 ze8

ze1: TCFEL2_sim.log
	$(eval NGSFIGS=$(shell grep _FIG $(basename $<).log | sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGPS=$(addsuffix .ps, $(NGSFIGS)))
	$(eval NGSTABS=$(shell grep _TAB $(basename $<).log | sed 's/_TAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(NGSTABS), sed -n '/^$i_TAB/,/^$i_END/{p;/^$i_END/q}' $(basename $<).log | grep -v $i_TAB | grep -v $i_END | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab.tex;)
	$(foreach i, $(NGSFIGPS), ps2pdf $i;)
	touch $(basename $<).tex
	cp $(basename $<).log $(basename $<).tex
	$(eval OBJTABS=$(shell grep _TAB $(basename $<).tex | sed 's/_TAB/\_TAB/g' | sed ':a;N;$!ba;s/\n/\par/g'))

ze2: TCFELab22.log
	$(eval NGSFIGS=$(shell grep _FIG $(basename $<).log | sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGPS=$(addsuffix .ps, $(NGSFIGS)))
	$(eval NGSTABS=$(shell grep _TAB $(basename $<).log | sed 's/_TAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(NGSTABS), sed -n '/^$i_TAB/,/^$i_END/{p;/^$i_END/q}' $(basename $<).log | grep -v $i_TAB | grep -v $i_END | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab.tex;)
	$(foreach i, $(NGSFIGPS), ps2pdf $i;)
	touch $(basename $<).tex
	cp $(basename $<).log $(basename $<).tex
	$(eval OBJTABS=$(shell grep _TAB $(basename $<).tex | sed 's/_TAB/\_TAB/g' | sed ':a;N;$!ba;s/\n/\par/g'))

ze3: TCFELab23.log
	$(eval NGSFIGS=$(shell grep _FIG $(basename $<).log | sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGPS=$(addsuffix .ps, $(NGSFIGS)))
	$(eval NGSTABS=$(shell grep _TAB $(basename $<).log | sed 's/_TAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(NGSTABS), sed -n '/^$i_TAB/,/^$i_END/{p;/^$i_END/q}' $(basename $<).log | grep -v $i_TAB | grep -v $i_END | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab.tex;)
	$(foreach i, $(NGSFIGPS), ps2pdf $i;)
	touch $(basename $<).tex
	cp $(basename $<).log $(basename $<).tex
	$(eval OBJTABS=$(shell grep _TAB $(basename $<).tex | sed 's/_TAB/\_TAB/g' | sed ':a;N;$!ba;s/\n/\par/g'))

ze4: TCFELab23.log
	$(eval NGSFIGS=$(shell grep _FIG $(basename $<).log | sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGPS=$(addsuffix .ps, $(NGSFIGS)))
	$(eval NGSTABS=$(shell grep _TAB $(basename $<).log | sed 's/_TAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(NGSTABS), sed -n '/^$i_TAB/,/^$i_END/{p;/^$i_END/q}' $(basename $<).log | grep -v $i_TAB | grep -v $i_END | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab.tex;)
	$(foreach i, $(NGSFIGPS), ps2pdf $i;)
	touch $(basename $<).tex
	cp $(basename $<).log $(basename $<).tex
	$(eval OBJTABS=$(shell grep _TAB $(basename $<).tex | sed 's/_TAB/\_TAB/g' | sed ':a;N;$!ba;s/\n/\par/g'))

ze5: TCFELab24.log
	$(eval NGSFIGS=$(shell grep _FIG $(basename $<).log | sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGPS=$(addsuffix .ps, $(NGSFIGS)))
	$(eval NGSTABS=$(shell grep _TAB $(basename $<).log | sed 's/_TAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(NGSTABS), sed -n '/^$i_TAB/,/^$i_END/{p;/^$i_END/q}' $(basename $<).log | grep -v $i_TAB | grep -v $i_END | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab.tex;)
	$(foreach i, $(NGSFIGPS), ps2pdf $i;)
	touch $(basename $<).tex
	cp $(basename $<).log $(basename $<).tex
	$(eval OBJTABS=$(shell grep _TAB $(basename $<).tex | sed 's/_TAB/\_TAB/g' | sed ':a;N;$!ba;s/\n/\par/g'))

ze6: TCFELab25.log
	$(eval NGSFIGS=$(shell grep _FIG $(basename $<).log | sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGPS=$(addsuffix .ps, $(NGSFIGS)))
	$(eval NGSTABS=$(shell grep _TAB $(basename $<).log | sed 's/_TAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(NGSTABS), sed -n '/^$i_TAB/,/^$i_END/{p;/^$i_END/q}' $(basename $<).log | grep -v $i_TAB | grep -v $i_END | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab.tex;)
	$(foreach i, $(NGSFIGPS), ps2pdf $i;)
	touch $(basename $<).tex
	cp $(basename $<).log $(basename $<).tex
	$(eval OBJTABS=$(shell grep _TAB $(basename $<).tex | sed 's/_TAB/\_TAB/g' | sed ':a;N;$!ba;s/\n/\par/g'))

ze7: TCFELab26.log
	$(eval NGSFIGS=$(shell grep _FIG $(basename $<).log | sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGPS=$(addsuffix .ps, $(NGSFIGS)))
	$(eval NGSTABS=$(shell grep _TAB $(basename $<).log | sed 's/_TAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(NGSTABS), sed -n '/^$i_TAB/,/^$i_END/{p;/^$i_END/q}' $(basename $<).log | grep -v $i_TAB | grep -v $i_END | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab.tex;)
	$(foreach i, $(NGSFIGPS), ps2pdf $i;)
	touch $(basename $<).tex
	cp $(basename $<).log $(basename $<).tex
	$(eval OBJTABS=$(shell grep _TAB $(basename $<).tex | sed 's/_TAB/\_TAB/g' | sed ':a;N;$!ba;s/\n/\par/g'))

ze8: TCFELab27.log
	$(eval NGSFIGS=$(shell grep _FIG $(basename $<).log | sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGPS=$(addsuffix .ps, $(NGSFIGS)))
	$(eval NGSTABS=$(shell grep _TAB $(basename $<).log | sed 's/_TAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(NGSTABS), sed -n '/^$i_TAB/,/^$i_END/{p;/^$i_END/q}' $(basename $<).log | grep -v $i_TAB | grep -v $i_END | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab.tex;)
	$(foreach i, $(NGSFIGPS), ps2pdf $i;)
	touch $(basename $<).tex
	cp $(basename $<).log $(basename $<).tex
	$(eval OBJTABS=$(shell grep _TAB $(basename $<).tex | sed 's/_TAB/\_TAB/g' | sed ':a;N;$!ba;s/\n/\par/g'))


TCFEL2_sim.log: ../sim/TCFEL2_sim.net
	ngspice -b $< -o $@

TCFELab22.log: ../sim/TCFELab22.net
	ngspice -b $< -o $@

TCFELab23.log: ../sim/TCFELab23.net
	ngspice -b $< -o $@

TCFELab24.log: ../sim/TCFELab24.net
	ngspice -b $< -o $@

TCFELab25.log: ../sim/TCFELab25.net
	ngspice -b $< -o $@

TCFELab26.log: ../sim/TCFELab26.net
	ngspice -b $< -o $@

TCFELab27.log: ../sim/TCFELab27.net
	ngspice -b $< -o $@

clean:
	@rm -f $(LOG) *.tex *.pdf *.ps *~

.PHONY: all clean
