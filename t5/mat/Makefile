FILES := $(wildcard *.m)
TRASH := $(wildcard *.log)
EPS := $(wildcard *.eps)

log:
	$(foreach var, $(FILES), octave $(var) > octave_$(basename $(var)).log;)

clean:
	@rm -f $(TRASH) octave-workspace *.eps *~ *.png

.PHONY: clean
