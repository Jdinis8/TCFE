FILES := $(wildcard *.m)
TRASH := $(wildcard *.log)

log:
	$(foreach var, $(FILES), octave $(var) > octave_$(basename $(var)).log)

clean:
	@rm -f $(TRASH) octave-workspace *.eps *~

.PHONY: clean
