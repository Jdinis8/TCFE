%.exe:
	@echo Olá
	make -C $(basename $@)

clean:
	make -C t0 clean
	make -C t1 clean