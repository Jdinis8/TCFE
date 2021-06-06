%.exe:
	@echo Olá
	make -C $(basename $@)

clean:
	make -C t0 clean
	make -C t1 clean
	make -C t2 clean
	make -C t3 clean
	make -C t4 clean
	make -C t5 clean