
all:
	$(MAKE) -s -C arrays/
	$(MAKE) -s -C container/
	$(MAKE) -s -C diffdir/
	$(MAKE) -s -C directory/
	$(MAKE) -s -C factorial/
	$(MAKE) -s -C fit/
	$(MAKE) -s -C gnuplot/
	$(MAKE) -s -C header/
	$(MAKE) -s -C jdtools/
	$(MAKE) -s -C linux/
	$(MAKE) -s -C list/
	$(MAKE) -s -C matrix-lapack/
	$(MAKE) -s -C parseur/
	$(MAKE) -s -C prime/
	$(MAKE) -s -C pso/
	$(MAKE) -s -C pstricks2pdf/
	$(MAKE) -s -C rand/
	$(MAKE) -s -C read-write/
	$(MAKE) -s -C rewrite/
	$(MAKE) -s -C rst/
	$(MAKE) -s -C run/
	$(MAKE) -s -C time/
	$(MAKE) -s -C young/

clean:$(SUBDIRS)
	$(MAKE) -s -C arrays clean
	$(MAKE) -s -C container clean
	$(MAKE) -s -C diffdir clean
	$(MAKE) -s -C directory clean
	$(MAKE) -s -C factorial clean
	$(MAKE) -s -C fit clean
	$(MAKE) -s -C gnuplot clean
	$(MAKE) -s -C header clean
	$(MAKE) -s -C jdtools clean
	$(MAKE) -s -C linux clean
	$(MAKE) -s -C list clean
	$(MAKE) -s -C matrix-lapack clean
	$(MAKE) -s -C parseur clean
	$(MAKE) -s -C prime clean
	$(MAKE) -s -C pso clean
	$(MAKE) -s -C pstricks2pdf clean
	$(MAKE) -s -C rand clean
	$(MAKE) -s -C read-write clean
	$(MAKE) -s -C rewrite clean
	$(MAKE) -s -C rst clean
	$(MAKE) -s -C run clean
	$(MAKE) -s -C time clean
	$(MAKE) -s -C young clean

ref:
	@echo Create the documentation
	@doxygen doxygen/Doxyfile
	@firefox doxygen/html/annotated.html &

