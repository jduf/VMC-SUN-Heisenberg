all:
	make -f */makefile

clean:
	rm */*.o */*.gch */test

ref:
	doxygen Doxyfile
	firefox html/files.html
