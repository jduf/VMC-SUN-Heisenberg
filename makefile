clean:
	rm */*.o */*.gch

ref:
	doxygen Doxyfile
	firefox html/files.html
