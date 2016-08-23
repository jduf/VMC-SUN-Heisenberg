DIRECTORY=$(shell find . -maxdepth 1 -type d -not -name "doxygen" -not -name ".*")

all:$(DIRECTORY)

$(DIRECTORY):
	@echo =========  Compile $@ =========  
	@$(MAKE) -j 12 -i -s -C $@
	@$(MAKE) -s -i -C $@ clean

.PHONY: $(DIRECTORY)

ref:
	@echo Create the documentation
	@doxygen doxygen/Doxyfile > doxygen/log 2> doxygen/err
	@firefox doxygen/html/annotated.html &

