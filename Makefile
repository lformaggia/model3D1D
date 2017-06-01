

DEFAULT_GOAL= all

.PHONY= all clean distclean library doc pdf MANwork

all:
	make -C MixedFormulation
	make -C CurvedFormulation
	make -C GraphGenerator

MANwork:
	make -C MixedFormulation
	make -C CurvedFormulation

library:
	make library -C MixedFormulation
	make library -C CurvedFormulation
	make library -C GraphGenerator

doc:
	make doc -C MixedFormulation
	make doc -C CurvedFormulation
	make doc -C GraphGenerator

pdf:
	make pdf -C MixedFormulation
	make pdf -C CurvedFormulation
	make pdf -C GraphGenerator

clean: 
	make clean -C MixedFormulation
	make clean -C CurvedFormulation
	make clean -C GraphGenerator

distclean: 
	make distclean -C MixedFormulation
	make distclean -C CurvedFormulation
	make distclean -C GraphGenerator
