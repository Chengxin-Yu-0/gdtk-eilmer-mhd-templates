LMR_DIR := lmrsim

help:
	@echo "Main targets available are:"
	@echo "   prep-MODE (substitute MODE with 'nom', 'mhd-dipole', or 'mhd-femm')"
	@echo "   run"
	@echo "   post"
	@echo "   clean"

gas:
	cp ${DGD_REPO}/examples/kinetics/air-chemistry-1T/air-7sp-1T.inp .
	cp ${DGD_REPO}/examples/kinetics/air-chemistry-1T/GuptaEtAl-air-reactions.lua .
	sed -i '/^SUBSET_SELECTION =/s/5/7/' GuptaEtAl-air-reactions.lua
	lmr prep-gas -i air-7sp-1T.inp -o air-7sp-1T.gas
	lmr prep-chem -g air-7sp-1T.gas -i GuptaEtAl-air-reactions.lua -o air-7sp-9r-1T.chem

grid:
	lmr prep-grid -j grid.lua

sim-nom:
	lmr prep-sim -j job-nom.lua

sim-dipole:
	lmr prep-sim -j job-dipole.lua

sim-femm:
	lmr prep-sim -j job-femm.lua

mhd:
	./mhd_prep

prep-nom: gas grid sim-nom

prep-mhd-dipole: gas grid sim-dipole

prep-mhd-femm: gas grid sim-femm mhd

run:
	lmr run

post:
	lmr snapshot2vtk --all
	paraview lmrsim/vtk/fluid.pvd
	
clean:
	-rm air-7sp-1T.gas
	-rm air-7sp-1T.inp
	-rm air-7sp-9r-1T.chem
	-rm GuptaEtAl-air-reactions.lua
	-rm prepgridlog.dat
	-rm prepsimlog.dat
	-rm field_table.lua
	-rm debugging_output.txt
	-rm -r lmrsim

