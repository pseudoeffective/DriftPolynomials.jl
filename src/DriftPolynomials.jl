################################################################################
# DriftPolynomials.jl
#
# An algebraic combinatorics package for Julia.
#
# Copyright (C) 2025 Dave Anderson, pseudoeffective.github.io
################################################################################


module DriftPolynomials

################################################################################
# Import
################################################################################


# AbstractAlgebra
import AbstractAlgebra:
	base_ring, gen, gens, parent_type, nvars, polynomial_ring, MPolyBuildCtx, 
	
	push_term!, finish

# Nemo
import Nemo:
	ZZ, QQ, libflint, Ring, MPolyRing, MPolyRingElem, ZZMPolyRing, ZZMPolyRingElem, 
	
	QQMPolyRing, QQMPolyRingElem, evaluate, vars, coefficients


# BumplessPipeDreams
import BumplessPipeDreams:
	 BPD, all_bpds, is_asm, draw_bpd, _draw_bpd_plots

# Printf 
import Printf:
	@sprintf


################################################################################
# Export
################################################################################

export
	ZZ, QQ, PolyRing, ZZMPolyRing, ZZMPolyRingElem, QQMPolyRing, QQMPolyRingElem,
	
	base_ring, polynomial_ring, gen, gens, nvars, vars, coefficients, evaluate, 

	BPD, Drift, is_asm, drift_class, nw_reset, se_reset, 
	
	can_drift, can_undrift, drift, undrift,
	
	partition2drift, bpd2drift, drift2rkmtx, rkmtx2asm, 

	isintegrable, random_drift, empty_drift, 

	drift_poly, markconfig, xy_ring,

	draw_drift, draw_bpd


################################################################################
# source files
################################################################################

include("drifts.jl")
include("xy_ring.jl")
include("drift_polys.jl") # must be included after `drifts.jl`
include("split.jl")
include("draw_drift.jl")


end # module DriftPolynomials
