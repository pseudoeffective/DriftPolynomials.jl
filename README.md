# DriftPolynomials

*A Julia package for computing drift polynomials*



This package provides functions for computing *drift polynomials*, which are certain generalizations of [Schur polynomials](https://en.wikipedia.org/wiki/Schur_polynomial) studied in [Anderson-Fulton]() *to appear*.

To use the package, first install Julia.  If you are new to this language, I recommend Ulrich Thiel's page on [JuLie](https://ulthiel.github.io/JuLie.jl/dev/), which in addition to providing software, has a nice introduction for the algebraically-minded mathematician.

To install the DriftPolynomials package do this:
```julia-repl
julia> using Pkg
julia> Pkg.add(url="https://github.com/pseudoeffective/DriftPolynomials.jl")
```

To start using the package, type
```julia-repl
julia> using DriftPolynomials
```

The basic objects are *drift configurations* (of type `Drift`).


