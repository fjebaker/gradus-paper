# Gradus.jl: spacetime-agnostic general relativistic ray-tracing for X-ray spectral modelling

We introduce Gradus.jl, an open-source and publicly available general
relativistic ray-tracing toolkit for spectral modelling in arbitrary
spacetimes. Our software is written in the Julia programming language, making
use of forward-mode automatic differentiation for computing the Christoffel
symbols during geodesic integration, and for propagating derivatives through
the entire ray-tracer. Relevant numerical methods are detailed, and our models
are validated using a number of tests and comparisons to other codes. The
differentiability is used to optimally calculate Cunningham transfer functions
-- used to efficiently pre-compute relativistic effects in spectral models. A
method is described for calculating such transfer functions for disc with
non-zero vertical height, including the treatment of self-obscuration. An
extension of the transfer function formalism that includes timing information
is described, and used to calculate high-resolution reverberation lag spectra
for a lamppost corona. The lag-frequency and lag-energy spectra for a
Shakura-Sunyaev accretion disc with various lamppost heights and Eddington
ratios are calculated, and the general impact of disc thickness in reflection
models is discussed.

## Building the paper

    make

To build the diffs from before the reviewer comments


    make diff

To bundle for distribution and final submission


    make bundle

## $\LaTeX$ document

This uses the [MNRAS style files](https://www.ctan.org/tex-archive/macros/latex/contrib/mnras). Further
information about preparing MNRAS papers in $\LaTeX$ can be found in the
[Instructions to Authors](https://academic.oup.com/mnras/pages/general_instructions).
