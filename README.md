# MRIRecoTutorial
Tutorial on the Julia Package MRIReco.jl

> **Warning**
> This tutorial is outdated and uses an old version of Julia and MRIReco.jl. We have transfered all content to the [MRIReco.jl documentation](https://magneticresonanceimaging.github.io/MRIReco.jl/latest/).

## Getting Started

This tutorial can be used and run in two ways:
* You can directly use Pluto notebook you are reading right now. This needs some preparation discussed next.
* You can also copy/past into a Julia REPL. This is in particular a good way if you are reading an HTML or PDF version of this tutorial.
* Keep in mind that this notebook will run all cells, which needs some time since we also showcase some more computationally intensive methods.

## Preparation

* Go to [https://julialang.org/downloads/](https://julialang.org/downloads/) and download the most recent stable release (version 1.6 when this notebook was prepared).
* Checkout [https://github.com/MagneticResonanceImaging/MRIRecoTutorial](https://github.com/MagneticResonanceImaging/MRIRecoTutorial).
* Start Julia from within a terminal within the folder you checked out the repository.
* We need to install some packages. With ] we can enter the Pkg mode and enter
  ```julia
  (@v1.6) pkg> activate .
  ```
  which will install all necessary packages.
* To run this notebook use the command `import Pluto; Pluto.run()` and open the notebook in the appearing file browser.
* In the box *Open from file:* enter the filename of the notebook (MRIRecoTutorial.jl) and press enter.
* We are using [Pluto.jl](https://github.com/fonsp/Pluto.jl), which is similar to Jupyter but has some advantages, in particular the notebook is stored as plain Julia code that you can access with a simple text editor.
* Instead of running the notebook in the browser you can also run it in the REPL
  ```julia
  julia> include("MRIRecoTutorial.jl")
  ```
  Note that this, however will only show the plot of the last reconstruction since this notebook is tailored for usage in Pluto.
