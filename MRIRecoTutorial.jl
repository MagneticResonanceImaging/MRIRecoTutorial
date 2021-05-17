### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ 96d4f44e-49d1-4398-b7f8-4daa7290b287
begin
	using PlutoUI
	PlutoUI.TableOfContents()
end

# ╔═╡ 915d5707-daf4-4127-abf8-68a5528d9209
begin
	using ImageUtils
	using Plots

	p = shepp_logan(256)
	
	heatmap(reverse(p,dims=1), c=:viridis)
end

# ╔═╡ 9a1941ff-3d6c-40b3-b47b-254e859c5b6c
begin
	using FFTW

	pFT = fftshift( fft(p) )
	heatmap(reverse(log.(abs.(pFT).+1.1),dims=1), c=:viridis)
end

# ╔═╡ 71d11567-6098-4a0e-b1ad-5d0319d80ecf
begin
	using Wavelets

	pW = dwt(p, wavelet(WT.db2))
	heatmap(reverse(log.(abs.(pW).+1.1),dims=1), c=:viridis)
end

# ╔═╡ cd73e41b-283e-4fe8-84c8-ee86ce7968ff
begin
	using SparseArrays

	function myfunc(a::AbstractMatrix{T}) where T<:Number
  		# generic fallback
	end

	function myfunc(a::StridedMatrix{T}) where T<:Number
  		# call BLAS routines
	end

	function myfunc(a::SparseMatrixCSC{T}) where T<:Number
 		# exploit sparsity
	end	
end

# ╔═╡ 1257912e-7eea-495b-9668-6b952cbd1b3b
begin
	using MRIReco
	
	# File loading
	f = BrukerFile("data/brukerfileCart")
	raw = RawAcquisitionData(f);
end

# ╔═╡ f12894c1-578d-47b4-8c71-5ccf31a5b84c
using LinearOperators, DelimitedFiles, LinearAlgebra

# ╔═╡ 2a2d075a-accb-11eb-36a3-5bacf8e2498d
md"""
# MRIReco.jl
#### Julia package for MRI reconstruction

by Tobias Knopp and Mirco Grosser
"""

# ╔═╡ 82453453-f2bb-423e-a7bb-099d1ed9dade
md"""
## Getting Started

This tutorial can be used and run in two ways:
* You can directly use Pluto notebook you are reading right now. This needs some preparation discussed next.
* You can also copy/past into a Julia REPL. This is in particular a good way if you are reading an HTML or PDF version of this tutorial.
* Keep in mind that this notebook will run all cells, which needs some time since we also showcase some more computationally intensive methods.
"""

# ╔═╡ 12dcb593-838f-4cf3-81e0-ca035dee56b5
md"""
### Preparation

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
"""

# ╔═╡ c8b39aa6-d85a-4e18-822a-5dccc3b7e5c7
md"""
## Introduction to Julia

* High level programming language with Matlab like syntax
* Provides a powerful REPL for interactive use
* Strongly influenced by C++ and Python
* Generates efficient machine code:
  no need to write performance critical things in C/C++

In more detail:
* Julia uses JIT-compilation to generate efficient machine code
* Julia has a sophisticated type system $\rightarrow$ allows the compiler to identify statically typed code fragments $\rightarrow$ efficient machine code
* Dynamic program style is fine where convenience is more important than performance
* Package manager allows efficient modularization
"""

# ╔═╡ b6a623f0-499c-4bb7-81c6-679bd491452f
md"""
In a nutshell (somewhat biased...): 

**Julia = (C++ $\cup$ Python $\cup$ Matlab) $\cap$ OnlyTheNiceThings**
"""

# ╔═╡ 056d2ff9-cc18-41bb-b35f-07936844cc60
md"""
### Julia Tooling

Commonly Julia is developed in

* Packages: Several source files / modules grouped together
* Scripts: Usually load packages and perform a dedicated operation

Let us explore some of the available functionality
"""

# ╔═╡ 56c864db-9b71-4d9d-81ee-3df9c839bfb4
md"""
Let's do some image processing
"""

# ╔═╡ 65ce8d0b-d4b7-4ce9-8954-803f5628d536
let
	pfilt = imfilter(p, Kernel.gaussian(3))
	heatmap(reverse(pfilt,dims=1).-0.1, c=:viridis)
end

# ╔═╡ ae8a1021-f85f-40a0-ad35-e571e974d3d2
md"Let's switch into frequency space"

# ╔═╡ bca47852-7812-4b3b-8f6a-b90dd727f31f
md"Or let's look at the wavelet transform"

# ╔═╡ f42ade9a-b8bc-499a-93f5-2cc077218e5b


# ╔═╡ 3ef6ac29-bcc6-4648-8c0c-e758ddd8f818
md"""
### Multiple Dispatch

Julia uses multiple dispatch to decide, which method is called
"""

# ╔═╡ e7d573d2-1bac-40c5-b7e6-2abdf94806c7
md"""
### How low level is Julia?

You can directly call C code with zero overhead
"""

# ╔═╡ 4adcb409-4632-4ddb-9e47-75e8d2dc4d83
t = ccall(:clock, Int32, ())

# ╔═╡ aeb4be94-bb33-44b0-9532-2d752fe2739b
md"""
One can look at the assembly code
"""

# ╔═╡ c1cea40b-6de6-4a83-a937-ece75dfca1a8
double(x) = 2*x

# ╔═╡ e1a3e63e-e3ce-4291-b54b-41db9d099c1d
with_terminal() do
  @code_native double(2)
end

# ╔═╡ 41e3e109-d1ad-4f51-835b-b4fa0f2c51ba
with_terminal() do
  @code_native double(2.0)
end

# ╔═╡ 7fac2740-1036-4cf6-88c9-9049f9ad7e41
md"""
## MRIReco: Architecture

* Julia package intended to provide a powerful MRI reconstruction toolbox
* In development (not as stable / feature complete as Gadgetron / BART)
* Focus on reconstruction researchers providing a hackable package
* Modular concept using existing goodies from Julialand
"""

# ╔═╡ 30660c21-4cce-4ed9-8b73-f68e60cbec92
md"""
### Dependency Graph
"""

# ╔═╡ f0cbdf18-6870-48be-8b42-a17305d6d7c1
md"The following shows a dependency graph of MRIReco and its sister package MPIReco (used for Magnetic Particle Imaging Reconstruction). This shows the philosophy to share common functionality into dedicated packages."

# ╔═╡ bfa0e160-a18f-4856-86c2-28485fa2e5e6
PlutoUI.LocalResource("./img/JuliaImagingPackages.png", :width => 700)

# ╔═╡ 2b7d6b3a-c4b9-4fb5-be8f-bc9ad4bd7723
md"""
### Data Types and Flow

Basic overview of the data flow during reconstruction
"""

# ╔═╡ 583111a4-2a79-4fe7-9122-294456c7478a
PlutoUI.LocalResource("./img/overview.svg", :width => 600)

# ╔═╡ f600f9b5-b759-4a54-abd9-cb3f14095583
md"## MRIReco: Examples"

# ╔═╡ bca790f3-da86-4be4-a97a-335741b259b5
md"""
### Example 1: Reconstruction Pipeline

Lets first get some data to perform a simple image reconstruction
"""

# ╔═╡ c3164c4a-8553-4188-aa15-2a7d9e9da106
md"""
The following code shows all the layers how to perform an image reconstruction and store the result:
```julia
using MRIReco

# File loading
f = BrukerFile("data/brukerfileCart")
raw = RawAcquisitionData(f)
acq = AcquisitionData(raw)

# Reconstruction
params = Dict{Symbol, Any}()
params[:reco] = "direct"
params[:reconSize] = (acq.encodingSize[1],acq.encodingSize[2])
Ireco = reconstruction(acq, params)

# Image handling
saveImage("reco.nii", Ireco, true);

```

Let's go through this step by step
"""

# ╔═╡ 984c46f2-e62c-4dac-97a6-eb9895d148a4
md"""
The `RawAcquisitionData` has a similar design as an ISMRMRD data file.
It has two members:
* `raw.params` is a dict that holds general parameter of the file
* `raw.profiles` is a vector with all measured profiles
"""

# ╔═╡ e9b2124f-5cb9-4705-80df-4e2f1b85d846
raw.params

# ╔═╡ 6b62ccaa-5ae8-43bf-a620-989af8754c6c
length(raw.profiles)

# ╔═╡ c948a18c-43c7-4bdc-a51b-591abcb57497
raw.profiles[1]

# ╔═╡ b572b96f-2bc8-4fc4-9113-d0a052c2bce4
md"""
#### Reconstruction

Next we will reconstruct the data. To do so we first convert the raw acquisition data
object into an acquisition data object. This will reorder the data in such a way, that
the reconstruction knows how to deal with it.
"""

# ╔═╡ 470557df-b938-4312-ba85-88e9de4a8564
acq = AcquisitionData(raw);

# ╔═╡ 78614fe5-e980-45fc-873a-105c8dba80e4
md" Now we can perform the reconstruction. It can be controlled using a dictionary of parameters"

# ╔═╡ 46ff4011-29c1-4c94-b794-0f655ec67daa
begin
	local params = Dict{Symbol, Any}()
	params[:reco] = "direct"
	params[:reconSize] = (acq.encodingSize[1],acq.encodingSize[2]);

	Ireco = reconstruction(acq, params);
end

# ╔═╡ 1e9452da-7710-4f3d-a389-2e0bc221acc7
md"The reconstructed image is stored as a 5-dimensional `AxisArray`. The latter is simply an array with named axis. 

Upon inspection of the reconstructed data we see that the first three dimensions correspond to the three spatial dimensions, whereas the fourth dimension encodes different contrasts and the fifth dimension encodes different receive channels."

# ╔═╡ 1f762afd-12bf-43f1-9ea8-1aba7a0a41b0
md"The image data can be stored as a NIfTI file by calling"

# ╔═╡ dc0545e8-e309-480c-8f6e-3a517b0c7c4f
saveImage("reco.nii", Ireco, true);

# ╔═╡ 5f1d779b-7796-4103-b512-b06c8dfc8776
md"To load this image back just call"

# ╔═╡ 52f2408a-1cde-451f-8422-da5f4791e668
Ireco2 = loadImage("reco.nii")

# ╔═╡ 35c98353-b1fd-4916-b143-e27a23174716
md"Lets plot the reconstructed data. We will be using Plots.jl for this"

# ╔═╡ da39843a-d344-4509-b7c7-13ef9ed37a65
heatmap(reverse( circshift(abs.(Ireco[:,:,10,1,1]),(0,5)),dims=1), c=:viridis)

# ╔═╡ 44b4863b-cd5f-4bb1-80ed-fe11458adf84
md"Bruker does provide reconstruction data as well:"

# ╔═╡ a7c32c4b-0299-4bd6-b82f-f5ffb2ed9cbc
begin
	b = BrukerFile("data/brukerfileCart")
	Ibruker = recoData(b)

	heatmap(reverse( Ibruker[:,:,10,1,1]',dims=1), c=:viridis)	
end

# ╔═╡ 617c5a93-36d8-424c-ab54-fa261be07439
md"""
#### Conversion

As we have loaded a `RawAcquistionData` from file, it is also possible to save it.
Currently, it is only possible to save into the ISMRMRD data format, which seems
to be the most sensible choice anyway.
"""

# ╔═╡ 023e6072-b02d-43d2-a958-711e989eaed7
begin
	fout = ISMRMRDFile("outputfile.h5")

	save(fout, raw);
end

# ╔═╡ 079bb385-df0a-47d1-ac9c-b352ae94f07f
md"We may load the data to check that it is still the same:"

# ╔═╡ fa80020c-1056-400d-a753-81993ace124e
begin
	rawCopy = RawAcquisitionData(fout)

	rawCopy.params
end

# ╔═╡ ed50e3d3-a228-4c51-a081-e98db113feb1
md"""
### Example 2: Simulation & Reconstruction

The following code shows a simple MRI simulation and reconstruction considering offresonance effects
"""

# ╔═╡ 107c4794-95a8-4e85-b789-2139091dc5de
begin
	N = 128

	I = MRIReco.shepp_logan(N)
	I = circularShutterFreq!(I,1)
	cmap = 1im*quadraticFieldmap(N,N,125*2pi)

	# simulation parameters
	local params = Dict{Symbol, Any}()
	params[:simulation] = "fast"
	params[:trajName] = "Spiral"
	params[:numProfiles] = 10
	params[:numSamplingPerProfile] = div(N*N,8)
	params[:windings] = 8
	params[:AQ] = 10.0e-3
	params[:correctionMap] = cmap[:,:,1]

	# do simulation
	acqData = simulation(I, params)
end

# ╔═╡ fc564eb0-ee21-43ea-8c65-f07677428ef2
begin
	tr = trajectory(acqData)
	nodes = kspaceNodes(tr)
	plot(nodes[1,1:2048],nodes[2,1:2048], label="trajectory")
end

# ╔═╡ 3ccafef9-5cb4-488f-a697-1c950399ee7c
begin
	# reco parameters
	local params = Dict{Symbol, Any}()
	params[:reco] = "standard"
	params[:regularization] = "L2"
	params[:iterations] = 3
	params[:solver] = "cgnr"
	params[:reconSize] = (N,N)

	INoCorr = reconstruction(acqData, params)
 
	params[:correctionMap] = cmap
	ICorr = reconstruction(acqData, params)

	heatmap(reverse( abs.(INoCorr[:,:,1,1,1]),dims=1), c=:viridis)
end

# ╔═╡ 6500d0af-d0f4-4dd7-a84a-e70b3424ee2f
heatmap(reverse( abs.(ICorr[:,:,1,1,1]),dims=1), c=:viridis)

# ╔═╡ 1c5c2818-fdc5-4214-bd3a-7c9a35e27153
md"""
### Example 3: Iterative Reconstruction
Iterative reconstructions can be performed using a high-level interface. All the relevant reconstruction parameters are collected in a `Dict{Symbol,Any}`. The reconstruction methods use those to formulate the image reconstruction problem and solve it.

To illustrate that, lets consider the following dataset, which contains one slice of a 3dFSE scan obtained from [mridata.org](mridata.org).
"""

# ╔═╡ f331e51f-a1e0-4959-8bc1-ef3b62d2d787
begin
	# load fully sampled data
	fKnee = ISMRMRDFile("data/knee_3dFSE_slice170.h5")
	acqDataKnee = AcquisitionData(fKnee);
end

# ╔═╡ f71ca727-1c64-41fd-a928-33eea4328598
md"For reference, let us reconstruct the fully sampled data"

# ╔═╡ 34f7fc8e-fd22-46e3-984a-940d504ce235
begin
	# reconstruct
	local params = Dict{Symbol, Any}()
	params[:reco] = "direct"
	params[:reconSize] = (320,320) # this size is also contained in acqData.encodingSize

	img = reconstruction(acqDataKnee, params)
	img = sqrt.(sum(img.^2,dims=5))

	# show image
	heatmap(reverse( abs.(img[:,:,1,1,1]),dims=1), c=:viridis)
end

# ╔═╡ c66c94b3-e5a2-4b87-a943-e69cb49092bd
md"To simulate an undersampled reconstruction, we retrospectively undersample the data using a Poisson disk pattern."

# ╔═╡ 96ed1677-3042-4aca-bb03-1719ebc0cba8
begin
	# sampling
	redFac= 4.0
	acqDataSub = sample_kspace(acqDataKnee,redFac,"poisson",calsize=30,profiles=false);

	# show sampling pattern
	msk = zeros(acqDataSub.encodingSize[1],acqDataSub.encodingSize[2])
	msk[acqDataSub.subsampleIndices[1]] .= 1

	heatmap(reverse( msk,dims=1), c=:grays)
end

# ╔═╡ 9592634e-ea3d-42eb-b013-cea1ff98601e
md"Estimate the coil sensitivities using ESPIRiT and reconstruct using SENSE"

# ╔═╡ 9aec589f-3615-49f6-9a92-c87afb7568a9
begin
	# coil sensitivities
	smaps = espirit(acqDataSub,(6,6),30,eigThresh_1=0.035,eigThresh_2=0.98)

	# SENSE reconstruction
	local params = Dict{Symbol, Any}()
	params[:reco] = "multiCoil"
	params[:reconSize] = (320,320)
	params[:senseMaps] = smaps

	params[:solver] = "cgnr"
	params[:regularization] = "L2"
	params[:λ] = 1.e1
	params[:iterations] = 5

	img_l2 = reconstruction(acqDataSub, params)

	# show image
	heatmap(reverse( abs.(img_l2[:,:,1,1,1]),dims=1), c=:viridis)
end

# ╔═╡ 936a7154-f601-4f06-ba41-e684e40677f3
md"Using ℓ¹-Wavelet-regularization recquires us to change some parameters."

# ╔═╡ d21a0949-6256-4325-84ca-3f78ab382a51
begin	
	local params = Dict{Symbol, Any}()
	params[:reco] = "multiCoil"
	params[:reconSize] = (320,320)
	params[:senseMaps] = smaps

	params[:solver] = "admm"
	params[:regularization] = "L1"
	params[:sparseTrafo] = "Wavelet"
	params[:λ] = 2.e1
	params[:iterations] = 50
	params[:ρ] = 1.0
	params[:absTol] = 1.e-4
	params[:relTol] = 1.e-2
	params[:tolInner] = 1.e-2

	img_wavelet = reconstruction(acqDataSub, params)

# show image
heatmap(reverse( abs.(img_wavelet[:,:,1,1,1]),dims=1), c=:viridis)
end

# ╔═╡ 63d00efb-c10e-4c8a-8160-d6156a6fcfd2
md"""
### Example 4: Extending MRIReco

Reconstructions parameters, such as sparsifying transforms or the form of regularization can be incorporated in two ways:
* specify name of of an already implemented transform or regularizer
* directly pass the `Regularization` or `LinearOperator` object to be used

The second option yields a simple way to
* modify reconstructions algorithms
* implement and test new reconstruction methods

#### Learned Sparsifying Transform

Let us implement a simple data-driven sparsifying transform (see: S. Ravishankar and Y. Bresler, IEEE Trans. Med. Imaging, 30 (5), 2011). 
* Learn a dictionary from a reference image (e.g. adjacent slice) using KSVD  $\color{green}\checkmark$
* Implement sparsifying transform which analyses the input image in terms of the dictionary
"""

# ╔═╡ 873a78e3-f4cd-425a-8507-0d24cef8b771
PlutoUI.LocalResource("img/dictTrafo.png", :width => 500)

# ╔═╡ 484b7d08-6af4-4795-b7d6-e741ff3fb503
function analyzeImage(x::Vector{T}, D::Matrix{T}, xsize::NTuple{2,Int64},
		              psize::NTuple{2,Int64}; t0::Int64=size(D,2),
		              tol=1e-3, nmax::Int64=size(D,1)) where T
	nx,ny = xsize
	px,py = psize
	x = reshape(x,nx,ny)
	x_pad = repeat(x,2,2)[1:nx+px-1,1:ny+py-1] # pad image using periodic boundary conditions
	α = zeros(T,size(D,2),nx,ny)
	patch = zeros(T,px*py)
	for j=1:ny, i=1:nx
		patch[:] .= vec(x_pad[i:i+px-1,j:j+py-1])
		norm(patch)==0 && continue
		# matchingpursuit is contained in Wavelets.jl
		α[:,i,j] .= matchingpursuit(patch, x->D*x, x->transpose(D)*x, tol, nmax)
	end
	return vec(α)
end

# ╔═╡ 316f47e9-58a8-45b7-89d4-375a8c0024a2
function synthesizeImage(α::Vector{T}, D::Matrix{T}, xsize::NTuple{2,Int64},
		                 psize::NTuple{2,Int64}) where T
	nx,ny = xsize
	px,py = psize
	x = zeros(T,nx+px-1,ny+py-1)
	α = reshape(α,:,nx,ny)
	for j=1:ny, i=1:nx
		x[i:i+px-1,j:j+py-1] .+= reshape(D*α[:,i,j],px,py)
	end
	return vec(x[1:nx,1:ny])/(px*py)
end

# ╔═╡ 41b190f7-f9f6-414b-b260-3c5316c08619
function dictOp(D::Matrix{T},xsize::NTuple{2,Int64},psize::NTuple{2,Int64},
		        tol::Float64=1.e-3,nmax::Int64=size(D,1)) where T
	produ = x->analyzeImage(x,D,xsize,psize,tol=tol,nmax=nmax)
	ctprodu = x->synthesizeImage(x,D,xsize,psize)
	return LinearOperator(prod(xsize)*size(D,2),prod(xsize),false,false
			, produ
			, nothing
			, ctprodu )
end

# ╔═╡ 38b1f34c-c71a-4f8f-a672-27e745898deb
md"Next, let us Load the Dictionary, build the sparsifying transform and perform reconstruct the knee data"

# ╔═╡ e85a04b1-7a0d-4b51-a18c-48e1cce6a09e
begin
	# load dict
	D = ComplexF64.(readdlm("data/kneeDict165.tsv"))
	
	local params = Dict{Symbol, Any}()
	params[:reco] = "multiCoil"
	params[:reconSize] = (320,320)
	params[:senseMaps] = smaps

	params[:solver] = "admm"
	params[:regularization] = "L1"
	params[:sparseTrafo] = dictOp(D,(320,320),(6,6),1.e-3,12)
	params[:λ] = 2.e1
	params[:iterations] = 30
	params[:ρ] = 1.0
	params[:absTol] = 1.e-4
	params[:relTol] = 1.e-2
	params[:tolInner] = 1.e-2

	img_dict = reconstruction(acqDataSub, params)

	# show image
	heatmap(reverse( abs.(img_dict[:,:,1,1,1]),dims=1), c=:viridis)
end

# ╔═╡ 2514b8f8-2492-41b1-8c85-0debe400af12
md"""
#### Developing MRIReco.jl

Changing MRIReco.jl is pretty easy. Go into the package manager and develop the package

```julia
] dev MRIReco
```
After developing the package, the code is located in

   `~/.julia/dev/MRIReco` 
   
and can be directly modified. 

Note that Julia requires a restart, if you change code. This limitation
can be circumvented by using the `Revise.jl` package (look at the Github page).

"""

# ╔═╡ f78b55aa-7f5b-4412-a293-0eccac11dd91
md"## Summary
MRIReco is a MRI reconstruction package written in julia, which aims to:
* be fast
* accessible
* easily extendable

The package is still relatively young. Thus, while being usable, further updates and extensions can be expected.

To further improve the package, we welcome feedback, bug reports and ideas for new features to be added.

For further details on *MRIReco.jl*, have a look at:
* The packages github repository:  [https://github.com/MagneticResonanceImaging/MRIRecoTutorial](https://github.com/MagneticResonanceImaging/MRIRecoTutorial)
* The accompanying paper: [https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.28792](https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.28792)
"

# ╔═╡ Cell order:
# ╟─2a2d075a-accb-11eb-36a3-5bacf8e2498d
# ╟─96d4f44e-49d1-4398-b7f8-4daa7290b287
# ╟─82453453-f2bb-423e-a7bb-099d1ed9dade
# ╟─12dcb593-838f-4cf3-81e0-ca035dee56b5
# ╟─c8b39aa6-d85a-4e18-822a-5dccc3b7e5c7
# ╟─b6a623f0-499c-4bb7-81c6-679bd491452f
# ╟─056d2ff9-cc18-41bb-b35f-07936844cc60
# ╠═915d5707-daf4-4127-abf8-68a5528d9209
# ╟─56c864db-9b71-4d9d-81ee-3df9c839bfb4
# ╠═65ce8d0b-d4b7-4ce9-8954-803f5628d536
# ╟─ae8a1021-f85f-40a0-ad35-e571e974d3d2
# ╠═9a1941ff-3d6c-40b3-b47b-254e859c5b6c
# ╟─bca47852-7812-4b3b-8f6a-b90dd727f31f
# ╠═71d11567-6098-4a0e-b1ad-5d0319d80ecf
# ╠═f42ade9a-b8bc-499a-93f5-2cc077218e5b
# ╟─3ef6ac29-bcc6-4648-8c0c-e758ddd8f818
# ╠═cd73e41b-283e-4fe8-84c8-ee86ce7968ff
# ╟─e7d573d2-1bac-40c5-b7e6-2abdf94806c7
# ╠═4adcb409-4632-4ddb-9e47-75e8d2dc4d83
# ╟─aeb4be94-bb33-44b0-9532-2d752fe2739b
# ╠═c1cea40b-6de6-4a83-a937-ece75dfca1a8
# ╠═e1a3e63e-e3ce-4291-b54b-41db9d099c1d
# ╠═41e3e109-d1ad-4f51-835b-b4fa0f2c51ba
# ╟─7fac2740-1036-4cf6-88c9-9049f9ad7e41
# ╟─30660c21-4cce-4ed9-8b73-f68e60cbec92
# ╟─f0cbdf18-6870-48be-8b42-a17305d6d7c1
# ╟─bfa0e160-a18f-4856-86c2-28485fa2e5e6
# ╟─2b7d6b3a-c4b9-4fb5-be8f-bc9ad4bd7723
# ╟─583111a4-2a79-4fe7-9122-294456c7478a
# ╟─f600f9b5-b759-4a54-abd9-cb3f14095583
# ╟─bca790f3-da86-4be4-a97a-335741b259b5
# ╟─c3164c4a-8553-4188-aa15-2a7d9e9da106
# ╠═1257912e-7eea-495b-9668-6b952cbd1b3b
# ╟─984c46f2-e62c-4dac-97a6-eb9895d148a4
# ╠═e9b2124f-5cb9-4705-80df-4e2f1b85d846
# ╠═6b62ccaa-5ae8-43bf-a620-989af8754c6c
# ╠═c948a18c-43c7-4bdc-a51b-591abcb57497
# ╟─b572b96f-2bc8-4fc4-9113-d0a052c2bce4
# ╠═470557df-b938-4312-ba85-88e9de4a8564
# ╟─78614fe5-e980-45fc-873a-105c8dba80e4
# ╠═46ff4011-29c1-4c94-b794-0f655ec67daa
# ╟─1e9452da-7710-4f3d-a389-2e0bc221acc7
# ╟─1f762afd-12bf-43f1-9ea8-1aba7a0a41b0
# ╠═dc0545e8-e309-480c-8f6e-3a517b0c7c4f
# ╟─5f1d779b-7796-4103-b512-b06c8dfc8776
# ╠═52f2408a-1cde-451f-8422-da5f4791e668
# ╟─35c98353-b1fd-4916-b143-e27a23174716
# ╠═da39843a-d344-4509-b7c7-13ef9ed37a65
# ╟─44b4863b-cd5f-4bb1-80ed-fe11458adf84
# ╠═a7c32c4b-0299-4bd6-b82f-f5ffb2ed9cbc
# ╟─617c5a93-36d8-424c-ab54-fa261be07439
# ╠═023e6072-b02d-43d2-a958-711e989eaed7
# ╟─079bb385-df0a-47d1-ac9c-b352ae94f07f
# ╠═fa80020c-1056-400d-a753-81993ace124e
# ╟─ed50e3d3-a228-4c51-a081-e98db113feb1
# ╠═107c4794-95a8-4e85-b789-2139091dc5de
# ╠═fc564eb0-ee21-43ea-8c65-f07677428ef2
# ╠═3ccafef9-5cb4-488f-a697-1c950399ee7c
# ╠═6500d0af-d0f4-4dd7-a84a-e70b3424ee2f
# ╟─1c5c2818-fdc5-4214-bd3a-7c9a35e27153
# ╠═f331e51f-a1e0-4959-8bc1-ef3b62d2d787
# ╟─f71ca727-1c64-41fd-a928-33eea4328598
# ╠═34f7fc8e-fd22-46e3-984a-940d504ce235
# ╟─c66c94b3-e5a2-4b87-a943-e69cb49092bd
# ╠═96ed1677-3042-4aca-bb03-1719ebc0cba8
# ╟─9592634e-ea3d-42eb-b013-cea1ff98601e
# ╠═9aec589f-3615-49f6-9a92-c87afb7568a9
# ╟─936a7154-f601-4f06-ba41-e684e40677f3
# ╠═d21a0949-6256-4325-84ca-3f78ab382a51
# ╟─63d00efb-c10e-4c8a-8160-d6156a6fcfd2
# ╟─873a78e3-f4cd-425a-8507-0d24cef8b771
# ╠═f12894c1-578d-47b4-8c71-5ccf31a5b84c
# ╠═484b7d08-6af4-4795-b7d6-e741ff3fb503
# ╠═316f47e9-58a8-45b7-89d4-375a8c0024a2
# ╠═41b190f7-f9f6-414b-b260-3c5316c08619
# ╟─38b1f34c-c71a-4f8f-a672-27e745898deb
# ╠═e85a04b1-7a0d-4b51-a18c-48e1cce6a09e
# ╟─2514b8f8-2492-41b1-8c85-0debe400af12
# ╟─f78b55aa-7f5b-4412-a293-0eccac11dd91
