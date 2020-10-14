<p align="center"><img width="35%" src="docs/umapxnet_big.png" /></p>

**UMapx.NET** is a dynamic link library for digital signal processing and analysis developed in C#.

# About
### Contains ready-made software tools:
* color spaces and their transformations,
* real and complex algebra,
* statistical distributions,
* special math functions,
* digital response filters,
* discrete orthogonal transforms and more.

### Suitable for a wide range of tasks:
* symbolic and graphical visualization of data,
* functional, vector and matrix analysis,
* interpolation, approximation and optimization of functions,
* numerical differentiation and integration,
* solving equations,
* matrix factorization.

### Includes special math toolboxes:
* **Wavelet Toolbox**. Provides wide functionality for the study of discrete and continuous wavelets. The toolbox also includes algorithms for discrete one-dimensional and two-dimensional wavelet transforms of real and complex signals.
* **Window Toolbox**. Includes a set of tools for synthesizing and orthogonalizing window functions. It implements discrete short-time Fourier transforms and Weyl-Heisenberg transforms ([**Gabor analysis**](https://github.com/asiryan/Weyl-Heisenberg-Bases-Toolbox)) for real and complex signals.
* **Image Processing Toolbox**. Contains efficient algorithms for processing, correcting and analyzing 32-bit raster images.

# Installation
Download from [**release**](release) folder and add **UMapx.dll** to your project or use [**nuget**](https://www.nuget.org/packages/UMapx/) package manager.  
```c#
using UMapx.Analysis;
using UMapx.Colorspace;
using UMapx.Core;
using UMapx.Decomposition;
using UMapx.Distribution;
using UMapx.Imaging;
using UMapx.Response;
using UMapx.Transform;
using UMapx.Wavelet;
using UMapx.Window;
```    

# License
**MIT**  

# References
A full list of references is given in a separate [**file**](docs/references.pdf).  
Imaging test data [**folder**](https://yadi.sk/d/rix2T-ARtGm4Dg).  
