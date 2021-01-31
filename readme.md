<p align="center"><img width="25%" src="docs/umapxnet_big.png" /></p>

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
* **Window Toolbox**. Includes a set of tools to synthesizing and orthogonalizing window functions. It implements discrete short-time Fourier transforms and Weyl-Heisenberg transforms ([**Gabor analysis**](https://github.com/asiryan/Weyl-Heisenberg-Toolbox)) for real and complex signals.
* **Image Processing Toolbox**. Contains efficient algorithms to processing, correcting and analyzing 32-bit raster images.
* **Video Processing Toolbox**. Includes a set of tools to video streaming and processing.

# Installation
You can build **UMapx.NET** from sources or install to your own project using nuget package manager.
| Version | Specification | OS | Platform | Download | Package |
|-------------|-------------|-------------|-------------|--------------|--------------|
| 4.0.3.3 | .NET Standard 2.0 | Windows | AnyCPU | [Release](https://github.com/asiryan/UMapx.NET/releases/) | [NuGet](https://www.nuget.org/packages/UMapx/) |

# Namespaces
```c#
using UMapx.Analysis;
using UMapx.Colorspace;
using UMapx.Core;
using UMapx.Decomposition;
using UMapx.Distribution;
using UMapx.Imaging;
using UMapx.Response;
using UMapx.Transform;
using UMapx.Video;
using UMapx.Wavelet;
using UMapx.Window;
```

# Relation to other frameworks
**UMapx.NET** is based on several separate frameworks and toolboxes ([AForge.NET](https://github.com/andrewkirillov/AForge.NET), [Accord.NET](https://github.com/accord-net/framework)) and some functions have been ported from other programming languages (MATLAB, Fortran, C++). At the same time, the purpose of this generalization was the best declarative understanding of algorithms and tools of digital signal processing, as well as optimization and performance improvement.  
**UMapx.NET** is faster than AForge.NET or Accord.NET in signal processing tasks and contains a larger set of functions for matrix analysis, linear algebra and functional analysis tasks.  

# License
**MIT**  

# References
A full list of references is given in a separate [**file**](docs/references.pdf).  
