<p align="center"><img width="25%" src="docs/umapxnet_big.png" /></p>
<p align="center"> Cross-platform .NET library for digital signal processing </p>    
<p align="center"><i> Every journey begins in the mind... </i></p>    

# UMapx
### Contains ready-made math tools:
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

### Includes special toolboxes:
* **Wavelet Toolbox**. Provides wide functionality for the study of discrete and continuous wavelets. The toolbox also includes algorithms for discrete one-dimensional and two-dimensional wavelet transforms of real and complex signals.
* **Window Toolbox**. Includes a set of tools for synthesizing and orthogonalizing window functions. It implements discrete short-time Fourier and Weyl–Heisenberg transforms ([Gabor analysis](https://github.com/asiryan/Weyl-Heisenberg-Toolbox)) for real and complex signals.
* **Image Processing Toolbox**. Contains efficient algorithms for processing, correcting, and analyzing 32-bit images.
* **Video Processing Toolbox**. Includes a set of tools for video streaming and processing.

# Supported types
**UMapx** supports only
* 32 bit types - `float`, `ComplexF`, etc (compatible with [System.Numerics](https://docs.microsoft.com/ru-ru/dotnet/api/system.numerics?view=netframework-4.8), [NAudio](https://github.com/naudio/NAudio) and other libraries),
* 32 bit image `BitmapData` format - `Format32bppArgb` (compatible with [AForge.NET](https://github.com/andrewkirillov/AForge.NET), [Accord.NET](https://github.com/accord-net/framework/) and so on),
* 24 bit video `BitmapData` format - `Format24bppRgb` (32 bit format not recommended).

# Installation
You can build **UMapx** from sources or install to your own project using nuget package manager.
| Specification | OS | Platform | Download | Package |
|-------------|-------------|-------------|--------------|--------------|
| .NET Standard 2.0 | Cross-platform | AnyCPU | [Release](https://github.com/asiryan/UMapx.NET/releases/) | [NuGet](https://www.nuget.org/packages/UMapx/) |

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
using UMapx.Visualization;
using UMapx.Wavelet;
using UMapx.Window;
```

# Examples of usage
* [Local Laplacian filters](https://github.com/asiryan/Local-Laplacian-filters) - NET Framework desktop application for HDR imaging.
* [Portrait mode effect](https://github.com/asiryan/Portrait-mode) - High quality implementation of the portrait mode effect using Neural Networks.
* [FaceONNX](https://github.com/FaceONNX/FaceONNX) - Face analytics library based on deep neural networks and ONNX runtime.

# Relation to other frameworks
**UMapx** builds on several existing frameworks (AForge.NET, Accord.NET, ALGLIB, etc.). Some functions have been ported from other programming languages, toolboxes, and libraries (Fortran, MATLAB, C++, Python). The goal of this generalization is to provide a declarative understanding of digital signal processing algorithms and to improve optimization and performance. **UMapx** is faster than AForge.NET and Accord.NET for common signal-processing tasks and includes a larger set of functions for matrix analysis, linear algebra, and functional analysis.

# License
**MIT**  

# References
A full list of references is given in a separate [file](docs/references.pdf).  
