using System.Numerics;
using UMapx.Core;
using UMapx.Transform;
using Xunit;
using static UMapx.Tests.NumericAssert;

namespace UMapx.Tests;

[Trait("Category", "Transform")]
public class TransformAuditTests
{
    internal static readonly string[] Names={"Cosine","FastCosine","Sine","FastSine","Chebyshev","FastChebyshev","Hartley","FastHartley","HartleyDirect","FastHartleyDirect","WalshHadamard","FastWalshHadamard","Fourier","FastFourier","Laplace","FastLaplace","Delta"};
    internal static ITransform Create(string name,Direction direction=Direction.Vertical,bool normalized=true)=>name switch
    {
        "Cosine"=>new CosineTransform(direction),"FastCosine"=>new FastCosineTransform(direction),
        "Sine"=>new SineTransform(direction),"FastSine"=>new FastSineTransform(direction),
        "Chebyshev"=>new ChebyshevTransform(direction),"FastChebyshev"=>new FastChebyshevTransform(direction),
        "Hartley"=>new HartleyTransform(normalized,SpectrumType.Fourier,direction),"FastHartley"=>new FastHartleyTransform(normalized,SpectrumType.Fourier,direction),
        "HartleyDirect"=>new HartleyTransform(normalized,SpectrumType.Hartley,direction),"FastHartleyDirect"=>new FastHartleyTransform(normalized,SpectrumType.Hartley,direction),
        "WalshHadamard"=>new WalshHadamardTransform(normalized,direction),"FastWalshHadamard"=>new FastWalshHadamardTransform(normalized,direction),
        "Fourier"=>new FourierTransform(normalized,direction),"FastFourier"=>new FastFourierTransform(normalized,direction),
        "Laplace"=>new LaplaceTransform(.07f,normalized,direction),"FastLaplace"=>new FastLaplaceTransform(.07f,normalized,direction),
        _=>new DeltaTransform(direction)
    };

    private static bool ComplexOnly(string name)=>name.Contains("Fourier")||name.Contains("Laplace");
    private static bool Scaled(string name)=>name.Contains("Hartley")||name.Contains("Walsh")||ComplexOnly(name);
    private static Complex Weight(string name,int k,int j,int n,bool normalized)
    {
        if(name.Contains("Cosine"))return Math.Cos(Math.PI*(j+.5)*k/n)*Math.Sqrt((k==0?1:2.0)/n);
        if(name.Contains("Sine"))return Math.Sin(Math.PI*(j+1)*(k+1)/(n+1))*Math.Sqrt(2.0/(n+1));
        if(name.Contains("Chebyshev"))return n==1?1:Math.Cos(Math.PI*k*j/(n-1))*Math.Sqrt(2.0/(n-1))*(k==0||k==n-1?1/Math.Sqrt(2):1)*(j==0||j==n-1?1/Math.Sqrt(2):1);
        double scale=normalized?Math.Sqrt(n):1;
        if(name.Contains("Hartley"))return (Math.Cos(2*Math.PI*j*k/n)+Math.Sin(2*Math.PI*j*k/n))/scale;
        if(name.Contains("Walsh"))return (BitOperations.PopCount((uint)(j&k))%2==0?1:-1)/scale;
        if(ComplexOnly(name))return Complex.Exp(-Complex.ImaginaryOne*2*Math.PI*j*k/n)*(name.Contains("Laplace")?Math.Exp(-.07f*j):1)/scale;
        return j==k?1:j==k-1?-1:0;
    }

    public static IEnumerable<object[]> VectorCases()
    {
        foreach(string name in Names)foreach(int n in new[]{1,2,3,4,5,8,17})
        {
            if(name.Contains("Walsh")&&(n&(n-1))!=0)continue;
            foreach(bool normalized in new[]{false,true})foreach(bool complex in new[]{false,true})
                if(complex||!ComplexOnly(name))yield return new object[]{name,n,normalized,complex};
        }
    }

    [Theory, MemberData(nameof(VectorCases))]
    public void OneDimensionalTransformsAgreeWithIndependentBasisFunctions(string name,int n,bool normalized,bool complex)
    {
        var x=Enumerable.Range(0,n).Select(i=>new Complex32((float)Math.Sin(i*.71)+.3f,complex?(float)Math.Cos(i*.37):0)).ToArray();
        var d=Create(name,normalized:normalized);
        var actual=complex?d.Forward(x):d.Forward(x.Select(z=>z.Real).ToArray()).Select(v=>new Complex32(v,0)).ToArray();
        Assert.Equal(n,actual.Length);
        for(int k=0;k<n;k++)
        {
            Complex expected=0;for(int j=0;j<n;j++)expected+=Weight(name,k,j,n,normalized)*(Complex)x[j];
            Close(expected,actual[k],3e-4);
        }
        var restored=complex?d.Backward(actual):d.Backward(actual.Select(z=>z.Real).ToArray()).Select(v=>new Complex32(v,0)).ToArray();
        double gain=!normalized&&Scaled(name)?n:1;
        for(int i=0;i<n;i++)Close((Complex)x[i]*gain,restored[i],.001);
    }

    public static IEnumerable<object[]> MatrixCases()
    {
        foreach(string name in Names)foreach(var direction in new[]{Direction.Horizontal,Direction.Vertical,Direction.Both})
            foreach(bool complex in new[]{false,true})if(complex||!ComplexOnly(name))yield return new object[]{name,direction,complex};
    }

    [Theory, MemberData(nameof(MatrixCases))]
    public void MatrixTransformsUseTheRequestedAxesAndPreserveComplexComponents(string name,Direction direction,bool complex)
    {
        int m=name.Contains("Walsh")?4:3,n=name.Contains("Walsh")?8:5;
        var a=new Complex32[m,n];var real=new float[m,n];
        for(int i=0;i<m;i++)for(int j=0;j<n;j++){real[i,j]=(float)Math.Sin(i+j*.3);a[i,j]=new(real[i,j],complex?(float)Math.Cos(i*.7-j*.4):0);}
        var d=Create(name,direction);var actual=complex?d.Forward(a):ToComplex(d.Forward(real));
        Assert.Equal(m,actual.GetLength(0));Assert.Equal(n,actual.GetLength(1));
        for(int i=0;i<m;i++)for(int j=0;j<n;j++)
        {
            Complex expected=0;
            // Complex matrix transforms use the two-sided basis convention U*A*V^H.
            for(int r=0;r<m;r++)for(int c=0;c<n;c++)
                expected+=(Complex)a[r,c]*(direction==Direction.Horizontal?(r==i?1:0):Weight(name,i,r,m,true))*(direction==Direction.Vertical?(c==j?1:0):Complex.Conjugate(Weight(name,j,c,n,true)));
            Close(expected,actual[i,j],.001);
        }
        var restored=complex?d.Backward(actual):ToComplex(d.Backward(ToReal(actual)));
        for(int i=0;i<m;i++)for(int j=0;j<n;j++)Close((Complex)a[i,j],restored[i,j],.002);
    }

    internal static Complex32[,] ToComplex(float[,] a)
    {var r=new Complex32[a.GetLength(0),a.GetLength(1)];for(int i=0;i<r.GetLength(0);i++)for(int j=0;j<r.GetLength(1);j++)r[i,j]=a[i,j];return r;}
    private static float[,] ToReal(Complex32[,] a)
    {var r=new float[a.GetLength(0),a.GetLength(1)];for(int i=0;i<r.GetLength(0);i++)for(int j=0;j<r.GetLength(1);j++)r[i,j]=a[i,j].Real;return r;}

    [Theory] [InlineData(3)] [InlineData(4)] [InlineData(5)] [InlineData(8)] [InlineData(17)]
    public void HilbertTransformHasTheCorrectPhaseAndRemovesDcAndNyquist(int n)
    {
        var x=Enumerable.Range(0,n).Select(j=>new Complex32(2+(float)Math.Cos(2*Math.PI*j/n),.25f*(float)Math.Sin(2*Math.PI*j/n))).ToArray();
        foreach(ITransform d in new ITransform[]{new HilbertTransform(),new FastHilbertTransform()})
        {
            var actual=d.Forward(x);var restored=d.Backward(actual);
            for(int j=0;j<n;j++)
            {
                Close(new Complex(Math.Sin(2*Math.PI*j/n),-.25*Math.Cos(2*Math.PI*j/n)),actual[j],1e-4);
                Close((Complex)x[j]-2,restored[j],1e-4);
            }
        }
    }

    [Theory] [InlineData(4,4)] [InlineData(8,12)] [InlineData(7,9)] [InlineData(16,16)]
    public void LaplacianPyramidsReconstructMatricesIncludingOddDimensions(int m,int n)
    {
        var x=NumericAssert.Matrix(m,n);var d=new LaplacianPyramidTransform(3,2);
        Close(x,d.Backward(d.Forward(x)),.001f);
        var z=ToComplex(x);for(int i=0;i<m;i++)for(int j=0;j<n;j++)z[i,j].Imag=(i-j)*.1f;
        var restored=d.Backward(d.Forward(z));
        for(int i=0;i<m;i++)for(int j=0;j<n;j++)Close((Complex)z[i,j],restored[i,j],.001);
    }

    [Theory] [InlineData(4,false)] [InlineData(8,false)] [InlineData(9,false)] [InlineData(4,true)] [InlineData(8,true)] [InlineData(9,true)]
    public void LaplacianPyramidsReconstructVectors(int n,bool complex)
    {
        var x=Enumerable.Range(0,n).Select(i=>(float)Math.Sin(i*.7)).ToArray();var d=new LaplacianPyramidTransform(3,2);
        if(!complex)Close(x,d.Backward(d.Forward(x)),.001f);
        else
        {
            var z=x.Select((v,i)=>new Complex32(v,.1f*i)).ToArray();var r=d.Backward(d.Forward(z));
            for(int i=0;i<n;i++)Close((Complex)z[i],r[i],.001);
        }
    }

    [Theory] [InlineData(2)] [InlineData(3)] [InlineData(4)]
    public void GaussianPyramidPreservesConstantSignalsAtEveryScale(int levels)
    {
        var d=new GaussianPyramidTransform(levels,2);var a=new float[16,24];for(int i=0;i<16;i++)for(int j=0;j<24;j++)a[i,j]=.375f;
        var x=Enumerable.Repeat(.375f,32).ToArray();var z=x.Select(v=>new Complex32(v,.125f)).ToArray();var c=ToComplex(a);
        var p=d.Forward(a);Assert.Equal(levels,p.Length);foreach(var v in p)Assert.All(v.Cast<float>(),t=>Close(.375,t,1e-5));
        var q=d.Forward(x);foreach(var v in q)Assert.All(v,t=>Close(.375,t,1e-5));
        var r=d.Forward(z);foreach(var v in r)Assert.All(v,t=>Close(new Complex(.375,.125),t,1e-5));
        var s=d.Forward(c);foreach(var v in s)Assert.All(v.Cast<Complex32>(),t=>Close(new Complex(.375,0),t,1e-5));
        Assert.Throws<NotSupportedException>(()=>d.Backward(p));Assert.Throws<NotSupportedException>(()=>d.Backward(q));
        Assert.Throws<NotSupportedException>(()=>d.Backward(r));Assert.Throws<NotSupportedException>(()=>d.Backward(s));
    }

    [Theory] [InlineData(ThresholdMode.Abs)] [InlineData(ThresholdMode.Under)] [InlineData(ThresholdMode.Over)]
    public void ThresholdFiltersRespectEqualityAndSignedComponents(ThresholdMode mode)
    {
        var input=new[]{-2f,-1,-.5f,0,.5f,1,2};float threshold=1;
        float F(float v)=>mode switch{ThresholdMode.Abs=>Math.Abs(v)<threshold?0:v,ThresholdMode.Under=>v<threshold?0:v,_=>v>threshold?0:v};
        var d=new ThresholdFilter(threshold,mode);var x=(float[])input.Clone();d.Apply(x);Close(input.Select(F).ToArray(),x);
        var a=new float[2,input.Length];var c=new Complex32[2,input.Length];var z=new Complex32[input.Length];
        for(int i=0;i<2;i++)for(int j=0;j<input.Length;j++){a[i,j]=input[j];c[i,j]=new(input[j],-.5f*input[j]);z[j]=c[i,j];}
        var expected=z.Select(v=>mode==ThresholdMode.Abs?(v.Abs<threshold?Complex32.Zero:v):new Complex32(F(v.Real),F(v.Imag))).ToArray();
        d.Apply(a);d.Apply(c);d.Apply(z);
        for(int j=0;j<input.Length;j++){Close((Complex)expected[j],z[j]);for(int i=0;i<2;i++){Close(F(input[j]),a[i,j]);Close((Complex)expected[j],c[i,j]);}}
    }

    public static IEnumerable<object[]> FilterCases()
    {
        foreach(string name in new[]{"Guided","Bilateral","BilateralGrid","Domain","LocalLaplacian","Laplacian"})
            foreach(bool matrix in new[]{false,true})foreach(bool complex in new[]{false,true})yield return new object[]{name,matrix,complex};
    }

    [Theory, MemberData(nameof(FilterCases))]
    public void SmoothingAndDetailFiltersPreserveConstants(string name,bool matrix,bool complex)
    {
        IFilter d=name switch{"Guided"=>new GuidedFilter(2),"Bilateral"=>new BilateralFilter(2,.1f,8),"BilateralGrid"=>new BilateralGridFilter(2,.1f),"Domain"=>new DomainTransformFilter(2,.1f),"LocalLaplacian"=>new LocalLaplacianFilter(2,.1f,5,3),_=>new LaplacianPyramidFilter(new LaplacianPyramidTransform(3,2))};
        var x=Enumerable.Repeat(.375f,16).ToArray();var z=x.Select(v=>new Complex32(v,.125f)).ToArray();
        var a=new float[8,12];var c=new Complex32[8,12];for(int i=0;i<8;i++)for(int j=0;j<12;j++){a[i,j]=.375f;c[i,j]=new(.375f,.125f);}
        if(name=="LocalLaplacian"&&complex)
        {if(matrix)Assert.Throws<NotSupportedException>(()=>d.Apply(c));else Assert.Throws<NotSupportedException>(()=>d.Apply(z));return;}
        if(matrix&&complex){d.Apply(c);Assert.All(c.Cast<Complex32>(),v=>Close(new Complex(.375,.125),v,.001));}
        else if(matrix){d.Apply(a);Assert.All(a.Cast<float>(),v=>Close(.375,v,.001));}
        else if(complex){d.Apply(z);Assert.All(z,v=>Close(new Complex(.375,.125),v,.001));}
        else{d.Apply(x);Assert.All(x,v=>Close(.375,v,.001));}

        if(name is "Guided" or "Bilateral" or "Domain")
        {
            // A small interior impulse must spread and lose height under smoothing.
            // Constant preservation alone would also accept an identity implementation.
            int center=matrix?4*12+6:8;
            x[8]+=.05f;z[8]+=new Complex32(.05f,.025f);
            a[4,6]+=.05f;c[4,6]+=new Complex32(.05f,.025f);
            Complex32[] response;
            if(matrix&&complex){d.Apply(c);response=c.Cast<Complex32>().ToArray();}
            else if(matrix){d.Apply(a);response=a.Cast<float>().Select(v=>new Complex32(v,0)).ToArray();}
            else if(complex){d.Apply(z);response=z;}
            else{d.Apply(x);response=x.Select(v=>new Complex32(v,0)).ToArray();}
            Assert.All(response,v=>
            {
                Assert.True(float.IsFinite(v.Real)&&float.IsFinite(v.Imag));
                Assert.InRange(v.Real,.375f-2e-6f,.425f+2e-6f);
                if(complex)Assert.InRange(v.Imag,.125f-2e-6f,.15f+2e-6f);
            });
            Assert.InRange(response[center].Real,.375f+1e-5f,.425f-1e-5f);
            Assert.Contains(Enumerable.Range(0,response.Length),i=>i!=center&&response[i].Real>.375f+1e-5f);
            if(complex)
            {
                Assert.InRange(response[center].Imag,.125f+1e-5f,.15f-1e-5f);
                Assert.Contains(Enumerable.Range(0,response.Length),i=>i!=center&&response[i].Imag>.125f+1e-5f);
            }
        }
    }
}
