using System.Numerics;
using UMapx.Core;
using UMapx.Transform;
using UMapx.Window;
using Xunit;
using static UMapx.Tests.NumericAssert;

namespace UMapx.Tests;

[Trait("Category", "WindowTransform")]
public class FramedTransformAuditTests
{
    public static IEnumerable<object[]> WeylMatrixCases(){foreach(string kind in new[]{"Complex","RealFourier","RealHartley"})foreach(var direction in Enum.GetValues<Direction>())yield return new object[]{kind,direction};}
    [Theory] [MemberData(nameof(WeylMatrixCases))]
    public void RectangularWeylTransformsAgreeAndInvertInEveryDirection(string kind,Direction direction)
    {
        const int n=16,m=24,shift=4;var window=new Hamming(shift);var x=Matrix(n,m);var z=new Complex32[n,m];for(int i=0;i<n;i++)for(int j=0;j<m;j++)z[i,j]=new Complex32(x[i,j],(float)Math.Cos(.3*i+.2*j));
        var spectrum=kind=="RealHartley"?SpectrumType.Hartley:SpectrumType.Fourier;
        ITransform direct=kind=="Complex"?new WeylHeisenbergTransform(window,shift,direction):new RealWeylHeisenbergTransform(window,shift,spectrum,direction);
        ITransform fast=kind=="Complex"?new FastWeylHeisenbergTransform(window,shift,direction):new FastRealWeylHeisenbergTransform(window,shift,spectrum,direction);
        var a=direct.Forward(z);var b=fast.Forward(z);Assert.Equal(a.GetLength(0),b.GetLength(0));Assert.Equal(a.GetLength(1),b.GetLength(1));
        for(int i=0;i<a.GetLength(0);i++)for(int j=0;j<a.GetLength(1);j++)Close((Complex)a[i,j],b[i,j],.001);
        foreach(var d in new[]{direct,fast}){var r=d.Backward(d.Forward(z));for(int i=0;i<n;i++)for(int j=0;j<m;j++)Close((Complex)z[i,j],r[i,j],.001);if(kind!="Complex")Close(x,d.Backward(d.Forward(x)),.001);}
    }
    [Theory] [InlineData(3,false)] [InlineData(4,false)] [InlineData(8,false)] [InlineData(3,true)] [InlineData(4,true)] [InlineData(8,true)]
    public void ShortTimeFourierTransformsMatchWindowedBlockDfts(int frame,bool fast)
    {
        var window=new Hamming(frame);var x=Enumerable.Range(0,frame*3).Select(i=>new Complex32((float)Math.Sin(i*.31),(float)Math.Cos(i*.27))).ToArray();
        ITransform d=fast?new FastShortTimeFourierTransform(window):new ShortTimeFourierTransform(window);var actual=d.Forward(x);
        for(int block=0;block<3;block++)for(int k=0;k<frame;k++)
        {
            // The implementation regularizes each window coefficient by adding 0.001.
            Complex expected=0;for(int j=0;j<frame;j++)expected+=(Complex)x[block*frame+j]*(.001+.53836-.46164*Math.Cos(2*Math.PI*j/(frame-1)))*Complex.Exp(-2*Math.PI*Complex.ImaginaryOne*j*k/frame)/Math.Sqrt(frame);
            Close(expected,actual[block*frame+k],1e-4);
        }
        var restored=d.Backward(actual);for(int i=0;i<x.Length;i++)Close((Complex)x[i],restored[i],1e-4);
        Assert.Throws<ArgumentException>(()=>d.Forward(new Complex32[frame+1]));
    }

    [Theory] [InlineData(12,3,false)] [InlineData(16,4,false)] [InlineData(32,8,false)] [InlineData(12,3,true)] [InlineData(16,4,true)] [InlineData(32,8,true)]
    public void ZakTransformMatchesItsFiniteSumAndInverse(int n,int m,bool fast)
    {
        var x=Enumerable.Range(0,n).Select(i=>new Complex32((float)Math.Sin(i*.31),(float)Math.Cos(i*.27))).ToArray();
        IZakTransform d=fast?new FastZakTransform(m):new ZakTransform(m);var actual=d.Forward(x);var real=d.Forward(x.Select(z=>z.Real).ToArray());int l=n/m;
        Assert.Equal(n,actual.GetLength(0));Assert.Equal(l,actual.GetLength(1));
        for(int row=0;row<n;row++)for(int k=0;k<l;k++)
        {
            Complex expected=0,expectedReal=0;
            for(int j=0;j<l;j++){int index=((row-m*j)%n+n)%n;var phase=Complex.Exp(-Complex.ImaginaryOne*2*Math.PI*j*k/l);expected+=(Complex)x[index]*phase;expectedReal+=x[index].Real*phase;}
            Close(expected,actual[row,k],1e-4);Close(expectedReal,real[row,k],1e-4);
        }
        var restored=d.Backward(actual);for(int i=0;i<n;i++)Close((Complex)x[i],restored[i],1e-4);
        Assert.Throws<ArgumentException>(()=>d.Forward(new float[n+1]));
    }

    public static IEnumerable<object[]> WeylCases()
    {
        foreach(var size in new[]{(16,4),(32,8),(48,8)})foreach(string kind in new[]{"ComplexCompare","ComplexDirect","ComplexFast","RealFourierCompare","RealHartleyCompare","RealFourierDirect","RealFourierFast","RealHartleyDirect","RealHartleyFast"})yield return new object[]{size.Item1,size.Item2,kind};
    }
    [Theory] [MemberData(nameof(WeylCases))]
    public void WeylHeisenbergFastAndDirectBasesAgreeAndReconstruct(int n,int m,string kind)
    {
        // A full-period Hamming window has a sparse spectrum and can produce a deficient frame.
        var window=new Hamming(m);var x=Enumerable.Range(0,n).Select(i=>(float)Math.Sin(i*.31)).ToArray();var z=x.Select((v,i)=>new Complex32(v,(float)Math.Cos(i*.19))).ToArray();
        bool real=kind.StartsWith("Real");var spectrum=kind.Contains("Hartley")?SpectrumType.Hartley:SpectrumType.Fourier;
        ITransform direct=real?new RealWeylHeisenbergTransform(window,m,spectrum):new WeylHeisenbergTransform(window,m),fast=real?new FastRealWeylHeisenbergTransform(window,m,spectrum):new FastWeylHeisenbergTransform(window,m);
        if(!kind.EndsWith("Compare"))
        {
            var d=kind.EndsWith("Fast")?fast:direct;var r=d.Backward(d.Forward(z));for(int i=0;i<n;i++)Close((Complex)z[i],r[i],.002);
            if(real)Close(x,d.Backward(d.Forward(x)),.002);return;
        }
        var expected=direct.Forward(z);var actual=fast.Forward(z);Assert.Equal(expected.Length,actual.Length);
        for(int i=0;i<actual.Length;i++)Close((Complex)expected[i],actual[i],.001);
        if(real)Close(direct.Forward(x),fast.Forward(x),.001);
    }
}
