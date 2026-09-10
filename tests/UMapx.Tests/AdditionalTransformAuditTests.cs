using System.Numerics;
using System.Text.Json;
using UMapx.Core;
using UMapx.Distribution;
using UMapx.Transform;
using UMapx.Wavelet;
using Xunit;
using static UMapx.Tests.NumericAssert;

namespace UMapx.Tests;

[Trait("Category","Transform")]
public class AdditionalTransformAuditTests
{
    public static IEnumerable<object[]> HankelCases()
    {
        using var stream=typeof(AdditionalTransformAuditTests).Assembly.GetManifestResourceStream("UMapx.Tests.Data.hankel.json");using var data=JsonDocument.Parse(stream!);
        foreach(var item in data.RootElement.EnumerateArray())yield return new object[]{item.GetProperty("order").GetInt32(),item.GetProperty("size").GetInt32(),item.GetProperty("values").GetRawText()};
    }
    [Theory] [MemberData(nameof(HankelCases))]
    public void HankelMatrixMatchesHighPrecisionBesselZeros(int order,int size,string values)
    {
        using var data=JsonDocument.Parse(values);var actual=HankelTransform.Matrix(size,order);
        for(int i=0;i<size;i++)for(int j=0;j<size;j++)Close(data.RootElement[i][j].GetDouble(),actual[i,j],5e-4,5e-4);
        var x=Enumerable.Range(0,size).Select(i=>(float)Math.Sin(.7*i)).ToArray();var d=new HankelTransform(order);var result=d.Forward(x);
        for(int i=0;i<size;i++){double sum=0;for(int j=0;j<size;j++)sum+=data.RootElement[i][j].GetDouble()*x[j];Close(sum,result[i],.001,.001);}
    }
    [Theory] [InlineData(false,false)] [InlineData(false,true)] [InlineData(true,false)] [InlineData(true,true)]
    public void WaveletFiltersWithZeroFactorPreserveTheirInput(bool complex,bool matrix)
    {
        var type=matrix?(complex?typeof(Complex32[,]):typeof(float[,])):(complex?typeof(Complex32[]):typeof(float[]));var input=(Array)MatrixAuditTests.Operand(type,1,16,16);
        // Bilateral grids quantize intensities on [0,1]. Keep magnitudes inside that domain.
        for(int y=0;y<(matrix?16:1);y++)for(int x=0;x<16;x++){var v=MatrixAuditTests.Value(input,y,x)/8;object sample=complex?(object)(Complex32)v:(float)v.Real;if(matrix)input.SetValue(sample,y,x);else input.SetValue(sample,x);}
        foreach(IFilter filter in new IFilter[]{new WaveletFilter(new WaveletDecomposition(WaveletPacket.D4,2),0),new EdgeAvoidingWaveletFilter(new EdgeAvoidingWaveletDecomposition(levels:2),0)})
        {
            var actual=(Array)input.Clone();typeof(IFilter).GetMethod("Apply",new[]{type})!.Invoke(filter,new object[]{actual});var expected=input.Cast<object>().ToArray();var result=actual.Cast<object>().ToArray();
            for(int i=0;i<expected.Length;i++)MatrixAuditTests.Check(MatrixAuditTests.Value(expected[i]),result[i],2e-4);
        }
    }
    [Theory] [InlineData(false,false)] [InlineData(false,true)] [InlineData(true,false)] [InlineData(true,true)]
    public void MultichannelWrappersPreserveChannelSeparationAndInvertTransforms(bool complex,bool matrix)
    {
        var element=matrix?(complex?typeof(Complex32[,]):typeof(float[,])):(complex?typeof(Complex32[]):typeof(float[]));var input=Array.CreateInstance(element,3);
        for(int c=0;c<3;c++)input.SetValue(MatrixAuditTests.Operand(element,c,8,16),c);
        var d=new MultidimensionalTransform(new CosineTransform());var p=new MultidimensionalPyramidTransform(new WaveletDecomposition(WaveletPacket.D4,2));
        foreach(object transform in new object[]{d,p})
        {
            var forward=(Array)transform.GetType().GetMethod("Forward",new[]{input.GetType()})!.Invoke(transform,new object[]{input})!;
            var restored=(Array)transform.GetType().GetMethod("Backward",new[]{forward.GetType()})!.Invoke(transform,new object[]{forward})!;Assert.Equal(3,restored.Length);
            for(int c=0;c<3;c++){var a=((Array)input.GetValue(c)!).Cast<object>().ToArray();var b=((Array)restored.GetValue(c)!).Cast<object>().ToArray();Assert.Equal(a.Length,b.Length);for(int i=0;i<a.Length;i++)MatrixAuditTests.Check(MatrixAuditTests.Value(a[i]),b[i],2e-4);}
        }
        var filter=new MultidimensionalFilter(new ThresholdFilter(100));filter.GetType().GetMethod("Apply",new[]{input.GetType()})!.Invoke(filter,new object[]{input});
        for(int c=0;c<3;c++)Assert.All(((Array)input.GetValue(c)!).Cast<object>(),v=>MatrixAuditTests.Check(0,v));
    }
    [Theory] [InlineData(-2f)] [InlineData(-.5f)] [InlineData(.5f)] [InlineData(2f)]
    public void TimeFrequencyKernelsMatchTheirFourierPair(float tau)
    {
        var choi=new ChoiWilliams(.3f);var cone=new ConeShape(.3f);
        foreach(float eta in new[]{-2f,-.1f,0f,.5f,2f})
        {
            Close(Math.Exp(-.3f*eta*eta*tau*tau),choi.Function(eta,tau));double v=Math.PI*eta*tau;Close((v==0?1:Math.Sin(v)/v)*Math.Exp(-2*Math.PI*.3f*tau*tau),cone.Function(eta,tau));
        }
        // Inverse Fourier transform of sinc(pi*eta*tau): an even rectangular pulse of width |tau|.
        foreach(float t in new[]{0f,.1f,3f})Close(Math.Abs(t)<Math.Abs(tau)/2?Math.Exp(-2*Math.PI*.3f*tau*tau)/Math.Abs(tau):0,cone.Distribution(t,tau));
    }
}
