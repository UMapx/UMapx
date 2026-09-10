using System.Numerics;
using UMapx.Core;
using Xunit;
using static UMapx.Tests.NumericAssert;

namespace UMapx.Tests;

[Trait("Category","Matrix")]
public class MatrixFilterAuditTests
{
    public static IEnumerable<object[]> ConvolutionCases(){foreach(bool a in new[]{false,true})foreach(bool b in new[]{false,true})foreach(bool normalized in new[]{false,true})foreach(string direction in new[]{"Matrix","Horizontal","Vertical","Both"})yield return new object[]{a,b,normalized,direction};}
    static Complex[,] Correlate(Complex[,] a,Complex[,] kernel,bool normalized)
    {
        int h=a.GetLength(0),w=a.GetLength(1),kh=kernel.GetLength(0),kw=kernel.GetLength(1);var r=new Complex[h,w];
        for(int y=0;y<h;y++)for(int x=0;x<w;x++){Complex sum=0,weight=0;for(int dy=0;dy<kh;dy++)for(int dx=0;dx<kw;dx++){int yy=y+dy-kh/2,xx=x+dx-kw/2;if(yy<0||yy>=h||xx<0||xx>=w)continue;sum+=a[yy,xx]*kernel[dy,dx];weight+=kernel[dy,dx];}r[y,x]=normalized?sum/weight:sum;}return r;
    }
    [Theory] [MemberData(nameof(ConvolutionCases))]
    public void MixedMatrixConvolutionOverloadsMatchIndependentNeighborhoodSums(bool complexData,bool complexKernel,bool normalized,string direction)
    {
        var a=(Array)MatrixAuditTests.Operand(complexData?typeof(Complex32[,]):typeof(float[,]),1,7,9);var k=(Array)MatrixAuditTests.Operand(direction=="Matrix"?(complexKernel?typeof(Complex32[,]):typeof(float[,])):(complexKernel?typeof(Complex32[]):typeof(float[])),2,3,3);
        object[] args=direction=="Matrix"?new object[]{a,k,normalized}:new object[]{a,k,Enum.Parse<Direction>(direction),normalized};var actual=(Array)MatrixAuditTests.Invoke(typeof(Matrice).GetMethod("Conv",args.Select(v=>v.GetType()).ToArray())!,args);
        var expected=new Complex[7,9];for(int y=0;y<7;y++)for(int x=0;x<9;x++)expected[y,x]=MatrixAuditTests.Value(a,y,x);
        if(direction=="Matrix"){var kernel=new Complex[3,3];for(int y=0;y<3;y++)for(int x=0;x<3;x++)kernel[y,x]=MatrixAuditTests.Value(k,y,x);expected=Correlate(expected,kernel,normalized);}
        else{foreach(bool vertical in direction=="Both"?new[]{false,true}:new[]{direction=="Vertical"}){var kernel=new Complex[vertical?3:1,vertical?1:3];for(int i=0;i<3;i++)kernel[vertical?i:0,vertical?0:i]=MatrixAuditTests.Value(k,0,i);expected=Correlate(expected,kernel,normalized);}}
        for(int y=0;y<7;y++)for(int x=0;x<9;x++)MatrixAuditTests.Check(expected[y,x],actual.GetValue(y,x)!,3e-4);
    }
    public static IEnumerable<object[]> StructureCases(){foreach(bool complex in new[]{false,true})foreach(string name in new[]{"Vander","Toeplitz","Hankeli","Hankel","Circulant","Symmetric","Companion","Diag"})yield return new object[]{complex,name};}
    [Theory] [MemberData(nameof(StructureCases))]
    public void StructuredMatricesMatchTheirEntryDefinitions(bool complex,string name)
    {
        var v=MatrixAuditTests.Operand(complex?typeof(Complex32[]):typeof(float[]),1,1,6);var actual=(Array)MatrixAuditTests.Invoke(typeof(Matrice).GetMethod(name,new[]{v.GetType()})!,v);int n=name=="Hankel"?3:6;Assert.Equal(n,actual.GetLength(0));Assert.Equal(n,actual.GetLength(1));
        for(int i=0;i<n;i++)for(int j=0;j<n;j++)
        {
            Complex expected=name switch{"Vander"=>Complex.Pow(MatrixAuditTests.Value(v,0,i),j),"Toeplitz" or "Symmetric"=>MatrixAuditTests.Value(v,0,Math.Abs(i-j)),"Hankeli"=>i+j<6?MatrixAuditTests.Value(v,0,i+j):0,"Hankel"=>MatrixAuditTests.Value(v,0,i+j),"Circulant"=>MatrixAuditTests.Value(v,0,(j-i+6)%6),"Diag"=>i==j?MatrixAuditTests.Value(v,0,i):0,_=>j==5?-MatrixAuditTests.Value(v,0,i):i==j+1?1:0};MatrixAuditTests.Check(expected,actual.GetValue(i,j)!);
        }
    }
    [Theory] [InlineData(1)] [InlineData(3)] [InlineData(5)] [InlineData(9)]
    public void MagicSquaresHaveTheCorrectSetAndEveryLineHasTheSameSum(int n)
    {
        var a=Matrice.Magic(n);Assert.Equal(Enumerable.Range(1,n*n).Select(x=>(float)x),a.Cast<float>().OrderBy(x=>x));double expected=n*(n*n+1)/2.0;
        for(int i=0;i<n;i++){Close(expected,Enumerable.Range(0,n).Sum(j=>(double)a[i,j]));Close(expected,Enumerable.Range(0,n).Sum(j=>(double)a[j,i]));}Close(expected,Enumerable.Range(0,n).Sum(i=>(double)a[i,i]));Close(expected,Enumerable.Range(0,n).Sum(i=>(double)a[i,n-i-1]));
    }
    public static IEnumerable<object[]> RotationCases(){foreach(bool complex in new[]{false,true})foreach(var mode in Enum.GetValues<InterpolationMode>())foreach(float angle in new[]{0f,180f})yield return new object[]{complex,mode,angle};}
    [Theory] [MemberData(nameof(RotationCases))]
    public void MatrixRotationsAtExactHalfTurnsMatchIndexReversal(bool complex,InterpolationMode mode,float angle)
    {
        var a=(Array)MatrixAuditTests.Operand(complex?typeof(Complex32[,]):typeof(float[,]),1,5,7);var result=(Array)MatrixAuditTests.Invoke(typeof(Matrice).GetMethod("Rotate",new[]{a.GetType(),typeof(float),typeof(InterpolationMode)})!,a,angle,mode);
        for(int y=0;y<5;y++)for(int x=0;x<7;x++)MatrixAuditTests.Check(MatrixAuditTests.Value(a,angle==0?y:4-y,angle==0?x:6-x),result.GetValue(y,x)!,2e-4);
    }
}

[CollectionDefinition("SIMD audit",DisableParallelization=true)]
public class SimdAuditCollection { }

[Collection("SIMD audit")]
[Trait("Category","Matrix")]
public class SimdAuditTests
{
    [Theory] [InlineData(false,false)] [InlineData(true,false)] [InlineData(false,true)] [InlineData(true,true)]
    public void SimdProductsMatchIndependentAccumulationIncludingTailElements(bool complexA,bool complexB)
    {
        var a=(Array)MatrixAuditTests.Operand(complexA?typeof(Complex32[,]):typeof(float[,]),1,13,37);var b=(Array)MatrixAuditTests.Operand(complexB?typeof(Complex32[,]):typeof(float[,]),2,37,19);bool original=Globals.SIMD;
        try{Globals.SIMD=true;var actual=(Array)MatrixAuditTests.Invoke(typeof(Matrice).GetMethod("Dot",new[]{a.GetType(),b.GetType()})!,a,b);for(int y=0;y<13;y++)for(int x=0;x<19;x++){Complex expected=0;for(int k=0;k<37;k++)expected+=MatrixAuditTests.Value(a,y,k)*MatrixAuditTests.Value(b,k,x);MatrixAuditTests.Check(expected,actual.GetValue(y,x)!,1e-4);}}
        finally{Globals.SIMD=original;}
    }
}
