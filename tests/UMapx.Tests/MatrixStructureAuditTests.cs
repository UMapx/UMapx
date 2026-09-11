using System.Numerics;
using UMapx.Core;
using UMapx.Transform;
using Xunit;
using static UMapx.Tests.NumericAssert;
using static UMapx.Tests.MatrixAuditTests;

namespace UMapx.Tests;

[Trait("Category","Matrix")]
public class MatrixStructureAuditTests
{
    public static IEnumerable<object[]> ShapeCases()
    {
        foreach(bool complex in new[]{false,true})foreach(var shape in new[]{(3,5),(5,3),(4,4)})foreach(var direction in Enum.GetValues<Direction>())yield return new object[]{complex,shape.Item1,shape.Item2,direction};
    }
    [Theory] [MemberData(nameof(ShapeCases))]
    public void SwappingRowsAndColumnsPermutesEveryEntry(bool complex,int rows,int cols,Direction direction)
    {
        var a=(Array)Operand(complex?typeof(Complex32[,]):typeof(float[,]),1,rows,cols);var original=(Array)a.Clone();
        Invoke(typeof(Matrice).GetMethod("Swap",new[]{a.GetType(),typeof(int),typeof(int),typeof(Direction)})!,a,0,1,direction);
        for(int i=0;i<rows;i++)for(int j=0;j<cols;j++)
        {int ii=direction!=Direction.Vertical&&i<2?1-i:i,jj=direction!=Direction.Horizontal&&j<2?1-j:j;Check(Value(original,ii,jj),a.GetValue(i,j)!);}
    }
    [Theory] [MemberData(nameof(ShapeCases))]
    public void MatrixSlicesAndFiniteDifferencesRespectTheSelectedAxes(bool complex,int rows,int cols,Direction direction)
    {
        var a=(Array)Operand(complex?typeof(Complex32[,]):typeof(float[,]),1,rows,cols);
        var slice=(Array)Invoke(typeof(Matrice).GetMethod("Remove",new[]{a.GetType(),typeof(int),typeof(int),typeof(Direction)})!,a,1,2,direction);
        int sr=direction==Direction.Vertical?rows:2,sc=direction==Direction.Horizontal?cols:2;
        Assert.Equal(sr,slice.GetLength(0));Assert.Equal(sc,slice.GetLength(1));
        for(int i=0;i<sr;i++)for(int j=0;j<sc;j++)Check(Value(a,i+(direction==Direction.Vertical?0:1),j+(direction==Direction.Horizontal?0:1)),slice.GetValue(i,j)!);
        // Diff uses Horizontal for adjacent columns, Vertical for adjacent rows.
        var d=(Array)Invoke(typeof(Matrice).GetMethod("Diff",new[]{a.GetType(),typeof(int),typeof(Direction),typeof(bool)})!,a,1,direction,false);
        int dr=rows-(direction==Direction.Horizontal?0:1),dc=cols-(direction==Direction.Vertical?0:1);
        Assert.Equal(dr,d.GetLength(0));Assert.Equal(dc,d.GetLength(1));
        for(int i=0;i<dr;i++)for(int j=0;j<dc;j++)
        {
            Complex expected=direction==Direction.Horizontal?Value(a,i,j+1)-Value(a,i,j):direction==Direction.Vertical?Value(a,i+1,j)-Value(a,i,j):Value(a,i+1,j+1)-Value(a,i+1,j)-Value(a,i,j+1)+Value(a,i,j);
            Check(expected,d.GetValue(i,j)!);
        }
    }
    public static IEnumerable<object[]> MeshCases(){foreach(bool x in new[]{false,true})foreach(bool y in new[]{false,true})foreach(var size in new[]{(3,5),(5,3)})yield return new object[]{x,y,size.Item1,size.Item2};}
    [Theory] [MemberData(nameof(MeshCases))]
    public void MeshEvaluationUsesBothIndependentCoordinates(bool complexX,bool complexY,int nx,int ny)
    {
        var x=Operand(complexX?typeof(Complex32[]):typeof(float[]),1,1,nx);var y=Operand(complexY?typeof(Complex32[]):typeof(float[]),2,1,ny);
        Delegate function=complexX||complexY?(Delegate)new IMeshComplex32((a,b)=>a+2*b):new IMeshFloat((a,b)=>a+2*b);
        var result=(Array)Invoke(typeof(Matrice).GetMethod("Compute",new[]{x.GetType(),y.GetType(),function.GetType()})!,x,y,function);
        Assert.Equal(nx,result.GetLength(0));Assert.Equal(ny,result.GetLength(1));
        for(int i=0;i<nx;i++)for(int j=0;j<ny;j++)Check(Value(x,0,i)+2*Value(y,0,j),result.GetValue(i,j)!);
    }
    [Theory] [InlineData(false)] [InlineData(true)]
    public void ElementMappingAndRowColumnReplacementPreserveCoordinates(bool complex)
    {
        if(complex)
        {
            var a=(Complex32[,])Operand(typeof(Complex32[,]),1);var mapped=a.Compute(z=>z*z+1);
            for(int i=0;i<3;i++)for(int j=0;j<5;j++)Close(Value(a,i,j)*Value(a,i,j)+1,mapped[i,j]);
            var row=a.GetRow(1);var col=a.GetCol(2);for(int j=0;j<5;j++)Close(Value(a,1,j),row[j]);for(int i=0;i<3;i++)Close(Value(a,i,2),col[i]);
            a=a.SetRow(row,0);a=a.SetCol(col,4);for(int j=0;j<4;j++)Close((Complex)row[j],a[0,j]);for(int i=0;i<3;i++)Close((Complex)col[i],a[i,4]);
            var v=row.Compute(z=>z*z+1);for(int i=0;i<v.Length;i++)Close((Complex)row[i]*(Complex)row[i]+1,v[i]);
        }
        else
        {
            var a=(float[,])Operand(typeof(float[,]),1);var mapped=a.Compute(z=>z*z+1);
            for(int i=0;i<3;i++)for(int j=0;j<5;j++)Close(a[i,j]*a[i,j]+1,mapped[i,j]);
            var row=a.GetRow(1);var col=a.GetCol(2);for(int j=0;j<5;j++)Close(a[1,j],row[j]);for(int i=0;i<3;i++)Close(a[i,2],col[i]);
            a=a.SetRow(row,0);a=a.SetCol(col,4);for(int j=0;j<4;j++)Close(row[j],a[0,j]);for(int i=0;i<3;i++)Close(col[i],a[i,4]);
            var v=row.Compute(z=>z*z+1);for(int i=0;i<v.Length;i++)Close(row[i]*row[i]+1,v[i]);
        }
    }
    [Theory] [InlineData(false)] [InlineData(true)]
    public void ArrayExtensionPreservesTheCenteredOriginalAndConstantBorders(bool complex)
    {
        var a=(Array)Operand(complex?typeof(Complex32[,]):typeof(float[,]),1,3,5);
        var extended=(Array)Invoke(typeof(Matrice).GetMethod("Extend",new[]{a.GetType(),typeof(int),typeof(int)})!,a,7,9);
        for(int i=0;i<3;i++)for(int j=0;j<5;j++)Check(Value(a,i,j),extended.GetValue(i+2,j+2)!);
        var v=Array.CreateInstance(complex?typeof(Complex32):typeof(float),5);for(int i=0;i<5;i++)v.SetValue(complex?(object)new Complex32(2,3):2f,i);
        var padded=(Array)Invoke(typeof(Matrice).GetMethod("Extend",new[]{v.GetType(),typeof(int)})!,v,9);foreach(var value in padded)Check(Value(v),value!);
    }
    public static IEnumerable<object[]> StructuredCases(){foreach(string name in new[]{"Exchange","Lehmer","Hilbert","GCD","Stirling"})foreach(int n in new[]{1,3,8})yield return new object[]{name,n};}
    [Theory] [MemberData(nameof(StructuredCases))]
    public void ClassicalMatricesMatchTheirIntegerOrRationalDefinitions(string name,int n)
    {
        foreach(bool second in name=="Stirling"?new[]{false,true}:new[]{false})
        {
            var actual=name=="Stirling"?Matrice.Stirling(n,second):(float[,])Invoke(typeof(Matrice).GetMethod(name,new[]{typeof(int)})!,n);
            var stirling=new long[n,n];stirling[0,0]=1;for(int i=1;i<n;i++)for(int j=1;j<=i;j++)stirling[i,j]=stirling[i-1,j-1]+(second?j:i-1)*stirling[i-1,j];
            for(int i=0;i<n;i++)for(int j=0;j<n;j++)Close(name switch{"Exchange"=>i+j==n-1?1:0,"Lehmer"=>(double)(Math.Min(i,j)+1)/(Math.Max(i,j)+1),"Hilbert"=>1.0/(i+j+1),"GCD"=>(double)BigInteger.GreatestCommonDivisor(i+1,j+1),_=>stirling[i,j]},actual[i,j]);
        }
    }
    [Theory] [InlineData(false)] [InlineData(true)]
    public void MatrixPredicatesRecognizeConstructedExamples(bool complex)
    {
        if(complex)
        {
            Complex32[,] hermitian={{2,new(1,3)},{new(1,-3),4}},skew={{new(0,2),new(1,3)},{new(-1,3),new(0,4)}};
            Assert.True(hermitian.IsSymmetric());Assert.False(hermitian.IsSkewSymmetric());Assert.True(skew.IsSkewSymmetric());Assert.False(skew.IsSymmetric());
            Assert.True(new Complex32[3,5].IsDiagonal());Assert.False(hermitian.IsDiagonal());Assert.True(new Complex32[1,3].IsVector());Assert.False(hermitian.IsVector());Assert.True(hermitian.IsSquare());Assert.False(new Complex32[3,5].IsSquare());
            Assert.True(hermitian.IsEquals((Complex32[,])hermitian.Clone()));Assert.False(hermitian.IsEquals(skew));
        }
        else
        {
            float[,] symmetric={{2,3},{3,4}},skew={{0,3},{-3,0}};
            Assert.True(symmetric.IsSymmetric());Assert.False(symmetric.IsSkewSymmetric());Assert.True(skew.IsSkewSymmetric());Assert.False(skew.IsSymmetric());
            Assert.True(new float[3,5].IsDiagonal());Assert.False(symmetric.IsDiagonal());Assert.True(new float[1,3].IsVector());Assert.False(symmetric.IsVector());Assert.True(symmetric.IsSquare());Assert.False(new float[3,5].IsSquare());
            Assert.True(symmetric.IsEquals((float[,])symmetric.Clone()));Assert.False(symmetric.IsEquals(skew));Assert.True(symmetric.IsNonNegative());Assert.False(skew.IsNonNegative());
        }
    }
    [Theory] [InlineData(false)] [InlineData(true)]
    public void FrequencyMasksPreserveTheRequestedBins(bool complex)
    {
        var matrix=(Array)Operand(complex?typeof(Complex32[,]):typeof(float[,]),1,7,9);var expected=(Array)matrix.Clone();var vector=(Array)Operand(complex?typeof(Complex32[]):typeof(float[]),2,1,9);var vcopy=(Array)vector.Clone();
        var filter=new FrequencyFilter(-2,2);if(complex)filter.Apply((Complex32[])vector);else filter.Apply((float[])vector);
        for(int i=0;i<9;i++)Check(i>=2&&i<=6?Value(vcopy,0,i):Complex.Zero,vector.GetValue(i)!);
        filter.FrequencyRange=new RangeInt(1,2);if(complex)filter.Apply((Complex32[,])matrix);else filter.Apply((float[,])matrix);
        // Matrix radii are quantized to integer bins by this API.
        for(int i=0;i<7;i++)for(int j=0;j<9;j++){int squared=(i-3)*(i-3)+(j-4)*(j-4);Check(squared>=1&&squared<9?Value(expected,i,j):Complex.Zero,matrix.GetValue(i,j)!);}
    }
}
