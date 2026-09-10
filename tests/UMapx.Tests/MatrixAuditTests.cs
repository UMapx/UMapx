using System.Numerics;
using System.Reflection;
using UMapx.Core;
using Xunit;
using static UMapx.Tests.NumericAssert;

namespace UMapx.Tests;

[Trait("Category", "Matrix")]
public class MatrixAuditTests
{
    static readonly Type[] NumericTypes = { typeof(float), typeof(Complex32), typeof(float[]), typeof(Complex32[]), typeof(float[,]), typeof(Complex32[,]) };
    static readonly MethodInfo[] Arithmetic = typeof(Matrice).GetMethods(BindingFlags.Public|BindingFlags.Static)
        .Where(m => new[]{"Add","Sub","Mul","Div","Pow"}.Contains(m.Name) && m.GetParameters().Length==2
            && m.GetParameters().All(p=>NumericTypes.Contains(p.ParameterType))).OrderBy(m=>m.ToString()).ToArray();

    public static IEnumerable<object[]> ArithmeticCases() => Arithmetic.Select(m=>new object[]{m.ToString()!});

    internal static object Operand(Type type,int seed,int rows=3,int columns=5)
    {
        object Scalar(Type t,int i,int j) => t==typeof(float)?(object)(float)(.7+.17*i+.11*j+.23*seed):new Complex32((float)(.7+.17*i+.11*j+.23*seed),(float)(.2+.07*i-.03*j));
        if(!type.IsArray)return Scalar(type,0,0);
        var e=type.GetElementType()!;
        var result=type.GetArrayRank()==1?Array.CreateInstance(e,columns):Array.CreateInstance(e,rows,columns);
        for(int i=0;i<(result.Rank==1?1:rows);i++)for(int j=0;j<columns;j++)
            if(result.Rank==1)result.SetValue(Scalar(e,0,j),j);else result.SetValue(Scalar(e,i,j),i,j);
        return result;
    }
    internal static Complex Value(object v,int i=0,int j=0)
    {
        if(v is Array a)v=a.Rank==1?a.GetValue(j)!:a.GetValue(i,j)!;
        return v is Complex32 z?(Complex)z:new Complex(Convert.ToDouble(v),0);
    }
    internal static void Check(Complex expected,object actual,double tolerance=2e-5)
    {
        if(actual is Complex32 z)Close(expected,z,tolerance,tolerance);
        else {Close(0,expected.Imaginary,tolerance);Close(expected.Real,Convert.ToDouble(actual),tolerance,tolerance);}
    }
    internal static object Invoke(MethodInfo method, params object[] args)
    {
        try{return method.Invoke(null,args)!;}
        catch(TargetInvocationException e) when(e.InnerException!=null){System.Runtime.ExceptionServices.ExceptionDispatchInfo.Capture(e.InnerException).Throw();throw;}
    }

    [Theory] [MemberData(nameof(ArithmeticCases))]
    public void EveryElementwiseArithmeticOverloadMatchesScalarArithmetic(string signature)
    {
        var m=Arithmetic.Single(m=>m.ToString()==signature);var p=m.GetParameters();var a=Operand(p[0].ParameterType,1);var b=Operand(p[1].ParameterType,2);
        var actual=(Array)Invoke(m,a,b);int rows=actual.Rank==1?1:actual.GetLength(0),cols=actual.GetLength(actual.Rank-1);
        Assert.Equal(5,cols);if(actual.Rank==2)Assert.Equal(3,rows);
        for(int i=0;i<rows;i++)for(int j=0;j<cols;j++)
        {
            Complex x=Value(a,i,j),y=Value(b,i,j);
            Complex expected=m.Name switch {"Add"=>x+y,"Sub"=>x-y,"Mul"=>x*y,"Div"=>x/y,"Pow"=>Complex.Pow(x,y),_=>throw new Exception()};
            Check(expected,actual.Rank==1?actual.GetValue(j)!:actual.GetValue(i,j)!);
        }
    }

    public static IEnumerable<object[]> ProductCases()
    {
        foreach(string op in new[]{"Dot","Kronecker"})foreach(bool ac in new[]{false,true})foreach(bool bc in new[]{false,true})yield return new object[]{op,ac,bc};
    }
    [Theory] [MemberData(nameof(ProductCases))]
    public void MixedMatrixProductsMatchIndependentComplexArithmetic(string operation,bool complexA,bool complexB)
    {
        var a=Operand(complexA?typeof(Complex32[,]):typeof(float[,]),1,3,5);
        var b=Operand(complexB?typeof(Complex32[,]):typeof(float[,]),2,5,2);
        var actual=(Array)Invoke(typeof(Matrice).GetMethod(operation,new[]{a.GetType(),b.GetType()})!,a,b);
        int rows=operation=="Dot"?3:15,cols=operation=="Dot"?2:10;
        Assert.Equal(rows,actual.GetLength(0));Assert.Equal(cols,actual.GetLength(1));
        for(int i=0;i<rows;i++)for(int j=0;j<cols;j++)
        {
            Complex expected=0;
            if(operation=="Dot")for(int k=0;k<5;k++)expected+=Value(a,i,k)*Value(b,k,j);
            else expected=Value(a,i/5,j/2)*Value(b,i%5,j%2);
            Check(expected,actual.GetValue(i,j)!,1e-4);
        }
    }

    public static IEnumerable<object[]> DiagonalCases()
    {
        foreach(bool left in new[]{false,true})foreach(bool inverse in new[]{false,true})foreach(bool mc in new[]{false,true})foreach(bool vc in new[]{false,true})yield return new object[]{left,inverse,mc,vc};
    }
    [Theory] [MemberData(nameof(DiagonalCases))]
    public void DiagonalMultiplicationScalesTheDocumentedRowsOrColumns(bool left,bool inverse,bool complexMatrix,bool complexVector)
    {
        var a=Operand(complexMatrix?typeof(Complex32[,]):typeof(float[,]),1,3,5);
        var v=Operand(complexVector?typeof(Complex32[]):typeof(float[]),2,1,left?3:5);
        var args=left?new[]{v,a,(object)inverse}:new[]{a,v,(object)inverse};
        var actual=(Array)Invoke(typeof(Matrice).GetMethod("Dot",args.Select(a=>a.GetType()).ToArray())!,args);
        for(int i=0;i<3;i++)for(int j=0;j<5;j++)Check(inverse?Value(a,i,j)/Value(v,0,left?i:j):Value(a,i,j)*Value(v,0,left?i:j),actual.GetValue(i,j)!);
    }

    [Theory] [InlineData(2)] [InlineData(5)] [InlineData(17)]
    public void RealStatisticsMatchSampleDefinitions(int n)
    {
        var x=Enumerable.Range(0,n).Select(i=>(float)(Math.Sin(i*.7)+i*.2)).ToArray();double mean=x.Average(v=>(double)v),variance=x.Sum(v=>Math.Pow(v-mean,2))/(n-1);
        Close(x.Sum(v=>(double)v),Matrice.Sum(x));Close(mean,Matrice.Mean(x));Close(variance,Matrice.Var(x));Close(Math.Sqrt(variance),Matrice.StnDev(x));Close(variance,Matrice.Cov(x));
        var sorted=x.OrderBy(v=>v).ToArray();Close(sorted[0],Matrice.Min(x));Close(sorted[^1],Matrice.Max(x));
        var a=new float[n,3];for(int i=0;i<n;i++)for(int j=0;j<3;j++)a[i,j]=x[i]*(j+1)+j;
        var means=a.Mean();var vars=a.Var();var deviations=a.StnDev();var sums=a.Sum();var cov=a.Cov();
        for(int j=0;j<3;j++){Close(mean*(j+1)+j,means[j]);Close(variance*(j+1)*(j+1),vars[j]);Close(Math.Sqrt(variance)*(j+1),deviations[j]);Close((mean*(j+1)+j)*n,sums[j]);for(int k=0;k<3;k++)Close(variance*(j+1)*(k+1),cov[j,k]);}
        var p=Enumerable.Range(1,n).Select(i=>(float)(2.0*i/(n*(n+1)))).ToArray();Close(-p.Sum(v=>v*Math.Log2(v)),Matrice.Entropy(p));
        var normalized=x.Normalized();for(int i=0;i<n;i++)Close((x[i]-sorted[0])/(sorted[^1]-sorted[0]),normalized[i]);
    }

    [Theory] [InlineData("Variance")] [InlineData("CovarianceVector")] [InlineData("CovarianceMatrix")] [InlineData("Norm")]
    public void ComplexStatisticsRespectHermitianInnerProducts(string operation)
    {
        var x=new[]{new Complex32(1,2),new Complex32(-2,1),new Complex32(3,-4),new Complex32(2,3)};
        Complex mean=x.Aggregate(Complex.Zero,(s,v)=>s+(Complex)v)/x.Length;
        double variance=x.Sum(v=>Complex.Abs((Complex)v-mean)*Complex.Abs((Complex)v-mean))/(x.Length-1);
        if(operation=="Variance")Check(variance,x.Var());
        else if(operation=="CovarianceVector")Check(variance,x.Cov());
        else if(operation=="Norm")Check(Math.Sqrt(x.Sum(v=>Complex.Abs((Complex)v)*Complex.Abs((Complex)v))),x.Abs());
        else {var a=new Complex32[x.Length,2];for(int i=0;i<x.Length;i++){a[i,0]=x[i];a[i,1]=new Complex32(0,1)*x[i];}var c=a.Cov();Check(variance,c[0,0]);Check(variance,c[1,1]);Check(Complex.ImaginaryOne*variance,c[0,1]);Check(-Complex.ImaginaryOne*variance,c[1,0]);}
    }

    public static IEnumerable<object[]> SpatialCases()
    {
        foreach(bool complex in new[]{false,true})foreach(string op in new[]{"Transpose","FlipHorizontal","FlipVertical","FlipBoth","Shift","Crop","Merge","Reshape","DiffHorizontal","DiffVertical","DiffBoth"})foreach(var size in new[]{(3,5),(5,3),(4,4)})yield return new object[]{op,complex,size.Item1,size.Item2};
    }
    [Theory] [MemberData(nameof(SpatialCases))]
    public void RectangularArrayOperationsMatchIndexDefinitions(string operation,bool complex,int rows,int cols)
    {
        var a=(Array)Operand(complex?typeof(Complex32[,]):typeof(float[,]),1,rows,cols);Array result;
        object Call(string name,params object[] tail){var args=new[]{(object)a}.Concat(tail).ToArray();return Invoke(typeof(Matrice).GetMethod(name,args.Select(x=>x.GetType()).ToArray())!,args);}
        if(operation=="Reshape")
        {
            result=(Array)Call("Reshape",rows*cols);for(int i=0;i<rows;i++)for(int j=0;j<cols;j++)Check(Value(a,i,j),result.GetValue(j*rows+i)!);
            var back=(Array)Invoke(typeof(Matrice).GetMethod("Reshape",new[]{result.GetType(),typeof(int)})!,result,rows);
            for(int i=0;i<rows;i++)for(int j=0;j<cols;j++)Check(Value(a,i,j),back.GetValue(i,j)!);return;
        }
        if(operation=="Transpose")result=(Array)Call("Transpose");
        else if(operation.StartsWith("Flip"))result=(Array)Call("Flip",Enum.Parse<Direction>(operation[4..]));
        else if(operation.StartsWith("Diff"))result=(Array)Call("Diff",1,Enum.Parse<Direction>(operation[4..]),false);
        else if(operation=="Shift")result=(Array)Call("Shift",1,-2);
        else if(operation=="Crop")result=(Array)Call("Crop",1,1,2,2,true);
        else
        {
            // A constant patch isolates placement from the independent bicubic resampling defect.
            var patch=Array.CreateInstance(a.GetType().GetElementType()!,2,2);
            for(int i=0;i<2;i++)for(int j=0;j<2;j++)patch.SetValue(complex?(object)new Complex32(2,.5f):2f,i,j);
            result=(Array)Call("Merge",patch,1,1,2,2);
        }
        for(int i=0;i<result.GetLength(0);i++)for(int j=0;j<result.GetLength(1);j++)
        {
            Complex expected=operation switch
            {
                "Transpose"=>Value(a,j,i),"FlipHorizontal"=>Value(a,i,cols-1-j),"FlipVertical"=>Value(a,rows-1-i,j),"FlipBoth"=>Value(a,rows-1-i,cols-1-j),
                "Shift"=>Value(a,(i-1+rows)%rows,(j+2)%cols),"Crop"=>Value(a,i+1,j+1),
                "Merge"=>i>=1&&i<3&&j>=1&&j<3?new Complex(2,complex?.5:0):Value(a,i,j),
                "DiffHorizontal"=>Value(a,i,j+1)-Value(a,i,j),"DiffVertical"=>Value(a,i+1,j)-Value(a,i,j),
                "DiffBoth"=>Value(a,i+1,j+1)-Value(a,i+1,j)-Value(a,i,j+1)+Value(a,i,j),_=>throw new Exception()
            };Check(expected,result.GetValue(i,j)!);
        }
    }

    [Theory] [InlineData(false,false)] [InlineData(false,true)] [InlineData(true,false)] [InlineData(true,true)]
    public void InverseSolveAndDeterminantSatisfyIndependentEquations(bool complex,bool pivot)
    {
        float[,] real=pivot?new float[,]{{0,2,1},{2,3,-1},{1,1,4}}:new float[,]{{4,2,1},{2,3,-1},{1,1,4}};
        var a=complex?(Array)real.ToComplex():(Array)real;
        if(complex){var c=(Complex32[,])a;c[0,2]+=new Complex32(0,.25f);c[2,0]-=new Complex32(0,.5f);}
        var inverse=(Array)Invoke(typeof(Matrice).GetMethod("Invert",new[]{a.GetType()})!,a);
        for(int i=0;i<3;i++)for(int j=0;j<3;j++){Complex sum=0;for(int k=0;k<3;k++)sum+=Value(a,i,k)*Value(inverse,k,j);Check(i==j?1:0,(Complex32)sum,2e-4);}
        var b=Operand(complex?typeof(Complex32[]):typeof(float[]),2,1,3);
        var solution=(Array)Invoke(typeof(Matrice).GetMethod("Solve",new[]{a.GetType(),b.GetType()})!,a,b);
        for(int i=0;i<3;i++){Complex sum=0;for(int j=0;j<3;j++)sum+=Value(a,i,j)*Value(solution,0,j);Check(Value(b,0,i),(Complex32)sum,2e-4);}
        Complex determinant=Value(a,0,0)*(Value(a,1,1)*Value(a,2,2)-Value(a,1,2)*Value(a,2,1))-Value(a,0,1)*(Value(a,1,0)*Value(a,2,2)-Value(a,1,2)*Value(a,2,0))+Value(a,0,2)*(Value(a,1,0)*Value(a,2,1)-Value(a,1,1)*Value(a,2,0));
        Check(determinant,Invoke(typeof(Matrice).GetMethod("Det",new[]{a.GetType()})!,a),2e-4);
    }

    [Theory] [InlineData(false,false)] [InlineData(false,true)] [InlineData(true,false)] [InlineData(true,true)]
    public void FilteringMatchesCenteredCorrelationWithTruncatedSupport(bool complex,bool normalize)
    {
        // Conv uses correlation orientation; this audit records that API convention.
        var v=Operand(complex?typeof(Complex32[]):typeof(float[]),1,1,9);var k=new float[]{1,2,4};
        var result=(Array)Invoke(typeof(Matrice).GetMethod("Conv",new[]{v.GetType(),typeof(float[]),typeof(bool)})!,v,k,normalize);
        for(int i=0;i<9;i++){Complex sum=0;double weight=0;for(int j=0;j<3;j++)if(i+j-1>=0&&i+j-1<9){sum+=Value(v,0,i+j-1)*k[j];weight+=k[j];}Check(normalize?sum/weight:sum,result.GetValue(i)!);}
    }

    [Theory] [InlineData(MorphologyMode.Median)] [InlineData(MorphologyMode.Erosion)] [InlineData(MorphologyMode.Dilatation)]
    public void MorphologyMatchesSortingOfReplicatedEdgeNeighborhoods(MorphologyMode mode)
    {
        // Public Matrice arguments are window sizes; the internal filter takes half sizes.
        var a=Matrix(5,7);var actual=a.Morph(3,5,mode);
        for(int i=0;i<5;i++)for(int j=0;j<7;j++)
        {
            var v=new List<float>();for(int di=-1;di<=1;di++)for(int dj=-2;dj<=2;dj++)v.Add(a[Math.Clamp(i+di,0,4),Math.Clamp(j+dj,0,6)]);v.Sort();
            Close(mode==MorphologyMode.Erosion?v[0]:mode==MorphologyMode.Dilatation?v[^1]:v[v.Count/2],actual[i,j]);
        }
        var row=Enumerable.Range(0,7).Select(j=>a[2,j]).ToArray();var filtered=row.Morph(5,mode);
        for(int j=0;j<7;j++){var v=Enumerable.Range(-2,5).Select(d=>row[Math.Clamp(j+d,0,6)]).OrderBy(x=>x).ToArray();Close(mode==MorphologyMode.Erosion?v[0]:mode==MorphologyMode.Dilatation?v[^1]:v[2],filtered[j]);}
    }

    [Theory] [InlineData(false)] [InlineData(true)]
    public void LocalMeanHasTheCorrectInteriorImpulseResponse(bool complex)
    {
        var x=new float[]{0,0,0,3,0,0,0,0,0};var v=complex?(object)x.ToComplex():x;
        var result=(Array)Invoke(typeof(Matrice).GetMethod("Mean",new[]{v.GetType(),typeof(int)})!,v,3);
        for(int i=1;i<x.Length-1;i++)Check((x[i-1]+x[i]+x[i+1])/3.0,result.GetValue(i)!);
    }

    [Theory] [InlineData(MorphologyMode.Median)] [InlineData(MorphologyMode.Dilatation)]
    public void MorphologyUsesZeroBasedRanksInAThreeSampleWindow(MorphologyMode mode)
    {
        float[] input={1,2,3,4,5};var actual=input.Morph(3,mode);Close(mode==MorphologyMode.Median?3:4,actual[2]);
    }
}
