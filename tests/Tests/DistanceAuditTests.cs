using System.Numerics;
using UMapx.Core;
using UMapx.Distance;
using Xunit;
using static UMapx.Tests.NumericAssert;

namespace UMapx.Tests;

[Trait("Category", "Distance")]
public class DistanceAuditTests
{
    private static readonly string[] NumericNames = { "Euclidean", "SquareEuclidean", "Manhattan", "Chebyshev", "Minkowski", "Canberra", "BrayCurtis", "Hellinger", "Angular", "Cosine" };
    private static readonly string[] BooleanNames = { "Dice", "Jaccard", "Kulczynski", "RusselRao", "SokalMichener", "SokalSneath", "Yule" };

    public static IEnumerable<object[]> NumericCases() => NumericNames.SelectMany(name => new[] { 1, 2, 5, 17 }.Select(n => new object[] { name, n }));
    public static IEnumerable<object[]> BooleanCases() => BooleanNames.SelectMany(name => new[] { 0, 1, 2 }.Select(seed => new object[] { name, seed }));

    private static DistanceBase Create(string name) => name switch
    {
        "Minkowski" => new Minkowski(3.5f),
        "Cosine" => new UMapx.Distance.Cosine(),
        _ => (DistanceBase)Activator.CreateInstance(typeof(IDistance).Assembly.GetType("UMapx.Distance." + name)!)!
    };

    private static double NumericReference(string name, Complex[] p, Complex[] q, bool complex)
    {
        var delta = p.Zip(q, (a, b) => Complex.Abs(a-b)).ToArray();
        double normProduct = Math.Sqrt(p.Sum(z => z.Magnitude*z.Magnitude)*q.Sum(z => z.Magnitude*z.Magnitude));
        Complex inner = p.Zip(q, (a,b) => a*Complex.Conjugate(b)).Aggregate(Complex.Zero, (a,b) => a+b);
        double cosine = normProduct == 0 ? 0 : (complex ? inner.Magnitude : inner.Real)/normProduct;
        return name switch
        {
            "Euclidean" => Math.Sqrt(delta.Sum(x => x*x)),
            "SquareEuclidean" => delta.Sum(x => x*x),
            "Manhattan" => delta.Sum(),
            "Chebyshev" => delta.Max(),
            "Minkowski" => Math.Pow(delta.Sum(x => Math.Pow(x,3.5)),1/3.5),
            "Canberra" => p.Zip(q,(a,b) => a.Magnitude+b.Magnitude == 0 ? 0 : Complex.Abs(a-b)/(a.Magnitude+b.Magnitude)).Sum(),
            "BrayCurtis" => delta.Sum()/p.Zip(q,(a,b) => Complex.Abs(a+b)).Sum(),
            "Hellinger" => Math.Sqrt(p.Zip(q,(a,b) => Math.Pow(Complex.Abs(Complex.Sqrt(a)-Complex.Sqrt(b)),2)).Sum()/2),
            "Cosine" => 1-cosine,
            "Angular" => Math.Acos(Math.Clamp(cosine,-1,1)),
            _ => throw new ArgumentOutOfRangeException(nameof(name))
        };
    }

    [Theory, MemberData(nameof(NumericCases))]
    public void NumericDistancesAgreeWithDefinitions(string name, int n)
    {
        var random = new Random(8801+n);
        var p = Enumerable.Range(0,n).Select(_ => (float)(.1+random.NextDouble())).ToArray();
        var q = Enumerable.Range(0,n).Select(_ => (float)(.1+random.NextDouble())).ToArray();
        var distance = Create(name);
        double expected = NumericReference(name,p.Select(x => new Complex(x,0)).ToArray(),q.Select(x => new Complex(x,0)).ToArray(),false);
        Close(expected,distance.Compute(p,q),5e-4);
        var pm = new float[2,n]; var qm = new float[2,n];
        for(int j=0;j<n;j++) { pm[0,j]=p[j]; qm[0,j]=q[j]; pm[1,j]=q[j]; qm[1,j]=p[j]; }
        foreach(float value in distance.Compute(pm,qm)) Close(expected,value,5e-4);
    }

    [Theory, MemberData(nameof(NumericCases))]
    public void ComplexDistancesAgreeWithHermitianDefinitions(string name, int n)
    {
        var random = new Random(9923+n);
        var p = Enumerable.Range(0,n).Select(_ => new Complex32((float)(.1+random.NextDouble()),(float)(2*random.NextDouble()-1))).ToArray();
        var q = Enumerable.Range(0,n).Select(_ => new Complex32((float)(.1+random.NextDouble()),(float)(2*random.NextDouble()-1))).ToArray();
        double expected = NumericReference(name,p.Select(z=>(Complex)z).ToArray(),q.Select(z=>(Complex)z).ToArray(),true);
        var distance = Create(name);
        Close(new Complex(expected,0),distance.Compute(p,q),5e-4);
        var pm = new Complex32[2,n]; var qm = new Complex32[2,n];
        for(int j=0;j<n;j++) { pm[0,j]=p[j]; qm[0,j]=q[j]; pm[1,j]=q[j]; qm[1,j]=p[j]; }
        foreach(var value in distance.Compute(pm,qm)) Close(new Complex(expected,0),value,5e-4);
    }

    [Theory, MemberData(nameof(BooleanCases))]
    public void BooleanDistancesUseContingencyCounts(string name, int seed)
    {
        var pairs = new List<(bool A,bool B)> { (true,true), (true,false), (false,true), (false,false) };
        var random = new Random(341+seed);
        pairs.AddRange(Enumerable.Range(0,11+seed).Select(_ => (random.Next(2)==1,random.Next(2)==1)));
        double tt=pairs.Count(x=>x.A&&x.B), tf=pairs.Count(x=>x.A&&!x.B), ft=pairs.Count(x=>!x.A&&x.B), ff=pairs.Count(x=>!x.A&&!x.B);
        double mismatch=tf+ft;
        double expected=name switch
        {
            "Dice" => mismatch/(2*tt+mismatch),
            "Jaccard" => mismatch/(tt+mismatch),
            "Kulczynski" => 1-.5*(tt/(tt+tf)+tt/(tt+ft)),
            "RusselRao" => 1-tt/pairs.Count,
            "SokalMichener" => mismatch/pairs.Count,
            "SokalSneath" => 2*mismatch/(tt+2*mismatch),
            "Yule" => 2*tf*ft/(tt*ff+tf*ft),
            _ => throw new ArgumentOutOfRangeException(nameof(name))
        };
        var p=pairs.Select(x=>x.A?1f:0f).ToArray(); var q=pairs.Select(x=>x.B?1f:0f).ToArray();
        var distance=Create(name);
        Close(expected,distance.Compute(p,q));
        Close(new Complex(expected,0),distance.Compute(p.Select(x=>(Complex32)x).ToArray(),q.Select(x=>(Complex32)x).ToArray()));
    }

    [Fact]
    public void DefinedZeroVectorConventionsArePreserved()
    {
        var zero=new float[3]; var unit=new[]{1f,0,0};
        Close(1,new UMapx.Distance.Cosine().Compute(zero,unit));
        Close(0,new UMapx.Distance.Cosine(true).Compute(zero,unit));
        Close(Math.PI/2,new Angular().Compute(zero,unit));
        Close(0,new Canberra().Compute(zero,zero));
        Close(0,new Jaccard().Compute(zero,zero));
        Assert.Throws<ArgumentException>(()=>new Hellinger().Compute(new[]{-1f},new[]{1f}));
        Assert.Throws<ArgumentException>(()=>new Minkowski(.5f));
    }
}
