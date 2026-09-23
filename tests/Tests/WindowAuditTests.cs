using UMapx.Window;
using Xunit;
using static UMapx.Tests.NumericAssert;

namespace UMapx.Tests;

[Trait("Category", "Window")]
public class WindowAuditTests
{
    private static readonly string[] Names = { "BartlettHann", "Blackman", "BlackmanHarris", "BlackmanNuttall", "Confined", "Cosine", "FlatTop", "Gabor", "Hamming", "Hann", "Kaiser", "Lanczos", "Normal", "Nuttall", "Parzen", "Planck", "Sine", "Tukey", "Welch" };
    public static IEnumerable<object[]> Cases() => Names.SelectMany(name => new[] { 2, 3, 8, 9, 32 }.Select(n => new object[] { name, n }));

    private static WindowBase Create(string name, int n) => name switch
    {
        "Normal" => new Normal(n,.7f,2),
        "Gabor" => new Gabor(n,.8f),
        "Kaiser" => new Kaiser(n,2),
        "Planck" => new Planck(n,.25f),
        "Tukey" => new Tukey(n,.5f),
        "Confined" => new Confined(n,.14f*n),
        _ => (WindowBase)Activator.CreateInstance(typeof(WindowBase).Assembly.GetType("UMapx.Window."+name)!, n)!
    };

    private static double I0(double x)
    {
        double sum=1,term=1;
        for(int k=1;k<100;k++) { term*=x*x/(4*k*k); sum+=term; if(term<1e-16*sum)break; }
        return sum;
    }

    private static double Reference(string name, int index, int n)
    {
        double t=(double)index/(n-1), centered=index-(n-1)/2.0;
        double Series(params double[] coefficients) => coefficients.Select((a,k)=>a*Math.Cos(2*Math.PI*k*t)).Sum();
        double ConfinedValue()
        {
            double sigma=(float)(.14f*n);
            double G(double x)=>Math.Exp(-Math.Pow((x-(n-1)/2.0)/(2*sigma),2));
            return G(index)-G(-.5)*(G(index+n)+G(index-n))/(G(-.5+n)+G(-.5-n));
        }
        double PlanckValue()
        {
            if(index==0||index==n-1)return 0;
            double u=Math.Min(t,1-t),a=.25;
            return u>=a?1:1/(1+Math.Exp(a/u+a/(u-a)));
        }
        return name switch
        {
            "BartlettHann" => .62-.48*Math.Abs(t-.5)-.38*Math.Cos(2*Math.PI*t),
            "Blackman" => Series(.42,-.5,.08),
            "BlackmanHarris" => Series(.35875,-.48829,.14128,-.01168),
            "BlackmanNuttall" => Series(.3635819,-.4891775,.1365995,-.0106411),
            "Confined" => ConfinedValue(),
            "Cosine" or "Sine" => Math.Sin(Math.PI*t),
            "FlatTop" => Series(1,-1.93,1.29,-.388,.028),
            "Gabor" => Math.Exp(-Math.Pow(2*Math.PI*centered/n/(double).8f,2)),
            "Hamming" => Series(.53836,-.46164),
            "Hann" => .5-.5*Math.Cos(2*Math.PI*t),
            "Kaiser" => I0(2*Math.PI*Math.Sqrt(Math.Max(0,1-Math.Pow(2*t-1,2))))/I0(2*Math.PI),
            "Lanczos" => t==.5?1:Math.Sin(Math.PI*(2*t-1))/(Math.PI*(2*t-1)),
            "Normal" => Math.Exp(-Math.Pow(centered/((double).7f*(n-1)/2),2)),
            "Nuttall" => Series(.355768,-.487396,.144232,-.012604),
            "Parzen" => Math.Abs(centered)<=n/4.0?1-6*Math.Pow(2*centered/n,2)*(1-2*Math.Abs(centered)/n):2*Math.Pow(1-2*Math.Abs(centered)/n,3),
            "Planck" => PlanckValue(),
            "Tukey" => Math.Min(t,1-t)>=.25?1:.5*(1+Math.Cos(Math.PI*(4*Math.Min(t,1-t)-1))),
            "Welch" => 1-Math.Pow(2*t-1,2),
            _ => throw new ArgumentOutOfRangeException(nameof(name))
        };
    }

    [Theory, MemberData(nameof(Cases))]
    public void SamplesAgreeWithIndependentWindowFormulas(string name, int n)
    {
        var window=Create(name,n); var actual=window.GetWindow();
        Assert.Equal(n,actual.Length);
        for(int i=0;i<n;i++) Close(Reference(name,i,n),actual[i],5e-5,5e-5);
    }

    [Theory, MemberData(nameof(Cases))]
    public void ExplicitFrameSizeDoesNotDependOnStoredFrameSize(string name, int n)
    {
        // This checks the explicit-size overload against an independently sized
        // object, including the centered coordinate conventions of each window.
        var window=Create(name,n+4); var target=Create(name,n);
        if(window is Confined oldConfined && target is Confined targetConfined)oldConfined.Sigma=targetConfined.Sigma;
        var expected=target.GetWindow();var actual=window.GetWindow(n);Assert.Equal(expected.Length,actual.Length);
        // Nonfinite mathematical values are caught by SamplesAgreeWithIndependentWindowFormulas.
        for(int i=0;i<expected.Length;i++)if(float.IsNaN(expected[i]))Assert.True(float.IsNaN(actual[i]));else Close(expected[i],actual[i],1e-5);
    }

    [Theory, MemberData(nameof(Cases))]
    public void VectorEvaluationAgreesWithScalarEvaluation(string name, int n)
    {
        var window=Create(name,n);
        float[] x={.125f,.5f,1f};
        Close(x.Select(v=>window.Function(v)).ToArray(),window.Function(x));
        Close(x.Select(v=>window.Function(v,n)).ToArray(),window.Function(x,n));
        Assert.Throws<ArgumentOutOfRangeException>(()=>window.FrameSize=0);
    }
}
