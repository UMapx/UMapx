using System.Text.Json;
using UMapx.Distribution;
using Xunit;
using static UMapx.Tests.NumericAssert;

namespace UMapx.Tests;

[Trait("Category","Distribution")]
public class DistributionShapeAuditTests
{
    public static IEnumerable<object[]> Cases()
    {
        using var stream=typeof(DistributionShapeAuditTests).Assembly.GetManifestResourceStream("UMapx.Tests.Data.distribution-catalog.json");using var document=JsonDocument.Parse(stream!);
        foreach(var law in document.RootElement.EnumerateArray())foreach(string member in new[]{"Mode","Median"})
            yield return new object[]{law.GetProperty("name").GetString()!,member,law.GetRawText()};
    }
    [Theory] [MemberData(nameof(Cases))]
    public void ModesAndMediansSatisfyTheirProbabilityDefinitions(string name,string member,string json)
    {
        using var document=JsonDocument.Parse(json);var law=document.RootElement;var type=typeof(IDistribution).Assembly.GetType("UMapx.Distribution."+name)!;
        var signature=law.GetProperty("signature").EnumerateArray().Select(v=>v.GetString()=="i"?typeof(int):typeof(float)).ToArray();
        var args=law.GetProperty("constructor").EnumerateArray().Select((v,i)=>signature[i]==typeof(int)?(object)(int)v.GetSingle():v.GetSingle()).ToArray();
        var d=(IDistribution)type.GetConstructor(signature)!.Invoke(args);
        if(law.GetProperty("unsupported").EnumerateArray().Any(v=>v.GetString()==member))
        {
            // An explicit unsupported getter is an API limitation, not a successful numerical calculation.
            var exception=Record.Exception(()=>type.GetProperty(member)!.GetValue(d));
            if(exception!=null)Assert.IsType<NotSupportedException>(exception.InnerException??exception);
            else Assert.True(float.IsNaN(Convert.ToSingle(type.GetProperty(member)!.GetValue(d))));
            return;
        }
        double Evaluate(string method,double x)=>Convert.ToDouble(type.GetMethod(method,new[]{typeof(float)})!.Invoke(d,new object[]{(float)x}));
        bool discrete=law.GetProperty("discrete").GetBoolean();
        if(member=="Median")
        {
            float median=d.Median;Assert.True(float.IsFinite(median),$"Expected a finite median; actual {median}.");
            if(law.GetProperty("unsupported").EnumerateArray().Any(v=>v.GetString()=="Distribution"))return;
            double right=Evaluate("Distribution",median),left=discrete?Evaluate("Distribution",Math.Ceiling(median)-1):right;
            Assert.True(left<=.501&&right>=.499,$"Median {median:G9}: F(m-)={left:G9}, F(m)={right:G9}.");
        }
        else
        {
            float[] modes=d.Mode;
            foreach(float mode in modes)
            {
                if(float.IsNaN(mode))
                {
                    // These APIs explicitly represent a flat mode set with NaN.
                    Assert.Contains(name,new[]{"Uniform","UniformDiscrete","Trapezoidal"});continue;
                }
                Assert.InRange(mode,d.Support.Min,d.Support.Max);
                if(!discrete&&(mode==d.Support.Min||mode==d.Support.Max))continue; // Endpoint density conventions vary.
                double center=Evaluate("Function",mode);Assert.False(double.IsNaN(center));
                double step=discrete?1:.02*(1+Math.Abs(mode));
                foreach(double x in new[]{mode-step,mode+step})if(x>d.Support.Min&&x<d.Support.Max)
                {double neighbor=Evaluate("Function",x);Assert.True(center+2e-5+Math.Abs(center)*2e-5>=neighbor,$"Mode {mode:G9}: density {center:G9} < density {neighbor:G9} at {x:G9}.");}
            }
        }
    }
    [Theory] [InlineData(0f,1f,1f)] [InlineData(.5f,2f,.7f)]
    public void BirnbaumSaundersModeIsThePositiveStationaryPoint(float location,float scale,float shape)
    {
        // Differentiate log f(x): t^3+(1+g^2)t^2+(3g^2-1)t-1=0, t=(x-location)/scale.
        double lo=0,hi=1,g2=(double)shape*shape;
        for(int i=0;i<80;i++){double t=(lo+hi)/2;double value=t*t*t+(1+g2)*t*t+(3*g2-1)*t-1;if(value>0)hi=t;else lo=t;}
        var mode=new BirnbaumSaunders(location,scale,shape).Mode;Assert.Single(mode);Close(location+scale*(lo+hi)/2,mode[0],3e-5);
    }

    [Fact]
    public void FoldedNormalWithLocationSmallerThanScaleHasModeZero()
    {
        // A positive stationary point would satisfy x=mu*tanh(mu*x/sigma^2).
        // tanh(t)<t excludes such a point when |mu|<=sigma.
        var modes=new FoldedNormal(.75f,1.25f).Mode;Assert.Single(modes);Close(0,modes[0]);
    }
}
