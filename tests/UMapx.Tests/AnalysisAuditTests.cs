using System.Numerics;
using UMapx.Analysis;
using UMapx.Core;
using Xunit;
using static UMapx.Tests.NumericAssert;

namespace UMapx.Tests;

[Trait("Category", "Analysis")]
public class AnalysisAuditTests
{
    public static IEnumerable<object[]> IntegrationCases()
    {
        foreach(var method in Enum.GetValues<IntegrationMethod>())
            foreach(bool reverse in new[]{false,true})
                foreach(bool complex in new[]{false,true})
                    yield return new object[]{method,reverse,complex};
    }

    [Theory, MemberData(nameof(IntegrationCases))]
    public void QuadratureAgreesWithAnalyticAntiderivatives(IntegrationMethod method,bool reverse,bool complex)
    {
        int degree=method==IntegrationMethod.Rectangle?0:method is IntegrationMethod.Simpson or IntegrationMethod.Romberg?3:1;
        int panels=method==IntegrationMethod.Romberg?64:12;
        int argument=method==IntegrationMethod.Romberg?7:panels;
        var q=new Integration(method); Assert.Equal(method,q.MethodType);
        Complex a=complex?new(.2,.3):new(.2,0), b=complex?new(1.1,.5):new(1.1,0);
        if(reverse)(a,b)=(b,a);
        // Round endpoints before deriving the reference so all overloads receive identical inputs.
        a=(Complex)(Complex32)a; b=(Complex)(Complex32)b;
        Complex expected=(Complex.Pow(b,degree+1)-Complex.Pow(a,degree+1))/(degree+1);
        bool sampledPanels=method is IntegrationMethod.Rectangle or IntegrationMethod.Midpoint;
        int count=sampledPanels?panels:panels+1;
        double offset=method==IntegrationMethod.Midpoint?.5:0;
        var samples=Enumerable.Range(0,count).Select(i=>Complex.Pow(a+(b-a)*(i+offset)/panels,degree)).ToArray();
        int sampledArgument=method==IntegrationMethod.Romberg?argument:count;
        if(complex)
        {
            Close(expected,q.Compute((IComplex32)(z=>(Complex32)Complex.Pow(z,degree)),(Complex32)a,(Complex32)b,argument),3e-5);
            Close(expected,q.Compute(samples.Select(z=>(Complex32)z).ToArray(),(Complex32)a,(Complex32)b,sampledArgument),3e-5);
        }
        else
        {
            Close(expected.Real,q.Compute((IFloat)(x=>(float)Math.Pow(x,degree)),(float)a.Real,(float)b.Real,argument),3e-5);
            Close(expected.Real,q.Compute(samples.Select(z=>(float)z.Real).ToArray(),(float)a.Real,(float)b.Real,sampledArgument),3e-5);
        }
    }

    public static IEnumerable<object[]> DifferentiationCases()
    {
        foreach(int points in new[]{2,3,4})foreach(int order in Enumerable.Range(0,points+1))
            foreach(bool complex in new[]{false,true})yield return new object[]{points,order,complex};
    }

    [Theory, MemberData(nameof(DifferentiationCases))]
    public void FiniteDifferencesDifferentiatePolynomials(int points,int order,bool complex)
    {
        var d=new Differentiation(points); Assert.Equal(points,d.Points);
        double factor=Enumerable.Range(points-order+1,order).Aggregate(1.0,(a,b)=>a*b);
        Complex x=complex?new(.5,.25):new(.5,0),h=complex?new(.25,.125):new(.25,0);
        Complex expected=factor*Complex.Pow(x,points-order);
        if(complex)
        {
            Close(expected,d.Compute((IComplex32)(z=>(Complex32)Complex.Pow(z,points)),(Complex32)x,(Complex32)h,order),.005);
            var samples=Enumerable.Range(0,points+1).Select(i=>(Complex32)Complex.Pow(x+i*.25,points)).ToArray();
            Close(expected,d.Compute(samples,0,.25f,order),.005);
        }
        else
        {
            Close(expected.Real,d.Compute((IFloat)(z=>(float)Math.Pow(z,points)),.5f,.25f,order),.005);
            var samples=Enumerable.Range(0,points+1).Select(i=>(float)Math.Pow(.5+i*.25,points)).ToArray();
            Close(expected.Real,d.Compute(samples,0,.25f,order),.005);
        }
    }

    [Theory] [InlineData(2)] [InlineData(3)] [InlineData(4)] [InlineData(5)]
    public void FiniteDifferenceWeightsSatisfyPolynomialMomentEquations(int n)
    {
        var c=Differentiation.GetCoefficients(n);
        for(int order=0;order<n;order++)for(int degree=0;degree<n;degree++)
        {
            double value=Enumerable.Range(0,n).Sum(j=>c[order,j]*Math.Pow(j,degree));
            double expected=order==degree?Enumerable.Range(1,degree).Aggregate(1.0,(a,b)=>a*b):0;
            Close(expected,value,.002);
        }
    }

    [Theory] [InlineData(NonlinearMethod.Bisection)] [InlineData(NonlinearMethod.Chord)] [InlineData(NonlinearMethod.Secant)] [InlineData(NonlinearMethod.FalsePosition)]
    public void RootSolversFindSimpleRootsAndHonorEndpointRoots(NonlinearMethod method)
    {
        var solver=new Nonlinear(1e-6f,method);
        Close(Math.Sqrt(2),solver.Compute((IFloat)(x=>x*x-2),1,2),2e-5);
        Close(1,solver.Compute((IFloat)(x=>x-1),1,2),2e-5);
        Close(2,solver.Compute((IFloat)(x=>x-2),1,2),2e-5);
        if(method is NonlinearMethod.Bisection or NonlinearMethod.FalsePosition)
        {
            Assert.Throws<ArgumentException>(()=>solver.Compute((IFloat)(x=>x*x+1),1,2));
            Assert.Throws<NotSupportedException>(()=>solver.Compute((IComplex32)(x=>x),new Complex32(0,0),new Complex32(1,1)));
        }
        else
        {
            var root=solver.Compute((IComplex32)(z=>z*z-new Complex32(0,2)),new Complex32(1,.5f),new Complex32(2,2));
            Close(new Complex(1,1),root,2e-5);
        }
    }

    [Theory] [InlineData(false)] [InlineData(true)]
    public void GoldenSectionFindsTheArgumentOfTheExtremum(bool maximum)
    {
        var search=new Optimization(1e-5f);
        IFloat f=x=>(maximum?-1:1)*(x-.375f)*(x-.375f);
        Close(.375,search.Compute(f,-2,3,maximum),2e-4);
        search.Eps=1e-4f; Close(1e-4,search.Eps);
    }

    public static IEnumerable<object[]> OdeCases()
    {
        foreach(var method in Enum.GetValues<DifferentialMethod>())foreach(bool irregular in new[]{false,true})
            foreach(bool complex in new[]{false,true})yield return new object[]{method,irregular,complex};
    }

    [Theory, MemberData(nameof(OdeCases))]
    public void OdeSolversFollowTheAnalyticExponentialOnUniformAndNonuniformGrids(DifferentialMethod method,bool irregular,bool complex)
    {
        var x=Enumerable.Range(0,65).Select(i=>(float)(irregular?Math.Pow(i/64.0,1.4):i/64.0)).ToArray();
        var solver=new Differential(method);
        double tolerance=method==DifferentialMethod.Euler?.04:method==DifferentialMethod.RungeKutta2?.001:2e-5;
        if(complex)
        {
            var nodes=x.Select(t=>new Complex32(t,.2f*t)).ToArray();
            var y=solver.Compute((IMeshComplex32)((t,v)=>v),nodes,new Complex32(1,0));
            Assert.Equal(x.Length-1,y.Length);
            for(int i=0;i<y.Length;i++)Close(Complex.Exp(nodes[i+1]),y[i],tolerance);
        }
        else
        {
            var y=solver.Compute((IMeshFloat)((t,v)=>v),x,1);
            Assert.Equal(x.Length-1,y.Length);
            for(int i=0;i<y.Length;i++)Close(Math.Exp(x[i+1]),y[i],tolerance);
        }
    }

    [Theory] [InlineData(2)] [InlineData(3)] [InlineData(4)]
    public void AdamsBashforthIntegratesPolynomialRightHandSides(int order)
    {
        var c=Differential.GetCoefficients(order);
        for(int p=0;p<order;p++)Close(1.0/(p+1),Enumerable.Range(0,order).Sum(j=>c[j]*Math.Pow(-j,p)),2e-5);
        var x=Enumerable.Range(0,33).Select(i=>i/32f).ToArray();
        var d=new Differential();
        var y=d.Compute((IMeshFloat)((t,v)=>2*t),x,1,order);
        var z=d.Compute((IMeshComplex32)((t,v)=>2*t),x.Select(t=>new Complex32(t,.25f*t)).ToArray(),new Complex32(1,0),order);
        for(int i=0;i<y.Length;i++)
        {
            Close(1+x[i+1]*x[i+1],y[i],3e-5);
            Close(1+Complex.Pow(new Complex(x[i+1],.25*x[i+1]),2),z[i],3e-5);
        }
    }

    public static IEnumerable<object[]> InterpolationCases()
    {
        foreach(var method in Enum.GetValues<InterpolationMethod>())foreach(float query in new[]{-1f,0,.3f,1,2.4f,3})
            foreach(bool complex in new[]{false,true})yield return new object[]{method,query,complex};
    }

    [Theory, MemberData(nameof(InterpolationCases))]
    public void InterpolationReproducesPolynomialsAtNodesAndBetweenThem(InterpolationMethod method,float query,bool complex)
    {
        var x=new[]{-1f,0,1,3}; int degree=method==InterpolationMethod.Linear?1:3;
        var d=new Interpolation(method);
        if(complex)
        {
            var nodes=x.Select(t=>new Complex32(t,.125f*t)).ToArray();
            var samples=nodes.Select(t=>(Complex32)(1+Complex.Pow(t,degree))).ToArray();
            var q=new Complex32(query,.125f*query);
            if(method==InterpolationMethod.Linear)
            {
                Assert.Throws<NotSupportedException>(()=>d.Compute(nodes,samples,q));
                return;
            }
            Close(1+Complex.Pow(q,degree),d.Compute(nodes,samples,q),2e-4);
        }
        else
        {
            var samples=x.Select(t=>(float)(1+Math.Pow(t,degree))).ToArray();
            Close(1+Math.Pow(query,degree),d.Compute(x,samples,query),2e-4);
        }
    }

    [Theory] [InlineData(ApproximationMethod.Polynomial)] [InlineData(ApproximationMethod.Logarithmic)] [InlineData(ApproximationMethod.Exponential)] [InlineData(ApproximationMethod.Power)]
    public void LeastSquaresRecoversExactModelFamilies(ApproximationMethod method)
    {
        var x=new[]{.5f,.75f,1,1.25f,1.5f,2};
        float F(float t)=>method switch
        {
            ApproximationMethod.Polynomial=>1+.75f*t,
            ApproximationMethod.Logarithmic=>1+.75f*(float)Math.Log(t),
            ApproximationMethod.Exponential=>(float)Math.Exp(1+.75*t),
            _=>(float)Math.Exp(1+.75*Math.Log(t))
        };
        var y=x.Select(F).ToArray(); var a=new Approximation(1,method);
        Close(y,a.Compute(x,y),1e-4); Close(y,a.Compute(x,y,out var c),1e-4);
        Close(new[]{1f,.75f},c,1e-4);
        Close(y,a.Compute(x,y,out c,out var fit),1e-4); Close(1,fit,2e-5);
        Close(y,a.Compute(x,y,out c,out fit,out var equation),1e-4); Assert.False(string.IsNullOrWhiteSpace(equation));
        var cx=x.Select(t=>new Complex32(t,0)).ToArray(); var cy=y.Select(t=>new Complex32(t,0)).ToArray();
        var results=new[]{a.Compute(cx,cy),a.Compute(cx,cy,out var cc),a.Compute(cx,cy,out cc,out var cf),a.Compute(cx,cy,out cc,out cf,out var ce)};
        foreach(var result in results)for(int i=0;i<x.Length;i++)Close((Complex)cy[i],result[i],2e-4);
        Close(new Complex(1,0),cc[0],2e-4); Close(new Complex(.75,0),cc[1],2e-4); Close(new Complex(1,0),cf,2e-4);
        Assert.False(string.IsNullOrWhiteSpace(ce));
    }

    [Theory] [InlineData(1)] [InlineData(2)] [InlineData(3)] [InlineData(4)]
    public void RootsSatisfyPolynomialsAndRecoverTheirCoefficients(int degree)
    {
        var expected=Enumerable.Range(1,degree).Select(i=>new Complex32(i,0)).ToArray();
        var coefficients=new double[]{1};
        foreach(var r in expected)
        {
            var next=new double[coefficients.Length+1];
            for(int i=0;i<coefficients.Length;i++){next[i]+=coefficients[i];next[i+1]-=r.Real*coefficients[i];}
            coefficients=next;
        }
        var solver=new Roots(1e-7f); var p=coefficients.Select(t=>(float)t).ToArray();
        Close(p,solver.Compute(expected),1e-4);
        foreach(var padded in new[]{p,new[]{0f,0f}.Concat(p).ToArray()})
        {
            var actual=solver.Compute(padded).OrderBy(z=>z.Real).ToArray(); Assert.Equal(degree,actual.Length);
            for(int i=0;i<degree;i++)Close((Complex)expected[i],actual[i],.002);
        }
    }

    [Theory] [InlineData(1,1)] [InlineData(1,3)] [InlineData(2,2)] [InlineData(3,1)] [InlineData(4,2)]
    public void PadeCoefficientsMatchTheTaylorSeriesThroughTheRequestedOrder(int m,int n)
    {
        var t=new float[m+n+1]; t[0]=1;for(int i=1;i<t.Length;i++)t[i]=t[i-1]/i;
        var d=new Pade(m,n);var (p,q)=d.Compute(t);
        void Check(float[] numerator,float[] denominator)
        {
            for(int k=0;k<t.Length;k++)
            {
                double coefficient=Enumerable.Range(0,Math.Min(k,denominator.Length-1)+1).Sum(j=>denominator[j]*t[k-j]);
                Close(k<numerator.Length?numerator[k]:0,coefficient,2e-5);
            }
        }
        Check(p,q); var (cp,cq)=d.Compute(t.Select(v=>new Complex32(v,0)).ToArray());
        Check(cp.Select(z=>z.Real).ToArray(),cq.Select(z=>z.Real).ToArray());
        foreach(float x in new[]{-.1f,0,.1f})
        {
            double ratio=p.Select((v,i)=>v*Math.Pow(x,i)).Sum()/q.Select((v,i)=>v*Math.Pow(x,i)).Sum();
            Close(ratio,d.Compute(x,p,q)); Close(new Complex(ratio,0),d.Compute(new Complex32(x,0),p,q));
            Close(new Complex(ratio,0),d.Compute(x,cp,cq)); Close(new Complex(ratio,0),d.Compute(new Complex32(x,0),cp,cq));
        }
        Assert.False(string.IsNullOrWhiteSpace(d.Equation(p,q))); Assert.False(string.IsNullOrWhiteSpace(d.Equation(cp,cq)));
    }
}
