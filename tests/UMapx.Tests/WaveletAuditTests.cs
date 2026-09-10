using System.Numerics;
using System.Reflection;
using UMapx.Core;
using UMapx.Wavelet;
using Xunit;
using static UMapx.Tests.NumericAssert;

namespace UMapx.Tests;

[Trait("Category", "Wavelet")]
public class WaveletAuditTests
{
    private static WaveletPacket Bank(string name)=>(WaveletPacket)typeof(WaveletPacket).GetProperty(name,BindingFlags.Static|BindingFlags.Public)!.GetValue(null)!;
    public static IEnumerable<object[]> Banks()=>typeof(WaveletPacket).GetProperties(BindingFlags.Static|BindingFlags.Public).Where(p=>p.PropertyType==typeof(WaveletPacket)).Select(p=>new object[]{p.Name});

    [Theory, MemberData(nameof(Banks))]
    public void EveryBankAgreesWithIndependentPeriodicFilterConvolution(string name)
    {
        var bank=Bank(name);var clone=bank.Clone();Assert.NotSame(bank.LowPass,clone.LowPass);Close(bank.LowPass,clone.LowPass);
        int n=64;var x=Enumerable.Range(0,n).Select(i=>new Complex32((float)Math.Sin(i*.37),(float)Math.Cos(i*.71))).ToArray();
        foreach(bool normalized in new[]{false,true})
        {
            var d=new WaveletDecomposition(bank,1,normalized);var actual=d.Forward(x);var real=d.Forward(x.Select(z=>z.Real).ToArray());
            foreach(int band in new[]{0,1})
            {
                float[] h=band==0?bank.LowPass:bank.HighPass;
                for(int i=0;i<n/2;i++)
                {
                    Complex sum=0;
                    for(int tap=0;tap<h.Length;tap++)sum+=h[tap]*(Complex)x[((2*i+tap-(h.Length/2-1))%n+n)%n];
                    if(normalized)sum/=Math.Sqrt(2);
                    Close(sum,actual[band][i],2e-4);Close(sum.Real,real[band][i],2e-4);
                }
            }
            // Check synthesis independently even for explicitly nonorthogonal Legendre banks.
            var restored=d.Backward(actual);var restoredReal=d.Backward(real);
            for(int i=0;i<n;i++)
            {
                Complex expected=0;double expectedReal=0;
                foreach(int band in new[]{0,1})
                {
                    float[] h=band==0?bank.ILowPass:bank.IHighPass;
                    for(int tap=0;tap<h.Length;tap++)
                    {
                        int index=((i+tap-(h.Length/2-1))%n+n)%n;
                        if(index%2==1){expected+=h[tap]*(Complex)actual[band][index/2];expectedReal+=h[tap]*real[band][index/2];}
                    }
                }
                if(normalized){expected*=Math.Sqrt(2);expectedReal*=Math.Sqrt(2);}
                Close(expected,restored[i],3e-4);Close(expectedReal,restoredReal[i],3e-4);
            }
        }
    }

    public static IEnumerable<object[]> InvertibleBanks()
    {
        foreach(var row in Banks())
        {
            string name=(string)row[0];
            if(name.StartsWith('L')&&name!="L1")continue; // Explicitly documented as nonorthogonal.
            if(name is "Fbsp103" or "Fbsp105")continue; // Prototype filters without a supplied dual bank.
            foreach(int levels in new[]{1,3})yield return new object[]{name,levels};
        }
    }

    [Theory, MemberData(nameof(InvertibleBanks))]
    public void ReconstructionBanksRecoverImpulsesAndDeterministicSignals(string name,int levels)
    {
        int n=256;var x=Enumerable.Range(0,n).Select(i=>(float)(.25*Math.Sin(i*.37)+.15*Math.Cos(i*.9))).ToArray();x[17]+=1;
        var d=new WaveletDecomposition(Bank(name),levels);
        // BL and Meyer are finite approximations to infinite filters; retain an explicit approximation budget.
        double tolerance=name.StartsWith("BL")||name=="Meyer"?.003:5e-4;
        Close(x,d.Backward(d.Forward(x)),tolerance);
    }

    [Theory] [InlineData(1)] [InlineData(2)] [InlineData(3)]
    public void WaveletPackingAndMatrixDecompositionPreserveAllComponents(int levels)
    {
        var d=new WaveletDecomposition(WaveletPacket.D4,levels);var t=new WaveletTransform(d);
        var x=Enumerable.Range(0,64).Select(i=>(float)Math.Sin(i*.31)).ToArray();var z=x.Select((v,i)=>new Complex32(v,(float)Math.Cos(i*.2))).ToArray();
        Close(x,t.Backward(t.Forward(x)),2e-4);var rz=t.Backward(t.Forward(z));for(int i=0;i<x.Length;i++)Close((Complex)z[i],rz[i],3e-4);
        var a=NumericAssert.Matrix(16,24);var c=TransformAuditTests.ToComplex(a);for(int i=0;i<16;i++)for(int j=0;j<24;j++)c[i,j].Imag=(i-j)*.01f;
        Close(a,d.Backward(d.Forward(a)),3e-4);Close(a,t.Backward(t.Forward(a)),3e-4);
        foreach(var r in new[]{d.Backward(d.Forward(c)),t.Backward(t.Forward(c))})for(int i=0;i<16;i++)for(int j=0;j<24;j++)Close((Complex)c[i,j],r[i,j],5e-4);
    }

    [Theory] [InlineData(1)] [InlineData(3)]
    public void EdgeAvoidingWaveletResidualsReconstructWithoutLosingPhase(int levels)
    {
        var d=new EdgeAvoidingWaveletDecomposition(2,.1f,levels);var x=Enumerable.Range(0,16).Select(i=>(float)(.5+.4*Math.Sin(i*.31))).ToArray();
        var z=x.Select((v,i)=>new Complex32(v,.125f)).ToArray();var a=new float[8,12];var c=new Complex32[8,12];
        for(int i=0;i<8;i++)for(int j=0;j<12;j++){a[i,j]=x[(i+j)%16];c[i,j]=new(a[i,j],.125f);}
        Close(x,d.Backward(d.Forward(x)),1e-4);Close(a,d.Backward(d.Forward(a)),1e-4);
        var rz=d.Backward(d.Forward(z));var rc=d.Backward(d.Forward(c));
        for(int i=0;i<z.Length;i++)Close((Complex)z[i],rz[i],1e-4);for(int i=0;i<8;i++)for(int j=0;j<12;j++)Close((Complex)c[i,j],rc[i,j],1e-4);
    }

    public static IEnumerable<object[]> GaussianCases()
    {foreach(int n in Enumerable.Range(1,8))foreach(float x in new[]{-2f,-.75f,0,.25f,1,3})yield return new object[]{n,x};}
    private static Complex GaussianDerivative(int n,Complex x,Complex carrier)
    {
        // Differentiate the polynomial recursively, independently of the tabulated expressions.
        Complex[] p={1};
        for(int order=0;order<n;order++)
        {
            var q=new Complex[p.Length+1];
            for(int j=0;j<p.Length;j++){if(j>0)q[j-1]+=j*p[j];q[j]-=carrier*p[j];q[j+1]-=2*p[j];}
            p=q;
        }
        Complex value=0;for(int j=p.Length-1;j>=0;j--)value=value*x+p[j];
        return value*Complex.Exp(-x*x-carrier*x);
    }

    [Theory, MemberData(nameof(GaussianCases))]
    public void GaussianWaveletsHaveTheRequestedDerivativeShape(int n,float x)
    {
        double oddFactorial=Enumerable.Range(1,n).Aggregate(1.0,(a,k)=>a*(2*k-1));
        double scale=Math.Pow(2/Math.PI,.25)/Math.Sqrt(oddFactorial);
        Close((GaussianDerivative(n,x,0)*scale).Real,new Gaussian(n).Wavelet(x),2e-5);
        var complex=new ComplexGaussian(n);
        // Complex Gaussian normalization conventions differ; test the derivative shape up to its constant multiplier.
        Complex anchor=complex.Wavelet(.25f),normalizer=anchor/GaussianDerivative(n,.25,Complex.ImaginaryOne);
        Close(GaussianDerivative(n,x,Complex.ImaginaryOne)*normalizer,complex.Wavelet(x),1e-4);
        Assert.Throws<NotSupportedException>(()=>new Gaussian(n).Scaling(x));Assert.Throws<NotSupportedException>(()=>complex.Scaling(x));
    }

    [Theory] [InlineData(-2f)] [InlineData(-.75f)] [InlineData(0f)] [InlineData(.25f)] [InlineData(.5f)] [InlineData(1f)] [InlineData(3f)]
    public void ContinuousWaveletsAgreeWithTheirAnalyticDefinitions(float x)
    {
        double t=x,pi=Math.PI;var haar=new Haar();
        Close(t>=0&&t<1?1:0,haar.Scaling(x));Close(t>=0&&t<.5?1:t>=.5&&t<1?-1:0,haar.Wavelet(x));
        Close(2/(Math.Sqrt(3)*Math.Pow(pi,.25))*(1-t*t)*Math.Exp(-t*t/2),new MexicanHat().Wavelet(x));
        Close(Math.Exp(-t*t/2)*(Math.Cos(5*t)-Math.Exp(-12.5))/Math.Pow(pi,.25),new Morlet(5).Wavelet(x));
        Close(2/Math.Sqrt(5)*Math.Pow(pi,-.25)*new Complex(1-t*t,t)*Math.Exp(-t*t/2),new HermitianHat().Wavelet(x));
        Close(Math.Pow(pi*.5,-.25)*(Complex.Exp(2*pi*Complex.ImaginaryOne*t)-Math.Exp(-pi*pi))*Math.Exp(-t*t),new ComplexMorlet().Wavelet(x));
        var g=new Gabor(.25f,1.5f,2);var gv=Complex.Exp(-(t-.25)*(t-.25)/4-Complex.ImaginaryOne*1.5*(t-.25));
        Close(gv,g.Wavelet(x));Close(gv.Real,g.WaveletReal(x));
        double Sinc(double v)=>v==0?1:Math.Sin(pi*v)/(pi*v);
        // Retain the library's fbsp parameterization with the unnormalized cardinal sine.
        Close(Math.Sqrt(2)*Math.Pow(t==0?1:Math.Sin(t/8)/(t/8),3)*Complex.Exp(4*pi*Complex.ImaginaryOne*t),new Fbsp(3,2,2).Wavelet(x));
        Close(Sinc(t/2)*Math.Cos(3*pi*t/2),new Shannon().Wavelet(x));
        foreach(int n in new[]{1,2,5})
        {
            double factorial=Enumerable.Range(1,n).Aggregate(1.0,(a,b)=>a*b);
            Close(t<0?0:(t-n)*Math.Pow(t,n-1)*Math.Exp(-t)/factorial,new Poisson(n).Wavelet(x));
        }
        foreach(int n in new[]{1,2,3})
        {
            double polynomial=n==1?Math.Sqrt(2)*t:n==2?2/Math.Sqrt(3)*(1-t*t):2*Math.Sqrt(30)/15*(t*t*t-3*t);
            Close(new Complex(polynomial*Math.Pow(pi,-.25)*Math.Exp(-t*t/2),0),new Hermitian(n).Wavelet(x));
        }
    }

    [Theory] [InlineData(true,.5f)] [InlineData(false,.75f)] [InlineData(false,-.75f)]
    public void MeyerWaveletAndScalingHaveFiniteRemovableSingularities(bool wavelet,float x)
    {
        var d=new Meyer();Close(wavelet?4/Math.PI:2/(3*Math.PI),wavelet?d.Wavelet(x):d.Scaling(x),1e-5);
    }
}
