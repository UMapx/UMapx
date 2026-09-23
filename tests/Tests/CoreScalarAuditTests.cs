using System.Numerics;
using UMapx.Core;
using Xunit;
using static UMapx.Tests.NumericAssert;

namespace UMapx.Tests;

[Trait("Category", "Core")]
public class CoreScalarAuditTests
{
    private static readonly Dictionary<string,Func<double,double>> RealFunctions = new()
    {
        ["Pow"]=x=>x*x,["Exp"]=Math.Exp,["Log"]=Math.Log,["Log10"]=Math.Log10,["Log2"]=Math.Log2,
        ["Sqrt"]=Math.Sqrt,["Abs"]=Math.Abs,["Floor"]=Math.Floor,["Ceil"]=Math.Ceiling,["Round"]=Math.Round,
        ["Cos"]=Math.Cos,["Sin"]=Math.Sin,["Tan"]=Math.Tan,["Ctan"]=x=>1/Math.Tan(x),
        ["Sec"]=x=>1/Math.Cos(x),["Cosc"]=x=>1/Math.Sin(x),["Asin"]=Math.Asin,["Acos"]=Math.Acos,
        ["Atan"]=Math.Atan,["Actan"]=x=>Math.PI/2-Math.Atan(x),["Asec"]=x=>Math.Acos(1/x),["Acosc"]=x=>Math.Asin(1/x),
        ["Sinh"]=Math.Sinh,["Cosh"]=Math.Cosh,["Tanh"]=Math.Tanh,["Ctanh"]=x=>1/Math.Tanh(x),
        ["Sech"]=x=>1/Math.Cosh(x),["Cosch"]=x=>1/Math.Sinh(x),["Asinh"]=Math.Asinh,["Acosh"]=Math.Acosh,
        ["Atanh"]=Math.Atanh,["Actanh"]=x=>Math.Atanh(1/x),["Asech"]=x=>Math.Acosh(1/x),["Acosch"]=x=>Math.Asinh(1/x)
    };

    public static IEnumerable<object[]> RealCases()
    {
        foreach(var pair in RealFunctions)
            foreach(float x in new[]{-10000f,-100,-5,-2,-.9f,-.1f,.1f,.5f,.9f,1,2,5,100,10000})
            {
                double expected=pair.Value(x);
                if(double.IsFinite(expected)&&Math.Abs(expected)<=float.MaxValue)
                    yield return new object[]{pair.Key,x,expected};
            }
    }

    [Theory, MemberData(nameof(RealCases))]
    public void RealElementaryFunctionsAgreeWithDoublePrecision(string name,float x,double expected)
    {
        var method=typeof(Maths).GetMethod(name,new[]{typeof(float)})!;
        Close(expected,Convert.ToDouble(method.Invoke(null,new object[]{x})),3e-6,3e-5);
    }

    private static Complex Asinh(Complex z) => z.Real<0?-Asinh(-z):Complex.Log(z+Complex.Sqrt(z*z+1));
    private static Complex Acosh(Complex z) => Complex.Log(z+Complex.Sqrt(z-1)*Complex.Sqrt(z+1));
    private static readonly Dictionary<string,Func<Complex,Complex>> ComplexFunctions = new()
    {
        ["Exp"]=Complex.Exp,["Log"]=Complex.Log,["Log10"]=Complex.Log10,["Log2"]=z=>Complex.Log(z)/Math.Log(2),
        ["Sqrt"]=Complex.Sqrt,["Sin"]=Complex.Sin,["Cos"]=Complex.Cos,["Tan"]=Complex.Tan,
        ["Ctan"]=z=>1/Complex.Tan(z),["Sec"]=z=>1/Complex.Cos(z),["Cosc"]=z=>1/Complex.Sin(z),
        ["Asin"]=Complex.Asin,["Acos"]=Complex.Acos,["Atan"]=Complex.Atan,["Actan"]=z=>Complex.Atan(1/z),
        ["Asec"]=z=>Complex.Acos(1/z),["Acosc"]=z=>Complex.Asin(1/z),
        ["Sinh"]=Complex.Sinh,["Cosh"]=Complex.Cosh,["Tanh"]=Complex.Tanh,
        ["Ctanh"]=z=>1/Complex.Tanh(z),["Sech"]=z=>1/Complex.Cosh(z),["Cosch"]=z=>1/Complex.Sinh(z),
        ["Asinh"]=Asinh,["Acosh"]=Acosh,["Atanh"]=z=>(Complex.Log(1+z)-Complex.Log(1-z))/2,
        ["Actanh"]=z=>(Complex.Log(1+1/z)-Complex.Log(1-1/z))/2,["Asech"]=z=>Acosh(1/z),["Acosch"]=z=>Asinh(1/z)
    };

    public static IEnumerable<object[]> ComplexCases()
    {
        foreach(string name in ComplexFunctions.Keys)
            foreach(var z in new[]{new Complex32(.2f,.4f),new Complex32(-2,.5f),new Complex32(-2,-.5f),new Complex32(4,1),new Complex32(1,20),new Complex32(-1,-20)})
                yield return new object[]{name,z.Real,z.Imag};
    }

    [Theory, MemberData(nameof(ComplexCases))]
    public void ComplexElementaryFunctionsRespectPrincipalValues(string name,float real,float imaginary)
    {
        var z=new Complex32(real,imaginary);
        var method=typeof(Maths).GetMethod(name,new[]{typeof(Complex32)})!;
        Close(ComplexFunctions[name](z),(Complex32)method.Invoke(null,new object[]{z})!,5e-5,8e-5);
    }

    [Theory] [InlineData(7)] [InlineData(53)] [InlineData(731)]
    public void ComplexArithmeticAndConversionsAgreeWithSystemNumerics(int seed)
    {
        var random=new Random(seed);
        for(int i=0;i<40;i++)
        {
            var a=new Complex32((float)(random.NextDouble()*4-2),(float)(random.NextDouble()*4-2));
            var b=new Complex32((float)(random.NextDouble()+.1),(float)(random.NextDouble()+.1));
            float s=.25f+(float)random.NextDouble(); Complex x=a,y=b;
            Close(x+y,a+b); Close(x-y,a-b); Close(x*y,a*b); Close(x/y,a/b);
            Close(x+s,a+s); Close(s+x,s+a); Close(x-s,a-s); Close(s-x,s-a);
            Close(x*s,a*s); Close(s*x,s*a); Close(x/s,a/s); Close(s/x,s/a);
            Close(-x,-a); Close(x,+a); Close(Complex.Conjugate(x),a.Conjugate);
            Close(x.Magnitude,a.Abs); Close(x.Magnitude*x.Magnitude,a.Abs2);
            Close(x.Phase,Maths.Angle(a)); Close(x.Magnitude,Maths.Abs(a));
            Close(x,Maths.FromPolar((float)x.Magnitude,(float)x.Phase));
            Close(x,Complex32.FromPolarCoordinates((float)x.Magnitude,(float)x.Phase));
            Close(Complex.Pow(x,s),Maths.Pow(a,s)); Close(Complex.Pow(x,y),Maths.Pow(a,b));
            Close(Complex.Pow(s,y),Maths.Pow(s,b)); Close(Complex.Pow(x,1/s),Maths.Sqrt(a,s));
            Close(Complex.Pow(x,1/y),Maths.Sqrt(a,b)); Close(Complex.Log(x)/Math.Log(s),Maths.Log(a,s));
            Close(new Complex(Math.Round(a.Real),Math.Round(a.Imag)),Maths.Round(a));
            Close(new Complex(Math.Round(a.Real,2),Math.Round(a.Imag,2)),Maths.Round(a,2));
            Assert.Equal(a,a.Clone()); Assert.True(a==a.Clone()); Assert.False(a!=a.Clone());
            Assert.False(Complex32.IsNaN(a)); Assert.False(Complex32.IsInfinity(a));
            Assert.Equal(a.GetHashCode(),a.Clone().GetHashCode());
        }
        Assert.True(Complex32.IsNaN(new Complex32(float.NaN,0)));
        Assert.True(Complex32.IsInfinity(new Complex32(0,float.PositiveInfinity)));
    }

    private static Quaternion Q(Quaternion32 value)=>new(value.X,value.Y,value.Z,value.W);
    private static Quaternion32 Q(Quaternion value)=>new(value.X,value.Y,value.Z,value.W);
    private static void EqualQuaternion(Quaternion expected,Quaternion32 actual)
    { Close(expected.X,actual.X); Close(expected.Y,actual.Y); Close(expected.Z,actual.Z); Close(expected.W,actual.W); }

    [Theory] [InlineData(0f)] [InlineData(.2f)] [InlineData(.75f)] [InlineData(1f)]
    public void QuaternionAlgebraAndInterpolationAgreeWithSystemNumerics(float amount)
    {
        var a=Quaternion32.FromYPR(.2f,-.4f,.7f); var b=Quaternion32.FromYPR(-.6f,.3f,-.8f);
        EqualQuaternion(Quaternion.CreateFromYawPitchRoll(.2f,-.4f,.7f),a);
        EqualQuaternion(Q(a)+Q(b),a+b); EqualQuaternion(Q(a)-Q(b),a-b);
        EqualQuaternion(Q(a)*Q(b),a*b); EqualQuaternion(Q(a)/Q(b),a/b);
        EqualQuaternion(-Q(a),-a); EqualQuaternion(Q(a)*2,a*2); EqualQuaternion(Q(a)*.5f,a/2);
        EqualQuaternion(Quaternion.Conjugate(Q(a)),a.Conjugate);
        EqualQuaternion(Quaternion.Inverse(Q(a)),a.Inverse);
        EqualQuaternion(Quaternion.Normalize(Q(a)*2),(a*2).Normalize);
        EqualQuaternion(Quaternion.Concatenate(Q(a),Q(b)),Quaternion32.Concatenate(a,b));
        EqualQuaternion(Quaternion.Lerp(Q(a),Q(b),amount),Quaternion32.Lerp(a,b,amount));
        EqualQuaternion(Quaternion.Slerp(Q(a),Q(b),amount),Quaternion32.Slerp(a,b,amount));
        Close(Quaternion.Dot(Q(a),Q(b)),Quaternion32.Dot(a,b)); Close(Q(a).Length(),a.Abs); Close(Q(a).LengthSquared(),a.SquaredAbs);
        Assert.True(Quaternion32.Identity.IsIdentity); Assert.Equal(a,a.Clone());
        Assert.True(a==a.Clone()); Assert.False(a!=a.Clone()); Assert.Equal(a.GetHashCode(),a.Clone().GetHashCode());
    }

    [Fact]
    public void ScalarRangesAndRoundingFollowTheirDefinitions()
    {
        foreach(float x in new[]{-300f,-1.25f,0,1.25f,100,300})
        {
            Close(Math.Clamp(x,-2,2),Maths.Range(x,-2,2));
            Assert.Equal(x>=-2&&x<=2,Maths.IsRange(x,-2,2));
            Close((x+2)/4,Maths.Normalize(x,-2,2));
            Assert.Equal((byte)Math.Clamp((int)x,0,255),Maths.Byte(x));
            Assert.Equal((sbyte)Math.Clamp((int)x,-128,127),Maths.sByte(x));
            Close(Math.Round(x,1),Maths.Round(x,1));
        }
        foreach(int x in new[]{-300,-1,0,1,100,300})
        {
            Assert.Equal(Math.Clamp(x,-2,2),Maths.Range(x,-2,2));
            Assert.Equal(x>=-2&&x<=2,Maths.IsRange(x,-2,2));
            Assert.Equal((byte)Math.Clamp(x,0,255),Maths.Byte(x));
            Assert.Equal((sbyte)Math.Clamp(x,-128,127),Maths.sByte(x));
        }
        foreach(float x in new[]{.25f,1,4,16})
        {
            Close(Math.Log(x,3),Maths.Log(x,3)); Close(Math.Pow(x,1.0/3),Maths.Sqrt(x,3)); Close(x*x,Maths.Pow(x));
        }
        foreach(float a in new[]{-2f,.5f,3})foreach(float b in new[]{-3f,.25f,2})
        {
            Close(Math.Max(a,b),Maths.Max(a,b)); Close(Math.Min(a,b),Maths.Min(a,b));
            Close(Math.Max(a,Math.Max(b,1)),Maths.Max(a,b,1)); Close(Math.Min(a,Math.Min(b,1)),Maths.Min(a,b,1));
            Close(Math.Atan2(a,b),Maths.Atan2(a,b)); Close(Math.Sqrt((double)a*a+(double)b*b),Maths.Hypotenuse(a,b));
        }
    }
}
