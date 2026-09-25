using System.Numerics;
using System.Reflection;
using UMapx.Core;
using Xunit;
using static UMapx.Tests.NumericAssert;

namespace UMapx.Tests;

internal static class MatrixTestSupport
{
    internal static object Operand(Type type, int seed, int rows = 3, int columns = 5)
    {
        object Scalar(Type t, int i, int j) => t == typeof(float) ? (object)(float)(.7 + .17 * i + .11 * j + .23 * seed) : new Complex32((float)(.7 + .17 * i + .11 * j + .23 * seed), (float)(.2 + .07 * i - .03 * j));
        if (!type.IsArray)
            return Scalar(type, 0, 0);
        var e = type.GetElementType()!;
        var result = type.GetArrayRank() == 1 ? Array.CreateInstance(e, columns) : Array.CreateInstance(e, rows, columns);
        for (int i = 0; i < (result.Rank == 1 ? 1 : rows); i++)
            for (int j = 0; j < columns; j++)
                if (result.Rank == 1)
                    result.SetValue(Scalar(e, 0, j), j);
                else
                    result.SetValue(Scalar(e, i, j), i, j);
        return result;
    }

    internal static Complex Value(object v, int i = 0, int j = 0)
    {
        if (v is Array a)
            v = a.Rank == 1 ? a.GetValue(j)! : a.GetValue(i, j)!;
        return v is Complex32 z ? (Complex)z : new Complex(Convert.ToDouble(v), 0);
    }

    internal static void Check(Complex expected, object actual, double tolerance = 2e-5)
    {
        if (actual is Complex32 z)
            Close(expected, z, tolerance, tolerance);
        else
        {
            Close(0, expected.Imaginary, tolerance);
            Close(expected.Real, Convert.ToDouble(actual), tolerance, tolerance);
        }
    }

    internal static object Invoke(MethodInfo method, params object[] args)
    {
        try
        {
            return method.Invoke(null, args)!;
        }
        catch (TargetInvocationException e) when (e.InnerException != null)
        {
            System.Runtime.ExceptionServices.ExceptionDispatchInfo.Capture(e.InnerException).Throw();
            throw;
        }
    }

    internal static object Call(string name, params object[] args) => Invoke(typeof(Matrice).GetMethod(name, args.Select(x => x.GetType()).ToArray())!, args);
    internal static void RealComplex(double expected, Complex32 actual)
    {
        Assert.Equal(0, actual.Imag);
        if (float.IsInfinity((float)expected))
            Assert.Equal((float)expected, actual.Real);
        else
            Close(expected, actual.Real, 4 * (double)float.Epsilon, 2e-5);
    }

    internal static Complex32[,] ToComplex(float[,] a)
    {
        var r = new Complex32[a.GetLength(0), a.GetLength(1)];
        for (int i = 0; i < r.GetLength(0); i++)
            for (int j = 0; j < r.GetLength(1); j++)
                r[i, j] = a[i, j];
        return r;
    }
}
