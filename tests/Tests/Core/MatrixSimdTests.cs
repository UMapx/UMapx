using System.Numerics;
using UMapx.Core;
using Xunit;

namespace UMapx.Tests;

[CollectionDefinition("SIMD settings", DisableParallelization = true)]
public class SimdCollection
{
}

[Collection("SIMD settings")]
[Trait("Category", "Matrix")]
public class MatrixSimdTests
{
    [Theory]
    [InlineData(false, false)]
    [InlineData(true, false)]
    [InlineData(false, true)]
    [InlineData(true, true)]
    public void SimdProductsMatchIndependentAccumulationIncludingTailElements(bool complexA, bool complexB)
    {
        var a = (Array)MatrixTestSupport.Operand(complexA ? typeof(Complex32[,]) : typeof(float[,]), 1, 13, 37);
        var b = (Array)MatrixTestSupport.Operand(complexB ? typeof(Complex32[,]) : typeof(float[,]), 2, 37, 19);
        bool original = Globals.SIMD;
        try
        {
            Globals.SIMD = true;
            var actual = (Array)MatrixTestSupport.Invoke(typeof(Matrice).GetMethod("Dot", new[] { a.GetType(), b.GetType() })!, a, b);
            for (int y = 0; y < 13; y++)
                for (int x = 0; x < 19; x++)
                {
                    Complex expected = 0;
                    for (int k = 0; k < 37; k++)
                        expected += MatrixTestSupport.Value(a, y, k) * MatrixTestSupport.Value(b, k, x);
                    MatrixTestSupport.Check(expected, actual.GetValue(y, x)!, 1e-4);
                }
        }
        finally
        {
            Globals.SIMD = original;
        }
    }
}
