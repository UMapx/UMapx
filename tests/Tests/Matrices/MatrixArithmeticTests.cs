using System.Numerics;
using System.Reflection;
using UMapx.Core;
using Xunit;
using static UMapx.Tests.MatrixTestSupport;

namespace UMapx.Tests;

[Trait("Category", "Matrix")]
public class MatrixArithmeticTests
{
    static readonly Type[] NumericTypes =
    {
        typeof(float),
        typeof(Complex32),
        typeof(float[]),
        typeof(Complex32[]),
        typeof(float[, ]),
        typeof(Complex32[, ])
    };
    static readonly MethodInfo[] Arithmetic = typeof(Matrice).GetMethods(BindingFlags.Public | BindingFlags.Static).Where(m => new[] { "Add", "Sub", "Mul", "Div", "Pow" }.Contains(m.Name) && m.GetParameters().Length == 2 && m.GetParameters().All(p => NumericTypes.Contains(p.ParameterType))).OrderBy(m => m.ToString()).ToArray();
    public static IEnumerable<object[]> ArithmeticCases() => Arithmetic.Select(m => new object[] { m.ToString()! });
    [Theory]
    [MemberData(nameof(ArithmeticCases))]
    public void EveryElementwiseArithmeticOverloadMatchesScalarArithmetic(string signature)
    {
        var m = Arithmetic.Single(m => m.ToString() == signature);
        var p = m.GetParameters();
        var a = Operand(p[0].ParameterType, 1);
        var b = Operand(p[1].ParameterType, 2);
        var actual = (Array)Invoke(m, a, b);
        int rows = actual.Rank == 1 ? 1 : actual.GetLength(0), cols = actual.GetLength(actual.Rank - 1);
        Assert.Equal(5, cols);
        if (actual.Rank == 2)
            Assert.Equal(3, rows);
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
            {
                Complex x = Value(a, i, j), y = Value(b, i, j);
                Complex expected = m.Name switch
                {
                    "Add" => x + y,
                    "Sub" => x - y,
                    "Mul" => x * y,
                    "Div" => x / y,
                    "Pow" => Complex.Pow(x, y),
                    _ => throw new Exception()
                };
                Check(expected, actual.Rank == 1 ? actual.GetValue(j)! : actual.GetValue(i, j)!);
            }
    }

    public static IEnumerable<object[]> ProductCases()
    {
        foreach (string op in new[]
        {
            "Dot",
            "Kronecker"
        }

        )
            foreach (bool ac in new[]
            {
                false,
                true
            }

            )
                foreach (bool bc in new[]
                {
                    false,
                    true
                }

                )
                    yield return new object[]
                    {
                        op,
                        ac,
                        bc
                    };
    }

    [Theory]
    [MemberData(nameof(ProductCases))]
    public void MixedMatrixProductsMatchIndependentComplexArithmetic(string operation, bool complexA, bool complexB)
    {
        var a = Operand(complexA ? typeof(Complex32[,]) : typeof(float[,]), 1, 3, 5);
        var b = Operand(complexB ? typeof(Complex32[,]) : typeof(float[,]), 2, 5, 2);
        var actual = (Array)Invoke(typeof(Matrice).GetMethod(operation, new[] { a.GetType(), b.GetType() })!, a, b);
        int rows = operation == "Dot" ? 3 : 15, cols = operation == "Dot" ? 2 : 10;
        Assert.Equal(rows, actual.GetLength(0));
        Assert.Equal(cols, actual.GetLength(1));
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
            {
                Complex expected = 0;
                if (operation == "Dot")
                    for (int k = 0; k < 5; k++)
                        expected += Value(a, i, k) * Value(b, k, j);
                else
                    expected = Value(a, i / 5, j / 2) * Value(b, i % 5, j % 2);
                Check(expected, actual.GetValue(i, j)!, 1e-4);
            }
    }

    public static IEnumerable<object[]> DiagonalCases()
    {
        foreach (bool left in new[]
        {
            false,
            true
        }

        )
            foreach (bool inverse in new[]
            {
                false,
                true
            }

            )
                foreach (bool mc in new[]
                {
                    false,
                    true
                }

                )
                    foreach (bool vc in new[]
                    {
                        false,
                        true
                    }

                    )
                        yield return new object[]
                        {
                            left,
                            inverse,
                            mc,
                            vc
                        };
    }

    [Theory]
    [MemberData(nameof(DiagonalCases))]
    public void DiagonalMultiplicationScalesTheDocumentedRowsOrColumns(bool left, bool inverse, bool complexMatrix, bool complexVector)
    {
        var a = Operand(complexMatrix ? typeof(Complex32[,]) : typeof(float[,]), 1, 3, 5);
        var v = Operand(complexVector ? typeof(Complex32[]) : typeof(float[]), 2, 1, left ? 3 : 5);
        var args = left ? new[]
        {
            v,
            a,
            (object)inverse
        }

        : new[]
        {
            a,
            v,
            (object)inverse
        };
        var actual = (Array)Invoke(typeof(Matrice).GetMethod("Dot", args.Select(a => a.GetType()).ToArray())!, args);
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 5; j++)
                Check(inverse ? Value(a, i, j) / Value(v, 0, left ? i : j) : Value(a, i, j) * Value(v, 0, left ? i : j), actual.GetValue(i, j)!);
    }

    [Theory]
    [InlineData(false, false)]
    [InlineData(false, true)]
    [InlineData(true, false)]
    [InlineData(true, true)]
    public void InverseSolveAndDeterminantSatisfyIndependentEquations(bool complex, bool pivot)
    {
        float[,] real = pivot ? new float[,]
        {
            {
                0,
                2,
                1
            },
            {
                2,
                3,
                -1
            },
            {
                1,
                1,
                4
            }
        }

        : new float[,]
        {
            {
                4,
                2,
                1
            },
            {
                2,
                3,
                -1
            },
            {
                1,
                1,
                4
            }
        };
        var a = complex ? (Array)real.ToComplex() : (Array)real;
        if (complex)
        {
            var c = (Complex32[,])a;
            c[0, 2] += new Complex32(0, .25f);
            c[2, 0] -= new Complex32(0, .5f);
        }

        var inverse = (Array)Invoke(typeof(Matrice).GetMethod("Invert", new[] { a.GetType() })!, a);
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                Complex sum = 0;
                for (int k = 0; k < 3; k++)
                    sum += Value(a, i, k) * Value(inverse, k, j);
                Check(i == j ? 1 : 0, (Complex32)sum, 2e-4);
            }

        var b = Operand(complex ? typeof(Complex32[]) : typeof(float[]), 2, 1, 3);
        var solution = (Array)Invoke(typeof(Matrice).GetMethod("Solve", new[] { a.GetType(), b.GetType() })!, a, b);
        for (int i = 0; i < 3; i++)
        {
            Complex sum = 0;
            for (int j = 0; j < 3; j++)
                sum += Value(a, i, j) * Value(solution, 0, j);
            Check(Value(b, 0, i), (Complex32)sum, 2e-4);
        }

        Complex determinant = Value(a, 0, 0) * (Value(a, 1, 1) * Value(a, 2, 2) - Value(a, 1, 2) * Value(a, 2, 1)) - Value(a, 0, 1) * (Value(a, 1, 0) * Value(a, 2, 2) - Value(a, 1, 2) * Value(a, 2, 0)) + Value(a, 0, 2) * (Value(a, 1, 0) * Value(a, 2, 1) - Value(a, 1, 1) * Value(a, 2, 0));
        Check(determinant, Invoke(typeof(Matrice).GetMethod("Det", new[] { a.GetType() })!, a), 2e-4);
    }

    public static IEnumerable<object[]> RectangularDiagonalCases()
    {
        foreach (var flags in DiagonalCases())
            foreach (var shape in new[]
            {
                (0, 3),
                (3, 0),
                (1, 7),
                (7, 1),
                (2, 4),
                (5, 3)
            }

            )
                yield return flags.Concat(new object[] { shape.Item1, shape.Item2 }).ToArray();
    }

    [Theory]
    [MemberData(nameof(RectangularDiagonalCases))]
    public void AllDiagonalProductsRespectRectangularDimensionsAndZeroPseudoinverses(bool left, bool inverse, bool complexMatrix, bool complexVector, int rows, int cols)
    {
        var a = (Array)Operand(complexMatrix ? typeof(Complex32[,]) : typeof(float[,]), 1, rows, cols);
        var v = (Array)Operand(complexVector ? typeof(Complex32[]) : typeof(float[]), 2, 1, left ? rows : cols);
        if (v.Length > 0)
            v.SetValue(complexVector ? (object)new Complex32(0, 0) : 0f, 0);
        var result = (Array)(left ? Call("Dot", v, a, inverse) : Call("Dot", a, v, inverse));
        Assert.Equal(rows, result.GetLength(0));
        Assert.Equal(cols, result.GetLength(1));
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
            {
                Complex diagonal = Value(v, 0, left ? i : j);
                Complex expected = diagonal == 0 ? 0 : inverse ? Value(a, i, j) / diagonal : Value(a, i, j) * diagonal;
                Check(expected, result.GetValue(i, j)!);
            }
    }
}
