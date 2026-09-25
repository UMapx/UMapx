using System.Numerics;
using UMapx.Core;
using Xunit;
using static UMapx.Tests.MatrixTestSupport;
using static UMapx.Tests.NumericAssert;

namespace UMapx.Tests;

[Trait("Category", "Matrix")]
public class MatrixStructureTests
{
    public static IEnumerable<object[]> ShapeCases()
    {
        foreach (bool complex in new[]
        {
            false,
            true
        }

        )
            foreach (var shape in new[]
            {
                (3, 5),
                (5, 3),
                (4, 4)
            }

            )
                foreach (var direction in Enum.GetValues<Direction>())
                    yield return new object[]
                    {
                        complex,
                        shape.Item1,
                        shape.Item2,
                        direction
                    };
    }

    [Theory]
    [MemberData(nameof(ShapeCases))]
    public void SwappingRowsAndColumnsPermutesEveryEntry(bool complex, int rows, int cols, Direction direction)
    {
        var a = (Array)Operand(complex ? typeof(Complex32[,]) : typeof(float[,]), 1, rows, cols);
        var original = (Array)a.Clone();
        Invoke(typeof(Matrice).GetMethod("Swap", new[] { a.GetType(), typeof(int), typeof(int), typeof(Direction) })!, a, 0, 1, direction);
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
            {
                int ii = direction != Direction.Vertical && i < 2 ? 1 - i : i, jj = direction != Direction.Horizontal && j < 2 ? 1 - j : j;
                Check(Value(original, ii, jj), a.GetValue(i, j)!);
            }
    }

    [Theory]
    [MemberData(nameof(ShapeCases))]
    public void MatrixSlicesAndFiniteDifferencesRespectTheSelectedAxes(bool complex, int rows, int cols, Direction direction)
    {
        var a = (Array)Operand(complex ? typeof(Complex32[,]) : typeof(float[,]), 1, rows, cols);
        var slice = (Array)Invoke(typeof(Matrice).GetMethod("Remove", new[] { a.GetType(), typeof(int), typeof(int), typeof(Direction) })!, a, 1, 2, direction);
        int sr = direction == Direction.Vertical ? rows : 2, sc = direction == Direction.Horizontal ? cols : 2;
        Assert.Equal(sr, slice.GetLength(0));
        Assert.Equal(sc, slice.GetLength(1));
        for (int i = 0; i < sr; i++)
            for (int j = 0; j < sc; j++)
                Check(Value(a, i + (direction == Direction.Vertical ? 0 : 1), j + (direction == Direction.Horizontal ? 0 : 1)), slice.GetValue(i, j)!);
        // Diff uses Horizontal for adjacent columns, Vertical for adjacent rows.
        var d = (Array)Invoke(typeof(Matrice).GetMethod("Diff", new[] { a.GetType(), typeof(int), typeof(Direction), typeof(bool) })!, a, 1, direction, false);
        int dr = rows - (direction == Direction.Horizontal ? 0 : 1), dc = cols - (direction == Direction.Vertical ? 0 : 1);
        Assert.Equal(dr, d.GetLength(0));
        Assert.Equal(dc, d.GetLength(1));
        for (int i = 0; i < dr; i++)
            for (int j = 0; j < dc; j++)
            {
                Complex expected = direction == Direction.Horizontal ? Value(a, i, j + 1) - Value(a, i, j) : direction == Direction.Vertical ? Value(a, i + 1, j) - Value(a, i, j) : Value(a, i + 1, j + 1) - Value(a, i + 1, j) - Value(a, i, j + 1) + Value(a, i, j);
                Check(expected, d.GetValue(i, j)!);
            }
    }

    public static IEnumerable<object[]> MeshCases()
    {
        foreach (bool x in new[]
        {
            false,
            true
        }

        )
            foreach (bool y in new[]
            {
                false,
                true
            }

            )
                foreach (var size in new[]
                {
                    (3, 5),
                    (5, 3)
                }

                )
                    yield return new object[]
                    {
                        x,
                        y,
                        size.Item1,
                        size.Item2
                    };
    }

    [Theory]
    [MemberData(nameof(MeshCases))]
    public void MeshEvaluationUsesBothIndependentCoordinates(bool complexX, bool complexY, int nx, int ny)
    {
        var x = Operand(complexX ? typeof(Complex32[]) : typeof(float[]), 1, 1, nx);
        var y = Operand(complexY ? typeof(Complex32[]) : typeof(float[]), 2, 1, ny);
        Delegate function = complexX || complexY ? (Delegate)new IMeshComplex32((a, b) => a + 2 * b) : new IMeshFloat((a, b) => a + 2 * b);
        var result = (Array)Invoke(typeof(Matrice).GetMethod("Compute", new[] { x.GetType(), y.GetType(), function.GetType() })!, x, y, function);
        Assert.Equal(nx, result.GetLength(0));
        Assert.Equal(ny, result.GetLength(1));
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                Check(Value(x, 0, i) + 2 * Value(y, 0, j), result.GetValue(i, j)!);
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void ElementMappingAndRowColumnReplacementPreserveCoordinates(bool complex)
    {
        if (complex)
        {
            var a = (Complex32[,])Operand(typeof(Complex32[,]), 1);
            var mapped = a.Compute(z => z * z + 1);
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 5; j++)
                    Close(Value(a, i, j) * Value(a, i, j) + 1, mapped[i, j]);
            var row = a.GetRow(1);
            var col = a.GetCol(2);
            for (int j = 0; j < 5; j++)
                Close(Value(a, 1, j), row[j]);
            for (int i = 0; i < 3; i++)
                Close(Value(a, i, 2), col[i]);
            a = a.SetRow(row, 0);
            a = a.SetCol(col, 4);
            for (int j = 0; j < 4; j++)
                Close((Complex)row[j], a[0, j]);
            for (int i = 0; i < 3; i++)
                Close((Complex)col[i], a[i, 4]);
            var v = row.Compute(z => z * z + 1);
            for (int i = 0; i < v.Length; i++)
                Close((Complex)row[i] * (Complex)row[i] + 1, v[i]);
        }
        else
        {
            var a = (float[,])Operand(typeof(float[,]), 1);
            var mapped = a.Compute(z => z * z + 1);
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 5; j++)
                    Close(a[i, j] * a[i, j] + 1, mapped[i, j]);
            var row = a.GetRow(1);
            var col = a.GetCol(2);
            for (int j = 0; j < 5; j++)
                Close(a[1, j], row[j]);
            for (int i = 0; i < 3; i++)
                Close(a[i, 2], col[i]);
            a = a.SetRow(row, 0);
            a = a.SetCol(col, 4);
            for (int j = 0; j < 4; j++)
                Close(row[j], a[0, j]);
            for (int i = 0; i < 3; i++)
                Close(col[i], a[i, 4]);
            var v = row.Compute(z => z * z + 1);
            for (int i = 0; i < v.Length; i++)
                Close(row[i] * row[i] + 1, v[i]);
        }
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void ArrayExtensionPreservesTheCenteredOriginalAndConstantBorders(bool complex)
    {
        var a = (Array)Operand(complex ? typeof(Complex32[,]) : typeof(float[,]), 1, 3, 5);
        var extended = (Array)Invoke(typeof(Matrice).GetMethod("Extend", new[] { a.GetType(), typeof(int), typeof(int) })!, a, 7, 9);
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 5; j++)
                Check(Value(a, i, j), extended.GetValue(i + 2, j + 2)!);
        var v = Array.CreateInstance(complex ? typeof(Complex32) : typeof(float), 5);
        for (int i = 0; i < 5; i++)
            v.SetValue(complex ? (object)new Complex32(2, 3) : 2f, i);
        var padded = (Array)Invoke(typeof(Matrice).GetMethod("Extend", new[] { v.GetType(), typeof(int) })!, v, 9);
        foreach (var value in padded)
            Check(Value(v), value!);
    }

    public static IEnumerable<object[]> StructuredCases()
    {
        foreach (string name in new[]
        {
            "Exchange",
            "Lehmer",
            "Hilbert",
            "GCD",
            "Stirling"
        }

        )
            foreach (int n in new[]
            {
                1,
                3,
                8
            }

            )
                yield return new object[]
                {
                    name,
                    n
                };
    }

    [Theory]
    [MemberData(nameof(StructuredCases))]
    public void ClassicalMatricesMatchTheirIntegerOrRationalDefinitions(string name, int n)
    {
        foreach (bool second in name == "Stirling" ? new[]
        {
            false,
            true
        }

        : new[]
        {
            false
        }

        )
        {
            var actual = name == "Stirling" ? Matrice.Stirling(n, second) : (float[,])Invoke(typeof(Matrice).GetMethod(name, new[] { typeof(int) })!, n);
            var stirling = new long[n, n];
            stirling[0, 0] = 1;
            for (int i = 1; i < n; i++)
                for (int j = 1; j <= i; j++)
                    stirling[i, j] = stirling[i - 1, j - 1] + (second ? j : i - 1) * stirling[i - 1, j];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    Close(name switch
                    {
                        "Exchange" => i + j == n - 1 ? 1 : 0,
                        "Lehmer" => (double)(Math.Min(i, j) + 1) / (Math.Max(i, j) + 1),
                        "Hilbert" => 1.0 / (i + j + 1),
                        "GCD" => (double)BigInteger.GreatestCommonDivisor(i + 1, j + 1),
                        _ => stirling[i, j]
                    }, actual[i, j]);
        }
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void MatrixPredicatesRecognizeConstructedExamples(bool complex)
    {
        if (complex)
        {
            Complex32[,] hermitian =
            {
                {
                    2,
                    new(1, 3)
                },
                {
                    new(1, -3),
                    4
                }
            }, skew =
            {
                {
                    new(0, 2),
                    new(1, 3)
                },
                {
                    new(-1, 3),
                    new(0, 4)
                }
            };
            Assert.True(hermitian.IsSymmetric());
            Assert.False(hermitian.IsSkewSymmetric());
            Assert.True(skew.IsSkewSymmetric());
            Assert.False(skew.IsSymmetric());
            Assert.True(new Complex32[3, 5].IsDiagonal());
            Assert.False(hermitian.IsDiagonal());
            Assert.True(new Complex32[1, 3].IsVector());
            Assert.False(hermitian.IsVector());
            Assert.True(hermitian.IsSquare());
            Assert.False(new Complex32[3, 5].IsSquare());
            Assert.True(hermitian.IsEquals((Complex32[,])hermitian.Clone()));
            Assert.False(hermitian.IsEquals(skew));
        }
        else
        {
            float[,] symmetric =
            {
                {
                    2,
                    3
                },
                {
                    3,
                    4
                }
            }, skew =
            {
                {
                    0,
                    3
                },
                {
                    -3,
                    0
                }
            };
            Assert.True(symmetric.IsSymmetric());
            Assert.False(symmetric.IsSkewSymmetric());
            Assert.True(skew.IsSkewSymmetric());
            Assert.False(skew.IsSymmetric());
            Assert.True(new float[3, 5].IsDiagonal());
            Assert.False(symmetric.IsDiagonal());
            Assert.True(new float[1, 3].IsVector());
            Assert.False(symmetric.IsVector());
            Assert.True(symmetric.IsSquare());
            Assert.False(new float[3, 5].IsSquare());
            Assert.True(symmetric.IsEquals((float[,])symmetric.Clone()));
            Assert.False(symmetric.IsEquals(skew));
            Assert.True(symmetric.IsNonNegative());
            Assert.False(skew.IsNonNegative());
        }
    }

    public static IEnumerable<object[]> SpatialCases()
    {
        foreach (bool complex in new[]
        {
            false,
            true
        }

        )
            foreach (string op in new[]
            {
                "Transpose",
                "FlipHorizontal",
                "FlipVertical",
                "FlipBoth",
                "Shift",
                "Crop",
                "Merge",
                "Reshape",
                "DiffHorizontal",
                "DiffVertical",
                "DiffBoth"
            }

            )
                foreach (var size in new[]
                {
                    (3, 5),
                    (5, 3),
                    (4, 4)
                }

                )
                    yield return new object[]
                    {
                        op,
                        complex,
                        size.Item1,
                        size.Item2
                    };
    }

    [Theory]
    [MemberData(nameof(SpatialCases))]
    public void RectangularArrayOperationsMatchIndexDefinitions(string operation, bool complex, int rows, int cols)
    {
        var a = (Array)Operand(complex ? typeof(Complex32[,]) : typeof(float[,]), 1, rows, cols);
        Array result;
        object Call(string name, params object[] tail)
        {
            var args = new[]
            {
                (object)a
            }.Concat(tail).ToArray();
            return Invoke(typeof(Matrice).GetMethod(name, args.Select(x => x.GetType()).ToArray())!, args);
        }

        if (operation == "Reshape")
        {
            result = (Array)Call("Reshape", rows * cols);
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                    Check(Value(a, i, j), result.GetValue(j * rows + i)!);
            var back = (Array)Invoke(typeof(Matrice).GetMethod("Reshape", new[] { result.GetType(), typeof(int) })!, result, rows);
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                    Check(Value(a, i, j), back.GetValue(i, j)!);
            return;
        }

        if (operation == "Transpose")
            result = (Array)Call("Transpose");
        else if (operation.StartsWith("Flip"))
            result = (Array)Call("Flip", Enum.Parse<Direction>(operation[4..]));
        else if (operation.StartsWith("Diff"))
            result = (Array)Call("Diff", 1, Enum.Parse<Direction>(operation[4..]), false);
        else if (operation == "Shift")
            result = (Array)Call("Shift", 1, -2);
        else if (operation == "Crop")
            result = (Array)Call("Crop", 1, 1, 2, 2, true);
        else
        {
            // A constant patch isolates placement from resampling.
            var patch = Array.CreateInstance(a.GetType().GetElementType()!, 2, 2);
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++)
                    patch.SetValue(complex ? (object)new Complex32(2, .5f) : 2f, i, j);
            result = (Array)Call("Merge", patch, 1, 1, 2, 2);
        }

        for (int i = 0; i < result.GetLength(0); i++)
            for (int j = 0; j < result.GetLength(1); j++)
            {
                Complex expected = operation switch
                {
                    "Transpose" => Value(a, j, i),
                    "FlipHorizontal" => Value(a, i, cols - 1 - j),
                    "FlipVertical" => Value(a, rows - 1 - i, j),
                    "FlipBoth" => Value(a, rows - 1 - i, cols - 1 - j),
                    "Shift" => Value(a, (i - 1 + rows) % rows, (j + 2) % cols),
                    "Crop" => Value(a, i + 1, j + 1),
                    "Merge" => i >= 1 && i < 3 && j >= 1 && j < 3 ? new Complex(2, complex ? .5 : 0) : Value(a, i, j),
                    "DiffHorizontal" => Value(a, i, j + 1) - Value(a, i, j),
                    "DiffVertical" => Value(a, i + 1, j) - Value(a, i, j),
                    "DiffBoth" => Value(a, i + 1, j + 1) - Value(a, i + 1, j) - Value(a, i, j + 1) + Value(a, i, j),
                    _ => throw new Exception()
                };
                Check(expected, result.GetValue(i, j)!);
            }
    }

    static int Wrap(long index, int length) => (int)((index % length + length) % length);
    public static IEnumerable<object[]> ShiftCases()
    {
        foreach (bool complex in new[]
        {
            false,
            true
        }

        )
            foreach (var shape in new[]
            {
                (0, 3),
                (3, 0),
                (1, 7),
                (7, 1),
                (2, 5),
                (5, 3)
            }

            )
                foreach (int shift in new[]
                {
                    int.MinValue,
                    -123456789,
                    -1,
                    0,
                    1,
                    int.MaxValue
                }

                )
                    yield return new object[]
                    {
                        complex,
                        shape.Item1,
                        shape.Item2,
                        shift
                    };
    }

    [Theory]
    [MemberData(nameof(ShiftCases))]
    public void ShiftsWrapBothAxesWithoutIntegerOverflow(bool complex, int rows, int cols, int shift)
    {
        var a = (Array)Operand(complex ? typeof(Complex32[,]) : typeof(float[,]), 1, rows, cols);
        var result = (Array)Call("Shift", a, shift, shift);
        Assert.Equal(rows, result.GetLength(0));
        Assert.Equal(cols, result.GetLength(1));
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                Check(Value(a, Wrap((long)i - shift, rows), Wrap((long)j - shift, cols)), result.GetValue(i, j)!);
        var v = (Array)Operand(complex ? typeof(Complex32[]) : typeof(float[]), 1, 1, cols);
        var shifted = (Array)Call("Shift", v, shift);
        for (int j = 0; j < cols; j++)
            Check(Value(v, 0, Wrap((long)j - shift, cols)), shifted.GetValue(j)!);
    }

    public static IEnumerable<object[]> MergeCases()
    {
        foreach (int kind in new[]
        {
            0,
            1,
            2
        }

        )
            foreach (bool vector in new[]
            {
                false,
                true
            }

            )
                foreach (var offset in new[]
                {
                    (-2, -1),
                    (0, 0),
                    (1, 2),
                    (3, 6),
                    (4, 7),
                    (-10, 2),
                    (int.MinValue, int.MinValue),
                    (int.MaxValue, int.MaxValue)
                }

                )
                    yield return new object[]
                    {
                        kind,
                        vector,
                        offset.Item1,
                        offset.Item2
                    };
    }

    [Theory]
    [MemberData(nameof(MergeCases))]
    public void MergeClipsThePatchIntersectionAndPreservesBothInputs(int kind, bool vector, int top, int left)
    {
        Type Element(bool complex) => vector ? complex ? typeof(Complex32[]) : typeof(float[]) : complex ? typeof(Complex32[,]) : typeof(float[,]);
        var a = (Array)Operand(Element(kind != 0), 1, 4, 7);
        var b = (Array)Operand(Element(kind == 2), 8, 3, 4);
        var original = (Array)a.Clone();
        var patch = (Array)b.Clone();
        var result = (Array)(vector ? Call("Merge", a, b, left, 4) : Call("Merge", a, b, top, left, 3, 4));
        Assert.NotSame(a, result);
        for (int i = 0; i < (vector ? 1 : 4); i++)
            for (int j = 0; j < 7; j++)
            {
                long pi = vector ? 0 : (long)i - top, pj = (long)j - left;
                Complex expected = pi >= 0 && pi < (vector ? 1 : 3) && pj >= 0 && pj < 4 ? Value(patch, (int)pi, (int)pj) : Value(original, i, j);
                Check(expected, vector ? result.GetValue(j)! : result.GetValue(i, j)!);
                Check(Value(original, i, j), vector ? a.GetValue(j)! : a.GetValue(i, j)!);
            }

        Assert.Equal(patch.Cast<object>(), b.Cast<object>());
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void EmptyPatchesAreIdentityInsertions(bool complex)
    {
        var a = (Array)Operand(complex ? typeof(Complex32[,]) : typeof(float[,]), 1);
        foreach (var size in new[]
        {
            (0, 2),
            (2, 0),
            (0, 0)
        }

        )
        {
            var b = Array.CreateInstance(a.GetType().GetElementType()!, size.Item1, size.Item2);
            var result = (Array)Call("Merge", a, b, 1, 1, size.Item1, size.Item2);
            Assert.Equal(a.Cast<object>(), result.Cast<object>());
        }
    }

    [Theory]
    [InlineData(false, false)]
    [InlineData(false, true)]
    [InlineData(true, false)]
    [InlineData(true, true)]
    public void NonseparableMeshesEvaluateEveryCartesianPair(bool complexX, bool complexY)
    {
        foreach (var size in new[]
        {
            (0, 4),
            (3, 0),
            (1, 7),
            (5, 3)
        }

        )
        {
            var x = (Array)Operand(complexX ? typeof(Complex32[]) : typeof(float[]), 2, 1, size.Item1);
            var y = (Array)Operand(complexY ? typeof(Complex32[]) : typeof(float[]), 6, 1, size.Item2);
            Delegate function = complexX || complexY ? (Delegate)new IMeshComplex32((a, b) => a * b + a) : new IMeshFloat((a, b) => a * b + a);
            var result = (Array)Call("Compute", x, y, function);
            Assert.Equal(size.Item1, result.GetLength(0));
            Assert.Equal(size.Item2, result.GetLength(1));
            for (int i = 0; i < size.Item1; i++)
                for (int j = 0; j < size.Item2; j++)
                    Check(Value(x, 0, i) * Value(y, 0, j) + Value(x, 0, i), result.GetValue(i, j)!);
        }
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void SwappingTheLastAxisIndexIsAnInvolutionForRectangularMatrices(bool complex)
    {
        foreach (var shape in new[]
        {
            (1, 7),
            (7, 1),
            (3, 5),
            (5, 3)
        }

        )
            foreach (Direction direction in Enum.GetValues<Direction>())
            {
                var a = (Array)Operand(complex ? typeof(Complex32[,]) : typeof(float[,]), 1, shape.Item1, shape.Item2);
                var original = (Array)a.Clone();
                int last = (direction == Direction.Horizontal ? shape.Item1 : direction == Direction.Vertical ? shape.Item2 : Math.Min(shape.Item1, shape.Item2)) - 1;
                Call("Swap", a, 0, last, direction);
                for (int i = 0; i < shape.Item1; i++)
                    for (int j = 0; j < shape.Item2; j++)
                    {
                        int row = direction != Direction.Vertical && (i == 0 || i == last) ? last - i : i;
                        int col = direction != Direction.Horizontal && (j == 0 || j == last) ? last - j : j;
                        Check(Value(original, row, col), a.GetValue(i, j)!);
                    }

                Call("Swap", a, 0, last, direction);
                Assert.Equal(original.Cast<object>(), a.Cast<object>());
            }
    }

    public static IEnumerable<object[]> StructureCases()
    {
        foreach (bool complex in new[]
        {
            false,
            true
        }

        )
            foreach (string name in new[]
            {
                "Vander",
                "Toeplitz",
                "Hankeli",
                "Hankel",
                "Circulant",
                "Symmetric",
                "Companion",
                "Diag"
            }

            )
                yield return new object[]
                {
                    complex,
                    name
                };
    }

    [Theory]
    [MemberData(nameof(StructureCases))]
    public void StructuredMatricesMatchTheirEntryDefinitions(bool complex, string name)
    {
        var v = MatrixTestSupport.Operand(complex ? typeof(Complex32[]) : typeof(float[]), 1, 1, 6);
        var actual = (Array)MatrixTestSupport.Invoke(typeof(Matrice).GetMethod(name, new[] { v.GetType() })!, v);
        int n = name == "Hankel" ? 3 : 6;
        Assert.Equal(n, actual.GetLength(0));
        Assert.Equal(n, actual.GetLength(1));
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
            {
                Complex expected = name switch
                {
                    "Vander" => Complex.Pow(MatrixTestSupport.Value(v, 0, i), j),
                    "Toeplitz" or "Symmetric" => MatrixTestSupport.Value(v, 0, Math.Abs(i - j)),
                    "Hankeli" => i + j < 6 ? MatrixTestSupport.Value(v, 0, i + j) : 0,
                    "Hankel" => MatrixTestSupport.Value(v, 0, i + j),
                    "Circulant" => MatrixTestSupport.Value(v, 0, (j - i + 6) % 6),
                    "Diag" => i == j ? MatrixTestSupport.Value(v, 0, i) : 0,
                    _ => j == 5 ? -MatrixTestSupport.Value(v, 0, i) : i == j + 1 ? 1 : 0
                };
                MatrixTestSupport.Check(expected, actual.GetValue(i, j)!);
            }
    }

    [Theory]
    [InlineData(1)]
    [InlineData(3)]
    [InlineData(5)]
    [InlineData(9)]
    public void MagicSquaresHaveTheCorrectSetAndEveryLineHasTheSameSum(int n)
    {
        var a = Matrice.Magic(n);
        Assert.Equal(Enumerable.Range(1, n * n).Select(x => (float)x), a.Cast<float>().OrderBy(x => x));
        double expected = n * (n * n + 1) / 2.0;
        for (int i = 0; i < n; i++)
        {
            Close(expected, Enumerable.Range(0, n).Sum(j => (double)a[i, j]));
            Close(expected, Enumerable.Range(0, n).Sum(j => (double)a[j, i]));
        }

        Close(expected, Enumerable.Range(0, n).Sum(i => (double)a[i, i]));
        Close(expected, Enumerable.Range(0, n).Sum(i => (double)a[i, n - i - 1]));
    }

    [Theory]
    [InlineData(1, 5)]
    [InlineData(5, 1)]
    [InlineData(3, 7)]
    public void JaggedArrayConversionsAndComponentsPreserveEveryEntry(int rows, int cols)
    {
        var real = Matrix(rows, cols);
        var jagged = real.ToJagged();
        Close(real, jagged.FromJagged());
        var copy = jagged.Copy();
        copy[0][0] += 5;
        Assert.NotEqual(copy[0][0], jagged[0][0]);
        var neg = jagged.Negate();
        var absolute = jagged.Abs();
        var complex = jagged.ToComplex();
        var z = new Complex32[rows, cols];
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
            {
                Close(-real[i, j], neg[i][j]);
                Close(Math.Abs(real[i, j]), absolute[i][j]);
                Close((System.Numerics.Complex)real[i, j], complex[i][j]);
                z[i, j] = new Complex32(real[i, j], i * .2f - j * .3f);
            }

        var zj = z.ToJagged();
        var restored = zj.FromJagged();
        var zn = zj.Negate();
        var zr = zj.Real();
        var zi = zj.Imag();
        var za = zj.Abs();
        var angle = zj.Angle();
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
            {
                var v = (System.Numerics.Complex)z[i, j];
                Close(v, restored[i, j]);
                Close(-v, zn[i][j]);
                Close(v.Real, zr[i][j]);
                Close(v.Imaginary, zi[i][j]);
                Close(v.Magnitude, za[i][j]);
                Close(v.Phase, angle[i][j]);
            }

        foreach (string name in new[]
        {
            "Zero",
            "One",
            "Eye"
        }

        )
        {
            var a = (float[][])typeof(Jagged).GetMethod(name)!.Invoke(null, new object[] { rows, cols })!;
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                    Close(name == "One" || name == "Eye" && i == j ? 1 : 0, a[i][j]);
        }
    }
}
