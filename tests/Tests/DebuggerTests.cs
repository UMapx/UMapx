using System.Globalization;
using System.Text.RegularExpressions;
using UMapx.Core;
using Xunit;

namespace UMapx.Tests;

[CollectionDefinition("Debugger console", DisableParallelization = true)]
public class DebuggerConsoleCollection { }

[Collection("Debugger console")]
[Trait("Category", "Core")]
public class DebuggerTests
{
    private static string Capture(Action action)
    {
        var previous = Console.Out;
        using var writer = new StringWriter(CultureInfo.InvariantCulture) { NewLine = "\n" };
        try
        {
            Console.SetOut(writer);
            action();
            return writer.ToString();
        }
        finally { Console.SetOut(previous); }
    }

    private static int ConsoleWidth()
    {
        try { return Math.Max(20, Console.WindowWidth); }
        catch { return 80; }
    }

    [Fact]
    public void PrintObjectWritesItsTextAndABlankLine()
    {
        Assert.Equal("sample\n\n", Capture(() => Debugger.Print((object)"sample")));
        Assert.Equal("\n\n", Capture(() => Debugger.Print((object)null!)));
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void EmptyAndNullVectorsHaveAnExplicitEmptyHeader(bool vertical)
    {
        Assert.Equal("Empty vector 0\n", Capture(() => Debugger.Print(Array.Empty<int>(), vertical)));
        Assert.Equal("Empty vector 0\n", Capture(() => Debugger.Print((int[])null!, vertical)));
    }

    [Fact]
    public void VectorDefaultsToHorizontalAndVerticalValuesAreRightAligned()
    {
        var values = new[] { 1, -20, 300 };
        Assert.Equal("Vector 3\nColumns 1 through 3\n  1   -20   300\n\n",
            Capture(() => Debugger.Print(values)));
        Assert.Equal("Vector 3\nColumns 1 through 1\n    1\n  -20\n  300\n\n",
            Capture(() => Debugger.Print(values, true)));
        Assert.Equal(new[] { 1, -20, 300 }, values);
    }

    [Fact]
    public void MatrixAndRectangularJaggedArrayAlignEachColumn()
    {
        var matrix = new[,] { { 1, -20 }, { 300, 4 } };
        var jagged = new[] { new[] { 1, -20 }, new[] { 300, 4 } };
        const string expected = "Matrix 2 x 2\nColumns 1 through 2\n    1   -20\n  300     4\n\n";
        Assert.Equal(expected, Capture(() => Debugger.Print(matrix)));
        Assert.Equal(expected, Capture(() => Debugger.Print(jagged)));
        Assert.Equal(new[] { 1, -20, 300, 4 }, matrix.Cast<int>());
        Assert.Equal(new[] { 1, -20, 300, 4 }, jagged.SelectMany(row => row));
    }

    [Theory]
    [InlineData(0, 0)]
    [InlineData(0, 3)]
    [InlineData(2, 0)]
    public void EmptyMatricesPreserveBothDimensions(int rows, int columns)
    {
        Assert.Equal($"Empty matrix {rows} x {columns}\n",
            Capture(() => Debugger.Print(new int[rows, columns])));
    }

    [Fact]
    public void NullMatricesHaveAnExplicitEmptyHeader()
    {
        Assert.Equal("Empty matrix 0 x 0\n", Capture(() => Debugger.Print((int[,])null!)));
        Assert.Equal("Empty matrix 0 x 0\n", Capture(() => Debugger.Print((int[][])null!)));
    }

    [Fact]
    public void EmptyJaggedArrayDoesNotAccessItsFirstRow()
    {
        Assert.Equal("Empty matrix 0 x 0\n", Capture(() => Debugger.Print(Array.Empty<int[]>())));
    }

    [Fact]
    public void JaggedArrayWithOnlyEmptyOrNullRowsPreservesRowCount()
    {
        var values = new[] { Array.Empty<int>(), null!, Array.Empty<int>() };
        Assert.Equal("Empty matrix 3 x 0\n", Capture(() => Debugger.Print(values)));
    }

    [Fact]
    public void JaggedArrayPrintsEveryValueAndPadsMissingCells()
    {
        var values = new[] { new[] { 1, 20 }, new[] { 300 }, new[] { 4, 5, 6 } };
        const string expected = "Matrix 3 x 3\nColumns 1 through 3\n    1   20    \n  300         \n    4    5   6\n\n";
        Assert.Equal(expected, Capture(() => Debugger.Print(values)));
        Assert.Equal(new[] { 2, 1, 3 }, values.Select(row => row.Length));
        Assert.Equal(new[] { 1, 20, 300, 4, 5, 6 }, values.SelectMany(row => row));
    }

    [Fact]
    public void JaggedArrayCanStartWithANullOrEmptyRow()
    {
        const string expected = "Matrix 2 x 2\nColumns 1 through 2\n       \n  1   2\n\n";
        Assert.Equal(expected, Capture(() => Debugger.Print(new[] { (int[])null!, new[] { 1, 2 } })));
        Assert.Equal(expected, Capture(() => Debugger.Print(new[] { Array.Empty<int>(), new[] { 1, 2 } })));
    }

    [Fact]
    public void NullVectorCellsAndNullTextArePrintedAsEmptyCells()
    {
        object?[] values = { null, new NullText(), "x" };
        Assert.Equal("Vector 3\nColumns 1 through 3\n        x\n\n",
            Capture(() => Debugger.Print(values)));
        Assert.Equal("Vector 3\nColumns 1 through 1\n   \n   \n  x\n\n",
            Capture(() => Debugger.Print(values, true)));
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void NullMatrixCellsAndNullTextArePrintedAsEmptyCells(bool jagged)
    {
        const string expected = "Matrix 2 x 2\nColumns 1 through 2\n      x\n  y    \n\n";
        var actual = Capture(() =>
        {
            if (jagged) Debugger.Print(new object?[][] { new object?[] { null, "x" }, new object?[] { "y", new NullText() } });
            else Debugger.Print(new object?[,] { { null, "x" }, { "y", new NullText() } });
        });
        Assert.Equal(expected, actual);
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void WideVectorsSplitIntoBlocksWithoutDroppingOversizedValues(bool oversized)
    {
        string longValue = new('x', ConsoleWidth() + (oversized ? 1 : -2));
        var values = new[] { longValue, "tail", "end" };
        Assert.Equal($"Vector 3\nColumns 1 through 1\n  {longValue}\n\nColumns 2 through 3\n  tail   end\n\n",
            Capture(() => Debugger.Print(values)));
    }

    [Theory]
    [InlineData(false, false)]
    [InlineData(false, true)]
    [InlineData(true, false)]
    [InlineData(true, true)]
    public void WideMatricesSplitIntoBlocksAndKeepColumnAlignment(bool jagged, bool oversized)
    {
        string longValue = new('x', ConsoleWidth() + (oversized ? 1 : -2));
        string paddedZ = new string(' ', longValue.Length - 1) + "z";
        string expected = $"Matrix 2 x 3\nColumns 1 through 1\n  {longValue}\n  {paddedZ}\n\n" +
            "Columns 2 through 3\n   x   yy\n  tt    q\n\n";
        var actual = Capture(() =>
        {
            if (jagged) Debugger.Print(new[] { new[] { longValue, "x", "yy" }, new[] { "z", "tt", "q" } });
            else Debugger.Print(new[,] { { longValue, "x", "yy" }, { "z", "tt", "q" } });
        });
        Assert.Equal(expected, actual);
    }

    [Theory]
    [InlineData("en-US", "1.5")]
    [InlineData("ru-RU", "1,5")]
    public void ArrayValuesUseCurrentCulture(string culture, string number)
    {
        var previous = CultureInfo.CurrentCulture;
        try
        {
            CultureInfo.CurrentCulture = CultureInfo.GetCultureInfo(culture);
            Assert.Equal($"Vector 1\nColumns 1 through 1\n  {number}\n\n",
                Capture(() => Debugger.Print(new[] { 1.5f })));
            Assert.Equal($"Matrix 1 x 1\nColumns 1 through 1\n  {number}\n\n",
                Capture(() => Debugger.Print(new[,] { { 1.5f } })));
            Assert.Equal($"Matrix 1 x 1\nColumns 1 through 1\n  {number}\n\n",
                Capture(() => Debugger.Print(new[] { new[] { 1.5f } })));
        }
        finally { CultureInfo.CurrentCulture = previous; }
    }

    [Fact]
    public void InfoOfNullProducesNoOutput()
    {
        Assert.Equal(string.Empty, Capture(() => Debugger.Info(null!)));
    }

    [Fact]
    public void InfoDefaultsToDeclaredInstanceMembersWithoutSignatures()
    {
        CheckInfo(Capture(() => Debugger.Info(new InfoSample())), false, false, false);
    }

    public static IEnumerable<object[]> InfoOptions()
    {
        foreach (bool inherited in new[] { false, true })
        foreach (bool statics in new[] { false, true })
        foreach (bool signatures in new[] { false, true })
            yield return new object[] { inherited, statics, signatures };
    }

    [Theory]
    [MemberData(nameof(InfoOptions))]
    public void InfoHonorsMemberSelectionAndSignatureOptions(bool inherited, bool statics, bool signatures)
    {
        CheckInfo(Capture(() => Debugger.Info(new InfoSample(), inherited, statics, signatures)),
            inherited, statics, signatures);
    }

    private static void CheckInfo(string output, bool inherited, bool statics, bool signatures)
    {
        var properties = new List<string> { "Value", "ThrowsWhenRead" };
        var methods = new List<string> { "Overload()", "Overload(Int32)", "Combine(Int32, String)" };
        if (statics)
        {
            properties.Add("StaticValue");
            methods.Add("StaticMethod()");
        }
        if (inherited)
        {
            properties.Add("InheritedValue");
            methods.AddRange(new[] { "InheritedMethod()", "GetType()", "ToString()", "Equals(Object)", "GetHashCode()" });
            if (statics)
            {
                properties.Add("InheritedStaticValue");
                methods.AddRange(new[] { "InheritedStaticMethod()", "Equals(Object, Object)", "ReferenceEquals(Object, Object)" });
            }
        }
        var expectedMethods = signatures ? methods : methods.Select(method => method[..method.IndexOf('(')]).Distinct();
        var lines = output.Split('\n');
        Assert.Equal(5, lines.Length);
        Assert.Equal(typeof(InfoSample).ToString(), lines[0]);
        Assert.StartsWith("Properties: ", lines[1]);
        Assert.StartsWith("Methods: ", lines[2]);
        Assert.Equal(string.Empty, lines[3]);
        Assert.Equal(string.Empty, lines[4]);
        Assert.Equal(properties.OrderBy(name => name), lines[1][12..].Split(", ").OrderBy(name => name));
        var actualMethods = Regex.Matches(lines[2][9..], @"\w+\([^)]*\)|\w+").Select(match => match.Value).ToArray();
        Assert.Equal("Methods: " + string.Join(", ", actualMethods), lines[2]);
        Assert.Equal(expectedMethods.OrderBy(name => name), actualMethods.OrderBy(name => name));
    }

    [Fact]
    public void TicAndTocReportElapsedMillisecondsWithoutResettingTheStart()
    {
        int beforeStart = Environment.TickCount;
        Debugger.Tic();
        int afterStart = Environment.TickCount;
        for (int i = 0; i < 2; i++)
        {
            int beforeOutput = Environment.TickCount;
            string output = Capture(Debugger.Toc);
            int afterOutput = Environment.TickCount;
            var match = Regex.Match(output, @"\AElapsed time is (-?\d+) milliseconds\n\n\z");
            Assert.True(match.Success, output);
            int elapsed = int.Parse(match.Groups[1].Value, CultureInfo.InvariantCulture);
            Assert.InRange(elapsed, unchecked(beforeOutput - afterStart), unchecked(afterOutput - beforeStart));
        }
    }

    private sealed class NullText
    {
        public override string ToString() => null!;
    }

    public class InfoBase
    {
        public int InheritedValue => 1;
        public static int InheritedStaticValue => 2;
        public void InheritedMethod() { }
        public static void InheritedStaticMethod() { }
    }

    public sealed class InfoSample : InfoBase
    {
        public int Value { get; set; }
        public int ThrowsWhenRead => throw new InvalidOperationException();
        public static int StaticValue => 3;
        public event Action Changed { add { } remove { } }
        public void Overload() { }
        public void Overload(int value) { }
        public void Combine(int count, string text) { }
        public static void StaticMethod() { }
        private void HiddenMethod() { }
        public static InfoSample operator +(InfoSample left, InfoSample right) => left;
    }
}
