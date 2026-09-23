using System.Numerics;
using System.Reflection;
using System.Text.Json;
using UMapx.Core;
using UMapx.Wavelet;
using Xunit;
using static UMapx.Tests.NumericAssert;

namespace UMapx.Tests;

[Trait("Category", "Wavelet")]
public class WaveletRepairTests
{
    public static IEnumerable<object[]> Banks() => WaveletAuditTests.Banks()
        .Where(row => ((string)row[0]).StartsWith("Bior") || ((string)row[0]).StartsWith("CDF"));

    private static WaveletPacket Bank(string name) => (WaveletPacket)typeof(WaveletPacket)
        .GetProperty(name, BindingFlags.Public | BindingFlags.Static)!.GetValue(null)!;

    public static IEnumerable<object[]> ReconstructionCases()
    {
        foreach (var row in Banks())
        foreach (bool normalized in new[] { false, true })
        foreach (int size in new[] { 8, 32, 128 })
            yield return new object[] { row[0], normalized, size };
    }

    [Theory, MemberData(nameof(ReconstructionCases))]
    public void DualBanksReconstructEveryImpulsePhaseAndComplexSignals(string name, bool normalized, int size)
    {
        var d = new WaveletDecomposition(Bank(name), 3, normalized);
        for (int phase = 0; phase < size; phase++)
        {
            var impulse = new float[size];
            impulse[phase] = 1;
            Close(impulse, d.Backward(d.Forward(impulse)), 8e-6);
        }
        var input = Enumerable.Range(0, size).Select(i => new Complex32(
            (float)Math.Sin(.31 * i), (float)Math.Cos(.17 * i))).ToArray();
        var restored = d.Backward(d.Forward(input));
        for (int i = 0; i < size; i++) Close((Complex)input[i], restored[i], 2e-5);
        var image = new Complex32[8, 16];
        for (int y = 0; y < 8; y++)
        for (int x = 0; x < 16; x++) image[y, x] = input[(y * 16 + x) % size];
        var result = d.Backward(d.Forward(image));
        for (int y = 0; y < 8; y++)
        for (int x = 0; x < 16; x++) Close((Complex)image[y, x], result[y, x], 4e-5);
    }

    [Theory, MemberData(nameof(Banks))]
    public void AnalysisWaveletsRejectConstantsAndHaveTheirAdvertisedMoments(string name)
    {
        var bank = Bank(name);
        Close(Math.Sqrt(2), bank.LowPass.Sum(v => (double)v), 2e-6, 0);
        Close(0, bank.HighPass.Sum(v => (double)v), 2e-6, 0);
        Close(0, bank.IHighPass.Sum(v => (double)v), 2e-6, 0);
        // CDF names give spline order / dual order; Bior names put the
        // analysis wavelet's number of vanishing moments first.
        int moments = name[4] - '0';
        if (name == "CDF97") moments = 4; // Nine/seven taps, four moments on each side.
        for (int power = 1; power < moments; power++)
        {
            double sum = 0, scale = 0;
            for (int i = 0; i < bank.HighPass.Length; i++)
            {
                double term = bank.HighPass[i] * Math.Pow(i - (bank.HighPass.Length - 1) / 2.0, power);
                sum += term;
                scale += Math.Abs(term);
            }
            Assert.True(Math.Abs(sum) <= 3e-7 * Math.Max(1, scale), $"{name}: moment {power} = {sum}");
        }
        var bands = new WaveletDecomposition(bank, 1).Forward(Enumerable.Repeat(.375f, 64).ToArray());
        Assert.All(bands[1], value => Close(0, value, 2e-6, 0));
    }

    public static IEnumerable<object[]> MeyerCases()
    {
        using var stream = typeof(WaveletRepairTests).Assembly.GetManifestResourceStream("UMapx.Tests.Data.meyer-hankel.json");
        using var data = JsonDocument.Parse(stream!);
        foreach (var row in data.RootElement.GetProperty("meyer").EnumerateArray())
            yield return new object[] { row.GetProperty("wavelet").GetBoolean(), row.GetProperty("x").GetSingle(), row.GetProperty("expected").GetDouble() };
    }

    [Theory, MemberData(nameof(MeyerCases))]
    public void MeyerAgreesWithIndependentFrequencyDomainQuadrature(bool wavelet, float x, double expected)
    {
        var function = new Meyer();
        Close(expected, wavelet ? function.Wavelet(x) : function.Scaling(x), 2e-7, 2e-7);
    }

    [Theory]
    [InlineData(8, 12, 3)] [InlineData(7, 16, 1)] [InlineData(16, 15, 2)]
    public void MatrixDecompositionRejectsDimensionsThatBecomeOddAtAnActiveLevel(int rows, int columns, int levels)
    {
        var d = new WaveletDecomposition(WaveletPacket.Bior13, levels);
        Assert.Throws<ArgumentException>(() => d.Forward(new float[rows, columns]));
        Assert.Throws<ArgumentException>(() => d.Forward(new Complex32[rows, columns]));
    }
}
