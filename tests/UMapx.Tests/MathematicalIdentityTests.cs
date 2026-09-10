using System.Numerics;
using UMapx.Analysis;
using UMapx.Core;
using UMapx.Decomposition;
using UMapx.Response;
using UMapx.Transform;
using UMapx.Wavelet;
using Xunit;
using static UMapx.Tests.NumericAssert;

namespace UMapx.Tests;

[Trait("Category", "Identity")]
public class MathematicalIdentityTests
{
    [Theory]
    [InlineData(1)] [InlineData(2)] [InlineData(3)] [InlineData(4)] [InlineData(5)]
    [InlineData(7)] [InlineData(8)] [InlineData(15)] [InlineData(16)] [InlineData(17)] [InlineData(31)]
    public void FftAgreesWithIndependentDftAndInverts(int n)
    {
        var random = new Random(731911 + n);
        var x = Enumerable.Range(0, n).Select(_ => new Complex32((float)(random.NextDouble() * 2 - 1), (float)(random.NextDouble() * 2 - 1))).ToArray();
        var fft = new FastFourierTransform();
        var spectrum = fft.Forward(x);
        var restored = fft.Backward(spectrum);
        for (int k = 0; k < n; k++)
        {
            Complex sum = 0;
            for (int j = 0; j < n; j++) sum += (Complex)x[j] * Complex.Exp(-Complex.ImaginaryOne * 2 * Math.PI * j * k / n);
            Close(sum / Math.Sqrt(n), spectrum[k], 1e-4);
            Close((Complex)x[k], restored[k], 1e-4);
        }
    }

    [Theory]
    [InlineData(Direction.Horizontal)] [InlineData(Direction.Vertical)] [InlineData(Direction.Both)]
    public void RectangularFftInvertsInEveryDirection(Direction direction)
    {
        var x = new Complex32[3, 5];
        for (int i = 0; i < 3; i++) for (int j = 0; j < 5; j++) x[i, j] = new Complex32(i + j * .3f, i * j - 1);
        var fft = new FastFourierTransform(true, direction);
        var actual = fft.Backward(fft.Forward(x));
        for (int i = 0; i < 3; i++) for (int j = 0; j < 5; j++) Close((Complex)x[i, j], actual[i, j], 1e-4);
    }

    [Theory]
    [InlineData(2)] [InlineData(4)] [InlineData(8)] [InlineData(16)]
    public void OrthogonalRealTransformsInvert(int n)
    {
        var x = Enumerable.Range(0, n).Select(i => (float)Math.Sin(i * .37)).ToArray();
        ITransform[] transforms = { new CosineTransform(), new FastCosineTransform(), new SineTransform(), new FastSineTransform(), new HartleyTransform(), new FastHartleyTransform(), new WalshHadamardTransform(), new FastWalshHadamardTransform(), new ChebyshevTransform(), new FastChebyshevTransform() };
        foreach (var transform in transforms) Close(x, transform.Backward(transform.Forward(x)));
    }

    [Theory]
    [InlineData(1, 1)] [InlineData(2, 2)] [InlineData(3, 3)] [InlineData(5, 3)] [InlineData(3, 5)] [InlineData(8, 8)]
    public void SvdReconstructsRectangularMatrices(int m, int n)
    {
        var a = Matrix(m, n); var d = new SVD(a);
        Close(a, Product(Product(d.U, Diagonal(d.S)), Transpose(d.V)));
        Assert.All(d.S, value => Assert.True(float.IsFinite(value) && value >= 0));
    }

    [Theory] [InlineData(2, 2)] [InlineData(5, 3)] [InlineData(8, 8)]
    public void QrReconstructsAndHasOrthonormalColumns(int m, int n)
    {
        var a = Matrix(m, n); var d = new QR(a);
        Close(a, Product(d.Q, d.R));
        Close(Diagonal(Enumerable.Repeat(1f, n).ToArray()), Product(Transpose(d.Q), d.Q));
    }

    [Theory] [InlineData(2)] [InlineData(3)] [InlineData(5)]
    public void SquareMatrixDecompositionsReconstruct(int n)
    {
        var a = Matrix(n, n);
        var lu = new LU(a); Close(a, Product(lu.L, lu.U));
        var ldu = new LDU(a); Close(a, Product(Product(ldu.L, Diagonal(ldu.D)), ldu.U));
        var lq = new LQ(a); Close(a, Product(lq.L, lq.Q));
        var ql = new QL(a); Close(a, Product(ql.Q, ql.L));
        var rq = new RQ(a); Close(a, Product(rq.R, rq.Q));
        var polar = new Polar(a); Close(a, Product(polar.U, polar.P));
        var h = new Hessenberg(a); Close(a, Product(Product(h.P, h.H), Transpose(h.P)));
        var schur = new Schur(a); Close(a, Product(Product(schur.Q, schur.T), Transpose(schur.Q)));
        var evd = new EVD(a); Close(Product(a, evd.V), Product(evd.V, evd.R));
    }

    [Theory] [InlineData(2)] [InlineData(3)] [InlineData(5)]
    public void PositiveDefiniteMatrixDecompositionsReconstruct(int n)
    {
        var r = Matrix(n, n); var a = Product(r, Transpose(r));
        for (int i = 0; i < n; i++) a[i, i] += 1;
        var chol = new Cholesky(a); Close(a, Product(chol.L, Transpose(chol.L)));
        var ldl = new LDL(a); Close(a, Product(Product(ldl.L, Diagonal(ldl.D)), ldl.U));
        var udl = new UDL(a); Close(a, Product(Product(udl.U, Diagonal(udl.D)), udl.L));
    }

    [Theory] [InlineData(2)] [InlineData(3)] [InlineData(5)] [InlineData(8)] [InlineData(9)]
    public void SimpsonIntegratesCubicsExactlyForEvenAndOddSubintervalCounts(int n)
    {
        var integrator = new Integration(IntegrationMethod.Simpson);
        Close(.25, integrator.Compute((IFloat)(x => x * x * x), 0, 1, n));
        var y = Enumerable.Range(0, n + 1).Select(i => (float)Math.Pow((double)i / n, 3)).ToArray();
        Close(.25, integrator.Compute(y, 0, 1, n + 1));
    }

    [Theory]
    [InlineData(InterpolationMethod.Lagrange)] [InlineData(InterpolationMethod.Newton)] [InlineData(InterpolationMethod.Barycentric)]
    public void PolynomialInterpolationReproducesQuadratic(InterpolationMethod method)
    {
        var x = new[] { 0f, 1f, 2f }; var y = new[] { 1f, 4f, 9f };
        Close(6.25, new Interpolation(method).Compute(x, y, 1.5f));
    }

    [Theory] [InlineData(1)] [InlineData(2)] [InlineData(3)]
    public void OrthonormalWaveletBanksReconstruct(int levels)
    {
        var x = Enumerable.Range(0, 64).Select(i => (float)Math.Sin(i * .37)).ToArray();
        foreach (var bank in new[] { WaveletPacket.Haar, WaveletPacket.D2, WaveletPacket.D4 })
        {
            var wavelet = new WaveletDecomposition(bank, levels);
            Close(x, wavelet.Backward(wavelet.Forward(x)));
        }
    }

    [Fact]
    public void FirImpulseResponseEqualsCoefficients()
    {
        var fir = new FIR(new[] { 1f, 2f, 3f });
        Close(new[] { 1f, 2f, 3f, 0f }, fir.Reaction(new[] { 1f, 0f, 0f, 0f }));
    }

    [Theory] [InlineData(2)] [InlineData(3)] [InlineData(7)]
    public void BasicDistancesAgreeWithDirectDefinitions(int n)
    {
        var p = Enumerable.Range(0, n).Select(i => (float)i / n).ToArray();
        var q = p.Reverse().ToArray();
        var deltas = p.Zip(q, (a, b) => Math.Abs((double)a - b)).ToArray();
        Close(Math.Sqrt(deltas.Sum(d => d * d)), new UMapx.Distance.Euclidean().Compute(p, q));
        Close(deltas.Sum(), new UMapx.Distance.Manhattan().Compute(p, q));
        Close(deltas.Max(), new UMapx.Distance.Chebyshev().Compute(p, q));
        Close(Math.Sqrt(deltas.Sum(d => d * d)), new UMapx.Distance.Minkowski(2).Compute(p, q));
    }
}
