using UMapx.Core;
using Xunit;
using static UMapx.Tests.NumericAssert;

namespace UMapx.Tests;

[Trait("Category", "Core")]
public class KernelTests
{
    public static IEnumerable<object[]> KernelCases()
    {
        foreach (string name in new[]
        {
            "Bicubic",
            "Gaussian",
            "Lanczos",
            "Uniform",
            "Triangular",
            "Trapezoid",
            "Epanechnikov",
            "Quartic",
            "Triweight",
            "Tricube",
            "Cosine",
            "Logistic",
            "Sigmoid",
            "Silverman"
        }

        )
            foreach (float x in new[]
            {
                -10f,
                -2f,
                -1f,
                -.75f,
                -.5f,
                0f,
                .25f,
                .5f,
                .75f,
                1f,
                2f,
                10f
            }

            )
                yield return new object[]
                {
                    name,
                    x
                };
    }

    [Theory]
    [MemberData(nameof(KernelCases))]
    public void KernelValuesMatchTheirIndependentDefinitions(string name, float x)
    {
        double a = Math.Abs(x), q = 1 - x * x, expected = name switch
        {
            "Bicubic" => a <= 1 ? 1.5 * a * a * a - 2.5 * a * a + 1 : a < 2 ? -.5 * a * a * a + 2.5 * a * a - 4 * a + 2 : 0,
            "Gaussian" => Math.Exp(-x * x / 2.0),
            "Lanczos" => x == 0 ? 1 : a >= 1 ? 0 : Math.Pow(Math.Sin(Math.PI * x) / (Math.PI * x), 2),
            "Uniform" => a <= 1 ? .5 : 0,
            "Triangular" => Math.Max(0, 1 - a),
            "Trapezoid" => a <= .5 ? 1 : a < 1 ? 2 * (1 - a) : 0,
            "Epanechnikov" => a <= 1 ? .75 * q : 0,
            "Quartic" => a <= 1 ? 15.0 / 16 * q * q : 0,
            "Triweight" => a <= 1 ? 35.0 / 32 * q * q * q : 0,
            "Tricube" => a <= 1 ? 70.0 / 81 * Math.Pow(1 - a * a * a, 3) : 0,
            "Cosine" => a <= 1 ? Math.PI / 4 * Math.Cos(Math.PI * x / 2) : 0,
            "Logistic" => 1 / (2 + 2 * Math.Cosh(x)),
            "Sigmoid" => 1 / (Math.PI * Math.Cosh(x)),
            "Silverman" => .5 * Math.Exp(-a / Math.Sqrt(2)) * Math.Sin(a / Math.Sqrt(2) + Math.PI / 4),
            _ => throw new Exception()
        };
        Close(expected, (float)typeof(Kernel).GetMethod(name, new[] { typeof(float) })!.Invoke(null, new object[] { x })!);
        if (name == "Gaussian")
            Close(Math.Exp(-x * x / (2 * 2.3f * 2.3f)), Kernel.Gaussian(x, 2.3f));
        if (name == "Lanczos")
            Close(x == 0 ? 1 : a >= 3 ? 0 : 3 * Math.Sin(Math.PI * x) * Math.Sin(Math.PI * x / 3) / (Math.PI * Math.PI * x * x), Kernel.Lanczos(x, 3));
    }

    [Theory]
    [InlineData(3, 5)]
    [InlineData(5, 3)]
    [InlineData(5, 5)]
    public void GaussianOperatorsHaveTheExpectedSymmetryAndWeights(int rows, int cols)
    {
        var gaussian = Operator.Gaussian(rows, cols, 1.2f, .8f);
        for (int y = 0; y < rows; y++)
            for (int x = 0; x < cols; x++)
                Close(gaussian[rows - 1 - y, cols - 1 - x], gaussian[y, x]);
        double center = gaussian[rows / 2, cols / 2];
        for (int y = 0; y < rows; y++)
            for (int x = 0; x < cols; x++)
                Close(Math.Exp(-.5 * (Math.Pow((y - rows / 2) / 1.2f, 2) + Math.Pow((x - cols / 2) / .8f, 2))), gaussian[y, x] / center, 1e-5);
    }

    [Theory]
    [InlineData("Roberts")]
    [InlineData("Prewitt")]
    [InlineData("Sobel")]
    [InlineData("Scharr")]
    [InlineData("Laplacian")]
    [InlineData("LaplacianDiagonal")]
    [InlineData("LaplacianInvert")]
    public void DerivativeOperatorsAnnihilateConstantImages(string name)
    {
        var kernel = (float[,])typeof(Operator).GetMethod(name, Type.EmptyTypes)!.Invoke(null, null)!;
        Close(0, kernel.Cast<float>().Sum());
    }
}
