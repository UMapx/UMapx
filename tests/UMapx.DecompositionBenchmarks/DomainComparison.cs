using System.Diagnostics;
using System.Linq.Expressions;
using System.Numerics;
using System.Reflection;
using System.Text.Json;

// Loads either assembly in its own process, with no compile-time dependency on UMapx.
internal static class DomainComparison
{
    public static void Run(string[] args)
    {
        if (args.Length != 6) throw new ArgumentException("Expected assembly label algorithm columns rows domain.");
        string label = args[1], name = args[2], domain = args[5], algorithm = name.Split('-')[0];
        int n = int.Parse(args[3]), rows = int.Parse(args[4]);
        try
        {
            bool real = domain == "real";
            if (!real && domain != "complex") throw new ArgumentException("Domain must be real or complex.");
            var assembly = Assembly.LoadFrom(Path.GetFullPath(args[0]));
            var complexType = assembly.GetType("UMapx.Core.Complex32", true)!;
            var scalar = real ? typeof(float) : complexType;
            var matrixType = scalar.MakeArrayType(2);
            var type = assembly.GetType("UMapx.Decomposition." + algorithm, true)!;
            var method = type.GetMethods().Single(m => m.Name == "Decompose" && m.GetParameters()[0].ParameterType == matrixType);
            var inputA = Expression.Parameter(typeof(object), "a"); var inputB = Expression.Parameter(typeof(object), "b");
            var inputs = method.GetParameters().Select((p, i) => p.ParameterType == matrixType ? (Expression)Expression.Convert(i == 0 ? inputA : inputB, matrixType)
                : p.ParameterType == typeof(int) ? Expression.Constant(algorithm is "Schur" or "QZ" ? 1000 : 100)
                : Expression.Constant(p.DefaultValue, p.ParameterType)).ToArray();
            var output = Expression.Variable(method.ReturnType, "result");
            var fields = method.ReturnType.GetFields().Where(f => f.Name.StartsWith("Item")).OrderBy(f => f.Name).ToArray();
            Expression[] factors = method.ReturnType.IsArray ? new Expression[] { Expression.Convert(output, typeof(object)) }
                : fields.Select(f => (Expression)Expression.Convert(Expression.Field(output, f), typeof(object))).ToArray();
            var run = Expression.Lambda<Func<object, object, object[]>>(Expression.Block(new[] { output },
                Expression.Assign(output, Expression.Call(method, inputs)), Expression.NewArrayInit(typeof(object), factors)), inputA, inputB).Compile();

            bool hermitian = name.EndsWith("SPD") || algorithm == "Householder";
            var a = Matrix(rows, n, 173, real, hermitian); var b = Matrix(rows, n, 917, real, false);
            if (algorithm is "QZ" or "GEVD") for (int i = 0; i < n; i++) b[i, i] += n;
            var ctor = complexType.GetConstructor(new[] { typeof(float), typeof(float) })!;
            object Encode(Complex[,] matrix)
            {
                var result = Array.CreateInstance(scalar, matrix.GetLength(0), matrix.GetLength(1));
                for (int i = 0; i < matrix.GetLength(0); i++) for (int j = 0; j < matrix.GetLength(1); j++)
                    result.SetValue(real ? (object)(float)matrix[i, j].Real : ctor.Invoke(new object[] { (float)matrix[i, j].Real, (float)matrix[i, j].Imaginary }), i, j);
                return result;
            }
            var re = complexType.GetField("Real")!; var im = complexType.GetField("Imag")!;
            Complex Decode(object x) => x is float f ? f : new Complex((float)re.GetValue(x)!, (float)im.GetValue(x)!);
            Complex[,] Work(object x)
            {
                var array = (Array)x;
                var result = new Complex[array.GetLength(0), array.GetLength(1)];
                for (int i = 0; i < result.GetLength(0); i++) for (int j = 0; j < result.GetLength(1); j++) result[i, j] = Decode(array.GetValue(i, j)!);
                return result;
            }
            Complex[] Values(object x) => ((Array)x).Cast<object>().Select(Decode).ToArray();
            object aa = Encode(a), bb = Encode(b);
            a = Work(aa); b = Work(bb);
            object[] result = run(aa, bb);
            double residual;
            double orthogonality = 0;
            void Orthogonal(object factor)
            {
                var q = Work(factor); var identity = Diagonal(Enumerable.Repeat(Complex.One, q.GetLength(1)).ToArray());
                orthogonality = Math.Max(orthogonality, Error(identity, Product(Adjoint(q), q)));
            }
            switch (algorithm)
            {
                case "SVD":
                    residual = Error(a, Product(Product(Work(result[0]), Diagonal(Values(result[1]))), Adjoint(Work(result[2]))));
                    Orthogonal(result[0]); Orthogonal(result[2]); break;
                case "Householder": case "Schur":
                    residual = Error(a, Product(Product(Work(result[0]), Work(result[1])), Adjoint(Work(result[0]))));
                    Orthogonal(result[0]); break;
                case "QR":
                    residual = Error(a, Product(Work(result[0]), Work(result[1]))); Orthogonal(result[0]); break;
                case "Polar": residual = Error(a, Product(Work(result[0]), Work(result[1]))); break;
                case "GSVD":
                    residual = Math.Max(Error(a, Product(Product(Work(result[0]), Diagonal(Values(result[1]))), Work(result[4]))),
                        Error(b, Product(Product(Work(result[2]), Diagonal(Values(result[3]))), Work(result[4]))));
                    Orthogonal(result[0]); Orthogonal(result[2]); break;
                case "QZ":
                    residual = Math.Max(Error(a, Product(Product(Work(result[0]), Work(result[1])), Adjoint(Work(result[3])))),
                        Error(b, Product(Product(Work(result[0]), Work(result[2])), Adjoint(Work(result[3])))));
                    Orthogonal(result[0]); Orthogonal(result[3]); break;
                case "EVD": case "GEVD":
                    var v = Work(result[0]); var alpha = Values(result[1]);
                    if (real)
                        for (int j = 0; j < n; j++)
                        {
                            if (alpha[j].Imaginary <= 0) continue;
                            for (int i = 0; i < n; i++) { var z = new Complex(v[i, j].Real, v[i, j + 1].Real); v[i, j] = z; v[i, j + 1] = Complex.Conjugate(z); }
                            j++;
                        }
                    var av = Product(a, v); var bv = algorithm == "EVD" ? v : Product(b, v);
                    var beta = algorithm == "EVD" ? Enumerable.Repeat(Complex.One, n).ToArray() : Values(result[2]);
                    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) { av[i, j] *= beta[j]; bv[i, j] *= alpha[j]; }
                    residual = Error(av, bv);
                    if (hermitian) Orthogonal(result[0]);
                    break;
                default: throw new ArgumentException("Unsupported domain benchmark: " + algorithm);
            }
            if (!double.IsFinite(residual) || residual > 1e-4 || !double.IsFinite(orthogonality) || orthogonality > 1e-4)
                throw new InvalidOperationException($"Validation failed: residual={residual}, orthogonality={orthogonality}");
            if (Error(a, Work(aa)) != 0 || Error(b, Work(bb)) != 0) throw new InvalidOperationException("Input changed.");
            var timer = Stopwatch.StartNew(); int warmCalls = 0;
            do { GC.KeepAlive(run(aa, bb)); warmCalls++; } while (timer.Elapsed.TotalMilliseconds < 250 || warmCalls < 3);
            int count = Math.Clamp((int)Math.Ceiling(80 / (timer.Elapsed.TotalMilliseconds / warmCalls)), 1, 100000);
            var samples = new double[7]; var allocations = new long[7];
            for (int sample = 0; sample < samples.Length; sample++)
            {
                GC.Collect(); GC.WaitForPendingFinalizers(); GC.Collect();
                long allocated = GC.GetTotalAllocatedBytes(true);
                timer.Restart();
                for (int i = 0; i < count; i++) GC.KeepAlive(run(aa, bb));
                timer.Stop();
                samples[sample] = timer.Elapsed.TotalMilliseconds / count;
                allocations[sample] = (GC.GetTotalAllocatedBytes(true) - allocated) / count;
            }
            Console.WriteLine(JsonSerializer.Serialize(new { label, name, domain, rows, n, milliseconds = samples.Order().ElementAt(3), bytes = allocations.Order().ElementAt(3), samples, count, residual, orthogonality, runtime = Environment.Version.ToString(), error = (string?)null }));
        }
        catch (Exception e)
        {
            Console.WriteLine(JsonSerializer.Serialize(new { label, name, domain, rows, n, error = e.ToString() }));
            Environment.ExitCode = 1;
        }
    }

    private static Complex[,] Matrix(int m, int n, int seed, bool real, bool hermitian)
    {
        var random = new Random(seed + 13 * m + n); var a = new Complex[m, n];
        for (int i = 0; i < m; i++) for (int j = 0; j < n; j++) a[i, j] = new Complex((float)(2 * random.NextDouble() - 1), real ? 0 : (float)(2 * random.NextDouble() - 1));
        if (hermitian) for (int i = 0; i < n; i++) for (int j = 0; j <= i; j++) { if (i == j) a[i, i] = n + 1; else a[j, i] = Complex.Conjugate(a[i, j]); }
        return a;
    }
    private static Complex[,] Diagonal(Complex[] d) { var a = new Complex[d.Length, d.Length]; for (int i = 0; i < d.Length; i++) a[i, i] = d[i]; return a; }
    private static Complex[,] Adjoint(Complex[,] a) { var b = new Complex[a.GetLength(1), a.GetLength(0)]; for (int i = 0; i < a.GetLength(0); i++) for (int j = 0; j < a.GetLength(1); j++) b[j, i] = Complex.Conjugate(a[i, j]); return b; }
    private static Complex[,] Product(Complex[,] a, Complex[,] b) { var c = new Complex[a.GetLength(0), b.GetLength(1)]; for (int i = 0; i < c.GetLength(0); i++) for (int j = 0; j < c.GetLength(1); j++) for (int k = 0; k < a.GetLength(1); k++) c[i, j] += a[i, k] * b[k, j]; return c; }
    private static double Error(Complex[,] a, Complex[,] b) { double error = 0, norm = 0; for (int i = 0; i < a.GetLength(0); i++) for (int j = 0; j < a.GetLength(1); j++) { error += Math.Pow(Complex.Abs(a[i, j] - b[i, j]), 2); norm += Math.Pow(Complex.Abs(a[i, j]), 2); } return norm == 0 ? Math.Sqrt(error) : Math.Sqrt(error / norm); }
}
