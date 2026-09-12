using System.Diagnostics;
using System.Linq.Expressions;
using System.Reflection;
using System.Text.Json;

if (args.Length < 4)
{
    Console.Error.WriteLine("Usage: UMapx.DecompositionBenchmarks <assembly.dll> <label> <algorithm> <columns> [rows]");
    Environment.ExitCode = 2;
    return;
}

// Load only one version per process. Reflection and expression compilation are outside timing.
string dll = Path.GetFullPath(args[0]), label = args[1], name = args[2];
int n = int.Parse(args[3]), rows = args.Length > 4 ? int.Parse(args[4]) : n;
var assembly = Assembly.LoadFrom(dll);
string algorithm = name.Split('-')[0];
try
{
    var run = Create(assembly, algorithm, n, name.EndsWith("Full"));
    var a = Matrix(rows, n, 173); var b = Matrix(rows, n, 917);
    if (name.EndsWith("SPD") || algorithm is "Cholesky" or "LDL" or "UDL" or "Householder" or "Lanczos" or "Diagonal" or "Power")
        for (int i = 0; i < n; i++) for (int j = 0; j <= i; j++) a[i, j] = a[j, i] = i == j ? n + 1 : a[i, j];
    if (algorithm is "LU" or "LDU") for (int i = 0; i < n; i++) a[i, i] += n;
    if (algorithm is "GEVD" or "QZ") for (int i = 0; i < n; i++) b[i, i] += n;
    if (algorithm == "NMF") for (int i = 0; i < rows; i++) for (int j = 0; j < n; j++) a[i, j] = Math.Abs(a[i, j]);
    var original = (float[,])a.Clone();
    object[] result = run(a, b);
    foreach (object factor in result) if (factor is Array array)
        foreach (object item in array) if (item is float f && !float.IsFinite(f)) throw new Exception("Nonfinite factor");
    if (!a.Cast<float>().SequenceEqual(original.Cast<float>())) throw new Exception("Input mutated");
    double? residual = Residual(algorithm, a, result);
    // Disable tiering via the launch environment. Warm for at least 250 ms and three calls.
    var warm = Stopwatch.StartNew(); int warmCalls = 0;
    do { GC.KeepAlive(run(a, b)); warmCalls++; } while (warm.Elapsed.TotalMilliseconds < 250 || warmCalls < 3);
    int count = Math.Clamp((int)Math.Ceiling(80 / (warm.Elapsed.TotalMilliseconds / warmCalls)), 1, 100000);
    var times = new double[7]; var allocations = new double[7];
    for (int sample = 0; sample < times.Length; sample++)
    {
        GC.Collect(); GC.WaitForPendingFinalizers(); GC.Collect();
        long bytes = GC.GetTotalAllocatedBytes(true);
        var timer = Stopwatch.StartNew();
        for (int i = 0; i < count; i++) GC.KeepAlive(run(a, b));
        timer.Stop();
        allocations[sample] = (GC.GetTotalAllocatedBytes(true) - bytes) / (double)count;
        times[sample] = timer.Elapsed.TotalMilliseconds / count;
    }
    var sorted = times.Order().ToArray();
    Console.WriteLine(JsonSerializer.Serialize(new { label, name, rows, n, milliseconds = sorted[3], min = sorted[0], max = sorted[6], bytes = allocations.Order().ElementAt(3), count, samples = times, residual, runtime = Environment.Version.ToString(), error = (string?)null }));
}
catch (Exception e)
{
    Console.WriteLine(JsonSerializer.Serialize(new { label, name, rows, n, error = e.ToString() }));
    Environment.ExitCode = 1;
}

static float[,] Matrix(int m, int n, int seed)
{
    var random = new Random(seed + 13 * m + n); var a = new float[m, n];
    for (int i = 0; i < m; i++) for (int j = 0; j < n; j++) a[i, j] = (float)(2 * random.NextDouble() - 1);
    return a;
}

static Func<float[,], float[,], object[]> Create(Assembly assembly, string name, int n, bool full)
{
    var type = assembly.GetType("UMapx.Decomposition." + name, true)!;
    bool modern = type.IsAbstract && type.IsSealed;
    MethodBase method = modern ? type.GetMethods().Single(m => m.Name == "Decompose" && m.GetParameters()[0].ParameterType == typeof(float[,]))
        : type.GetConstructors().Single(c => c.GetParameters()[0].ParameterType == typeof(float[,]));
    var a = Expression.Parameter(typeof(float[,]), "a"); var b = Expression.Parameter(typeof(float[,]), "b");
    var parameters = method.GetParameters();
    var inputs = parameters.Select((p, i) => p.ParameterType == typeof(float[,]) ? (Expression)(i == 0 ? a : b)
        : p.ParameterType == typeof(int) ? Expression.Constant(name == "NMF" && i == 1 ? Math.Min(n, 8) : 100)
        : p.ParameterType == typeof(float) ? Expression.Constant(1e-16f)
        : Expression.Constant(full)).ToArray();
    Expression call = modern ? Expression.Call((MethodInfo)method, inputs) : Expression.New((ConstructorInfo)method, inputs);
    var output = Expression.Variable(call.Type, "output");
    string[] members = name switch
    {
        "QR" or "GramSchmidt" => modern || name == "QR" ? ["Q", "R"] : ["Q"],
        "LQ" => ["L", "Q"], "QL" => ["Q", "L"], "RQ" => ["R", "Q"],
        "LU" => modern ? ["L", "U", "P"] : ["L", "U"],
        "LDU" => modern ? ["L", "D", "U", "P"] : ["L", "D", "U"],
        "LDL" => ["L", "D"], "UDL" => ["U", "D"], "Cholesky" => ["L"],
        "SVD" => ["U", "S", "V"], "Polar" => ["U", "P"],
        "EVD" => ["V", "D"], "GEVD" => ["V", "Alpha", "Beta"],
        "QZ" => ["Q", "S", "T", "Z"], "GSVD" => ["U1", "S1", "U2", "S2", "X"],
        "Schur" or "Lanczos" => ["Q", "T"], "Hessenberg" => ["P", "H"],
        "Arnoldi" => ["Q", "H"], "Householder" => ["H", "T"],
        "Bidiagonal" => ["U", "B", "V"], "Diagonal" => ["B", "D"],
        "NMF" => ["W", "H"], "Power" => ["V"], _ => throw new Exception(name)
    };
    var values = members.Select((member, i) => Expression.Convert(modern
        ? name == "Cholesky" ? (Expression)output : Expression.Field(output, "Item" + (i + 1))
        : Expression.Property(output, member), typeof(object)));
    return Expression.Lambda<Func<float[,], float[,], object[]>>(Expression.Block([output], Expression.Assign(output, call), Expression.NewArrayInit(typeof(object), values)), a, b).Compile();
}

static double? Residual(string name, float[,] a, object[] d)
{
    // Independent double-precision reconstruction, entirely outside timing.
    double[,] W(int index) => Convert((float[,])d[index]);
    double[,] D(int index)
    {
        var v = (float[])d[index]; var r = new double[v.Length, v.Length];
        for (int i = 0; i < v.Length; i++) r[i, i] = v[i]; return r;
    }
    double[,]? reconstructed = name switch
    {
        "QR" or "LQ" or "QL" or "RQ" or "Polar" => Mul(W(0), W(1)),
        "SVD" => Mul(Mul(W(0), D(1)), Transpose(W(2))),
        "Bidiagonal" => Mul(Mul(W(0), W(1)), Transpose(W(2))),
        "Cholesky" => Mul(W(0), Transpose(W(0))),
        "LDL" or "UDL" => Mul(Mul(W(0), D(1)), Transpose(W(0))),
        "Hessenberg" or "Schur" or "Arnoldi" or "Lanczos" or "Householder" => Mul(Mul(W(0), W(1)), Transpose(W(0))),
        "GSVD" => Mul(Mul(W(0), D(1)), W(4)),
        "QZ" => Mul(Mul(W(0), W(1)), Transpose(W(3))),
        _ => null
    };
    if (reconstructed == null) return null;
    double error = 0, norm = 0;
    for (int i = 0; i < a.GetLength(0); i++) for (int j = 0; j < a.GetLength(1); j++)
    { double delta = reconstructed[i, j] - a[i, j]; error += delta * delta; norm += (double)a[i, j] * a[i, j]; }
    double residual = Math.Sqrt(error / norm);
    return double.IsFinite(residual) ? residual : double.MaxValue;
}
static double[,] Convert(float[,] a)
{
    var r = new double[a.GetLength(0), a.GetLength(1)];
    for (int i = 0; i < r.GetLength(0); i++) for (int j = 0; j < r.GetLength(1); j++) r[i, j] = a[i, j]; return r;
}
static double[,] Transpose(double[,] a)
{
    var r = new double[a.GetLength(1), a.GetLength(0)];
    for (int i = 0; i < r.GetLength(0); i++) for (int j = 0; j < r.GetLength(1); j++) r[i, j] = a[j, i]; return r;
}
static double[,] Mul(double[,] a, double[,] b)
{
    if (a.GetLength(1) != b.GetLength(0)) throw new Exception("Factor dimensions mismatch");
    var r = new double[a.GetLength(0), b.GetLength(1)];
    for (int i = 0; i < r.GetLength(0); i++) for (int k = 0; k < a.GetLength(1); k++)
        for (int j = 0; j < r.GetLength(1); j++) r[i, j] += a[i, k] * b[k, j]; return r;
}
