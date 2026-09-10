using System.Globalization;
using System.Runtime.CompilerServices;

namespace UMapx.Tests;

internal static class TestCulture
{
    [ModuleInitializer]
    internal static void Initialize()
    {
        // Stable numeric diagnostics and English exception messages on every developer machine.
        CultureInfo.DefaultThreadCurrentCulture = CultureInfo.InvariantCulture;
        CultureInfo.DefaultThreadCurrentUICulture = CultureInfo.InvariantCulture;
        CultureInfo.CurrentCulture = CultureInfo.InvariantCulture;
        CultureInfo.CurrentUICulture = CultureInfo.InvariantCulture;
    }
}
