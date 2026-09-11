using System.Diagnostics;
using Xunit;

namespace UMapx.Tests;

internal static class AuditProcess
{
    internal static async Task<string> RunAsync(string operation,string argument)
    {
        var start=new ProcessStartInfo("dotnet") { UseShellExecute=false, CreateNoWindow=true, RedirectStandardOutput=true, RedirectStandardError=true };
        foreach(string arg in new[]{"exec","--runtimeconfig",Path.Combine(AppContext.BaseDirectory,"UMapx.Tests.runtimeconfig.json"),"--depsfile",Path.Combine(AppContext.BaseDirectory,"UMapx.Tests.deps.json"),typeof(AuditProbe.Program).Assembly.Location,operation,argument})
            start.ArgumentList.Add(arg);
        using var process=Process.Start(start)!;
        var stdout=process.StandardOutput.ReadToEndAsync();var stderr=process.StandardError.ReadToEndAsync();
        using var deadline=new CancellationTokenSource(TimeSpan.FromSeconds(5));
        try { await process.WaitForExitAsync(deadline.Token); }
        catch(OperationCanceledException)
        {
            process.Kill(entireProcessTree:true);await process.WaitForExitAsync();
            Assert.Fail($"{operation}({argument}) did not terminate within five seconds.");
        }
        Assert.True(process.ExitCode==0,await stderr);
        return (await stdout).Trim();
    }
}
