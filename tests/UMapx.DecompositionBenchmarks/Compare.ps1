param(
    [Parameter(Mandatory)][string]$PreviousAssembly,
    [Parameter(Mandatory)][string]$CurrentAssembly,
    [Parameter(Mandatory)][string]$OutputFile,
    [int[]]$Sizes = @(32, 128),
    [string[]]$Names = @('QR','LQ','QL','RQ','LU','LDU','Cholesky','LDL','UDL','SVD','EVD','EVD-SPD','Hessenberg','Householder','Bidiagonal','GramSchmidt','Arnoldi','Lanczos','Polar','GEVD','QZ','GSVD','Schur','Diagonal','Power','NMF'),
    [int]$Rows = 0
)
$ErrorActionPreference = 'Stop'
$runner = Join-Path $PSScriptRoot 'bin/Release/net8.0/UMapx.DecompositionBenchmarks.dll'
if (!(Test-Path -LiteralPath $runner)) { throw 'Build UMapx.DecompositionBenchmarks.csproj in Release first.' }
$previous = (Resolve-Path -LiteralPath $PreviousAssembly).Path
$current = (Resolve-Path -LiteralPath $CurrentAssembly).Path
$output = [IO.Path]::GetFullPath($OutputFile)
if (Test-Path -LiteralPath $output) { throw 'Use a new output file for each independent comparison.' }
$directory = [IO.Path]::GetDirectoryName($output)
[IO.Directory]::CreateDirectory($directory) | Out-Null
foreach ($n in $Sizes) {
    foreach ($name in $Names) {
        # Reverse order for alternating sizes to reduce systematic ordering bias.
        $versions = if ([Array]::IndexOf($Sizes, $n) % 2) { @('current','previous') } else { @('previous','current') }
        foreach ($version in $versions) {
            $m = if ($Rows -gt 0) { $Rows } else { $n }
            $library = if ($version -eq 'previous') { $previous } else { $current }
            $info = [Diagnostics.ProcessStartInfo]::new('dotnet')
            $info.UseShellExecute = $false
            $info.CreateNoWindow = $true
            $info.RedirectStandardOutput = $true
            $info.RedirectStandardError = $true
            $info.Environment['DOTNET_TieredCompilation'] = '0'
            $info.Environment['DOTNET_ReadyToRun'] = '0'
            foreach ($argument in @($runner,$library,$version,$name,"$n","$m")) { $info.ArgumentList.Add($argument) }
            $process = [Diagnostics.Process]::Start($info)
            try {
                $stdout = $process.StandardOutput.ReadToEndAsync()
                $stderr = $process.StandardError.ReadToEndAsync()
                if (!$process.WaitForExit(60000)) {
                    $process.Kill()
                    $line = @{label=$version;name=$name;rows=$m;n=$n;error='Process exceeded 60 seconds'} | ConvertTo-Json -Compress
                } else {
                    $line = $stdout.Result.Trim()
                    if (!$line) { $line = @{label=$version;name=$name;rows=$m;n=$n;error=$stderr.Result} | ConvertTo-Json -Compress }
                }
                Add-Content -LiteralPath $output -Value $line -Encoding utf8
                $result = $line | ConvertFrom-Json
                if ($result.error) { Write-Output "$version $name ${m}x${n}: $($result.error.Split([Environment]::NewLine)[0])" }
                else { Write-Output ('{0} {1} {2}x{3}: {4:N3} ms; {5:N0} B/op' -f $version,$name,$m,$n,$result.milliseconds,$result.bytes) }
            } finally { $process.Dispose() }
        }
    }
}
