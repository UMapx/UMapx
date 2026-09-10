param(
    [string]$ResultsDirectory = 'artifacts/math-audit/run',
    [switch]$NoRestore
)

$ErrorActionPreference = 'Stop'
$workspace = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
Push-Location $workspace
try {
    $arguments = @('test', 'sources/UMapx.sln', '-c', 'Release', '--nologo',
        '-p:GeneratePackageOnBuild=false', '-p:DebugType=portable', '-p:DebugSymbols=true',
        '--collect', 'XPlat Code Coverage', '--settings', 'tests/UMapx.Tests/coverage.runsettings',
        '--logger', 'trx;LogFileName=full-audit.trx', '--results-directory', $ResultsDirectory)
    if ($NoRestore) { $arguments += '--no-restore' }
    & dotnet @arguments
    $auditExit = $LASTEXITCODE
    $coverage = Get-ChildItem -LiteralPath $ResultsDirectory -Filter 'coverage.cobertura.xml' -Recurse |
        Sort-Object LastWriteTimeUtc -Descending | Select-Object -First 1
    $trx = Join-Path $ResultsDirectory 'full-audit.trx'
    if ((Test-Path -LiteralPath $trx) -and $null -ne $coverage) {
        & python -X utf8 tools/summarize_audit.py --trx $trx --coverage $coverage.FullName --output (Join-Path $ResultsDirectory 'summary')
        if ($LASTEXITCODE -ne 0) { throw 'Audit evidence export failed.' }
    }
    else { throw 'The test run did not produce TRX and coverage evidence.' }
    # Keep the failing exit code: known defects must not make CI appear green.
    exit $auditExit
}
finally { Pop-Location }
