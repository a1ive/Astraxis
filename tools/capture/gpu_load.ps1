<#
.SYNOPSIS
Launches astraxis.exe on a scene event and reports the average 3D-engine
utilization of that process (Windows GPU performance counters, no admin needed).

.DESCRIPTION
The swapchain presents with vsync, so the utilization is about the fraction of
each frame the GPU is busy: compare it before and after a rendering change on
the same machine, window size and driver. An empty settings file reproduces the
default settings. Sampling ends before the 60 s idle tour starts if
Settle + Samples stays below ~55 s.

.EXAMPLE
pwsh tools/capture/gpu_load.ps1 -Scene sgr_a -Event 1 -Driver vulkan
#>
param(
    [string]$Exe = (Join-Path $PSScriptRoot '..\..\build\win-msvc\Release\astraxis.exe'),
    [string]$Scene = 'sgr_a',
    [int]$Event = 1,
    [int]$Settle = 20,       # seconds before sampling (camera move, warm-up)
    [int]$Samples = 20,      # one-second counter samples
    [string]$Driver = ''     # SDL_GPU_DRIVER: direct3d12 or vulkan (empty: default)
)

$ErrorActionPreference = 'Stop'
$cfg = Join-Path ([System.IO.Path]::GetTempPath()) 'astraxis_gpu_load_empty.ini'
Set-Content -Path $cfg -Value '' -NoNewline
$exePath = (Resolve-Path $Exe).Path
$psi = New-Object System.Diagnostics.ProcessStartInfo
$psi.FileName = $exePath
$psi.Arguments = "--scene $Scene --event $Event --mode window --config `"$cfg`""
$psi.UseShellExecute = $false
$psi.WorkingDirectory = Split-Path $exePath
if ($Driver) { $psi.Environment['SDL_GPU_DRIVER'] = $Driver }
$p = [System.Diagnostics.Process]::Start($psi)
try {
    Start-Sleep -Seconds $Settle
    $values = @()
    for ($i = 0; $i -lt $Samples; $i++) {
        $s = (Get-Counter "\GPU Engine(pid_$($p.Id)_*engtype_3D)\Utilization Percentage" -ErrorAction SilentlyContinue).CounterSamples
        $values += ($s | Measure-Object CookedValue -Sum).Sum
    }
    $stats = $values | Measure-Object -Average -Minimum -Maximum
    '{0} event {1} {2}: 3D utilization avg {3:N1}% min {4:N1} max {5:N1} (n={6})' -f $Scene, $Event, $Driver, $stats.Average, $stats.Minimum, $stats.Maximum, $stats.Count
} finally {
    if (-not $p.HasExited) { Stop-Process -Id $p.Id -Force }
}
