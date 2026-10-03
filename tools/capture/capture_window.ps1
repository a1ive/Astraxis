<#
.SYNOPSIS
Launches astraxis.exe, optionally sends key presses, and saves the window's
client area to a PNG with PrintWindow (works for D3D12 and Vulkan swapchains,
also when the window cannot be brought to the foreground).

.DESCRIPTION
The capture goes through a plain GDI memory DC and Image.FromHbitmap. Do not
use Graphics.GetHdc() on a System.Drawing.Bitmap (or Graphics.CopyFromScreen,
which uses it internally): GDI+ pre-fills the HDC with the sentinel colour
RGB(13, 11, 12) and on ReleaseHdc() turns every pixel that still has that
value into transparent black, whatever the bitmap's pixel format. Dark
gradients such as a star's glow then show isolated "holes" along the contour
where the image happens to be exactly (13, 11, 12); they look white in a
viewer with a light background.

.EXAMPLE
pwsh tools/capture/capture_window.ps1 -Exe build/win-msvc/Release/astraxis.exe -Scene kepler223 -Out shot.png
pwsh tools/capture/capture_window.ps1 -Exe build/win-msvc/Release/astraxis.exe -Scene jupiter -Event 2 -Keys OLH -Out shot.png
#>
param(
    [Parameter(Mandatory = $true)][string]$Exe,
    [Parameter(Mandatory = $true)][string]$Out,
    [string]$Scene = "",
    [int]$Event = -1,
    [string]$Keys = "",          # letters pressed after startup, e.g. "OLH" hides orbits, labels and UI
    [int]$Width = 1280,          # client area size
    [int]$Height = 720,
    [double]$WaitSeconds = 5.0   # after the key presses, before the capture
)

$ErrorActionPreference = "Stop"
Add-Type -AssemblyName System.Drawing
Add-Type @"
using System;
using System.Runtime.InteropServices;
public static class AstraxisCapture {
    [StructLayout(LayoutKind.Sequential)] public struct RECT { public int Left, Top, Right, Bottom; }
    [DllImport("user32.dll")] static extern bool PrintWindow(IntPtr hwnd, IntPtr hdc, uint flags);
    [DllImport("user32.dll")] public static extern bool GetClientRect(IntPtr hwnd, out RECT rect);
    [DllImport("user32.dll")] public static extern bool AdjustWindowRect(ref RECT rect, uint style, bool menu);
    [DllImport("user32.dll")] public static extern bool SetWindowPos(IntPtr hwnd, IntPtr after, int x, int y, int cx, int cy, uint flags);
    [DllImport("user32.dll")] public static extern bool PostMessage(IntPtr hwnd, uint msg, IntPtr wparam, IntPtr lparam);
    [DllImport("user32.dll")] static extern IntPtr GetDC(IntPtr hwnd);
    [DllImport("user32.dll")] static extern int ReleaseDC(IntPtr hwnd, IntPtr dc);
    [DllImport("gdi32.dll")] static extern IntPtr CreateCompatibleDC(IntPtr dc);
    [DllImport("gdi32.dll")] static extern IntPtr CreateCompatibleBitmap(IntPtr dc, int width, int height);
    [DllImport("gdi32.dll")] static extern IntPtr SelectObject(IntPtr dc, IntPtr obj);
    [DllImport("gdi32.dll")] static extern bool DeleteDC(IntPtr dc);
    [DllImport("gdi32.dll")] public static extern bool DeleteObject(IntPtr obj);

    // Client area into a new HBITMAP (the caller deletes it).
    public static IntPtr Grab(IntPtr hwnd, int width, int height) {
        IntPtr screen = GetDC(IntPtr.Zero);
        IntPtr mem = CreateCompatibleDC(screen);
        IntPtr bitmap = CreateCompatibleBitmap(screen, width, height);
        IntPtr old = SelectObject(mem, bitmap);
        PrintWindow(hwnd, mem, 3); // PW_CLIENTONLY | PW_RENDERFULLCONTENT
        SelectObject(mem, old);
        DeleteDC(mem);
        ReleaseDC(IntPtr.Zero, screen);
        return bitmap;
    }
}
"@

$exePath = (Resolve-Path $Exe).Path
$arguments = @()
if ($Scene) { $arguments += @("--scene", $Scene) }
if ($Event -ge 0) { $arguments += @("--event", "$Event") }
$startArgs = @{ FilePath = $exePath; WorkingDirectory = (Split-Path $exePath); PassThru = $true }
if ($arguments.Count -gt 0) { $startArgs.ArgumentList = $arguments }
$process = Start-Process @startArgs
try {
    for ($i = 0; $i -lt 100 -and $process.MainWindowHandle -eq [IntPtr]::Zero; $i++) {
        Start-Sleep -Milliseconds 100
        $process.Refresh()
    }
    $hwnd = $process.MainWindowHandle
    if ($hwnd -eq [IntPtr]::Zero) { throw "astraxis window did not appear" }

    $rect = New-Object AstraxisCapture+RECT
    $rect.Right = $Width
    $rect.Bottom = $Height
    [void][AstraxisCapture]::AdjustWindowRect([ref]$rect, 0x00CF0000, $false) # WS_OVERLAPPEDWINDOW
    # SWP_NOZORDER | SWP_NOACTIVATE
    [void][AstraxisCapture]::SetWindowPos($hwnd, [IntPtr]::Zero, 50, 50, $rect.Right - $rect.Left, $rect.Bottom - $rect.Top, 0x0014)
    Start-Sleep -Seconds 2

    foreach ($key in $Keys.ToUpper().ToCharArray()) {
        $vk = [IntPtr][int]$key
        [void][AstraxisCapture]::PostMessage($hwnd, 0x0100, $vk, [IntPtr]1)          # WM_KEYDOWN
        Start-Sleep -Milliseconds 50
        [void][AstraxisCapture]::PostMessage($hwnd, 0x0101, $vk, [IntPtr]0xC0000001) # WM_KEYUP
        Start-Sleep -Milliseconds 150
    }
    Start-Sleep -Milliseconds ([int]($WaitSeconds * 1000))

    $client = New-Object AstraxisCapture+RECT
    [void][AstraxisCapture]::GetClientRect($hwnd, [ref]$client)
    $hbitmap = [AstraxisCapture]::Grab($hwnd, $client.Right, $client.Bottom)
    try {
        $image = [System.Drawing.Image]::FromHbitmap($hbitmap)
        $image.Save([System.IO.Path]::GetFullPath($Out), [System.Drawing.Imaging.ImageFormat]::Png)
        $image.Dispose()
    } finally {
        [void][AstraxisCapture]::DeleteObject($hbitmap)
    }
    "$($client.Right)x$($client.Bottom) -> $Out"
} finally {
    Stop-Process -Id $process.Id -ErrorAction SilentlyContinue
}
