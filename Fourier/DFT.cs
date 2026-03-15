namespace Fourier;

using System;
using System.Numerics;

public static class DFT
{
    /// <summary>
    /// returns n-th roots of unity, the points of a regular n-gon on the unit circle in the complex plane.
    /// (see https://algebrica.org/roots-of-unity/ )
    /// </summary>
    /// <param name="n">how many points</param>
    /// <param name="direction">1 for counter clockwise, -1 for clockwise</param>
    /// <returns></returns>
    private static Complex[] RootsOfUnity(int n, int direction)
        => [.. Enumerable.Range(0, n).Select(i => Complex.FromPolarCoordinates(1, direction * i * Math.Tau / n))];

    /// <summary>
    /// Calculates the (forward or reverse) complex discrete fourier transform.
    /// </summary>
    /// <param name="data">Complex signal or spectrum</param>
    /// <param name="forward">true to transform the data from the time domain to frequency domain, false for the inverse</param>
    /// <returns></returns>
    public static Complex[] DiscreteFourierTransform(Complex[] data, bool forward)
    {
        var n = data.Length;

        // precalculate zetas, clockwise roots of unity for forward DFT, counterclockwise for reverse DFT
        var zetas = RootsOfUnity(n, forward ? -1 : 1);
        var scaleFactor = forward ? 1.0 : 1.0 / n;

        return
            [.. Enumerable.Range(0, n)
                .Select(i =>
                    data.Select((sample, j) => sample * zetas[(i * j) % n])
                    .Aggregate(Complex.Zero, Complex.Add) * scaleFactor
                )
            ];
    }
}
