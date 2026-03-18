namespace Fourier;

using System;
using System.Numerics;

/// <summary>
/// Beautiful but slower recursive FFT implementation, it shows the basic principle of the divide and conquer strategy
/// by splitting the data in even and odd values, doing FFTs on those and then recombining the results. This one is slower
/// because it allocates temporary arrays for the evens, odds and zetas.
/// </summary>
public static class RecursiveFFTA
{
    private static Complex[] SemiRootsOfUnity(int n, int direction)
        => [.. Enumerable.Range(0, n / 2).Select(i => Complex.FromPolarCoordinates(1, direction * i * Math.Tau / n))];

    public static void FastFourierTransform(Complex[] data)
    {
        if (data.Length < 2) return;

        // take even data points using linq:
        var evens = data.Where((_, i) => i % 2 == 0).ToArray();
        var odds = data.Where((_, i) => i % 2 == 1).ToArray();

        FastFourierTransform(evens);
        FastFourierTransform(odds);

        var zetas = SemiRootsOfUnity(data.Length, -1);

        var halfLength = data.Length / 2;
        for (var i = 0; i < halfLength; i++)
        {
            data[i] = evens[i] + zetas[i] * odds[i];
            // -zetas[i] is the same as +zetas[i+halfLength], so that allows cutting the zeta table in half
            data[halfLength + i] = evens[i] - zetas[i] * odds[i];
        }
    }
}
