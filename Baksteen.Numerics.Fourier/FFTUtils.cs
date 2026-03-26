namespace Baksteen.Numerics.Fourier;

using System;
using System.Numerics;
using System.Runtime.CompilerServices;

public static class FFTUtils
{
    private static readonly Complex[] _rotations;

    static FFTUtils()
    {
        // contains Wn for powers of two n. Wn = e^(-i*2π/n) = e^(-iτ/n)
        _rotations = [.. Enumerable.Range(0, 32)
            .Select(lg2 => Complex.FromPolarCoordinates(1, -Math.Tau / Math.Pow(2.0, lg2)))];
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Complex GetRotation(int stage, bool isInverse)
    {
        return isInverse ? Complex.Conjugate(_rotations[stage]) : _rotations[stage];
    }

    public static void Scale(Span<Complex> data)
    {
        var scaleFactor = Math.ScaleB(1.0, -BitOperations.Log2((uint)data.Length));
        foreach (ref var c in data) { c *= scaleFactor; }
    }
}
