namespace Fourier;

using System;
using System.Numerics;
using System.Runtime.CompilerServices;

/// <summary>
/// Allocation free recursive FFT. It's easy to read because it uses SkipSpan to allow slicing data 
/// into even and odd elements.
/// </summary>
public static class RecursiveFFTD
{
    private static readonly Complex[] _rotations;

    static RecursiveFFTD()
    {
        _rotations = [.. Enumerable.Range(0, 32)
            .Select(lg2 => Complex.FromPolarCoordinates(1, -Math.Tau / Math.Pow(2.0, lg2)))];
    }

    public static void FastFourierTransform(Span<Complex> data)
    {
        FastFourierTransform(new SkipSpan<Complex>(data));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void Butterfly(ref Complex even, ref Complex odd, Complex w)
    {
        var odd_w = odd * w;
        (even, odd) = (even + odd_w, even - odd_w);
    }

    private static void FastFourierTransform(SkipSpan<Complex> data)
    {
        if (data.Length >= 2)
        {
            var evens = data.SliceEvens();
            var odds = data.SliceOdds();

            FastFourierTransform(evens);
            FastFourierTransform(odds);

            // for i=0, w0 is always 1+0i, so we can avoid a multiplication here
            (evens[0], odds[0]) = (evens[0] + odds[0], evens[0] - odds[0]);

            if (data.Length >= 4)
            {
                var rotationstep = _rotations[BitOperations.Log2((uint)data.Length)];
                var w = rotationstep;

                for (var i = 1; i < (data.Length >> 1); i++)
                {
                    Butterfly(ref evens[i], ref odds[i], w);
                    w *= rotationstep;
                }

                // above, we can't directly store the combined results into the correct 'data' location because
                // 'evens' and 'odds' alias 'data' so it would result in clobbering results. So I temporarily store
                // the results in evens and odds itself. Now we have to reorder those so that 'data' becomes the
                // concatenation of evens and odds.
                ReorderEvenOdd(data);
            }
        }
    }

#if SIMPLER_BUT_SLOWER
    private static void FastFourierTransformO(SkipSpan<Complex> data)
    {
        if (data.Length >= 2)
        {
            var evens = data.SliceEvens();
            var odds = data.SliceOdds();

            FastFourierTransformO(evens);
            FastFourierTransformO(odds);

            var rotationstep = _rotations[BitOperations.Log2((uint)data.Length)];

            var w = Complex.One;
            for (var i = 0; i < (data.Length >> 1); i++)
            {
                (evens[i], odds[i]) = (evens[i] + w * odds[i], evens[i] - w * odds[i]);
                w *= rotationstep;
            }

            ReorderEvenFirst(data);
        }
    }
#endif

    // Reorder in-place so all even-indexed elements come first, then all odd-indexed ones.
    // n must be a power of two
    // example:
    // e0 o0 e1 o1 e2 o2 e3 o3 becomes: e0 e1 e2 e3 o0 o1 o2 o3
    private static void ReorderEvenOdd<T>(SkipSpan<T> data)
    {
        int n = data.Length;
        if (n < 4) return;

        static int P(int i, int n) => ((i & 1) == 0) ? (i >> 1) : ((n >> 1) + ((i - 1) >> 1));

        for (int s = 0; s < n; s++)
        {
            // Find minimum in cycle; only process if s is the minimum
            int min = s, i = s;
            while ((i = P(i, n)) != s)
                if (i < min) { min = -1; break; }

            if (min == s)
            {
                T tmp = data[s];
                for (i = s; (i = P(i, n)) != s;)
                    (tmp, data[i]) = (data[i], tmp);
                data[s] = tmp;
            }
        }
    }
}
