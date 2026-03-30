namespace Baksteen.Numerics.Fourier;

using System;
using System.Numerics;

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

    private static void FastFourierTransform(SkipSpan<Complex> data)
    {
        if (data.Length < 2)
        {
            return;
        }
        else if (data.Length == 2)
        {
            // radix 2 FFT, 180 degrees
            (data[0], data[1]) = (data[0] + data[1], data[0] - data[1]);
        }
        else if (data.Length == 4)
        {
            // radix 4 FFT, 90 degrees
            // Pairwise sums and differences
            var s0 = data[0] + data[2];
            var s1 = data[1] + data[3];
            var d0 = data[0] - data[2];
            var d1 = data[1] - data[3];
            d1 = new Complex(d1.Imaginary, -d1.Real);                  // rotate 90 degrees clockwise
            data[0] = s0 + s1;                                         // X0 = s0+s1 = x0+x1+x2+x3
            data[2] = s0 - s1;                                         // X2 = s0-s1 = x0+x2-x1-x3 = x0-x1+x2-x3
            data[1] = d0 + d1;                                         // X1 = d0 - j*d1 (90 degrees clockwise)
            data[3] = d0 - d1;                                         // X3 = d0 + j*d1 (90 degrees counterclockwise)
        }
        else
        {
            var evens = data.SliceEvens();
            FastFourierTransform(evens);

            var odds = data.SliceOdds();
            FastFourierTransform(odds);

            var rotationstep = _rotations[BitOperations.Log2((uint)data.Length)];

            var w = Complex.One;
            for (var i = 0; i < (data.Length >> 1); i++)
            {
                Butterflies.Butterfly(ref evens[i], ref odds[i], w);
                w *= rotationstep;
            }

            // above, we can't directly store the combined results into the correct 'data' location because
            // 'evens' and 'odds' alias 'data' so it would result in clobbering results. So I temporarily store
            // the results in evens and odds itself. Now we have to reorder those so that 'data' becomes the
            // concatenation of evens and odds.
            ReorderEvenOdd2(data);
        }
    }

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

    public static void ReorderEvenOdd2<T>(SkipSpan<T> array)
    {
        if (array.Length < 2) return;

        // Helper to compute the original source index for target position i
        // (this defines the exact permutation that preserves relative order)
        static int Source(int i, int n)
        {
            var m = (n + 1) >> 1;   // ⌈n/2⌉ – size of the first half
            var j = i;
            do
            {
                j = j < m ? (j << 1) : ((j - m) << 1) + 1;
            }
            while (j < i);
            return j;
        }

        for (var i = 1; i < array.Length; i++)
        {
            var j = Source(i, array.Length);

            // Perform the swap (self-swap is harmless and never occurs here)
            if (j != i)
            {
                (array[i], array[j]) = (array[j], array[i]);
            }
        }
    }
}
