namespace Baksteen.Numerics.Fourier;

using System;
using System.Numerics;

public static class RecursiveFFTE
{
    public static void FastFourierTransform(Span<Complex> data, bool isInverse)
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
            var halfLength = data.Length >> 1;
            var evens = new Complex[halfLength];
            var odds = new Complex[halfLength];

            for (var i = 0; i < data.Length; i += 2)
            {
                evens[i >> 1] = data[i];
            }
            FastFourierTransform(evens, isInverse);

            for (var i = 0; i < data.Length; i += 2)
            {
                odds[i >> 1] = data[i + 1];
            }
            FastFourierTransform(odds, isInverse);

            var rotationstep = FFTUtils.GetRotation(BitOperations.Log2((uint)data.Length), isInverse);
            var w = Complex.One;

            for (var i = 0; i < halfLength; i++)
            {
                (evens[i], odds[i]) = Butterflies.ButterflyTuple(evens[i], odds[i], w);
                w *= rotationstep;
            }

            for (var i = 0; i < halfLength; i++)
            {
                (data[i], data[halfLength + i]) = (evens[i], odds[i]);
            }
        }
    }
}
