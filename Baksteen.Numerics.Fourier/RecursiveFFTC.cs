namespace Baksteen.Numerics.Fourier;

using System;
using System.Numerics;

/// <summary>
/// Same as <see cref="RecursiveFFTB"/> , but without precomputed rotation table.
/// </summary>
public static class RecursiveFFTC
{
    public static void FastFourierTransform(Complex[] data)
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
            // take even data points using linq:
            var evens = data.Where((_, i) => i % 2 == 0).ToArray();
            var odds = data.Where((_, i) => i % 2 == 1).ToArray();

            FastFourierTransform(evens);
            FastFourierTransform(odds);

            var rotationstep = Complex.FromPolarCoordinates(1, -Math.Tau / data.Length);
            var halfLength = data.Length / 2;

            data[0] = evens[0] + odds[0];
            data[halfLength] = evens[0] - odds[0];

            var w = rotationstep;
            for (var i = 1; i < halfLength; i++, w *= rotationstep)
            {
                data[i] = evens[i] + w * odds[i];
                data[halfLength + i] = evens[i] - w * odds[i];
            }
        }
    }
}
