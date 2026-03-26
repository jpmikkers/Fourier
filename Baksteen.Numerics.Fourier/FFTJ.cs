namespace Baksteen.Numerics.Fourier;

using System;
using System.Numerics;

public static class FFTJ
{
    public static void FastFourierTransform(Span<Complex> data, bool isInverse)
    {
        if (!BitOperations.IsPow2(data.Length))
        {
            throw new ArgumentException("fft not a power of two", nameof(data));
        }

        Reorder.Shuffle(data);

        var butterfliesPerPart = 1;             // a single butterfly does 2 angles, +w and -w (=w+pi radians)
        var nrOfParts = data.Length >> 1;       // so the first layer is len/2 parts of single butterflies
        var rotationLookupIndex = 1;
        const int maxstackbutterflies = 16;
        Span<Complex> tw = stackalloc Complex[maxstackbutterflies];

        if (nrOfParts > 0)
        {
            // no mul needed in the first combination layer, every w is 1+0i and -1+0i
            for (var p = 0; p < nrOfParts; p++)
            {
                Butterflies.Butterfly(ref data[p << 1], ref data[(p << 1) + 1]);
            }
            butterfliesPerPart <<= 1;
            nrOfParts >>= 1;
            rotationLookupIndex++;
        }

        while (nrOfParts > 0)
        {
            var wr = FFTUtils.GetRotation(rotationLookupIndex, isInverse);

            if (butterfliesPerPart <= maxstackbutterflies)
            {
                tw[0] = Complex.One;
                var w = wr;

                for (var a = 1; a <= (butterfliesPerPart >> 1); a++)
                {
                    tw[a] = w;
                    tw[butterfliesPerPart - a] = new Complex(-w.Real, w.Imaginary);
                    w *= wr;
                }

                for (var p = 0; p < nrOfParts; p++)
                {
                    var evenindex = p << rotationLookupIndex;
                    var oddindex = evenindex + butterfliesPerPart;

                    for (var a = 0; a < butterfliesPerPart; a++)
                    {
                        Butterflies.Butterfly(ref data[evenindex], ref data[oddindex], tw[a]);
                        evenindex++;
                        oddindex++;
                    }
                }
            }
            else
            {
                for (var p = 0; p < nrOfParts; p++)
                {
                    var w = Complex.One;
                    var evenindex = p << rotationLookupIndex;
                    var oddindex = evenindex + butterfliesPerPart;

                    for (var a = 0; a < butterfliesPerPart; a++)
                    {
                        Butterflies.Butterfly(ref data[evenindex], ref data[oddindex], w);
                        w *= wr;
                        evenindex++;
                        oddindex++;
                    }
                }
            }

            butterfliesPerPart <<= 1;
            nrOfParts >>= 1;
            rotationLookupIndex++;
        }

        if (isInverse) FFTUtils.Scale(data);
    }

}
