namespace Baksteen.Numerics.Fourier;

using System;
using System.Numerics;

public static class FFTK
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

            if (butterfliesPerPart > nrOfParts)
            //if (data.Length >= 2048)
            {
                for (var p = 0; p < nrOfParts; p++)
                {
                    var evenindex = p << rotationLookupIndex;
                    var oddindex = evenindex + butterfliesPerPart;
                    var w = Complex.One;

                    for (var a = 0; a < butterfliesPerPart; a++)
                    {
                        Butterflies.Butterfly2(ref data[evenindex], ref data[oddindex], w);
                        w *= wr;
                        evenindex++;
                        oddindex++;
                    }
                }
            }
            else
            {
                var w = Complex.One;
                var evenindex = 0;
                var oddindex = evenindex + butterfliesPerPart;

                for (var a = 0; a < butterfliesPerPart; a++)
                {
                    var rli = 0;

                    for (var p = 0; p < nrOfParts; p++)
                    {
                        Butterflies.Butterfly(ref data[evenindex + rli], ref data[oddindex + rli], w);
                        rli += (1 << rotationLookupIndex);
                    }

                    w *= wr;
                    evenindex++;
                    oddindex++;
                }
            }

            butterfliesPerPart <<= 1;
            nrOfParts >>= 1;
            rotationLookupIndex++;
        }

        if (isInverse) FFTUtils.Scale(data);
    }
}
