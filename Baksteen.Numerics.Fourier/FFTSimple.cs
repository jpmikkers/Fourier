namespace Baksteen.Numerics.Fourier;

using System;
using System.Numerics;

public static class FFTSimple
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

        while (nrOfParts > 0)
        {
            var wr = FFTUtils.GetRotation(rotationLookupIndex, isInverse);

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

            butterfliesPerPart <<= 1;
            nrOfParts >>= 1;
            rotationLookupIndex++;
        }

        if (isInverse) FFTUtils.Scale(data);
    }
}
