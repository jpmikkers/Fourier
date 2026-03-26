namespace Baksteen.Numerics.Fourier;

using System;
using System.Numerics;

public static class FFTL
{
    private static readonly Complex[] _rotations;

    static FFTL()
    {
        // contains 2π/2^n for n=0-31, so π/1, π/2, π/4, π/8 etc..
        _rotations = [.. Enumerable.Range(0, 32)
            .Select(lg2 => Complex.FromPolarCoordinates(1, -Math.Tau / Math.Pow(2.0, lg2)))];
    }

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

        if (nrOfParts > 0)
        {
            if (isInverse)
            {
                var ti = 0;
                for (var p = 0; p < nrOfParts; p++)
                {
                    Butterflies.Butterfly(ref data[ti], ref data[ti + 2]);
                    Butterflies.Butterfly90CCW(ref data[ti + 1], ref data[ti + 3]);
                    ti += (1 << rotationLookupIndex);
                }
            }
            else
            {
                var ti = 0;
                for (var p = 0; p < nrOfParts; p++)
                {
                    Butterflies.Butterfly(ref data[ti], ref data[ti + 2]);
                    Butterflies.Butterfly90CW(ref data[ti + 1], ref data[ti + 3]);
                    ti += (1 << rotationLookupIndex);
                }
            }

            butterfliesPerPart <<= 1;
            nrOfParts >>= 1;
            rotationLookupIndex++;
        }

        while (nrOfParts > 0)
        {
            var wr = isInverse ? Complex.Conjugate(_rotations[rotationLookupIndex]) : _rotations[rotationLookupIndex];

            if (butterfliesPerPart > nrOfParts)
            //if (data.Length >= 2048)
            {
                for (var p = 0; p < nrOfParts; p++)
                {
                    var evenindex = p << rotationLookupIndex;
                    var oddindex = evenindex + butterfliesPerPart;

                    Butterflies.Butterfly(ref data[evenindex], ref data[oddindex]);
                    var w = wr;
                    evenindex++;
                    oddindex++;

                    for (var a = 1; a < butterfliesPerPart; a++)
                    {
                        Butterflies.Butterfly(ref data[evenindex], ref data[oddindex], w);
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

        if (isInverse)
        {
            var scaleFactor = Math.ScaleB(1.0, -BitOperations.Log2((uint)data.Length));
            foreach (ref var c in data) { c *= scaleFactor; }
        }
    }
}
