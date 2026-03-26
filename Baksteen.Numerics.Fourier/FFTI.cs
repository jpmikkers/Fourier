namespace Baksteen.Numerics.Fourier;

using System;
using System.Numerics;

public static class FFTI
{
#if NEVER
    public static void FastFourierTransformO(Span<Complex> data, bool isInverse)
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
            var wr = isInverse ? Complex.Conjugate(_rotations[rotationLookupIndex]) : _rotations[rotationLookupIndex];
            var w = Complex.One;
            var evenindex = 0;
            var oddindex = evenindex + butterfliesPerPart;

            if (butterfliesPerPart <= maxstackbutterflies)
            {
                for (var a = 0; a < butterfliesPerPart; a++)
                {
                    tw[a] = w;
                    w *= wr;
                }

                for (var a = 0; a < butterfliesPerPart; a++)
                {
                    for (var p = 0; p < nrOfParts; p++)
                    {
                        Butterflies.Butterfly(ref data[evenindex + (p << rotationLookupIndex)], ref data[oddindex + (p << rotationLookupIndex)], tw[a]);
                    }
                    evenindex++;
                    oddindex++;
                }
            }
            else
            {
                for (var a = 0; a < butterfliesPerPart; a++)
                {
                    for (var p = 0; p < nrOfParts; p++)
                    {
                        Butterflies.Butterfly(ref data[evenindex + (p << rotationLookupIndex)], ref data[oddindex + (p << rotationLookupIndex)], w);
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
#endif

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
                var w = Complex.One;

                for (var a = 0; a < butterfliesPerPart; a++)
                {
                    tw[a] = w;
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
