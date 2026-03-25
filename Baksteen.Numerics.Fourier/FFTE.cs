namespace Baksteen.Numerics.Fourier;

using System;
using System.Numerics;

public static class FFTE
{
    private static readonly Complex[] _rotations;

    static FFTE()
    {
        // contains 1/2^n rad for n=0-31, so 1/1 rad, 1/2 rad, 1/4 rad, 1/8 rad etc..
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

        var anglesPerPart = 2;              // a single butterfly does 2 angles, +w and -w (=w+pi radians)
        var nrOfParts = data.Length >> 1;   // so the first layer is len/2 parts of single butterflies
        var rotationLookupIndex = 1;

        if (nrOfParts > 0)
        {
            // no mul needed in the first combination layer, every w is 1+0i and -1+0i
            for (var p = 0; p < nrOfParts; p++)
            {
                Butterflies.Butterfly(ref data[p << 1], ref data[(p << 1) + 1]);
            }
            anglesPerPart <<= 1;
            nrOfParts >>= 1;
            rotationLookupIndex++;
        }

        while (nrOfParts > 0)
        {
            var wr = isInverse ? Complex.Conjugate(_rotations[rotationLookupIndex]) : _rotations[rotationLookupIndex];

            //Console.WriteLine($"nrOfParts {nrOfParts} anglesPerPart {anglesPerPart} {rotationLookupIndex}");

            for (var p = 0; p < nrOfParts; p++)
            {
                var w = Complex.One;
                var evenindex = p << rotationLookupIndex;
                var oddindex = evenindex + (anglesPerPart >> 1);

                for (var a = 0; a < anglesPerPart; a += 2)
                {
                    Butterflies.Butterfly(ref data[evenindex], ref data[oddindex], w);
                    w *= wr;
                    evenindex++;
                    oddindex++;
                }

                //for (var a = 0; a < anglesPerPart; a += 2)
                //{
                //    var evenindex = (p << rotationLookupIndex) + (a >> 1);
                //    var oddindex = evenindex + (anglesPerPart >> 1);
                //    Console.WriteLine($"  part {p} {a}\tei {evenindex} eo {oddindex}");
                //    Butterfly(ref data[evenindex], ref data[oddindex], w);
                //    w *= wr;
                //}
            }

            anglesPerPart <<= 1;
            nrOfParts >>= 1;
            rotationLookupIndex++;
        }

        if (isInverse)
        {
            var scaleFactor = Math.ScaleB(1.0, -BitOperations.Log2((uint)data.Length));
            foreach (ref var c in data)
            {
                c *= scaleFactor;
            }
        }
    }
}
