namespace Fourier;

using System;
using System.Numerics;
using System.Runtime.CompilerServices;

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

        if (nrOfParts > 0)
        {
            // first combination layer is a special case, every w is 1+0i and -1+0i, no mul needed
            for (var p = 0; p < nrOfParts; p++)
            {
                Butterfly(ref data[p << 1], ref data[(p << 1) + 1]);
            }
            anglesPerPart <<= 1;
            nrOfParts >>= 1;
        }

        while (nrOfParts > 0)
        {
            var widxinc = data.Length / anglesPerPart;

            //var nrOfParts = data.Length / anglesPerPart;
            Console.WriteLine($"nrOfParts {nrOfParts} anglesPerPart {anglesPerPart}");

            for (var p = 0; p < nrOfParts; p++)
            {
                for (var a = 0; a < anglesPerPart; a += 2)
                {
                    var evenindex = p * anglesPerPart + (a >> 1);
                    var oddindex = evenindex + (anglesPerPart >> 1);
                    var we = (a >> 1) * data.Length / anglesPerPart;
                    var wo = we + (data.Length >> 1);
                    Console.WriteLine($"  part {p} {a}\tei {evenindex} eo {oddindex}\tangle w{we} w{wo}\twidxinc{widxinc} {nrOfParts}");
                    var w = Complex.FromPolarCoordinates(1.0, -we * Math.Tau / data.Length);
                    Butterfly(ref data[evenindex], ref data[oddindex], w);
                }
            }

            anglesPerPart <<= 1;
            nrOfParts >>= 1;
        }
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void Butterfly(ref Complex even, ref Complex odd, Complex w)
    {
        var odd_w = odd * w;
        (even, odd) = (even + odd_w, even - odd_w);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void Butterfly(ref Complex even, ref Complex odd)
    {
        (even, odd) = (even + odd, even - odd);
    }
}
