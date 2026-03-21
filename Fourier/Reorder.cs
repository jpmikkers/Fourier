namespace Fourier;

using System;
using System.Buffers.Binary;
using System.Numerics;

static class Reorder
{
    private static uint BitReverse(uint x, int resultbits)
    {
        x = ((x & 0xaaaaaaaa) >> 1) | ((x & 0x55555555) << 1);
        x = ((x & 0xcccccccc) >> 2) | ((x & 0x33333333) << 2);
        x = ((x & 0xf0f0f0f0) >> 4) | ((x & 0x0f0f0f0f) << 4);
        return BinaryPrimitives.ReverseEndianness(x) >> (32 - resultbits);
    }

    public static void Shuffle<T>(Span<T> samples)
    {
        var numBits = BitOperations.Log2((uint)samples.Length);
        for (var i = 1; i < samples.Length - 1; i++)
        {
            var j = (int)BitReverse((uint)i, numBits);
            if (i < j)
            {
                (samples[i], samples[j]) = (samples[j], samples[i]);
            }
        }
    }
}
