namespace Baksteen.Numerics.Fourier;

using System;
using System.Buffers.Binary;
using System.Numerics;
using System.Runtime.CompilerServices;

static class Reorder
{
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
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

#if NEVER
    public static void TestFastFourierTransform(Span<Complex> data, bool isInverse)
    {
        if (!BitOperations.IsPow2(data.Length))
        {
            throw new ArgumentException("fft not a power of two", nameof(data));
        }

        Reorder.Shuffle(data);

        if (data.Length == 2)
        {
            Butterfly(ref data[0], ref data[1], Complex.One);
        }
        else if (data.Length == 4)
        {
            Butterfly(ref data[0], ref data[1], Complex.One);
            Butterfly(ref data[2], ref data[3], Complex.One);

            Butterfly(ref data[0], ref data[2], Complex.One);
            Butterfly(ref data[1], ref data[3], _rotations[2]);
        }
        else if (data.Length == 8)
        {
            Butterfly(ref data[0], ref data[1], Complex.One);
            Butterfly(ref data[2], ref data[3], Complex.One);
            Butterfly(ref data[4], ref data[5], Complex.One);
            Butterfly(ref data[6], ref data[7], Complex.One);

            Butterfly(ref data[0], ref data[2], Complex.One);                                       // w0_8 w4_8
            Butterfly(ref data[1], ref data[3], _rotations[3] * _rotations[3]);                     // w2_8 w6_8

            Butterfly(ref data[4], ref data[6], Complex.One);                                       // w0_8 w4_8
            Butterfly(ref data[5], ref data[7], _rotations[3] * _rotations[3]);                     // w2_8 w6_8

            Butterfly(ref data[0], ref data[4], Complex.One);                                       // w0_8 w4_8
            Butterfly(ref data[1], ref data[5], _rotations[3]);                                     // w1_8 w5_8
            Butterfly(ref data[2], ref data[6], _rotations[3] * _rotations[3]);                     // w2_8 w6_8
            Butterfly(ref data[3], ref data[7], _rotations[3] * _rotations[3] * _rotations[3]);     // w3_8 w7_8
        }
        else if (data.Length == 16)
        {
            var w1_16 = _rotations[4];

            Butterfly(ref data[0], ref data[1], Complex.One);
            Butterfly(ref data[2], ref data[3], Complex.One);
            Butterfly(ref data[4], ref data[5], Complex.One);
            Butterfly(ref data[6], ref data[7], Complex.One);
            Butterfly(ref data[8], ref data[9], Complex.One);
            Butterfly(ref data[10], ref data[11], Complex.One);
            Butterfly(ref data[12], ref data[13], Complex.One);
            Butterfly(ref data[14], ref data[15], Complex.One);

            Butterfly(ref data[0], ref data[2], Complex.One);                                   // w0_16 w8_16
            Butterfly(ref data[1], ref data[3], Complex.Pow(w1_16, 4));                         // w4_16 w12_16

            Butterfly(ref data[4], ref data[6], Complex.One);                                   // w0_16 w8_16
            Butterfly(ref data[5], ref data[7], Complex.Pow(w1_16, 4));                         // w4_16 w12_16

            Butterfly(ref data[8], ref data[10], Complex.One);                                  // w0_16 w8_16
            Butterfly(ref data[9], ref data[11], Complex.Pow(w1_16, 4));                        // w4_16 w12_16 

            Butterfly(ref data[12], ref data[14], Complex.One);                                 // w0_16 w8_16
            Butterfly(ref data[13], ref data[15], Complex.Pow(w1_16, 4));                       // w4_16 w12_16                      // w2_8 w6_8

            Butterfly(ref data[0], ref data[4], Complex.One);                                   // w0_16 w8_16
            Butterfly(ref data[1], ref data[5], Complex.Pow(w1_16, 2));                         // w2_16 w10_16
            Butterfly(ref data[2], ref data[6], Complex.Pow(w1_16, 4));                         // w4_16 w12_16
            Butterfly(ref data[3], ref data[7], Complex.Pow(w1_16, 6));                         // w6_16 w14_16

            Butterfly(ref data[8], ref data[12], Complex.One);                                  // w0_16 w8_16
            Butterfly(ref data[9], ref data[13], Complex.Pow(w1_16, 2));                        // w2_16 w10_16 
            Butterfly(ref data[10], ref data[14], Complex.Pow(w1_16, 4));                       // w4_16 w12_16
            Butterfly(ref data[11], ref data[15], Complex.Pow(w1_16, 6));                       // w6_16 w14_16                      // w2_8 w6_8

            Butterfly(ref data[0], ref data[8], Complex.One);                                   // w0_16 w8_16
            Butterfly(ref data[1], ref data[9], Complex.Pow(w1_16, 1));                         // w1_16 w9_16
            Butterfly(ref data[2], ref data[10], Complex.Pow(w1_16, 2));                        // w2_16 w10_16
            Butterfly(ref data[3], ref data[11], Complex.Pow(w1_16, 3));                        // w3_16 w11_16
            Butterfly(ref data[4], ref data[12], Complex.Pow(w1_16, 4));                        // w4_16 w12_16
            Butterfly(ref data[5], ref data[13], Complex.Pow(w1_16, 5));                        // w5_16 w13_16
            Butterfly(ref data[6], ref data[14], Complex.Pow(w1_16, 6));                        // w6_16 w14_16
            Butterfly(ref data[7], ref data[15], Complex.Pow(w1_16, 7));                        // w7_16 w15_16
        }
    }
#endif

}
