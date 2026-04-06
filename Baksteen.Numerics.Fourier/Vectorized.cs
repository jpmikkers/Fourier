namespace Baksteen.Numerics.Fourier;

using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;

public class Vectorized
{
    private static readonly Vector256<double> _negateOddMask = Vector256.Create(0.0, 0.0, -0.0, -0.0);

    public static unsafe void AvxButterfliesOnesUnrolled(Complex* pevenodd, int nrOfparts)
    {
        var p = 0;
        for (; p <= (nrOfparts - 2); p += 2)
        {
            var even0odd0 = Avx.LoadAlignedVector256((double*)pevenodd);       // (e0r,e0i,o0r,o0i)
            var evenegodd0 = Avx.Xor(even0odd0, _negateOddMask);             // (e0r,e0i,-o0r,-o0i)
            var odd0even0 = Avx.Permute2x128(even0odd0, even0odd0, 0x01);   // (o0r,o0i, e0r,e0i)

            var even1odd1 = Avx.LoadAlignedVector256((double*)(pevenodd + 2));
            var evenegodd1 = Avx.Xor(even1odd1, _negateOddMask);
            var odd1even1 = Avx.Permute2x128(even1odd1, even1odd1, 0x01);

            Avx.StoreAligned((double*)pevenodd, Avx.Add(odd0even0, evenegodd0));// (e0r+o0r,e0i+o0i,e0r-o0r,e0i-o0i)
            Avx.StoreAligned((double*)(pevenodd + 2), Avx.Add(odd1even1, evenegodd1));// (e1r+o1r,e1i+o1i,e1r-o0r,e1i-o0i)
            pevenodd += 4;
        }

        for (; p < nrOfparts; p++)
        {
            var even0odd0 = Avx.LoadAlignedVector256((double*)pevenodd);       // (e0r,e0i,o0r,o0i)
            var evenegodd = Avx.Xor(even0odd0, _negateOddMask);             // (e0r,e0i,-o0r,-o0i)
            var odd0even0 = Avx.Permute2x128(even0odd0, even0odd0, 0x01);   // (o0r,o0i, e0r,e0i)
            Avx.StoreAligned((double*)pevenodd, Avx.Add(odd0even0, evenegodd));// (e0r+o0r,e0i+o0i,e0r-o0r,e0i-o0i)
            pevenodd += 2;
        }
    }

    public static unsafe void AvxButterfliesOnes(Complex* pevenodd, int nrOfparts)
    {
        for (var p = 0; p < nrOfparts; p++)
        {
            var even0odd0 = Avx.LoadAlignedVector256((double*)pevenodd);       // (e0r,e0i,o0r,o0i)
            var evenegodd = Avx.Xor(even0odd0, _negateOddMask);             // (e0r,e0i,-o0r,-o0i)
            var odd0even0 = Avx.Permute2x128(even0odd0, even0odd0, 0x01);   // (o0r,o0i, e0r,e0i)
            var result = Avx.Add(odd0even0, evenegodd);                     // (e0r-o0r,e0i-o0i,e0r+o0r,e0i+o0i)
            Avx.StoreAligned((double*)pevenodd, result);
            pevenodd += 2;
        }

        //for (var p = 0; p < nrOfparts; p += 2)
        //{
        //    var even0odd0 = Avx.LoadAlignedVector256((double*)peven);                   // (e0r,e0i,o0r,o0i)
        //    var even1odd1 = Avx.LoadAlignedVector256((double*)(peven + 2));             // (e1r,e1i,o1r,o1i)

        //    var evens = Avx.Permute2x128(even0odd0, even1odd1, 0x20);                   // (e0r e0i e1r e1i)
        //    var odds = Avx.Permute2x128(even0odd0, even1odd1, 0x31);                    // (o0r o0i o1r o1i)

        //    var ra = Avx.Add(evens, odds);                                              // (e0r+o0r, e0i+o0i, e1r+o1r, e1i+o1i)
        //    var rb = Avx.Subtract(evens, odds);                                         // (e0r-o0r, e0i-o0i, e1r-o1r, e1i-o1i)

        //    Avx.StoreAligned((double*)peven, Avx.Permute2x128(ra, rb, 0x20));           // (e0r+o0r, e0i+o0i, e0r-o0r, e0i-o0i)
        //    Avx.StoreAligned((double*)(peven + 2), Avx.Permute2x128(ra, rb, 0x31));     // (e1r+o1r, e1i+o1i, e1r-o1r, e1i-o1i)
        //    peven += 4;
        //}
    }

    public static unsafe void Sse2ButterfliesOnes(Complex* pevenodd, int nrOfparts)
    {
        for (var p = 0; p < nrOfparts; p++)
        {
            var even0 = Sse2.LoadAlignedVector128((double*)pevenodd);
            var odd0 = Sse2.LoadAlignedVector128((double*)(pevenodd + 1));
            var re0 = Sse2.Add(even0, odd0);
            var ro0 = Sse2.Subtract(even0, odd0);
            Sse2.StoreAligned((double*)pevenodd, re0);
            Sse2.StoreAligned((double*)(pevenodd + 1), ro0);
            pevenodd += 2;
        }
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static unsafe void Sse2ButterflyOne(Complex* peven, Complex* podd)
    {
        var even0 = Sse2.LoadAlignedVector128((double*)peven);
        var odd0 = Sse2.LoadAlignedVector128((double*)podd);
        var re0 = Sse2.Add(even0, odd0);
        var ro0 = Sse2.Subtract(even0, odd0);
        Sse2.StoreAligned((double*)peven, re0);
        Sse2.StoreAligned((double*)podd, ro0);
    }

    public static unsafe void AvxButterflies(Complex* peven, Complex* pw, int numbtf, int wstride)
    {
        for (var t = 0; t < numbtf; t += 2)
        {
            var evens = Avx.LoadAlignedVector256((double*)peven);   // Avx.LoadVector256((double*)peven);
            var odds = Avx.LoadAlignedVector256((double*)(peven + numbtf)); // Avx.LoadVector256((double*)(peven + numbtf));
            var w = Vector256.Create(Sse2.LoadAlignedVector128((double*)pw), Sse2.LoadAlignedVector128((double*)(pw + wstride)));
            var oddw = Vectorized.ComplexMulAvx(odds, w);
            pw += (wstride << 1);

            var re = Avx.Add(evens, oddw);
            Avx.StoreAligned((double*)peven, re);

            var ro = Avx.Subtract(evens, oddw);
            Avx.StoreAligned((double*)(peven + numbtf), ro);

            peven += 2;
        }
    }

    public static unsafe void AvxButterfliesConjugate(Complex* peven, Complex* pw, int numbtf, int wstride)
    {
        for (var t = 0; t < numbtf; t += 2)
        {
            var evens = Avx.LoadAlignedVector256((double*)peven);   // Avx.LoadVector256((double*)peven);
            var odds = Avx.LoadAlignedVector256((double*)(peven + numbtf)); // Avx.LoadVector256((double*)(peven + numbtf));
            var w = Vector256.Create(Sse2.LoadAlignedVector128((double*)pw), Sse2.LoadAlignedVector128((double*)(pw + wstride)));

            var oddw = Vectorized.ComplexMulAvxConjugate(odds, w);
            pw += (wstride << 1);

            var re = Avx.Add(evens, oddw);
            Avx.StoreAligned((double*)peven, re);

            var ro = Avx.Subtract(evens, oddw);
            Avx.StoreAligned((double*)(peven + numbtf), ro);

            peven += 2;
        }
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector256<double> ComplexMulAvx(Vector256<double> a, Vector256<double> b)
    {
        // Multiplication:  (a + bi)(c + di) = (ac - bd) + (bc + ad)i
        var aRe = Avx.Permute(a, 0b00_00_00_00);              // a0r a0r a1r a1r
        var aIm = Avx.Permute(a, 0b00_00_11_11);              // a0i a0i a1i a1i
        var aImb = Avx.Multiply(aIm, b);                        // a0i*b0r a0i*b0i a1i*b1r a1i*b1i ..
        return Fma.MultiplyAddSubtract(aRe, b, Avx.Permute(aImb, 0b00_00_01_01));   // (a0r*b0r - a0i*b0i) (a0r*b0i + a0i*b0r) ...
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector256<double> ComplexMulAvxConjugate(Vector256<double> a_b, Vector256<double> c_d)
    {
        // Multiplication:  (a + bi)(c + (-di)) = (ac + bd) + (bc - ad)i
        // nodig:  ac bc,  dus a_b * c_c
        // en dan: bd ad,  dus b_a * d_d

        var c_c = Avx.Permute(c_d, 0b00_00_00_00);
        var d_d = Avx.Permute(c_d, 0b00_00_11_11);
        //var ac_bc = Avx.Multiply(a_b, c_c);
        var b_a = Avx.Permute(a_b, 0b00_00_01_01);
        var bd_ad = Avx.Multiply(b_a, d_d);
        return Fma.MultiplySubtractAdd(a_b, c_c, bd_ad);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector128<double> ComplexMulSse2(Vector128<double> ab, Vector128<double> cd)
    {
        // Multiplication:  (a + bi)(c + di) = (ac -bd) + (bc + ad)i
        var c_c = Sse2.UnpackLow(cd, cd);        // var c_c = Sse3.MoveAndDuplicate(cd); is a bit slower?
        var d_d = Sse2.UnpackHigh(cd, cd);
        var b_a = Sse2.Shuffle(ab, ab, 0x01);
        //var ac_bc = Sse2.Multiply(ab, c_c);
        //var bd_ad = Sse2.Multiply(b_a, d_d);
        //return Sse3.AddSubtract(ac_bc, bd_ad);  // (ac - bd) (bc + ad)
        //var ac_bc = Sse2.Multiply(ab, c_c);
        //var bd_ad = Sse2.Multiply(b_a, d_d);
        return Sse3.AddSubtract(Sse2.Multiply(ab, c_c), Sse2.Multiply(b_a, d_d));  // (ac - bd) (bc + ad)
    }
}
