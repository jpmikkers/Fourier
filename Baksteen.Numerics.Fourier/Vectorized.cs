namespace Baksteen.Numerics.Fourier;

using System.Runtime.CompilerServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;

public class Vectorized
{
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Vector256<double> ComplexMulAvx(Vector256<double> a, Vector256<double> b)
    {
        var aRe = Avx.Permute(a, 0b00_00_00_00);              // a0r a0r a1r a1r
        var aIm = Avx.Permute(a, 0b00_00_11_11);              // a0i a0i a1i a1i
        var aImb = Avx.Multiply(aIm, b);                        // a0i*b0r a0i*b0i a1i*b1r a1i*b1i ..
        return Fma.MultiplyAddSubtract(aRe, b, Avx.Permute(aImb, 0b00_00_01_01));   // (a0r*b0r - a0i*b0i) (a0r*b0i + a0i*b0r) ...
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
