namespace Baksteen.Numerics.Fourier;

using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;

internal static class Butterflies
{
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static (Complex even, Complex odd) ButterflyTuple(Complex even, Complex odd, Complex w)
    {
        var odd_w = KaratsubaMultiply(w, odd);
        return (even + odd_w, even - odd_w);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static unsafe void ButterflyVectorized(ref Complex even, ref Complex odd, Vector128<double> w)
    {
        //var s_even = MemoryMarshal.Cast<Complex, double>(MemoryMarshal.CreateSpan(ref even, 1));
        //var s_odd = MemoryMarshal.Cast<Complex, double>(MemoryMarshal.CreateSpan(ref odd, 1));
        //var s_w = MemoryMarshal.Cast<Complex, double>(MemoryMarshal.CreateReadOnlySpan(ref w, 1));

        fixed (Complex* pce = &even)
        fixed (Complex* pco = &odd)
        {
            var pe = (double*)pce;
            var po = (double*)pco;

            var vec_even = Sse2.LoadVector128(pe);
            var vec_odd_w = Vectorized.ComplexMulSse2(Sse2.LoadVector128(po), w);

            var result_even = Sse2.Add(vec_even, vec_odd_w);
            var result_odd = Sse2.Subtract(vec_even, vec_odd_w);

            Sse2.Store(pe, result_even);
            Sse2.Store(po, result_odd);
            //result_even.Store(pe);
            //result_odd.Store(po);
        }
        //var odd_w = w * odd;
        //(even, odd) = (even + odd_w, even - odd_w);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static (Vector128<double> re, Vector128<double> ro) ButterflyVectorized(Vector128<double> even, Vector128<double> odd, Vector128<double> w)
    {
        var vec_odd_w = Vectorized.ComplexMulSse2(odd, w);
        return (Sse2.Add(even, vec_odd_w), Sse2.Subtract(even, vec_odd_w));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static void Butterfly(ref Complex even, ref Complex odd, Complex w)
    {
        var odd_w = w * odd;
        (even, odd) = (even + odd_w, even - odd_w);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static void Butterfly2(ref Complex even, ref Complex odd, Complex w)
    {
        var tmpe = even;
        var odd_w = w * odd;
        even = tmpe + odd_w;
        odd = tmpe - odd_w;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static void Butterfly(ref Complex even, ref Complex odd)
    {
        (even, odd) = (even + odd, even - odd);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static void Butterfly90CW(ref Complex even, ref Complex odd)
    {
        var odd_w = new Complex(odd.Imaginary, -odd.Real);  // 90 degree clockwise
        (even, odd) = (even + odd_w, even - odd_w);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static void Butterfly90CCW(ref Complex even, ref Complex odd)
    {
        var odd_w = new Complex(-odd.Imaginary, odd.Real);  // 90 degree counterclockwise
        (even, odd) = (even + odd_w, even - odd_w);
    }

    // Multiply using the 3‑multiply algorithm
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Complex KaratsubaMultiply(Complex x, Complex y)
    {
        var p = x.Real * y.Real;
        var q = x.Imaginary * y.Imaginary;
        var r = (x.Real + x.Imaginary) * (y.Real + y.Imaginary);
        return new Complex(p - q, r - p - q);
    }
}
