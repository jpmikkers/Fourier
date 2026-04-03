namespace Benchmark;

using Baksteen.Numerics.Fourier;
using BenchmarkDotNet.Attributes;
using System;
using System.Linq;
using System.Numerics;

// For more information on the VS BenchmarkDotNet Diagnosers see https://learn.microsoft.com/visualstudio/profiling/profiling-with-benchmark-dotnet
//[CPUUsageDiagnoser]
public class Benchmarks
{
    //private SHA256 sha256 = SHA256.Create();
    //private byte[] data;
    private AlignedMemoryManager<Complex> alignedMemoryManager;
    private const int fftsize = 4096;
    private const int seed = 42;
    private const int repeats = 100;
    private Complex[] data;
    private Memory<Complex> tmpdata;
    private Fft64 fft64;
    private FFTM fftm;
    private FFTSimpleVectorizedC fftsvc;
    private FFTSimpleVectorizedD fftsvd;
    private FFTSimpleVectorizedE fftsve;
    private FFTSimpleVectorizedF fftsvf;
    private FFTSimpleVectorizedG fftsvg;
    private FFTSimpleVectorizedH fftsvh;
    private FFTAvxVectorizedI fftsvi;
    private FFTSimpleBigLut fftsbl;

    [GlobalSetup]
    public void Setup()
    {
        var rnd = new Random(seed);
        data = [.. Enumerable.Range(0, fftsize).Select(_ => new Complex(rnd.NextDouble(), rnd.NextDouble()))];
        alignedMemoryManager = new AlignedMemoryManager<Complex>(fftsize, 32);
        tmpdata = alignedMemoryManager.Memory;
        fft64 = new Fft64(fftsize);
        fftm = new FFTM(fftsize);
        fftsvc = new FFTSimpleVectorizedC(fftsize);
        fftsvd = new FFTSimpleVectorizedD(fftsize);
        fftsve = new FFTSimpleVectorizedE(fftsize);
        fftsvf = new FFTSimpleVectorizedF(fftsize);
        fftsvg = new FFTSimpleVectorizedG(fftsize);
        fftsvh = new FFTSimpleVectorizedH(fftsize);
        fftsvi = new FFTAvxVectorizedI(fftsize);
        fftsbl = new FFTSimpleBigLut(fftsize);
    }

    private void ResetTmpData()
    {
        data.AsSpan().CopyTo(tmpdata.Span);
    }

    [Benchmark(Baseline = true)]
    public void BenchFFT64()
    {
        ResetTmpData();
        fft64.Direct(tmpdata.Span, false);
    }

    //[Benchmark]
    //public void BenchFFTK()
    //{
    //    Array.Copy(data, tmpdata, data.Length);
    //    FFTK.FastFourierTransform(tmpdata, false);
    //}

    //[Benchmark()]
    //public void BenchFFTL()
    //{
    //    ResetTmpData();
    //    FFTL.FastFourierTransform(tmpdata.Span, false);
    //}

    //[Benchmark()]
    //public void BenchFFTM()
    //{
    //    Array.Copy(data, tmpdata, data.Length);
    //    fftm.FastFourierTransform(tmpdata, false);
    //}

    //[Benchmark()]
    //public void BenchFFTSimple()
    //{
    //    Array.Copy(data, tmpdata, data.Length);
    //    FFTSimple.FastFourierTransform(data, false);
    //}

    //[Benchmark(Baseline = true)]
    //public void BenchFFTSimpleVectorizedC()
    //{
    //    ResetTmpData();
    //    fftsvc.FastFourierTransform(tmpdata.Span, false);
    //}

    //[Benchmark()]
    //public void BenchFFTSimpleVectorizedE()
    //{
    //    ResetTmpData();
    //    fftsve.FastFourierTransform(tmpdata.Span, false);
    //}

    //[Benchmark(Baseline = true)]
    [Benchmark()]
    public void BenchFFTSimpleVectorizedF()
    {
        ResetTmpData();
        fftsvf.FastFourierTransform(tmpdata.Span, false);
    }

    [Benchmark()]
    public void BenchFFTSimpleVectorizedG()
    {
        ResetTmpData();
        fftsvg.FastFourierTransform(tmpdata.Span, false);
    }

    //[Benchmark()]
    //public void BenchFFTSimpleVectorizedH()
    //{
    //    ResetTmpData();
    //    fftsvh.FastFourierTransform(tmpdata.Span, false);
    //}

    [Benchmark()]
    public void BenchFFTSimpleVectorizedI()
    {
        ResetTmpData();
        fftsvi.FastFourierTransform(tmpdata.Span, false);
    }

    //[Benchmark()]
    //public void BenchFFTSimpleBigLut()
    //{
    //    ResetTmpData();
    //    fftsbl.FastFourierTransform(tmpdata.Span, false);
    //}

    //[Benchmark()]
    //public void BenchFFTSimpleVectorizedD()
    //{
    //    Array.Copy(data, tmpdata, data.Length);
    //    fftsvd.FastFourierTransform(data, false);
    //}

    //[Benchmark()]
    //public void BenchRecursiveFFTD()
    //{
    //    Array.Copy(data, tmpdata, data.Length);
    //    RecursiveFFTD.FastFourierTransform(data);
    //    //fftm.FastFourierTransform(tmpdata, false);
    //}

    //[Benchmark()]
    //public void BenchRecursiveFFTE()
    //{
    //    Array.Copy(data, tmpdata, data.Length);
    //    RecursiveFFTE.FastFourierTransform(data, false);
    //    //fftm.FastFourierTransform(tmpdata, false);
    //}

    //[Benchmark()]
    //public void BenchFFTE()
    //{
    //    Array.Copy(data, tmpdata, data.Length);
    //    FFTE.FastFourierTransform(tmpdata, false);
    //}

    //[Benchmark()]
    //public void BenchFFTF()
    //{
    //    Array.Copy(data, tmpdata, data.Length);
    //    FFTF.FastFourierTransform(tmpdata, false);
    //}

    //[Benchmark()]
    //public void BenchFFTG()
    //{
    //    Array.Copy(data, tmpdata, data.Length);
    //    FFTG.FastFourierTransform(tmpdata, false);
    //}

    //[Benchmark()]
    //public void BenchFFTH()
    //{
    //    Array.Copy(data, tmpdata, data.Length);
    //    FFTH.FastFourierTransform(tmpdata, false);
    //}

    //[Benchmark()]
    //public void BenchFFTI()
    //{
    //    Array.Copy(data, tmpdata, data.Length);
    //    FFTI.FastFourierTransform(tmpdata, false);
    //}

    //[Benchmark()]
    //public void BenchFFTJ()
    //{
    //    Array.Copy(data, tmpdata, data.Length);
    //    FFTJ.FastFourierTransform(tmpdata, false);
    //}
}
