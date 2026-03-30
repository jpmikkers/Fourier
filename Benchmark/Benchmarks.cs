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
    private const int fftsize = 4096;
    private const int seed = 42;
    private const int repeats = 100;
    private Complex[] data;
    private Complex[] tmpdata;
    private Fft64 fft64;
    private FFTM fftm;

    [GlobalSetup]
    public void Setup()
    {
        var rnd = new Random(seed);
        data = [.. Enumerable.Range(0, fftsize).Select(_ => new Complex(rnd.NextDouble(), rnd.NextDouble()))];
        tmpdata = new Complex[fftsize];
        fft64 = new Fft64(fftsize);
        fftm = new FFTM(fftsize);
    }

    //[Benchmark(Baseline = true)]
    public void BenchFFT64()
    {
        Array.Copy(data, tmpdata, data.Length);
        fft64.Direct(tmpdata, false);
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
    //    Array.Copy(data, tmpdata, data.Length);
    //    FFTL.FastFourierTransform(tmpdata, false);
    //}

    //[Benchmark()]
    //public void BenchFFTM()
    //{
    //    Array.Copy(data, tmpdata, data.Length);
    //    fftm.FastFourierTransform(tmpdata, false);
    //}

    [Benchmark(Baseline = true)]
    public void BenchFFTSimple()
    {
        Array.Copy(data, tmpdata, data.Length);
        FFTSimple.FastFourierTransform(data, false);
    }

    [Benchmark()]
    public void BenchFFTSimpleVectorizedB()
    {
        Array.Copy(data, tmpdata, data.Length);
        FFTSimpleVectorizedB.FastFourierTransform(data, false);
    }

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
