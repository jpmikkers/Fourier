using Baksteen.Numerics.Fourier;
using ScottPlot;
using System.Numerics;

//Vector256<double> a = Vector256.Create(1.0, 2.0, 3.0, 4.0);
//Vector256<double> b = Vector256.Create(5.0, 6.0, 7.0, 8.0);

//var result1 = Vectorized.ComplexMulAvx(a, b);
//var result2 = new Complex(1, 2) * new Complex(5, 6);
//var result3 = new Complex(3, 4) * new Complex(7, 8);

//var result4 = Vectorized.ComplexMulSse2(Vector128.Create(1.0, 2.0), Vector128.Create(5.0, 6.0));
//var result5 = Vectorized.ComplexMulSse2(Vector128.Create(3.0, 4.0), Vector128.Create(7.0, 8.0));


////ArrayUtils.ReorderEvenOddRelative2(even);
////Console.WriteLine("n=6 result : " + string.Join(" ", blah));
//return;

var signal = new double[16];
for (var i = 0; i < signal.Length; i++)
{
    signal[i] = Random.Shared.NextDouble();
}

//for (var i = 0; i < signal.Length; i++)
//{
//    signal[i] =
//        Math.Sin(i * 20 * Math.Tau / signal.Length)
//        + Math.Sin(i * 40 * Math.Tau / signal.Length)
//        + Math.Sin(i * 80 * Math.Tau / signal.Length)
//        + Math.Sin(i * 160 * Math.Tau / signal.Length)
//        + Math.Sin(i * 320 * Math.Tau / signal.Length);
//}

var plot1 = new Plot();
plot1.Add.Signal(signal);
plot1.Title("Signal");
plot1.SavePng("signal.png", 1024, 768);

Complex[] complexSignal = [.. signal.Select(x => (Complex)x)];
var spectrum = DFT.DiscreteFourierTransform(complexSignal, forward: true);

var plot2 = new Plot();
plot2.Add.Signal(spectrum.Select(Complex.Abs).ToList());
plot2.Title("DFT Spectrum");
plot2.SavePng("spectrum dft.png", 1024, 768);

var spectrum_alt = complexSignal.ToArray();
//new Fft64(spectrum_alt.Length).Direct(spectrum_alt, isInverse: false);
//FFTL.FastFourierTransform(spectrum_alt, isInverse: false);
//new FFTM(spectrum_alt.Length).FastFourierTransform(spectrum_alt, isInverse: false);
new FFTSimpleVectorizedE(spectrum.Length).FastFourierTransform(spectrum_alt, false);
//RecursiveFFTE.FastFourierTransform(spectrum_alt, isInverse: false);

var plot2alt = new Plot();
plot2alt.Add.Signal(spectrum_alt.Select(Complex.Abs).ToList());
plot2alt.Title("FFT Spectrum");
plot2alt.SavePng("spectrum fft.png", 1024, 768);

FFTH.FastFourierTransform(spectrum_alt, isInverse: true);
var reconstructed = spectrum_alt.Select(x => x.Real).ToArray();
//var reconstructed = DFT.DiscreteFourierTransform(spectrum_alt, forward: false).Select(x => x.Real).ToArray();

var plot3 = new Plot();
plot3.Add.Signal(reconstructed);
plot3.Title("Reconstructed signal");
plot3.SavePng("reconstructed.png", 1024, 768);

