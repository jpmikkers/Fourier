using Baksteen.Numerics.Fourier;
using ScottPlot;
using System.Numerics;

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
new FFTM(spectrum_alt.Length).FastFourierTransform(spectrum_alt, isInverse: false);

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
