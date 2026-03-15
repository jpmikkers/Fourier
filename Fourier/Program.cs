using Fourier;
using ScottPlot;
using System.Numerics;

var signal = new double[1024];
for (var i = 0; i < signal.Length; i++)
{
    signal[i] =
        Math.Sin(i * 20 * Math.Tau / signal.Length)
        + Math.Sin(i * 40 * Math.Tau / signal.Length)
        + Math.Sin(i * 80 * Math.Tau / signal.Length)
        + Math.Sin(i * 160 * Math.Tau / signal.Length)
        + Math.Sin(i * 320 * Math.Tau / signal.Length);
}

var plot1 = new Plot();
plot1.Add.Signal(signal);
plot1.Title("Signal");
plot1.SavePng("signal.png", 1024, 768);

var spectrum = DFT.DiscreteFourierTransform([.. signal.Select(x => (Complex)x)], forward: true);

var plot2 = new Plot();
plot2.Add.Signal(spectrum.Select(Complex.Abs).ToList());
plot2.Title("Spectrum");
plot2.SavePng("spectrum.png", 1024, 768);

var reconstructed = DFT.DiscreteFourierTransform(spectrum, forward: false).Select(x => x.Real).ToArray();

var plot3 = new Plot();
plot3.Add.Signal(reconstructed);
plot3.Title("Reconstructed signal");
plot3.SavePng("reconstructed.png", 1024, 768);
