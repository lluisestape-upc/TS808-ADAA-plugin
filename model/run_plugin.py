import numpy as np
import scipy.io
import pedalboard

# Load the TS808ADAA VST3 plugin
from pathlib import Path
# Release build produced by Projucer + Visual Studio 2022 (see README)
PLUGIN = str(Path(__file__).resolve().parent.parent / "plugin" / "Builds" / "VisualStudio2022" / "x64" / "Release"
             / "VST3" / "TS808ADAA.vst3" / "Contents" / "x86_64-win" / "TS808ADAA.vst3")

plugin = pedalboard.load_plugin(PLUGIN)
plugin.input_db = 0.0
plugin.output_db = 0.0

for amp in ["0.0981", "0.9798"]:
    d = scipy.io.loadmat(f"plugin_ref_amp{amp}.mat")                    # Load the reference data for the current amplitude
    u = d["u_raw"].flatten().astype(np.float32)                         # Extract the raw input signal
    fs = float(d["fs"].item())                                          # Extract the sample rate       

    plugin.method = "K-method"
    yK = plugin.process(u, sample_rate=fs, reset=True)                      # Process the input signal using the K-method

    plugin.method = "K x2"
    yK2 = plugin.process(u, sample_rate=fs, reset=True)                     # Process the input signal using the K-method 2x

    plugin.method = "K x4"
    yK4 = plugin.process(u, sample_rate=fs, reset=True)                     # Process the input signal using the K-method 4x

    plugin.method = "ADAA1"
    yA = plugin.process(u, sample_rate=fs, reset=True)                      # Process the input signal using the ADAA1 method

    scipy.io.savemat(f"plugin_out_amp{amp}.mat", {"yK": yK, "yA": yA})      # Save the processed output for the current amplitude
    print(amp, "ok", yK.shape, yA.shape)                                    # Print the amplitude, status, and shapes of the processed outputs