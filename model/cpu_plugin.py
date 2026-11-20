import numpy as np
import scipy.io
import pedalboard
import time

# Load the TS808ADAA VST3 plugin
from pathlib import Path
# Release build produced by Projucer + Visual Studio 2022 (see README)
PLUGIN = str(Path(__file__).resolve().parent.parent / "plugin" / "Builds" / "VisualStudio2022" / "x64" / "Release"
             / "VST3" / "TS808ADAA.vst3" / "Contents" / "x86_64-win" / "TS808ADAA.vst3")

plugin = pedalboard.load_plugin(PLUGIN)
plugin.input_db = 0.0
plugin.output_db = 0.0

amp = 0.9798

d = scipy.io.loadmat(f"plugin_ref_amp{amp}.mat")                    # Load the reference data for the current amplitude
u = d["u_raw"].flatten().astype(np.float32)                         # Extract the raw input signal
fs = float(d["fs"].item())

u20 = np.tile(u, 10)
dur = len(u20) / fs

def medir(metodo, bypass, repeticiones=20):
    plugin.method = metodo
    plugin.bypass = bypass
    plugin.process(u20[:48000], fs, buffer_size=512, reset=True)
    tiempos = []
    for r in range(repeticiones):
        t0 = time.perf_counter()
        plugin.process(u20, fs, buffer_size=512, reset=True)
        t1 = time.perf_counter()
        tiempos.append(t1 - t0)
    return tiempos

tK = medir("K-method", False)
print(f"K-method : {min(tK)/len(u20)*1e9:6.1f} ns/muestra \t| CPU {min(tK)/dur*100:.3f} % \t| mediana {np.median(tK)/len(u20)*1e9:6.1f} ns")

tB = medir("K-method", True)
print(f"Bypass : {min(tB)/len(u20)*1e9:6.1f} ns/muestra \t| CPU {min(tB)/dur*100:.3f} % \t| mediana {np.median(tB)/len(u20)*1e9:6.1f} ns")

tK2 = medir("K x2", False)
print(f"K x2 : {min(tK2)/len(u20)*1e9:6.1f} ns/muestra \t| CPU {min(tK2)/dur*100:.3f} % \t| mediana {np.median(tK2)/len(u20)*1e9:6.1f} ns")

tK4 = medir("K x4", False)
print(f"K x4 : {min(tK4)/len(u20)*1e9:6.1f} ns/muestra \t| CPU {min(tK4)/dur*100:.3f} % \t| mediana {np.median(tK4)/len(u20)*1e9:6.1f} ns")

tA = medir("ADAA1", False)
print(f"ADAA1 : {min(tA)/len(u20)*1e9:6.1f} ns/muestra \t| CPU {min(tA)/dur*100:.3f} % \t| mediana {np.median(tA)/len(u20)*1e9:6.1f} ns")


    