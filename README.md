# TS808 ADAA plugin

Clipping stage of the Ibanez TS808 Tube Screamer, simulated with the K-method
(state-space, trapezoidal rule) and first-order antiderivative antialiasing (ADAA1),
implemented as a JUCE VST3 plugin. Part of my bachelor's thesis (TFG) at UPC.

The model and the ADAA1 formulation follow Martin Jedla's
[K_method_ADAA](https://github.com/jedlamartin/K_method_ADAA) (included as the
`jedla_ref` submodule, commit `dc32f76`). Code comments are in Spanish.

## Contents

| Folder | What it contains |
|---|---|
| `plugin/` | Projucer project and C++ source of the plugin |
| `model/` | MATLAB reference model, metrics and the Python scripts that drive the plugin |
| `ltspice/` | LTspice schematics of the full pedal, used as ground truth |
| `jedla_ref/` | Martin Jedla's original MATLAB code (submodule) |

## Plugin

Parameters: input gain, output gain and method:

- **K-method**: K-method at the host sample rate.
- **ADAA1**: K-method with ADAA1. The model is computed at 2/3 fs (rescaled
  sample rate), with input compensation and averaging of the parallel path.
  F(p[n-1]) is cached, so there is one table lookup per sample.
- **K x2 / K x4**: K-method with oversampling (`juce::dsp::Oversampling`,
  polyphase IIR half-band filters). Work in progress: the latency is not yet
  reported to the host.

Only the selected method is computed.

Build: open `plugin/TS808ADAA.jucer` in Projucer (modules include `juce_dsp`),
export to Visual Studio 2022 and build the VST3 target in Release.

## Validation against MATLAB

1. `model/ts_sine_real.m`: runs the reference model for a 1073 Hz sine
   (amplitude 0.0981 or 0.9798, set in the script) and saves
   `plugin_ref_amp<amp>.mat`. Run it once for each amplitude.
2. `model/run_plugin.py`: loads the Release VST3 with
   [pedalboard](https://github.com/spotify/pedalboard), processes the same input
   with each method and saves `plugin_out_amp<amp>.mat`.
3. `model/compare_plugin.m`: compares the plugin with MATLAB sample by sample
   and computes the metrics in `alias_metrics.m`.

`alias_metrics.m`: Hann window, FFT zero-padded to 10N, the bins within ±5 Hz of
each harmonic are set to zero, and the remaining power is summed over 0 to fs/2
and over 0 to 20 kHz (in dB). It also returns the fundamental amplitude H1.

## CPU measurement

`model/cpu_plugin.py` uses pedalboard and `time.perf_counter`. The plugin
processes the reference signal repeated 10 times, in 512-sample buffers. Each
method is run 20 times (after one warm-up run) and the minimum time is reported
in ns per sample, together with the bypass time as a baseline. Ratios between
methods are measured in the same run.

## LTspice

- `ts808_tube_screamer.asc`: full pedal (Rgain = 200k), transient analysis with
  a 1073 Hz sine at A = 0.1 and 1 V, and `.four` on the clipping node.
- `ts808_ac_sin_diodos.asc`: the same circuit without the diodes, AC analysis
  of the linear part.
