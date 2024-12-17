# Super-heterodyne Receiver Modeling

## Project Description

The objective of this project is to simulate the operation of a superheterodyne receiver using MATLAB. The receiver processes multiple AM-modulated signals, downconverts them to an intermediate frequency (IF), and reconstructs the original baseband signals. Key steps include signal modulation, filtering, frequency conversion, and final detection.

The simulation involves:

- Modulating input stereo signals (provided as WAV files) to form a Frequency Division Multiplexed (FDM) signal.
- Implementing the core blocks of a superheterodyne receiver, including RF Bandpass Filters, Mixers, IF Bandpass Filters, and Baseband Detection stages.

## Implementation

1. **Input Signals:**

   - Stereo (2-channel) audio signals are converted to mono by summing the two channels.
   - Signals are padded with zeros to equalize lengths and upsampled to achieve higher sampling frequencies for accurate processing.

2. **Modulation and FDM Signal Generation:**

   - Each signal is modulated using DSB-SC (Double-Sideband Suppressed Carrier) modulation with carrier frequencies incrementing by 50 kHz.
   - The modulated signals are summed to create the FDM signal.

3. **RF Stage:**

   - Designed as a tunable Bandpass Filter (BPF) to isolate a specific modulated signal centered at the desired carrier frequency.

4. **Oscillator and Mixing:**

   - The oscillator generates a carrier for mixing with the filtered signal to downconvert it to an intermediate frequency (IF).

5. **IF Stage:**

   - An IF Bandpass Filter isolates the desired signal at the intermediate frequency (25 kHz).

6. **Baseband Detection:**

   - The IF signal is mixed with a local carrier and passed through a Lowpass Filter (LPF) to extract the baseband signal.

7. **Demodulation:**

   - The recovered signal is downsampled and played to verify its quality.

## MATLAB Key Functions

- **RF and IF Bandpass Filters:**
  Designed using MATLAB's `fdesign.bandpass` and `design` functions with Butterworth filters for sharp transition bands.

- **Mixer:**
  Implemented as a simple pointwise multiplication of the signal with the oscillator signal.

- **Lowpass Filter:**
  Designed with `fdesign.lowpass` to remove high-frequency components and recover the baseband signal.

## Output Figures

- Plots of the signal at various stages, including:
  - FDM signal spectrum.
  - Filtered RF and IF signals.
  - Final baseband signal in both time and frequency domains.
- Audio playback of the reconstructed baseband signal to verify successful demodulation.


more details is in attached pdf file.

