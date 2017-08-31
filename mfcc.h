/* ETSI ES 201 108 V1.1.3 implementation
 *  Written in 2017 by Marc Burns (m4burns@uwaterloo.ca)
 *  Public domain
 */

#include <functional>
#include <array>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fftw3.h>

/* Table 4.1 */
template<int sampleRate>
struct SampleRateParameters { };

template<>
struct SampleRateParameters<8000> {
  static const int frameLength = 200, shiftInterval = 80, fftLength = 256;
};

template<>
struct SampleRateParameters<11000> {
  static const int frameLength = 256, shiftInterval = 110, fftLength = 256;
};

template<>
struct SampleRateParameters<16000> {
  static const int frameLength = 400, shiftInterval = 160, fftLength = 512;
};

template<int sampleRate = 8000>
class MFCCExtractor {
public:
  using params = SampleRateParameters<sampleRate>;

  constexpr static const int
    mfscCount = 23,
    mfccCount = 13,
    featureCount = mfccCount + 1,
    frameLength = params::frameLength,
    shiftInterval = params::shiftInterval,
    fftLength = params::fftLength,
    spectrumLength = fftLength / 2 + 1;

  constexpr static const double mfscStartFrequency = 64.0;

  using feature_vector_t = std::array<double, featureCount>;
  using handle_feature_vector_t = std::function<void(const feature_vector_t &)>;

  MFCCExtractor(handle_feature_vector_t handleFeatureVector)
    : handleFeatureVector(handleFeatureVector),
      last_s_of(0.0),
      last_s_in(0.0),
      remainderSamplesCount(0)
  {
    fillHammingWindow(hammingWindow);
    computeMelBins(cbins);

    fftInput = (double*) fftw_malloc(fftLength * sizeof(double));
    fftOutput = (fftw_complex*) fftw_malloc(spectrumLength * sizeof(fftw_complex));
    fft_plan = fftw_plan_dft_r2c_1d(fftLength, fftInput, fftOutput, FFTW_ESTIMATE);
  }

  ~MFCCExtractor()
  {
    fftw_destroy_plan(fft_plan);
    fftw_free(fftOutput);
    fftw_free(fftInput);
  }

  void processSamples(const double * samples, int count)
  {
    std::vector<double> mutableSamples(remainderSamplesCount + count);

    std::copy(remainderSamples.begin(), remainderSamples.begin() + remainderSamplesCount, mutableSamples.begin());
    std::copy(samples, samples + count, mutableSamples.begin() + remainderSamplesCount);

    offsetCompensation(&mutableSamples[remainderSamplesCount], count);

    remainderSamplesCount = frameSamples(&mutableSamples[0], mutableSamples.size());

    std::copy(mutableSamples.end() - remainderSamplesCount, mutableSamples.end(), remainderSamples.begin());
  }

private:
  const handle_feature_vector_t handleFeatureVector;

  double last_s_of, last_s_in;
  int remainderSamplesCount;

  using frame_t = std::array<double, frameLength>;
  using spectrum_frame_t = std::array<double, spectrumLength>;
  using mfsc_frame_t = std::array<double, mfscCount>;
  using cbins_t = std::array<int, mfscCount + 2>;

  frame_t hammingWindow;
  frame_t remainderSamples, frame;
  spectrum_frame_t powerSpectrum;
  mfsc_frame_t mfscFrame;
  cbins_t cbins;
  feature_vector_t feature_vec;

  fftw_plan fft_plan;
  double * fftInput;
  fftw_complex * fftOutput;

  const double expFloor = 2.0e-22;

  /* 4.2.3 */
  void offsetCompensation(double * samples, int count) {
    for(int i = 0; i < count; i++) {
      last_s_of = samples[i] - last_s_in + 0.999 * last_s_of;
      last_s_in = samples[i];
      samples[i] = last_s_of;
    }
  }

  /* 4.2.4 */
  /* returns count of unprocessed samples at end of array */
  int frameSamples(double * samples, int count) {
    int i = 0;
    for(; i + frameLength <= count; i += shiftInterval) {
      std::copy(samples + i, samples + i + frameLength, frame.begin());
      processFrame(frame);
    }
    return count + shiftInterval - i;
  }

  void processFrame(frame_t & frame) {
    double energy = logE(frame);
    preemphasis(frame);
    window(fftInput, frame);
    fft(powerSpectrum);
    melFilter(mfscFrame, powerSpectrum, cbins);
    dct(feature_vec, mfscFrame);
    feature_vec[featureCount - 1] = energy;
    handleFeatureVector(feature_vec);
  }

  /* 4.2.5 */
  double logE(frame_t & frame) {
    double sum_power = 0.0;
    for(int i = 0; i < frameLength; i++) {
      const double sample = frame[i];
      sum_power += sample * sample;
    }
    return std::log(std::max(sum_power, expFloor));
  }

  /* 4.2.6 */
  void preemphasis(frame_t & frame) {
    double last_s_of = frame[0];
    for(int i = 0; i < frameLength; i++) {
      const double save_s_of = frame[i];
      frame[i] -= 0.97 * last_s_of;
      last_s_of = save_s_of;
    }
  }

  /* 4.2.7 - precomputation */
  void fillHammingWindow(frame_t & frame) {
    for(int i = 0; i < frameLength; i++) {
      frame[i] = 0.54 - 0.46 * std::cos((2.0 * M_PI * (double)i) / (frameLength - 1));
    }
  }

  /* 4.2.7 */
  void window(double * output, const frame_t & frame) {
    for(int i = 0; i < frameLength; i++) {
      output[i] = frame[i] * hammingWindow[i];
    }
  }

  /* 4.2.8 */
  void fft(spectrum_frame_t & powerSpectrum) {
    for(int i = frameLength; i < fftLength; i++) {
      fftInput[i] = 0.0;
    }
    fftw_execute(fft_plan);
    for(int i = 0; i < spectrumLength; i++) {
      const fftw_complex & spec = fftOutput[i];
      powerSpectrum[i] = std::sqrt(spec[0] * spec[0] + spec[1] * spec[1]);
    }
  }

  /* 4.2.9 - precomputation */
  static double mel(double f) {
    return 2595.0 * std::log(1.0 + f / 700.0) / std::log(10.0);
  }
  static double mel_inv(double mel) {
    return 700.0 * (std::exp(mel / 2595.0) - 1.0);
  }
  void computeMelBins(cbins_t & cbins) {
    const double mel_start = mel(mfscStartFrequency);
    const double mel_nyquist = mel(sampleRate / 2);
    for(int i = 0; i < mfscCount + 2; i++) {
      double center_freq = mel_inv(mel_start + (mel_nyquist - mel_start) * (double)i / (double)(mfscCount + 1));
      cbins[i] = std::round(center_freq * (double)fftLength / (double)sampleRate);
    }
  }

  /* 4.2.9 */
  void melFilter(mfsc_frame_t & mfscFrame, const spectrum_frame_t & powerSpectrum, const cbins_t & cbins) {
    for(int i = 1; i < mfscCount + 1; i++) {
      double fbank = 0.0;
      for(int j = cbins[i-1]; j <= cbins[i]; j++) {
        double coef = (double)(j - cbins[i-1] + 1) / (double)(cbins[i] - cbins[i-1] + 1);
        fbank += coef * powerSpectrum[j];
      }
      for(int j = cbins[i] + 1; j <= cbins[i + 1]; j++) {
        double coef = 1.0 - (double)(j - cbins[i]) / (double)(cbins[i+1] - cbins[i] + 1);
        fbank += coef * powerSpectrum[j];
      }
      /* 4.2.10 */
      mfscFrame[i-1] = std::log(std::max(fbank, expFloor));
    }
  }

  /* 4.2.11 */
  void dct(feature_vector_t & output, mfsc_frame_t & input) {
    for(int i = 0; i < mfccCount; i++) {
      double C_i = 0.0;
      for(int j = 0; j < mfscCount; j++) {
        C_i += input[j] * std::cos(M_PI * (double)i * ((double)j + 0.5) / mfscCount);
      }
      output[i] = C_i;
    }
  }
};
