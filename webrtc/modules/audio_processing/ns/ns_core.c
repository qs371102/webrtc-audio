/*
 *  Copyright (c) 2012 The WebRTC project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "webrtc/base/checks.h"
#include "webrtc/common_audio/fft4g.h"
#include "webrtc/common_audio/signal_processing/include/signal_processing_library.h"
#include "webrtc/modules/audio_processing/ns/noise_suppression.h"
#include "webrtc/modules/audio_processing/ns/ns_core.h"
#include "webrtc/modules/audio_processing/ns/windows_private.h"

// Set Feature Extraction Parameters.
static void set_feature_extraction_parameters(NoiseSuppressionC* nssc) {
  // Bin size of histogram.
  nssc->featureExtractionParams.binSizeLrt = 0.1f;
  nssc->featureExtractionParams.binSizeSpecFlat = 0.05f;
  nssc->featureExtractionParams.binSizeSpecDiff = 0.1f;

  // Range of histogram over which LRT threshold is computed.
  nssc->featureExtractionParams.rangeAvgHistLrt = 1.f;

  // Scale parameters: multiply dominant peaks of the histograms by scale factor
  // to obtain thresholds for prior model.
  // For LRT and spectral difference.
  nssc->featureExtractionParams.factor1ModelPars = 1.2f;
  // For spectral_flatness: used when noise is flatter than speech.
  nssc->featureExtractionParams.factor2ModelPars = 0.9f;

  // Peak limit for spectral flatness (varies between 0 and 1).
  nssc->featureExtractionParams.thresPosSpecFlat = 0.6f;

  // Limit on spacing of two highest peaks in histogram: spacing determined by
  // bin size.
  nssc->featureExtractionParams.limitPeakSpacingSpecFlat =
      2 * nssc->featureExtractionParams.binSizeSpecFlat;
  nssc->featureExtractionParams.limitPeakSpacingSpecDiff =
      2 * nssc->featureExtractionParams.binSizeSpecDiff;

  // Limit on relevance of second peak.
  nssc->featureExtractionParams.limitPeakWeightsSpecFlat = 0.5f;
  nssc->featureExtractionParams.limitPeakWeightsSpecDiff = 0.5f;

  // Fluctuation limit of LRT feature.
  nssc->featureExtractionParams.thresFluctLrt = 0.05f;

  // Limit on the max and min values for the feature thresholds.
  nssc->featureExtractionParams.maxLrt = 1.f;
  nssc->featureExtractionParams.minLrt = 0.2f;

  nssc->featureExtractionParams.maxSpecFlat = 0.95f;
  nssc->featureExtractionParams.minSpecFlat = 0.1f;

  nssc->featureExtractionParams.maxSpecDiff = 1.f;
  nssc->featureExtractionParams.minSpecDiff = 0.16f;

  // Criteria of weight of histogram peak to accept/reject feature.
  nssc->featureExtractionParams.thresWeightSpecFlat =
      (int)(0.3 * (nssc->modelUpdatePars[1]));  // For spectral flatness.
  nssc->featureExtractionParams.thresWeightSpecDiff =
      (int)(0.3 * (nssc->modelUpdatePars[1]));  // For spectral difference.
}

// Initialize state.
int WebRtcNs_InitCore(NoiseSuppressionC* nssc, uint32_t fs) {
  int i;
  // Check for valid pointer.
  if (nssc == NULL) {
    return -1;
  }

  // Initialization of struct.
  if (fs == 8000 || fs == 16000 || fs == 32000 || fs == 48000) {
    nssc->fs = fs;
  } else {
    return -1;
  }
  nssc->windShift = 0;
  // We only support 10ms frames.
  if (fs == 8000) {
    nssc->blockLen = 80;
    nssc->anaLen = 128;
    nssc->window = kBlocks80w128;
  } else {
    nssc->blockLen = 160;
    nssc->anaLen = 256;
    nssc->window = kBlocks160w256;
  }
  nssc->magnLen = nssc->anaLen / 2 + 1;  // Number of frequency bins.

  // Initialize FFT work arrays.
  nssc->ip[0] = 0;  // Setting this triggers initialization.
  memset(nssc->dataBuf, 0, sizeof(float) * ANAL_BLOCKL_MAX);
  WebRtc_rdft(nssc->anaLen, 1, nssc->dataBuf, nssc->ip, nssc->wfft);

  memset(nssc->analyzeBuf, 0, sizeof(float) * ANAL_BLOCKL_MAX);
  memset(nssc->dataBuf, 0, sizeof(float) * ANAL_BLOCKL_MAX);
  memset(nssc->syntBuf, 0, sizeof(float) * ANAL_BLOCKL_MAX);

  // For HB processing.
  memset(nssc->dataBufHB,
         0,
         sizeof(float) * NUM_HIGH_BANDS_MAX * ANAL_BLOCKL_MAX);

  // For quantile noise estimation.
  memset(nssc->quantile, 0, sizeof(float) * HALF_ANAL_BLOCKL);
  for (i = 0; i < SIMULT * HALF_ANAL_BLOCKL; i++) {
    nssc->lquantile[i] = 8.f;
    nssc->density[i] = 0.3f;
  }

  for (i = 0; i < SIMULT; i++) {
    nssc->counter[i] =
        (int)floor((float)(END_STARTUP_LONG * (i + 1)) / (float)SIMULT);
  }

  nssc->updates = 0;

  // Wiener filter initialization.
  for (i = 0; i < HALF_ANAL_BLOCKL; i++) {
    nssc->smooth[i] = 1.f;
  }

  // Set the aggressiveness: default.
  nssc->aggrMode = 0;

  // Initialize variables for new method.
  nssc->priorSpeechProb = 0.5f;  // Prior prob for speech/noise.
  // Previous analyze mag spectrum.
  memset(nssc->magnPrevAnalyze, 0, sizeof(float) * HALF_ANAL_BLOCKL);
  // Previous process mag spectrum.
  memset(nssc->magnPrevProcess, 0, sizeof(float) * HALF_ANAL_BLOCKL);
  // Current noise-spectrum.
  memset(nssc->noise, 0, sizeof(float) * HALF_ANAL_BLOCKL);
  // Previous noise-spectrum.
  memset(nssc->noisePrev, 0, sizeof(float) * HALF_ANAL_BLOCKL);
  // Conservative noise spectrum estimate.
  memset(nssc->magnAvgPause, 0, sizeof(float) * HALF_ANAL_BLOCKL);
  // For estimation of HB in second pass.
  memset(nssc->speechProb, 0, sizeof(float) * HALF_ANAL_BLOCKL);
  // Initial average magnitude spectrum.
  memset(nssc->initMagnEst, 0, sizeof(float) * HALF_ANAL_BLOCKL);
  for (i = 0; i < HALF_ANAL_BLOCKL; i++) {
    // Smooth LR (same as threshold).
    nssc->logLrtTimeAvg[i] = LRT_FEATURE_THR;
  }

  // Feature quantities.
  // Spectral flatness (start on threshold).
  nssc->featureData[0] = SF_FEATURE_THR;
  nssc->featureData[1] = 0.f;  // Spectral entropy: not used in this version.
  nssc->featureData[2] = 0.f;  // Spectral variance: not used in this version.
  // Average LRT factor (start on threshold).
  nssc->featureData[3] = LRT_FEATURE_THR;
  // Spectral template diff (start on threshold).
  nssc->featureData[4] = SF_FEATURE_THR;
  nssc->featureData[5] = 0.f;  // Normalization for spectral difference.
  // Window time-average of input magnitude spectrum.
  nssc->featureData[6] = 0.f;

  // Histogram quantities: used to estimate/update thresholds for features.
  memset(nssc->histLrt, 0, sizeof(int) * HIST_PAR_EST);
  memset(nssc->histSpecFlat, 0, sizeof(int) * HIST_PAR_EST);
  memset(nssc->histSpecDiff, 0, sizeof(int) * HIST_PAR_EST);


  nssc->blockInd = -1;  // Frame counter.
  // Default threshold for LRT feature.
  nssc->priorModelPars[0] = LRT_FEATURE_THR;
  // Threshold for spectral flatness: determined on-line.
  nssc->priorModelPars[1] = 0.5f;
  // sgn_map par for spectral measure: 1 for flatness measure.
  nssc->priorModelPars[2] = 1.f;
  // Threshold for template-difference feature: determined on-line.
  nssc->priorModelPars[3] = 0.5f;
  // Default weighting parameter for LRT feature.
  nssc->priorModelPars[4] = 1.f;
  // Default weighting parameter for spectral flatness feature.
  nssc->priorModelPars[5] = 0.f;
  // Default weighting parameter for spectral difference feature.
  nssc->priorModelPars[6] = 0.f;

  // Update flag for parameters:
  // 0 no update, 1 = update once, 2 = update every window.
  nssc->modelUpdatePars[0] = 2;
  nssc->modelUpdatePars[1] = 500;  // Window for update.
  // Counter for update of conservative noise spectrum.
  nssc->modelUpdatePars[2] = 0;
  // Counter if the feature thresholds are updated during the sequence.
  nssc->modelUpdatePars[3] = nssc->modelUpdatePars[1];

  nssc->signalEnergy = 0.0;
  nssc->sumMagn = 0.0;
  nssc->whiteNoiseLevel = 0.0;
  nssc->pinkNoiseNumerator = 0.0;
  nssc->pinkNoiseExp = 0.0;

  set_feature_extraction_parameters(nssc);

  // Default mode.
  WebRtcNs_set_policy_core(nssc, 0);

  nssc->initFlag = 1;
  return 0;
}

// Estimate noise.
static void NoiseEstimation(NoiseSuppressionC* nssc,
                            float* magn,
                            float* noise) {
  size_t i, s, offset;
  float lmagn[HALF_ANAL_BLOCKL], delta;

  if (nssc->updates < END_STARTUP_LONG) {
    nssc->updates++;
  }

  for (i = 0; i < nssc->magnLen; i++) {
    lmagn[i] = (float)log(magn[i]);
  }

  // Loop over simultaneous estimates.
  for (s = 0; s < SIMULT; s++) {
    offset = s * nssc->magnLen;

    // newquantest(...)
    for (i = 0; i < nssc->magnLen; i++) {
      // Compute delta.
      if (nssc->density[offset + i] > 1.0) {
        delta = FACTOR * 1.f / nssc->density[offset + i];
      } else {
        delta = FACTOR;
      }

      // Update log quantile estimate.
      if (lmagn[i] > nssc->lquantile[offset + i]) {
        nssc->lquantile[offset + i] +=
            QUANTILE * delta / (float)(nssc->counter[s] + 1);
      } else {
        nssc->lquantile[offset + i] -=
            (1.f - QUANTILE) * delta / (float)(nssc->counter[s] + 1);
      }

      // Update density estimate.
      if (fabs(lmagn[i] - nssc->lquantile[offset + i]) < WIDTH) {
        nssc->density[offset + i] =
            ((float)nssc->counter[s] * nssc->density[offset + i] +
             1.f / (2.f * WIDTH)) /
            (float)(nssc->counter[s] + 1);
      }
    }  // End loop over magnitude spectrum.

    if (nssc->counter[s] >= END_STARTUP_LONG) {
      nssc->counter[s] = 0;
      if (nssc->updates >= END_STARTUP_LONG) {
        for (i = 0; i < nssc->magnLen; i++) {
          nssc->quantile[i] = (float)exp(nssc->lquantile[offset + i]);
        }
      }
    }

    nssc->counter[s]++;
  }  // End loop over simultaneous estimates.

  // Sequentially update the noise during startup.
  if (nssc->updates < END_STARTUP_LONG) {
    // Use the last "s" to get noise during startup that differ from zero.
    for (i = 0; i < nssc->magnLen; i++) {
      nssc->quantile[i] = (float)exp(nssc->lquantile[offset + i]);
    }
  }

  for (i = 0; i < nssc->magnLen; i++) {
    noise[i] = nssc->quantile[i];
  }
}

// Extract thresholds for feature parameters.
// Histograms are computed over some window size (given by
// nssc->modelUpdatePars[1]).
// Thresholds and weights are extracted every window.
// |flag| = 0 updates histogram only, |flag| = 1 computes the threshold/weights.
// Threshold and weights are returned in: nssc->priorModelPars.
static void FeatureParameterExtraction(NoiseSuppressionC* nssc, int flag) {
  int i, useFeatureSpecFlat, useFeatureSpecDiff, numHistLrt;
  int maxPeak1, maxPeak2;
  int weightPeak1SpecFlat, weightPeak2SpecFlat, weightPeak1SpecDiff,
      weightPeak2SpecDiff;

  float binMid, featureSum;
  float posPeak1SpecFlat, posPeak2SpecFlat, posPeak1SpecDiff, posPeak2SpecDiff;
  float fluctLrt, avgHistLrt, avgSquareHistLrt, avgHistLrtCompl;

  // 3 features: LRT, flatness, difference.
  // lrt_feature = nssc->featureData[3];
  // flat_feature = nssc->featureData[0];
  // diff_feature = nssc->featureData[4];

  // Update histograms.
  if (flag == 0) {
    // LRT
    if ((nssc->featureData[3] <
         HIST_PAR_EST * nssc->featureExtractionParams.binSizeLrt) &&
        (nssc->featureData[3] >= 0.0)) {
      i = (int)(nssc->featureData[3] /
                nssc->featureExtractionParams.binSizeLrt);
      nssc->histLrt[i]++;
    }
    // Spectral flatness.
    if ((nssc->featureData[0] <
         HIST_PAR_EST * nssc->featureExtractionParams.binSizeSpecFlat) &&
        (nssc->featureData[0] >= 0.0)) {
      i = (int)(nssc->featureData[0] /
                nssc->featureExtractionParams.binSizeSpecFlat);
      nssc->histSpecFlat[i]++;
    }
    // Spectral difference.
    if ((nssc->featureData[4] <
         HIST_PAR_EST * nssc->featureExtractionParams.binSizeSpecDiff) &&
        (nssc->featureData[4] >= 0.0)) {
      i = (int)(nssc->featureData[4] /
                nssc->featureExtractionParams.binSizeSpecDiff);
      nssc->histSpecDiff[i]++;
    }
  }

  // Extract parameters for speech/noise probability.
  if (flag == 1) {
    // LRT feature: compute the average over
    // nssc->featureExtractionParams.rangeAvgHistLrt.
    avgHistLrt = 0.0;
    avgHistLrtCompl = 0.0;
    avgSquareHistLrt = 0.0;
    numHistLrt = 0;
    for (i = 0; i < HIST_PAR_EST; i++) {
      binMid = ((float)i + 0.5f) * nssc->featureExtractionParams.binSizeLrt;
      if (binMid <= nssc->featureExtractionParams.rangeAvgHistLrt) {
        avgHistLrt += nssc->histLrt[i] * binMid;
        numHistLrt += nssc->histLrt[i];
      }
      avgSquareHistLrt += nssc->histLrt[i] * binMid * binMid;
      avgHistLrtCompl += nssc->histLrt[i] * binMid;
    }
    if (numHistLrt > 0) {
      avgHistLrt = avgHistLrt / ((float)numHistLrt);
    }
    avgHistLrtCompl = avgHistLrtCompl / ((float)nssc->modelUpdatePars[1]);
    avgSquareHistLrt = avgSquareHistLrt / ((float)nssc->modelUpdatePars[1]);
    fluctLrt = avgSquareHistLrt - avgHistLrt * avgHistLrtCompl;
    // Get threshold for LRT feature.
    if (fluctLrt < nssc->featureExtractionParams.thresFluctLrt) {
      // Very low fluctuation, so likely noise.
      nssc->priorModelPars[0] = nssc->featureExtractionParams.maxLrt;
    } else {
      nssc->priorModelPars[0] =
          nssc->featureExtractionParams.factor1ModelPars * avgHistLrt;
      // Check if value is within min/max range.
      if (nssc->priorModelPars[0] < nssc->featureExtractionParams.minLrt) {
        nssc->priorModelPars[0] = nssc->featureExtractionParams.minLrt;
      }
      if (nssc->priorModelPars[0] > nssc->featureExtractionParams.maxLrt) {
        nssc->priorModelPars[0] = nssc->featureExtractionParams.maxLrt;
      }
    }
    // Done with LRT feature.

    // For spectral flatness and spectral difference: compute the main peaks of
    // histogram.
    maxPeak1 = 0;
    maxPeak2 = 0;
    posPeak1SpecFlat = 0.0;
    posPeak2SpecFlat = 0.0;
    weightPeak1SpecFlat = 0;
    weightPeak2SpecFlat = 0;

    // Peaks for flatness.
    for (i = 0; i < HIST_PAR_EST; i++) {
      binMid =
          (i + 0.5f) * nssc->featureExtractionParams.binSizeSpecFlat;
      if (nssc->histSpecFlat[i] > maxPeak1) {
        // Found new "first" peak.
        maxPeak2 = maxPeak1;
        weightPeak2SpecFlat = weightPeak1SpecFlat;
        posPeak2SpecFlat = posPeak1SpecFlat;

        maxPeak1 = nssc->histSpecFlat[i];
        weightPeak1SpecFlat = nssc->histSpecFlat[i];
        posPeak1SpecFlat = binMid;
      } else if (nssc->histSpecFlat[i] > maxPeak2) {
        // Found new "second" peak.
        maxPeak2 = nssc->histSpecFlat[i];
        weightPeak2SpecFlat = nssc->histSpecFlat[i];
        posPeak2SpecFlat = binMid;
      }
    }

    // Compute two peaks for spectral difference.
    maxPeak1 = 0;
    maxPeak2 = 0;
    posPeak1SpecDiff = 0.0;
    posPeak2SpecDiff = 0.0;
    weightPeak1SpecDiff = 0;
    weightPeak2SpecDiff = 0;
    // Peaks for spectral difference.
    for (i = 0; i < HIST_PAR_EST; i++) {
      binMid =
          ((float)i + 0.5f) * nssc->featureExtractionParams.binSizeSpecDiff;
      if (nssc->histSpecDiff[i] > maxPeak1) {
        // Found new "first" peak.
        maxPeak2 = maxPeak1;
        weightPeak2SpecDiff = weightPeak1SpecDiff;
        posPeak2SpecDiff = posPeak1SpecDiff;

        maxPeak1 = nssc->histSpecDiff[i];
        weightPeak1SpecDiff = nssc->histSpecDiff[i];
        posPeak1SpecDiff = binMid;
      } else if (nssc->histSpecDiff[i] > maxPeak2) {
        // Found new "second" peak.
        maxPeak2 = nssc->histSpecDiff[i];
        weightPeak2SpecDiff = nssc->histSpecDiff[i];
        posPeak2SpecDiff = binMid;
      }
    }

    // For spectrum flatness feature.
    useFeatureSpecFlat = 1;
    // Merge the two peaks if they are close.
    if ((fabs(posPeak2SpecFlat - posPeak1SpecFlat) <
         nssc->featureExtractionParams.limitPeakSpacingSpecFlat) &&
        (weightPeak2SpecFlat >
         nssc->featureExtractionParams.limitPeakWeightsSpecFlat *
             weightPeak1SpecFlat)) {
      weightPeak1SpecFlat += weightPeak2SpecFlat;
      posPeak1SpecFlat = 0.5f * (posPeak1SpecFlat + posPeak2SpecFlat);
    }
    // Reject if weight of peaks is not large enough, or peak value too small.
    if (weightPeak1SpecFlat <
            nssc->featureExtractionParams.thresWeightSpecFlat ||
        posPeak1SpecFlat < nssc->featureExtractionParams.thresPosSpecFlat) {
      useFeatureSpecFlat = 0;
    }
    // If selected, get the threshold.
    if (useFeatureSpecFlat == 1) {
      // Compute the threshold.
      nssc->priorModelPars[1] =
          nssc->featureExtractionParams.factor2ModelPars * posPeak1SpecFlat;
      // Check if value is within min/max range.
      if (nssc->priorModelPars[1] < nssc->featureExtractionParams.minSpecFlat) {
        nssc->priorModelPars[1] = nssc->featureExtractionParams.minSpecFlat;
      }
      if (nssc->priorModelPars[1] > nssc->featureExtractionParams.maxSpecFlat) {
        nssc->priorModelPars[1] = nssc->featureExtractionParams.maxSpecFlat;
      }
    }
    // Done with flatness feature.

    // For template feature.
    useFeatureSpecDiff = 1;
    // Merge the two peaks if they are close.
    if ((fabs(posPeak2SpecDiff - posPeak1SpecDiff) <
         nssc->featureExtractionParams.limitPeakSpacingSpecDiff) &&
        (weightPeak2SpecDiff >
         nssc->featureExtractionParams.limitPeakWeightsSpecDiff *
             weightPeak1SpecDiff)) {
      weightPeak1SpecDiff += weightPeak2SpecDiff;
      posPeak1SpecDiff = 0.5f * (posPeak1SpecDiff + posPeak2SpecDiff);
    }
    // Get the threshold value.
    nssc->priorModelPars[3] =
        nssc->featureExtractionParams.factor1ModelPars * posPeak1SpecDiff;
    // Reject if weight of peaks is not large enough.
    if (weightPeak1SpecDiff <
        nssc->featureExtractionParams.thresWeightSpecDiff) {
      useFeatureSpecDiff = 0;
    }
    // Check if value is within min/max range.
    if (nssc->priorModelPars[3] < nssc->featureExtractionParams.minSpecDiff) {
      nssc->priorModelPars[3] = nssc->featureExtractionParams.minSpecDiff;
    }
    if (nssc->priorModelPars[3] > nssc->featureExtractionParams.maxSpecDiff) {
      nssc->priorModelPars[3] = nssc->featureExtractionParams.maxSpecDiff;
    }
    // Done with spectral difference feature.

    // Don't use template feature if fluctuation of LRT feature is very low:
    // most likely just noise state.
    if (fluctLrt < nssc->featureExtractionParams.thresFluctLrt) {
      useFeatureSpecDiff = 0;
    }

    // Select the weights between the features.
    // nssc->priorModelPars[4] is weight for LRT: always selected.
    // nssc->priorModelPars[5] is weight for spectral flatness.
    // nssc->priorModelPars[6] is weight for spectral difference.
    featureSum = (float)(1 + useFeatureSpecFlat + useFeatureSpecDiff);
    nssc->priorModelPars[4] = 1.f / featureSum;
    nssc->priorModelPars[5] = ((float)useFeatureSpecFlat) / featureSum;
    nssc->priorModelPars[6] = ((float)useFeatureSpecDiff) / featureSum;

    // Set hists to zero for next update.
    if (nssc->modelUpdatePars[0] >= 1) {
      for (i = 0; i < HIST_PAR_EST; i++) {
        nssc->histLrt[i] = 0;
        nssc->histSpecFlat[i] = 0;
        nssc->histSpecDiff[i] = 0;
      }
    }
  }  // End of flag == 1.
}

// Compute spectral flatness on input spectrum.
// |magnIn| is the magnitude spectrum.
// Spectral flatness is returned in nssc->featureData[0].
static void ComputeSpectralFlatness(NoiseSuppressionC* nssc,
                                    const float* magnIn) {
  size_t i;
  size_t shiftLP = 1;  // Option to remove first bin(s) from spectral measures.
  float avgSpectralFlatnessNum, avgSpectralFlatnessDen, spectralTmp;

  // Compute spectral measures.
  // For flatness.
  avgSpectralFlatnessNum = 0.0;
  avgSpectralFlatnessDen = nssc->sumMagn;
  for (i = 0; i < shiftLP; i++) {
    avgSpectralFlatnessDen -= magnIn[i];
  }
  // Compute log of ratio of the geometric to arithmetic mean: check for log(0)
  // case.
  for (i = shiftLP; i < nssc->magnLen; i++) {
    if (magnIn[i] > 0.0) {
      avgSpectralFlatnessNum += (float)log(magnIn[i]);
    } else {
      nssc->featureData[0] -= SPECT_FL_TAVG * nssc->featureData[0];
      return;
    }
  }
  // Normalize.
  avgSpectralFlatnessDen = avgSpectralFlatnessDen / nssc->magnLen;
  avgSpectralFlatnessNum = avgSpectralFlatnessNum / nssc->magnLen;

  // Ratio and inverse log: check for case of log(0).
  spectralTmp = (float)exp(avgSpectralFlatnessNum) / avgSpectralFlatnessDen;

  // Time-avg update of spectral flatness feature.
  nssc->featureData[0] += SPECT_FL_TAVG * (spectralTmp - nssc->featureData[0]);
  // Done with flatness feature.
}

// Compute prior and post SNR based on quantile noise estimation.
// Compute DD estimate of prior SNR.
// Inputs:
//   * |magn| is the signal magnitude spectrum estimate.
//   * |noise| is the magnitude noise spectrum estimate.
// Outputs:
//   * |snrLocPrior| is the computed prior SNR.
//   * |snrLocPost| is the computed post SNR.
static void ComputeSnr(const NoiseSuppressionC* nssc,
                       const float* magn,
                       const float* noise,
                       float* snrLocPrior,
                       float* snrLocPost) {
  size_t i;

  for (i = 0; i < nssc->magnLen; i++) {
    // Previous post SNR.
    // Previous estimate: based on previous frame with gain filter.
    float previousEstimateStsa = nssc->magnPrevAnalyze[i] /
        (nssc->noisePrev[i] + 0.0001f) * nssc->smooth[i];
    // Post SNR.
    snrLocPost[i] = 0.f;
    if (magn[i] > noise[i]) {
      snrLocPost[i] = magn[i] / (noise[i] + 0.0001f) - 1.f;
    }
    // DD estimate is sum of two terms: current estimate and previous estimate.
    // Directed decision update of snrPrior.
    snrLocPrior[i] =
        DD_PR_SNR * previousEstimateStsa + (1.f - DD_PR_SNR) * snrLocPost[i];
  }  // End of loop over frequencies.
}

// Compute the difference measure between input spectrum and a template/learned
// noise spectrum.
// |magnIn| is the input spectrum.
// The reference/template spectrum is nssc->magnAvgPause[i].
// Returns (normalized) spectral difference in nssc->featureData[4].
static void ComputeSpectralDifference(NoiseSuppressionC* nssc,
                                      const float* magnIn) {
  // avgDiffNormMagn = var(magnIn) - cov(magnIn, magnAvgPause)^2 /
  // var(magnAvgPause)
  size_t i;
  float avgPause, avgMagn, covMagnPause, varPause, varMagn, avgDiffNormMagn;

  avgPause = 0.0;
  avgMagn = nssc->sumMagn;
  // Compute average quantities.
  for (i = 0; i < nssc->magnLen; i++) {
    // Conservative smooth noise spectrum from pause frames.
    avgPause += nssc->magnAvgPause[i];
  }
  avgPause /= nssc->magnLen;
  avgMagn /= nssc->magnLen;

  covMagnPause = 0.0;
  varPause = 0.0;
  varMagn = 0.0;
  // Compute variance and covariance quantities.
  for (i = 0; i < nssc->magnLen; i++) {
    covMagnPause += (magnIn[i] - avgMagn) * (nssc->magnAvgPause[i] - avgPause);
    varPause +=
        (nssc->magnAvgPause[i] - avgPause) * (nssc->magnAvgPause[i] - avgPause);
    varMagn += (magnIn[i] - avgMagn) * (magnIn[i] - avgMagn);
  }
  covMagnPause /= nssc->magnLen;
  varPause /= nssc->magnLen;
  varMagn /= nssc->magnLen;
  // Update of average magnitude spectrum.
  nssc->featureData[6] += nssc->signalEnergy;

  avgDiffNormMagn =
      varMagn - (covMagnPause * covMagnPause) / (varPause + 0.0001f);
  // Normalize and compute time-avg update of difference feature.
  avgDiffNormMagn = (float)(avgDiffNormMagn / (nssc->featureData[5] + 0.0001f));
  nssc->featureData[4] +=
      SPECT_DIFF_TAVG * (avgDiffNormMagn - nssc->featureData[4]);
}

// Compute speech/noise probability.
// Speech/noise probability is returned in |probSpeechFinal|.
// |magn| is the input magnitude spectrum.
// |noise| is the noise spectrum.
// |snrLocPrior| is the prior SNR for each frequency.
// |snrLocPost| is the post SNR for each frequency.
static void SpeechNoiseProb(NoiseSuppressionC* nssc,
                            float* probSpeechFinal,
                            const float* snrLocPrior,
                            const float* snrLocPost) {
  size_t i;
  int sgnMap;
  float invLrt, gainPrior, indPrior;
  float logLrtTimeAvgKsum, besselTmp;
  float indicator0, indicator1, indicator2;
  float tmpFloat1, tmpFloat2;
  float weightIndPrior0, weightIndPrior1, weightIndPrior2;
  float threshPrior0, threshPrior1, threshPrior2;
  float widthPrior, widthPrior0, widthPrior1, widthPrior2;

  widthPrior0 = WIDTH_PR_MAP;
  // Width for pause region: lower range, so increase width in tanh map.
  widthPrior1 = 2.f * WIDTH_PR_MAP;
  widthPrior2 = 2.f * WIDTH_PR_MAP;  // For spectral-difference measure.

  // Threshold parameters for features.
  threshPrior0 = nssc->priorModelPars[0];
  threshPrior1 = nssc->priorModelPars[1];
  threshPrior2 = nssc->priorModelPars[3];

  // Sign for flatness feature.
  sgnMap = (int)(nssc->priorModelPars[2]);

  // Weight parameters for features.
  weightIndPrior0 = nssc->priorModelPars[4];
  weightIndPrior1 = nssc->priorModelPars[5];
  weightIndPrior2 = nssc->priorModelPars[6];

  // Compute feature based on average LR factor.
  // This is the average over all frequencies of the smooth log LRT.
  logLrtTimeAvgKsum = 0.0;
  for (i = 0; i < nssc->magnLen; i++) {
    tmpFloat1 = 1.f + 2.f * snrLocPrior[i];
    tmpFloat2 = 2.f * snrLocPrior[i] / (tmpFloat1 + 0.0001f);
    besselTmp = (snrLocPost[i] + 1.f) * tmpFloat2;
    nssc->logLrtTimeAvg[i] +=
        LRT_TAVG * (besselTmp - (float)log(tmpFloat1) - nssc->logLrtTimeAvg[i]);
    logLrtTimeAvgKsum += nssc->logLrtTimeAvg[i];
  }
  logLrtTimeAvgKsum = (float)logLrtTimeAvgKsum / (nssc->magnLen);
  nssc->featureData[3] = logLrtTimeAvgKsum;
  // Done with computation of LR factor.

  // Compute the indicator functions.
  // Average LRT feature.
  widthPrior = widthPrior0;
  // Use larger width in tanh map for pause regions.
  if (logLrtTimeAvgKsum < threshPrior0) {
    widthPrior = widthPrior1;
  }
  // Compute indicator function: sigmoid map.
  indicator0 =
      0.5f *
      ((float)tanh(widthPrior * (logLrtTimeAvgKsum - threshPrior0)) + 1.f);

  // Spectral flatness feature.
  tmpFloat1 = nssc->featureData[0];
  widthPrior = widthPrior0;
  // Use larger width in tanh map for pause regions.
  if (sgnMap == 1 && (tmpFloat1 > threshPrior1)) {
    widthPrior = widthPrior1;
  }
  if (sgnMap == -1 && (tmpFloat1 < threshPrior1)) {
    widthPrior = widthPrior1;
  }
  // Compute indicator function: sigmoid map.
  indicator1 =
      0.5f *
      ((float)tanh((float)sgnMap * widthPrior * (threshPrior1 - tmpFloat1)) +
       1.f);

  // For template spectrum-difference.
  tmpFloat1 = nssc->featureData[4];
  widthPrior = widthPrior0;
  // Use larger width in tanh map for pause regions.
  if (tmpFloat1 < threshPrior2) {
    widthPrior = widthPrior2;
  }
  // Compute indicator function: sigmoid map.
  indicator2 =
      0.5f * ((float)tanh(widthPrior * (tmpFloat1 - threshPrior2)) + 1.f);

  // Combine the indicator function with the feature weights.
  indPrior = weightIndPrior0 * indicator0 + weightIndPrior1 * indicator1 +
             weightIndPrior2 * indicator2;
  // Done with computing indicator function.

  // Compute the prior probability.
  nssc->priorSpeechProb += PRIOR_UPDATE * (indPrior - nssc->priorSpeechProb);
  // Make sure probabilities are within range: keep floor to 0.01.
  if (nssc->priorSpeechProb > 1.f) {
    nssc->priorSpeechProb = 1.f;
  }
  if (nssc->priorSpeechProb < 0.01f) {
    nssc->priorSpeechProb = 0.01f;
  }

  // Final speech probability: combine prior model with LR factor:.
  gainPrior = (1.f - nssc->priorSpeechProb) / (nssc->priorSpeechProb + 0.0001f);
  for (i = 0; i < nssc->magnLen; i++) {
    invLrt = (float)exp(-nssc->logLrtTimeAvg[i]);
    invLrt = (float)gainPrior * invLrt;
    probSpeechFinal[i] = 1.f / (1.f + invLrt);
  }
}

// Update the noise features.
// Inputs:
//   * |magn| is the signal magnitude spectrum estimate.
//   * |updateParsFlag| is an update flag for parameters.
static void FeatureUpdate(NoiseSuppressionC* nssc,
                          const float* magn,
                          int updateParsFlag) {
  // Compute spectral flatness on input spectrum.
  ComputeSpectralFlatness(nssc, magn);
  // Compute difference of input spectrum with learned/estimated noise spectrum.
  ComputeSpectralDifference(nssc, magn);
  // Compute histograms for parameter decisions (thresholds and weights for
  // features).
  // Parameters are extracted once every window time.
  // (=nssc->modelUpdatePars[1])
  if (updateParsFlag >= 1) {
    // Counter update.
    nssc->modelUpdatePars[3]--;
    // Update histogram.
    if (nssc->modelUpdatePars[3] > 0) {
      FeatureParameterExtraction(nssc, 0);
    }
    // Compute model parameters.
    if (nssc->modelUpdatePars[3] == 0) {
      FeatureParameterExtraction(nssc, 1);
      nssc->modelUpdatePars[3] = nssc->modelUpdatePars[1];
      // If wish to update only once, set flag to zero.
      if (updateParsFlag == 1) {
        nssc->modelUpdatePars[0] = 0;
      } else {
        // Update every window:
        // Get normalization for spectral difference for next window estimate.
        nssc->featureData[6] =
            nssc->featureData[6] / ((float)nssc->modelUpdatePars[1]);
        nssc->featureData[5] =
            0.5f * (nssc->featureData[6] + nssc->featureData[5]);
        nssc->featureData[6] = 0.f;
      }
    }
  }
}

// Update the noise estimate.
// Inputs:
//   * |magn| is the signal magnitude spectrum estimate.
//   * |snrLocPrior| is the prior SNR.
//   * |snrLocPost| is the post SNR.
// Output:
//   * |noise| is the updated noise magnitude spectrum estimate.
static void UpdateNoiseEstimate(NoiseSuppressionC* nssc,
                                const float* magn,
                                const float* snrLocPrior,
                                const float* snrLocPost,
                                float* noise) {
  size_t i;
  float probSpeech, probNonSpeech;
  // Time-avg parameter for noise update.
  float gammaNoiseTmp = NOISE_UPDATE;
  float gammaNoiseOld;
  float noiseUpdateTmp;

  for (i = 0; i < nssc->magnLen; i++) {
    probSpeech = nssc->speechProb[i];
    probNonSpeech = 1.f - probSpeech;
    // Temporary noise update:
    // Use it for speech frames if update value is less than previous.
    noiseUpdateTmp = gammaNoiseTmp * nssc->noisePrev[i] +
                     (1.f - gammaNoiseTmp) * (probNonSpeech * magn[i] +
                                              probSpeech * nssc->noisePrev[i]);
    // Time-constant based on speech/noise state.
    gammaNoiseOld = gammaNoiseTmp;
    gammaNoiseTmp = NOISE_UPDATE;
    // Increase gamma (i.e., less noise update) for frame likely to be speech.
    if (probSpeech > PROB_RANGE) {
      gammaNoiseTmp = SPEECH_UPDATE;
    }
    // Conservative noise update.
    if (probSpeech < PROB_RANGE) {
      nssc->magnAvgPause[i] += GAMMA_PAUSE * (magn[i] - nssc->magnAvgPause[i]);
    }
    // Noise update.
    if (gammaNoiseTmp == gammaNoiseOld) {
      noise[i] = noiseUpdateTmp;
    } else {
      noise[i] = gammaNoiseTmp * nssc->noisePrev[i] +
                 (1.f - gammaNoiseTmp) * (probNonSpeech * magn[i] +
                                          probSpeech * nssc->noisePrev[i]);
      // Allow for noise update downwards:
      // If noise update decreases the noise, it is safe, so allow it to
      // happen.
      if (noiseUpdateTmp < noise[i]) {
        noise[i] = noiseUpdateTmp;
      }
    }
  }  // End of freq loop.
}

// Updates |buffer| with a new |frame|.
// Inputs:
//   * |frame| is a new speech frame or NULL for setting to zero.
//   * |frame_length| is the length of the new frame.
//   * |buffer_length| is the length of the buffer.
// Output:
//   * |buffer| is the updated buffer.
static void UpdateBuffer(const float* frame,
                         size_t frame_length,
                         size_t buffer_length,
                         float* buffer) {
  RTC_DCHECK_LT(buffer_length, 2 * frame_length);

  memcpy(buffer,
         buffer + frame_length,
         sizeof(*buffer) * (buffer_length - frame_length));
  if (frame) {
    memcpy(buffer + buffer_length - frame_length,
           frame,
           sizeof(*buffer) * frame_length);
  } else {
    memset(buffer + buffer_length - frame_length,
           0,
           sizeof(*buffer) * frame_length);
  }
}

// Transforms the signal from time to frequency domain.
// Inputs:
//   * |time_data| is the signal in the time domain.
//   * |time_data_length| is the length of the analysis buffer.
//   * |magnitude_length| is the length of the spectrum magnitude, which equals
//     the length of both |real| and |imag| (time_data_length / 2 + 1).
// Outputs:
//   * |time_data| is the signal in the frequency domain.
//   * |real| is the real part of the frequency domain.
//   * |imag| is the imaginary part of the frequency domain.
//   * |magn| is the calculated signal magnitude in the frequency domain.
static void FFT(NoiseSuppressionC* nssc,
                float* time_data,
                size_t time_data_length,
                size_t magnitude_length,
                float* real,
                float* imag,
                float* magn) {
  size_t i;

  RTC_DCHECK_EQ(magnitude_length, time_data_length / 2 + 1);

  WebRtc_rdft(time_data_length, 1, time_data, nssc->ip, nssc->wfft);

  imag[0] = 0;
  real[0] = time_data[0];
  magn[0] = fabsf(real[0]) + 1.f;
  imag[magnitude_length - 1] = 0;
  real[magnitude_length - 1] = time_data[1];
  magn[magnitude_length - 1] = fabsf(real[magnitude_length - 1]) + 1.f;
  for (i = 1; i < magnitude_length - 1; ++i) {
    real[i] = time_data[2 * i];
    imag[i] = time_data[2 * i + 1];
    // Magnitude spectrum.
    magn[i] = sqrtf(real[i] * real[i] + imag[i] * imag[i]) + 1.f;
  }
}

// Transforms the signal from frequency to time domain.
// Inputs:
//   * |real| is the real part of the frequency domain.
//   * |imag| is the imaginary part of the frequency domain.
//   * |magnitude_length| is the length of the spectrum magnitude, which equals
//     the length of both |real| and |imag|.
//   * |time_data_length| is the length of the analysis buffer
//     (2 * (magnitude_length - 1)).
// Output:
//   * |time_data| is the signal in the time domain.
static void IFFT(NoiseSuppressionC* nssc,
                 const float* real,
                 const float* imag,
                 size_t magnitude_length,
                 size_t time_data_length,
                 float* time_data) {
  size_t i;

  RTC_DCHECK_EQ(time_data_length, 2 * (magnitude_length - 1));

  time_data[0] = real[0];
  time_data[1] = real[magnitude_length - 1];
  for (i = 1; i < magnitude_length - 1; ++i) {
    time_data[2 * i] = real[i];
    time_data[2 * i + 1] = imag[i];
  }
  WebRtc_rdft(time_data_length, -1, time_data, nssc->ip, nssc->wfft);

  for (i = 0; i < time_data_length; ++i) {
    time_data[i] *= 2.f / time_data_length;  // FFT scaling.
  }
}

// Calculates the energy of a buffer.
// Inputs:
//   * |buffer| is the buffer over which the energy is calculated.
//   * |length| is the length of the buffer.
// Returns the calculated energy.
static float Energy(const float* buffer, size_t length) {
  size_t i;
  float energy = 0.f;

  for (i = 0; i < length; ++i) {
    energy += buffer[i] * buffer[i];
  }

  return energy;
}

// Windows a buffer.
// Inputs:
//   * |window| is the window by which to multiply.
//   * |data| is the data without windowing.
//   * |length| is the length of the window and data.
// Output:
//   * |data_windowed| is the windowed data.
static void Windowing(const float* window,
                      const float* data,
                      size_t length,
                      float* data_windowed) {
  size_t i;

  for (i = 0; i < length; ++i) {
    data_windowed[i] = window[i] * data[i];
  }
}

// Estimate prior SNR decision-directed and compute DD based Wiener Filter.
// Input:
//   * |magn| is the signal magnitude spectrum estimate.
// Output:
//   * |theFilter| is the frequency response of the computed Wiener filter.
static void ComputeDdBasedWienerFilter(const NoiseSuppressionC* nssc,
                                       const float* magn,
                                       float* theFilter) {
  size_t i;
  float snrPrior, previousEstimateStsa, currentEstimateStsa;

  for (i = 0; i < nssc->magnLen; i++) {
    // Previous estimate: based on previous frame with gain filter.
    previousEstimateStsa = nssc->magnPrevProcess[i] /
                           (nssc->noisePrev[i] + 0.0001f) * nssc->smooth[i];
    // Post and prior SNR.
    currentEstimateStsa = 0.f;
    if (magn[i] > nssc->noise[i]) {
      currentEstimateStsa = magn[i] / (nssc->noise[i] + 0.0001f) - 1.f;
    }
    // DD estimate is sum of two terms: current estimate and previous estimate.
    // Directed decision update of |snrPrior|.
    snrPrior = DD_PR_SNR * previousEstimateStsa +
               (1.f - DD_PR_SNR) * currentEstimateStsa;
    // Gain filter.
    theFilter[i] = snrPrior / (nssc->overdrive + snrPrior);
  }  // End of loop over frequencies.
}

// Changes the aggressiveness of the noise suppression method.
// |mode| = 0 is mild (6dB), |mode| = 1 is medium (10dB) and |mode| = 2 is
// aggressive (15dB).
// Returns 0 on success and -1 otherwise.
int WebRtcNs_set_policy_core(NoiseSuppressionC* nssc, int mode) {
  // Allow for modes: 0, 1, 2, 3.
  if (mode < 0 || mode > 3) {
    return (-1);
  }

  nssc->aggrMode = mode;
  if (mode == 0) {
    nssc->overdrive = 1.f;
    nssc->denoiseBound = 0.5f;
    nssc->gainmap = 0;
  } else if (mode == 1) {
    // nssc->overdrive = 1.25f;
    nssc->overdrive = 1.f;
    nssc->denoiseBound = 0.25f;
    nssc->gainmap = 1;
  } else if (mode == 2) {
    // nssc->overdrive = 1.25f;
    nssc->overdrive = 1.1f;
    nssc->denoiseBound = 0.125f;
    nssc->gainmap = 1;
  } else if (mode == 3) {
    // nssc->overdrive = 1.3f;
    nssc->overdrive = 1.25f;
    nssc->denoiseBound = 0.09f;
    nssc->gainmap = 1;
  }
  return 0;
}

void WebRtcNs_AnalyzeCore(NoiseSuppressionC* nssc, const float* speechFrame) {
  size_t i;
  const size_t kStartBand = 5;  // Skip first frequency bins during estimation.
  int updateParsFlag;
  float energy;
  float signalEnergy = 0.f;
  float sumMagn = 0.f;
  float tmpFloat1, tmpFloat2, tmpFloat3;
  float winData[ANAL_BLOCKL_MAX];
  float magn[HALF_ANAL_BLOCKL], noise[HALF_ANAL_BLOCKL];
  float snrLocPost[HALF_ANAL_BLOCKL], snrLocPrior[HALF_ANAL_BLOCKL];
  float real[ANAL_BLOCKL_MAX], imag[HALF_ANAL_BLOCKL];
  // Variables during startup.
  float sum_log_i = 0.0;
  float sum_log_i_square = 0.0;
  float sum_log_magn = 0.0;
  float sum_log_i_log_magn = 0.0;
  float parametric_exp = 0.0;
  float parametric_num = 0.0;

  // Check that initiation has been done.
  RTC_DCHECK_EQ(1, nssc->initFlag);
  updateParsFlag = nssc->modelUpdatePars[0];

  // Update analysis buffer for L band.
  UpdateBuffer(speechFrame, nssc->blockLen, nssc->anaLen, nssc->analyzeBuf);

  Windowing(nssc->window, nssc->analyzeBuf, nssc->anaLen, winData);
  energy = Energy(winData, nssc->anaLen);
  if (energy == 0.0) {
    // We want to avoid updating statistics in this case:
    // Updating feature statistics when we have zeros only will cause
    // thresholds to move towards zero signal situations. This in turn has the
    // effect that once the signal is "turned on" (non-zero values) everything
    // will be treated as speech and there is no noise suppression effect.
    // Depending on the duration of the inactive signal it takes a
    // considerable amount of time for the system to learn what is noise and
    // what is speech.
    return;
  }

  nssc->blockInd++;  // Update the block index only when we process a block.

  FFT(nssc, winData, nssc->anaLen, nssc->magnLen, real, imag, magn);

  for (i = 0; i < nssc->magnLen; i++) {
    signalEnergy += real[i] * real[i] + imag[i] * imag[i];
    sumMagn += magn[i];
    if (nssc->blockInd < END_STARTUP_SHORT) {
      if (i >= kStartBand) {
        tmpFloat2 = logf((float)i);
        sum_log_i += tmpFloat2;
        sum_log_i_square += tmpFloat2 * tmpFloat2;
        tmpFloat1 = logf(magn[i]);
        sum_log_magn += tmpFloat1;
        sum_log_i_log_magn += tmpFloat2 * tmpFloat1;
      }
    }
  }
  signalEnergy /= nssc->magnLen;
  nssc->signalEnergy = signalEnergy;
  nssc->sumMagn = sumMagn;

  // Quantile noise estimate.
  NoiseEstimation(nssc, magn, noise);
  // Compute simplified noise model during startup.
  if (nssc->blockInd < END_STARTUP_SHORT) {
    // Estimate White noise.
    nssc->whiteNoiseLevel += sumMagn / nssc->magnLen * nssc->overdrive;
    // Estimate Pink noise parameters.
    tmpFloat1 = sum_log_i_square * (nssc->magnLen - kStartBand);
    tmpFloat1 -= (sum_log_i * sum_log_i);
    tmpFloat2 =
        (sum_log_i_square * sum_log_magn - sum_log_i * sum_log_i_log_magn);
    tmpFloat3 = tmpFloat2 / tmpFloat1;
    // Constrain the estimated spectrum to be positive.
    if (tmpFloat3 < 0.f) {
      tmpFloat3 = 0.f;
    }
    nssc->pinkNoiseNumerator += tmpFloat3;
    tmpFloat2 = (sum_log_i * sum_log_magn);
    tmpFloat2 -= (nssc->magnLen - kStartBand) * sum_log_i_log_magn;
    tmpFloat3 = tmpFloat2 / tmpFloat1;
    // Constrain the pink noise power to be in the interval [0, 1].
    if (tmpFloat3 < 0.f) {
      tmpFloat3 = 0.f;
    }
    if (tmpFloat3 > 1.f) {
      tmpFloat3 = 1.f;
    }
    nssc->pinkNoiseExp += tmpFloat3;

    // Calculate frequency independent parts of parametric noise estimate.
    if (nssc->pinkNoiseExp > 0.f) {
      // Use pink noise estimate.
      parametric_num =
          expf(nssc->pinkNoiseNumerator / (float)(nssc->blockInd + 1));
      parametric_num *= (float)(nssc->blockInd + 1);
      parametric_exp = nssc->pinkNoiseExp / (float)(nssc->blockInd + 1);
    }
    for (i = 0; i < nssc->magnLen; i++) {
      // Estimate the background noise using the white and pink noise
      // parameters.
      if (nssc->pinkNoiseExp == 0.f) {
        // Use white noise estimate.
        nssc->parametricNoise[i] = nssc->whiteNoiseLevel;
      } else {
        // Use pink noise estimate.
        float use_band = (float)(i < kStartBand ? kStartBand : i);
        nssc->parametricNoise[i] =
            parametric_num / powf(use_band, parametric_exp);
      }
      // Weight quantile noise with modeled noise.
      noise[i] *= (nssc->blockInd);
      tmpFloat2 =
          nssc->parametricNoise[i] * (END_STARTUP_SHORT - nssc->blockInd);
      noise[i] += (tmpFloat2 / (float)(nssc->blockInd + 1));
      noise[i] /= END_STARTUP_SHORT;
    }
  }
  // Compute average signal during END_STARTUP_LONG time:
  // used to normalize spectral difference measure.
  if (nssc->blockInd < END_STARTUP_LONG) {
    nssc->featureData[5] *= nssc->blockInd;
    nssc->featureData[5] += signalEnergy;
    nssc->featureData[5] /= (nssc->blockInd + 1);
  }

  // Post and prior SNR needed for SpeechNoiseProb.
  ComputeSnr(nssc, magn, noise, snrLocPrior, snrLocPost);

  FeatureUpdate(nssc, magn, updateParsFlag);
  SpeechNoiseProb(nssc, nssc->speechProb, snrLocPrior, snrLocPost);
  UpdateNoiseEstimate(nssc, magn, snrLocPrior, snrLocPost, noise);

  // Keep track of noise spectrum for next frame.
  memcpy(nssc->noise, noise, sizeof(*noise) * nssc->magnLen);
  memcpy(nssc->magnPrevAnalyze, magn, sizeof(*magn) * nssc->magnLen);
}

void WebRtcNs_ProcessCore(NoiseSuppressionC* nssc,
                          const float* const* speechFrame,
                          size_t num_bands,
                          float* const* outFrame) {
  // Main routine for noise reduction.
  int flagHB = 0;
  size_t i, j;

  float energy1, energy2, gain, factor, factor1, factor2;
  float fout[BLOCKL_MAX];
  float winData[ANAL_BLOCKL_MAX];
  float magn[HALF_ANAL_BLOCKL];
  float theFilter[HALF_ANAL_BLOCKL], theFilterTmp[HALF_ANAL_BLOCKL];
  float real[ANAL_BLOCKL_MAX], imag[HALF_ANAL_BLOCKL];

  // SWB variables.
  int deltaBweHB = 1;
  int deltaGainHB = 1;
  float decayBweHB = 1.0;
  float gainMapParHB = 1.0;
  float gainTimeDomainHB = 1.0;
  float avgProbSpeechHB, avgProbSpeechHBTmp, avgFilterGainHB, gainModHB;
  float sumMagnAnalyze, sumMagnProcess;

  // Check that initiation has been done.
  RTC_DCHECK_EQ(1, nssc->initFlag);
  RTC_DCHECK_LE(num_bands - 1, NUM_HIGH_BANDS_MAX);

  const float* const* speechFrameHB = NULL;
  float* const* outFrameHB = NULL;
  size_t num_high_bands = 0;
  if (num_bands > 1) {
    speechFrameHB = &speechFrame[1];
    outFrameHB = &outFrame[1];
    num_high_bands = num_bands - 1;
    flagHB = 1;
    // Range for averaging low band quantities for H band gain.
    deltaBweHB = (int)nssc->magnLen / 4;
    deltaGainHB = deltaBweHB;
  }

  // Update analysis buffer for L band.
  UpdateBuffer(speechFrame[0], nssc->blockLen, nssc->anaLen, nssc->dataBuf);

  if (flagHB == 1) {
    // Update analysis buffer for H bands.
    for (i = 0; i < num_high_bands; ++i) {
      UpdateBuffer(speechFrameHB[i],
                   nssc->blockLen,
                   nssc->anaLen,
                   nssc->dataBufHB[i]);
    }
  }

  Windowing(nssc->window, nssc->dataBuf, nssc->anaLen, winData);
  energy1 = Energy(winData, nssc->anaLen);
  if (energy1 == 0.0) {
    // Synthesize the special case of zero input.
    // Read out fully processed segment.
    for (i = nssc->windShift; i < nssc->blockLen + nssc->windShift; i++) {
      fout[i - nssc->windShift] = nssc->syntBuf[i];
    }
    // Update synthesis buffer.
    UpdateBuffer(NULL, nssc->blockLen, nssc->anaLen, nssc->syntBuf);

    for (i = 0; i < nssc->blockLen; ++i)
      outFrame[0][i] =
          WEBRTC_SPL_SAT(WEBRTC_SPL_WORD16_MAX, fout[i], WEBRTC_SPL_WORD16_MIN);

    // For time-domain gain of HB.
    if (flagHB == 1) {
      for (i = 0; i < num_high_bands; ++i) {
        for (j = 0; j < nssc->blockLen; ++j) {
          outFrameHB[i][j] = WEBRTC_SPL_SAT(WEBRTC_SPL_WORD16_MAX,
                                            nssc->dataBufHB[i][j],
                                            WEBRTC_SPL_WORD16_MIN);
        }
      }
    }

    return;
  }

  FFT(nssc, winData, nssc->anaLen, nssc->magnLen, real, imag, magn);

  if (nssc->blockInd < END_STARTUP_SHORT) {
    for (i = 0; i < nssc->magnLen; i++) {
      nssc->initMagnEst[i] += magn[i];
    }
  }

  ComputeDdBasedWienerFilter(nssc, magn, theFilter);

  for (i = 0; i < nssc->magnLen; i++) {
    // Flooring bottom.
    if (theFilter[i] < nssc->denoiseBound) {
      theFilter[i] = nssc->denoiseBound;
    }
    // Flooring top.
    if (theFilter[i] > 1.f) {
      theFilter[i] = 1.f;
    }
    if (nssc->blockInd < END_STARTUP_SHORT) {
      theFilterTmp[i] =
          (nssc->initMagnEst[i] - nssc->overdrive * nssc->parametricNoise[i]);
      theFilterTmp[i] /= (nssc->initMagnEst[i] + 0.0001f);
      // Flooring bottom.
      if (theFilterTmp[i] < nssc->denoiseBound) {
        theFilterTmp[i] = nssc->denoiseBound;
      }
      // Flooring top.
      if (theFilterTmp[i] > 1.f) {
        theFilterTmp[i] = 1.f;
      }
      // Weight the two suppression filters.
      theFilter[i] *= (nssc->blockInd);
      theFilterTmp[i] *= (END_STARTUP_SHORT - nssc->blockInd);
      theFilter[i] += theFilterTmp[i];
      theFilter[i] /= (END_STARTUP_SHORT);
    }

    nssc->smooth[i] = theFilter[i];
    real[i] *= nssc->smooth[i];
    imag[i] *= nssc->smooth[i];
  }
  // Keep track of |magn| spectrum for next frame.
  memcpy(nssc->magnPrevProcess, magn, sizeof(*magn) * nssc->magnLen);
  memcpy(nssc->noisePrev, nssc->noise, sizeof(nssc->noise[0]) * nssc->magnLen);
  // Back to time domain.
  IFFT(nssc, real, imag, nssc->magnLen, nssc->anaLen, winData);

  // Scale factor: only do it after END_STARTUP_LONG time.
  factor = 1.f;
  if (nssc->gainmap == 1 && nssc->blockInd > END_STARTUP_LONG) {
    factor1 = 1.f;
    factor2 = 1.f;

    energy2 = Energy(winData, nssc->anaLen);
    gain = (float)sqrt(energy2 / (energy1 + 1.f));

    // Scaling for new version.
    if (gain > B_LIM) {
      factor1 = 1.f + 1.3f * (gain - B_LIM);
      if (gain * factor1 > 1.f) {
        factor1 = 1.f / gain;
      }
    }
    if (gain < B_LIM) {
      // Don't reduce scale too much for pause regions:
      // attenuation here should be controlled by flooring.
      if (gain <= nssc->denoiseBound) {
        gain = nssc->denoiseBound;
      }
      factor2 = 1.f - 0.3f * (B_LIM - gain);
    }
    // Combine both scales with speech/noise prob:
    // note prior (priorSpeechProb) is not frequency dependent.
    factor = nssc->priorSpeechProb * factor1 +
             (1.f - nssc->priorSpeechProb) * factor2;
  }  // Out of nssc->gainmap == 1.

  Windowing(nssc->window, winData, nssc->anaLen, winData);

  // Synthesis.
  for (i = 0; i < nssc->anaLen; i++) {
    nssc->syntBuf[i] += factor * winData[i];
  }
  // Read out fully processed segment.
  for (i = nssc->windShift; i < nssc->blockLen + nssc->windShift; i++) {
    fout[i - nssc->windShift] = nssc->syntBuf[i];
  }
  // Update synthesis buffer.
  UpdateBuffer(NULL, nssc->blockLen, nssc->anaLen, nssc->syntBuf);

  for (i = 0; i < nssc->blockLen; ++i)
    outFrame[0][i] =
        WEBRTC_SPL_SAT(WEBRTC_SPL_WORD16_MAX, fout[i], WEBRTC_SPL_WORD16_MIN);

  // For time-domain gain of HB.
  if (flagHB == 1) {
    // Average speech prob from low band.
    // Average over second half (i.e., 4->8kHz) of frequencies spectrum.
    avgProbSpeechHB = 0.0;
    for (i = nssc->magnLen - deltaBweHB - 1; i < nssc->magnLen - 1; i++) {
      avgProbSpeechHB += nssc->speechProb[i];
    }
    avgProbSpeechHB = avgProbSpeechHB / ((float)deltaBweHB);
    // If the speech was suppressed by a component between Analyze and
    // Process, for example the AEC, then it should not be considered speech
    // for high band suppression purposes.
    sumMagnAnalyze = 0;
    sumMagnProcess = 0;
    for (i = 0; i < nssc->magnLen; ++i) {
      sumMagnAnalyze += nssc->magnPrevAnalyze[i];
      sumMagnProcess += nssc->magnPrevProcess[i];
    }
    avgProbSpeechHB *= sumMagnProcess / sumMagnAnalyze;
    // Average filter gain from low band.
    // Average over second half (i.e., 4->8kHz) of frequencies spectrum.
    avgFilterGainHB = 0.0;
    for (i = nssc->magnLen - deltaGainHB - 1; i < nssc->magnLen - 1; i++) {
      avgFilterGainHB += nssc->smooth[i];
    }
    avgFilterGainHB = avgFilterGainHB / ((float)(deltaGainHB));
    avgProbSpeechHBTmp = 2.f * avgProbSpeechHB - 1.f;
    // Gain based on speech probability.
    gainModHB = 0.5f * (1.f + (float)tanh(gainMapParHB * avgProbSpeechHBTmp));
    // Combine gain with low band gain.
    gainTimeDomainHB = 0.5f * gainModHB + 0.5f * avgFilterGainHB;
    if (avgProbSpeechHB >= 0.5f) {
      gainTimeDomainHB = 0.25f * gainModHB + 0.75f * avgFilterGainHB;
    }
    gainTimeDomainHB = gainTimeDomainHB * decayBweHB;
    // Make sure gain is within flooring range.
    // Flooring bottom.
    if (gainTimeDomainHB < nssc->denoiseBound) {
      gainTimeDomainHB = nssc->denoiseBound;
    }
    // Flooring top.
    if (gainTimeDomainHB > 1.f) {
      gainTimeDomainHB = 1.f;
    }
    // Apply gain.
    for (i = 0; i < num_high_bands; ++i) {
      for (j = 0; j < nssc->blockLen; j++) {
        outFrameHB[i][j] =
            WEBRTC_SPL_SAT(WEBRTC_SPL_WORD16_MAX,
                           gainTimeDomainHB * nssc->dataBufHB[i][j],
                           WEBRTC_SPL_WORD16_MIN);
      }
    }
  }  // End of H band gain computation.
}
