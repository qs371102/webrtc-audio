#include <iostream>

#include "audio_processing.h"

using namespace webrtc;

int main() {
    std::cout << "Hello, World!" << std::endl;
    // Usage example, omitting error checking:
    AudioProcessing* apm = AudioProcessing::Create();

    apm->high_pass_filter()->Enable(true);

    apm->echo_cancellation()->enable_drift_compensation(false);
    apm->echo_cancellation()->Enable(true);
    apm->gain_control()->set_analog_level_limits(0, 255);
    apm->gain_control()->set_mode(GainControl::kAdaptiveAnalog);
    apm->gain_control()->Enable(true);

    apm->voice_detection()->Enable(true);

    // Start a voice call...

    // ... Render frame arrives bound for the audio HAL ...
    // apm->ProcessReverseStream(render_frame);

    // ... Capture frame arrives from the audio HAL ...
    // Call required set_stream_ functions.
    // apm->set_stream_delay_ms(delay_ms);
    // apm->gain_control()->set_stream_analog_level(analog_level);

    // apm->ProcessStream(capture_frame);

    // Call required stream_ functions.
    int analog_level = apm->gain_control()->stream_analog_level();

    // Repeate render and capture processing for the duration of the call...
    // Start a new call...
    apm->Initialize();

    // Close the application...
    delete apm;
    return 0;
}