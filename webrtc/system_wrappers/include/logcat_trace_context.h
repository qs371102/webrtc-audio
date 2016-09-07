/*
 *  Copyright (c) 2013 The WebRTC project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef WEBRTC_SYSTEM_WRAPPERS_INCLUDE_LOGCAT_TRACE_CONTEXT_H_
#define WEBRTC_SYSTEM_WRAPPERS_INCLUDE_LOGCAT_TRACE_CONTEXT_H_

#include "webrtc/system_wrappers/include/trace.h"
#ifdef WEBRTC_ANDROID

namespace webrtc {

// Scoped helper class for directing Traces to Android's logcat facility.  While
// this object lives, Trace output will be sent to logcat.
class LogcatTraceContext : public webrtc::TraceCallback {
 public:
  LogcatTraceContext();
  ~LogcatTraceContext() override;

  // TraceCallback impl.
  void Print(TraceLevel level, const char* message, int length) override;
};

}  // namespace webrtc
#endif
#endif  // WEBRTC_SYSTEM_WRAPPERS_INCLUDE_LOGCAT_TRACE_CONTEXT_H_
