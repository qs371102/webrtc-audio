/*
 *  Copyright 2009 The WebRTC Project Authors. All rights reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef WEBRTC_BASE_WIN32SOCKETINIT_H_
#define WEBRTC_BASE_WIN32SOCKETINIT_H_
 #ifdef WEBRTC_WIN

namespace rtc {

void EnsureWinsockInit();

}  // namespace rtc
#endif
#endif  // WEBRTC_BASE_WIN32SOCKETINIT_H_
