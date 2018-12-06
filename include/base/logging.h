// Copyright 2018 Google Inc. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS-IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//

#ifndef EG_BASE_LOGGING_H_
#define EG_BASE_LOGGING_H_

#ifdef EG_USE_GLOG

#include <glog/logging.h>

// The names CHECK, etc. are too common and may conflict with other
// packages.  We use EG_CHECK to make it easier to switch to
// something other than GLOG for logging.

#define EG_LOG LOG
#define EG_LOG_IF LOG_IF
#define EG_DLOG_IF DLOG_IF

#define EG_CHECK CHECK
#define EG_CHECK_EQ CHECK_EQ
#define EG_CHECK_NE CHECK_NE
#define EG_CHECK_LT CHECK_LT
#define EG_CHECK_LE CHECK_LE
#define EG_CHECK_GT CHECK_GT
#define EG_CHECK_GE CHECK_GE

#define EG_DCHECK DCHECK
#define EG_DCHECK_EQ DCHECK_EQ
#define EG_DCHECK_NE DCHECK_NE
#define EG_DCHECK_LT DCHECK_LT
#define EG_DCHECK_LE DCHECK_LE
#define EG_DCHECK_GT DCHECK_GT
#define EG_DCHECK_GE DCHECK_GE

#define EG_VLOG VLOG
#define EG_VLOG_IS_ON VLOG_IS_ON

#else  // !defined(EG_USE_GLOG)

#include <iostream>

#include "base/log_severity.h"
#include "third_party/absl/base/attributes.h"
#include "third_party/absl/base/log_severity.h"

class EGLogMessage {
 public:
  EGLogMessage(const char* file, int line,
               absl::LogSeverity severity, std::ostream& stream)
    : severity_(severity), stream_(stream) {
    if (enabled()) {
      stream_ << file << ":" << line << " "
              << absl::LogSeverityName(severity) << " ";
    }
  }
  ~EGLogMessage() { if (enabled()) stream_ << std::endl; }

  std::ostream& stream() { return stream_; }

 private:
  bool enabled() const {
#ifdef ABSL_MIN_LOG_LEVEL
    return (static_cast<int>(severity_) >= ABSL_MIN_LOG_LEVEL ||
            severity_ >= absl::LogSeverity::kFatal);
#else
    return true;
#endif
  }

  absl::LogSeverity severity_;
  std::ostream& stream_;
};

// Same as EGLogMessage, but destructor is marked no-return to avoid
// "no return value warnings" in functions that return non-void.
class EGFatalLogMessage : public EGLogMessage {
 public:
  EGFatalLogMessage(const char* file, int line,
                    absl::LogSeverity severity, std::ostream& stream)
      ABSL_ATTRIBUTE_COLD
    : EGLogMessage(file, line, severity, stream) {}
  ABSL_ATTRIBUTE_NORETURN ~EGFatalLogMessage() { abort(); }
};

// Logging stream that does nothing.
struct EGNullStream {
  template <typename T>
  EGNullStream& operator<<(const T& v) { return *this; }
};

// Used to suppress "unused value" warnings.
struct EGLogMessageVoidify {
  // Must have precedence lower than << but higher than ?:.
  void operator&(std::ostream&) {}
};

#define EG_LOG_MESSAGE_(LogMessageClass, log_severity) \
    LogMessageClass(__FILE__, __LINE__, log_severity, std::cerr)
#define EG_LOG_INFO \
    EG_LOG_MESSAGE_(EGLogMessage, absl::LogSeverity::kInfo)
#define EG_LOG_WARNING \
    EG_LOG_MESSAGE_(EGLogMessage, absl::LogSeverity::kWarning)
#define EG_LOG_ERROR \
    EG_LOG_MESSAGE_(EGLogMessage, absl::LogSeverity::kError)
#define EG_LOG_FATAL \
    EG_LOG_MESSAGE_(EGFatalLogMessage, absl::LogSeverity::kFatal)
#ifndef NDEBUG
#define EG_LOG_DFATAL EG_LOG_FATAL
#else
#define EG_LOG_DFATAL EG_LOG_ERROR
#endif

#define EG_LOG(severity) EG_LOG_##severity.stream()

// Implementing this as if (...) {} else EG_LOG(...) will cause dangling else
// warnings when someone does if (...) EG_LOG_IF(...), so do this tricky
// thing instead.
#define EG_LOG_IF(severity, condition) \
    !(condition) ? (void)0 : EGLogMessageVoidify() & EG_LOG(severity)

#define EG_CHECK(condition) \
    EG_LOG_IF(FATAL, ABSL_PREDICT_FALSE(!(condition))) \
        << ("Check failed: " #condition " ")

#ifndef NDEBUG

#define EG_DLOG_IF EG_LOG_IF
#define EG_DCHECK EG_CHECK

#else  // defined(NDEBUG)

#define EG_DLOG_IF(severity, condition) \
    while (false && (condition)) EGNullStream()
#define EG_DCHECK(condition) \
    while (false && (condition)) EGNullStream()

#endif  // defined(NDEBUG)

#define EG_CHECK_OP(op, val1, val2) EG_CHECK((val1) op (val2))
#define EG_CHECK_EQ(val1, val2) EG_CHECK_OP(==, val1, val2)
#define EG_CHECK_NE(val1, val2) EG_CHECK_OP(!=, val1, val2)
#define EG_CHECK_LT(val1, val2) EG_CHECK_OP(<, val1, val2)
#define EG_CHECK_LE(val1, val2) EG_CHECK_OP(<=, val1, val2)
#define EG_CHECK_GT(val1, val2) EG_CHECK_OP(>, val1, val2)
#define EG_CHECK_GE(val1, val2) EG_CHECK_OP(>=, val1, val2)

#define EG_DCHECK_OP(op, val1, val2) EG_DCHECK((val1) op (val2))
#define EG_DCHECK_EQ(val1, val2) EG_DCHECK_OP(==, val1, val2)
#define EG_DCHECK_NE(val1, val2) EG_DCHECK_OP(!=, val1, val2)
#define EG_DCHECK_LT(val1, val2) EG_DCHECK_OP(<, val1, val2)
#define EG_DCHECK_LE(val1, val2) EG_DCHECK_OP(<=, val1, val2)
#define EG_DCHECK_GT(val1, val2) EG_DCHECK_OP(>, val1, val2)
#define EG_DCHECK_GE(val1, val2) EG_DCHECK_OP(>=, val1, val2)

// We don't support VLOG.
#define EG_VLOG(verbose_level) EGNullStream()
#define EG_VLOG_IS_ON(verbose_level) (false)

#endif  // !defined(EG_USE_GLOG)

#endif  // EG_BASE_LOGGING_H_
