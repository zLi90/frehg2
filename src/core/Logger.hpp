/// \file Logger.hpp
/// \brief Rank-aware logging with severity levels and a throwing fatal path.
///
/// All user-visible text produced by the libraries flows through this logger
/// (direct \c printf / \c std::cout use outside Logger.cpp and main.cpp is
/// forbidden and CI-checked, plan §11.3). Informational messages print on
/// rank 0 only; warnings and errors print on every rank that emits them,
/// prefixed with the rank id.

#ifndef FREHG_CORE_LOGGER_HPP
#define FREHG_CORE_LOGGER_HPP

#include <mpi.h>

#include <sstream>
#include <string>

namespace frehg::log {

/// Message severities, in increasing order of importance.
enum class Level {
  Debug,  ///< developer diagnostics, hidden by default
  Info,   ///< normal progress reporting (rank 0 only)
  Warn,   ///< recoverable anomalies
  Error   ///< serious failures (fatal errors are logged at this level)
};

/// Bind the logger to a communicator and minimum severity. Called once by
/// PetscSession; safe to call again in tests. Before init() the logger
/// behaves as a single-rank logger at Info level.
void init(MPI_Comm comm, Level minLevel = Level::Info);

/// Set the minimum severity that is actually printed.
void setLevel(Level minLevel);

/// \return the current minimum severity.
Level level();

/// Log a developer-diagnostic message (printed on every rank when enabled).
void debug(const std::string& message);

/// Log an informational message (printed on rank 0 only).
void info(const std::string& message);

/// Log a warning (printed on every rank that emits it).
void warn(const std::string& message);

/// Log an error (printed on every rank that emits it).
void error(const std::string& message);

/// Log \p message at Error severity, then throw frehg::FatalError.
///
/// This is the single fatal-error path of the code base: callers must not
/// continue after invalid input, failed library calls, or solver divergence
/// (plan §11.1 rule 5). \c main catches the exception and terminates the MPI
/// job; unit tests assert on it.
[[noreturn]] void fatal(const std::string& message);

/// Convenience for building messages from streamable parts:
/// \code log::info(log::msg() << "nx=" << nx); \endcode
class msg {
 public:
  /// Append any streamable value to the message.
  template <class T>
  msg& operator<<(const T& value) {
    stream_ << value;
    return *this;
  }
  /// Implicit conversion so msg can be passed where std::string is expected.
  operator std::string() const { return stream_.str(); }  // NOLINT(google-explicit-constructor)

 private:
  std::ostringstream stream_;
};

}  // namespace frehg::log

#endif  // FREHG_CORE_LOGGER_HPP
