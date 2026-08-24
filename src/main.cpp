/// \file main.cpp
/// \brief Frehg2 executable entry point.
///
/// `frehg <config.yaml>` validates the configuration and runs the
/// simulation; `frehg --validate <config.yaml>` runs the full schema +
/// cross-field validation only and reports every problem found. Extra
/// arguments after the configuration path go to PETSc (e.g. -fs_pc_type
/// jacobi for the strict rank-invariance mode of plan §8.2).

#include "core/Config.hpp"
#include "core/Logger.hpp"
#include "core/PetscSession.hpp"
#include "core/Timer.hpp"
#include "core/Types.hpp"
#include "driver/Simulation.hpp"

#include <mpi.h>

#include <iostream>
#include <string>
#include <vector>

namespace {

int printUsage() {
  std::cout << "frehg " << FREHG_VERSION << " (git " << FREHG_GIT_SHA << ")\n"
            << "usage:\n"
            << "  frehg <config.yaml> [petsc options]   run a simulation\n"
            << "  frehg --validate <config.yaml>        validate a configuration and exit\n";
  return 2;
}

int runValidate(const std::string& path) {
  const frehg::ValidationResult result = frehg::validateConfigFile(path);
  if (result.ok()) {
    std::cout << "VALID: " << path << "\n";
    const frehg::FrehgConfig config = frehg::loadConfig(path);
    std::cout << frehg::describeConfig(config) << "\n";
    return 0;
  }
  std::cout << "INVALID: " << path << " (" << result.errors.size() << " error"
            << (result.errors.size() == 1 ? "" : "s") << ")\n";
  for (const std::string& error : result.errors) {
    std::cout << "  - " << error << "\n";
  }
  return 1;
}

int runSimulation(const std::string& path) {
  const frehg::FrehgConfig config = frehg::loadConfig(path);
  frehg::log::info(frehg::log::msg() << "frehg " << FREHG_VERSION << " (git " << FREHG_GIT_SHA
                                     << ") running " << config.simulation.id);
  frehg::log::info(frehg::describeConfig(config));
  frehg::driver::Simulation simulation(MPI_COMM_WORLD, config);
  simulation.run();
  const std::string timers = frehg::Timer::report(MPI_COMM_WORLD);
  if (!timers.empty()) {
    frehg::log::info(timers);
  }
  frehg::log::info("run complete");
  return 0;
}

}  // namespace

/// Entry point: dispatches --validate or a run; fatal errors terminate the
/// MPI job.
int main(int argc, char** argv) {
  const std::vector<std::string> args(argv + 1, argv + argc);
  const bool validate = !args.empty() && args[0] == "--validate";
  if (args.empty() || (validate && args.size() != 2)) {
    return printUsage();
  }

  int status = 0;
  try {
    frehg::PetscSession session(argc, argv);
    status = validate ? runValidate(args[1]) : runSimulation(args[0]);
  } catch (const frehg::FatalError& e) {
    std::cerr << "frehg: fatal: " << e.what() << "\n";
    int mpiInitialized = 0;
    MPI_Initialized(&mpiInitialized);
    if (mpiInitialized != 0) {
      int size = 1;
      MPI_Comm_size(MPI_COMM_WORLD, &size);
      if (size > 1) {
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }
    status = 1;
  }
  return status;
}
