/// \file PetscSession.cpp
/// \brief Implementation of the RAII MPI/Kokkos/PETSc session.

#include "core/PetscSession.hpp"

#include "core/Logger.hpp"
#include "core/Types.hpp"

#include <Kokkos_Core.hpp>
#include <petscsys.h>

namespace frehg {

PetscSession::PetscSession(int& argc, char**& argv, const std::string& petscOptionsFile) {
  int mpiInitialized = 0;
  MPI_Initialized(&mpiInitialized);
  if (mpiInitialized == 0) {
    int provided = 0;
    if (MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided) != MPI_SUCCESS) {
      log::fatal("MPI_Init_thread failed");
    }
    ownsMpi_ = true;
  }
  comm_ = MPI_COMM_WORLD;
  MPI_Comm_rank(comm_, &rank_);
  MPI_Comm_size(comm_, &size_);
  log::init(comm_);

  if (!Kokkos::is_initialized()) {
    Kokkos::initialize(argc, argv);
  }

  PetscBool petscInitialized = PETSC_FALSE;
  PetscInitialized(&petscInitialized);
  if (petscInitialized == PETSC_FALSE) {
    if (PetscInitialize(&argc, &argv, nullptr, nullptr) != PETSC_SUCCESS) {
      log::fatal("PetscInitialize failed");
    }
  }

  if (!petscOptionsFile.empty()) {
    if (PetscOptionsInsertFile(comm_, nullptr, petscOptionsFile.c_str(), PETSC_TRUE) !=
        PETSC_SUCCESS) {
      log::fatal(log::msg() << "failed to read PETSc options file '" << petscOptionsFile << "'");
    }
    log::info(log::msg() << "loaded PETSc options file '" << petscOptionsFile << "'");
  }

  log::info(log::msg() << "frehg2 " << FREHG_VERSION << " (git " << FREHG_GIT_SHA << "), "
                       << size_ << " MPI rank" << (size_ > 1 ? "s" : "") << ", Kokkos backend "
                       << ExecSpace::name());
}

PetscSession::~PetscSession() {
  PetscBool petscFinalized = PETSC_FALSE;
  PetscFinalized(&petscFinalized);
  if (petscFinalized == PETSC_FALSE) {
    PetscFinalize();
  }
  if (Kokkos::is_initialized()) {
    Kokkos::finalize();
  }
  if (ownsMpi_) {
    int mpiFinalized = 0;
    MPI_Finalized(&mpiFinalized);
    if (mpiFinalized == 0) {
      MPI_Finalize();
    }
  }
}

}  // namespace frehg
