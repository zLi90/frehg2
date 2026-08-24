/// \file PetscSession.hpp
/// \brief RAII initialization of MPI, Kokkos, and PETSc (plan §10 P0).
///
/// Exactly one PetscSession is constructed at process start (in \c main or a
/// test main). Construction initializes MPI, then Kokkos, then PETSc;
/// destruction finalizes them in reverse order. No other component calls the
/// libraries' init/finalize routines.

#ifndef FREHG_CORE_PETSCSESSION_HPP
#define FREHG_CORE_PETSCSESSION_HPP

#include <mpi.h>

#include <string>

namespace frehg {

/// RAII owner of the MPI + Kokkos + PETSc runtime.
class PetscSession {
 public:
  /// Initialize the runtime.
  /// \param argc, argv command-line arguments, forwarded to Kokkos and PETSc.
  /// \param petscOptionsFile optional PETSc options file inserted into the
  ///        global options database (config key solver.petsc_options_file).
  PetscSession(int& argc, char**& argv, const std::string& petscOptionsFile = std::string());

  /// Finalize PETSc, Kokkos, and MPI in reverse initialization order.
  ~PetscSession();

  PetscSession(const PetscSession&) = delete;
  PetscSession& operator=(const PetscSession&) = delete;

  /// \return the world communicator for this run.
  MPI_Comm comm() const { return comm_; }

  /// \return this process's rank in comm().
  int rank() const { return rank_; }

  /// \return the number of ranks in comm().
  int size() const { return size_; }

 private:
  MPI_Comm comm_ = MPI_COMM_NULL;
  int rank_ = 0;
  int size_ = 1;
  bool ownsMpi_ = false;
};

}  // namespace frehg

#endif  // FREHG_CORE_PETSCSESSION_HPP
