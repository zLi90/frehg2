/// \file main.cpp
/// \brief Test entry point: one PetscSession wraps the whole run.

#include "core/PetscSession.hpp"

#include <gtest/gtest.h>

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  frehg::PetscSession session(argc, argv);
  return RUN_ALL_TESTS();
}
