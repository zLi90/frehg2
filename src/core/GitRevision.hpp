/// \file GitRevision.hpp
/// \brief Source revision of the running binary, for output-file provenance.
///
/// The definition is generated into the build tree by
/// cmake/FrehgGitRevision.cmake (v2 plan V2-A20 finding 5). The run record's
/// provenance.git_sha, the HDF5 /frehg2 git_sha attribute and the log header
/// all read it, so every output names the code that produced it.

#ifndef FREHG_CORE_GIT_REVISION_HPP
#define FREHG_CORE_GIT_REVISION_HPP

namespace frehg {

/// Git revision of the source tree this binary was built from.
///
/// \return "<sha12>" (the 12-hex-digit abbreviated commit SHA) for a clean
/// tree, "<sha12>-dirty" when tracked files differed from that commit (staged
/// or not; untracked files are not considered), or "unknown" when the build
/// had no git or the source tree is not a git checkout (a release tarball).
const char* gitRevision() noexcept;

}  // namespace frehg

#endif  // FREHG_CORE_GIT_REVISION_HPP
