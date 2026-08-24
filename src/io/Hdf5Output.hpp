/// \file Hdf5Output.hpp
/// \brief Single-file parallel HDF5 output per the plan §7 contract.
///
/// Layout (all datasets float64, times keyed by integer seconds):
/// \verbatim
/// /frehg2                      root attrs: version, git sha, config, created
/// /grid/{x_center, y_center, z_center, bottom, ktop, dz}
/// /surface/<var>/<t>           [NY*NX]      index j*NX + i
/// /groundwater/<var>/<t>       [NY*NX*NZ]   index (j*NX + i)*NZ + k
/// /groundwater/zcell/0         static (b3 script compatibility)
/// /transport/<var>/<t>
/// /monitor/<name>              extendable (time, value...) table
/// /checkpoint/<t>/...          full prognostic state; attrs t, step
/// \endverbatim
///
/// Inactive 3D cells are written as NaN. Every field dataset carries units,
/// time, and long_name attributes. Writes are collective
/// (H5Pset_dxpl_mpio); each rank writes its (j, i) block as a regular 1D
/// hyperslab (start j0*NX+i0 [*NZ], stride NX [*NZ], count nyLocal, block
/// nxLocal [*NZ]).

#ifndef FREHG_IO_HDF5OUTPUT_HPP
#define FREHG_IO_HDF5OUTPUT_HPP

#include "core/Grid.hpp"
#include "core/Types.hpp"

#include <hdf5.h>
#include <mpi.h>

#include <string>
#include <vector>

namespace frehg::io {

/// Descriptive attributes attached to every field dataset (plan §7).
struct VarMeta {
  std::string units;     ///< e.g. "m", "m s-1", "psu"
  std::string longName;  ///< e.g. "free-surface elevation"
};

/// One parallel HDF5 output file.
class Hdf5Output {
 public:
  /// Create (truncate) the output file collectively and write the root
  /// provenance attributes. Parent directories are created if needed.
  /// \param grid decomposition the file layout follows (kept by reference).
  /// \param filename output path.
  /// \param configText full configuration text embedded as an attribute.
  Hdf5Output(const Grid& grid, const std::string& filename, const std::string& configText);

  /// Open an existing file read/write (used by restart and tests).
  Hdf5Output(const Grid& grid, const std::string& filename);

  /// Closes the file.
  ~Hdf5Output();

  Hdf5Output(const Hdf5Output&) = delete;
  Hdf5Output& operator=(const Hdf5Output&) = delete;

  /// Write the static /grid group: cell-center coordinates, layer
  /// thicknesses, bottom elevation, ktop, and /groundwater/zcell/0 (per-cell
  /// centers tiled from the uniform layer table).
  /// \param bottomInterior bed elevation, host view sized (nyLocal, nxLocal).
  void writeGridMeta(const HostField2<real_t>& bottomInterior);

  /// \copydoc writeGridMeta
  /// \param zcellInterior per-cell center elevations sized
  ///        (nyLocal, nxLocal, nz) — the geometry-aware zcell the
  ///        groundwater mesh provides (partial cells, terrain following);
  ///        b3's analysis script reads the (j*NX+i)*NZ+k layout (plan §7).
  void writeGridMeta(const HostField2<real_t>& bottomInterior,
                     const HostField3<real_t>& zcellInterior);

  /// Write one 2D field snapshot to /\<group\>/\<var\>/\<t\>. \p field includes
  /// halos; only the interior is written. Cells of fully inactive columns
  /// (ktop == nz) are written as NaN.
  void writeField2(const std::string& group, const std::string& var, real_t t,
                   const Field2<real_t>& field, const VarMeta& meta);

  /// Write one 3D field snapshot to /\<group\>/\<var\>/\<t\>; cells above ktop are
  /// written as NaN.
  void writeField3(const std::string& group, const std::string& var, real_t t,
                   const Field3<real_t>& field, const VarMeta& meta);

  /// \name Distributed raw dataset access (shared with Checkpoint/Monitor)
  ///@{
  /// Collectively write a [NY*NX] dataset from per-rank interior blocks.
  void writeDistributed2(const std::string& datasetPath, const HostField2<real_t>& interior);
  /// Collectively write a [NY*NX*NZ] dataset from per-rank interior blocks.
  void writeDistributed3(const std::string& datasetPath, const HostField3<real_t>& interior);
  /// Collectively read this rank's block of a [NY*NX] dataset.
  void readDistributed2(const std::string& datasetPath, HostField2<real_t>& interior) const;
  /// Collectively read this rank's block of a [NY*NX*NZ] dataset.
  void readDistributed3(const std::string& datasetPath, HostField3<real_t>& interior) const;
  ///@}

  /// \name Small helpers used across the io component and by tests
  ///@{
  /// \return true if a dataset or group exists at \p path.
  bool exists(const std::string& path) const;
  /// Read an entire 1D dataset, replicated on every rank.
  std::vector<real_t> readAll1D(const std::string& path) const;
  /// Attach / read scalar and string attributes on an existing object.
  void writeStringAttribute(const std::string& objectPath, const std::string& name,
                            const std::string& value);
  /// \copydoc writeStringAttribute
  void writeDoubleAttribute(const std::string& objectPath, const std::string& name, double value);
  /// \copydoc writeStringAttribute
  void writeLongAttribute(const std::string& objectPath, const std::string& name, long value);
  /// Read back attributes written by the setters above.
  std::string readStringAttribute(const std::string& objectPath, const std::string& name) const;
  /// \copydoc readStringAttribute
  double readDoubleAttribute(const std::string& objectPath, const std::string& name) const;
  /// \copydoc readStringAttribute
  long readLongAttribute(const std::string& objectPath, const std::string& name) const;
  ///@}

  /// Ensure all groups along \p path exist (collective).
  void ensureGroup(const std::string& path);

  /// Flush the file to disk (collective).
  void flush();

  /// The §7 time key: str(int(round(t_seconds))).
  static std::string timeKey(real_t t);

  /// \return the underlying HDF5 file handle (io-component collaborators).
  hid_t fileId() const { return file_; }

  /// \return the grid this file is laid out for.
  const Grid& grid() const { return grid_; }

 private:
  void writeRootAttributes(const std::string& configText);
  hid_t createDataset1D(const std::string& datasetPath, hsize_t globalSize);

  const Grid& grid_;
  std::string filename_;
  hid_t file_ = H5I_INVALID_HID;
};

}  // namespace frehg::io

#endif  // FREHG_IO_HDF5OUTPUT_HPP
