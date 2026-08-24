/// \file GridDataReader.hpp
/// \brief Text readers for gridded input data (rasters, 3D fields, id maps).
///
/// Two formats are accepted, auto-detected by the first non-comment token:
///
///  1. **Headered raster**: `key value` header lines with keys ncols, nrows,
///     xllcorner, yllcorner, cellsize, nodata_value (any subset, ncols/nrows
///     required), followed by nrows * ncols whitespace-separated values.
///     Data rows map to grid rows in file order: the first data row is
///     j = 0. This matches how the legacy benchmark inputs were consumed and
///     is therefore the fidelity-preserving convention; it is *not* the
///     ESRI north-up convention. nodata values are returned as NaN.
///
///  2. **Flat list**: bare whitespace-separated values, `#` comments allowed.
///     2D fields are ordered j*nx + i, 3D fields (j*nx + i)*nz + k — the
///     same flattening as the legacy ASCII output and the §7 HDF5 contract.
///
/// The value count must match the expected extents exactly; mismatches are
/// fatal with the observed and expected counts in the message.

#ifndef FREHG_IO_GRIDDATAREADER_HPP
#define FREHG_IO_GRIDDATAREADER_HPP

#include "core/Types.hpp"

#include <string>
#include <vector>

namespace frehg::io {

/// Read a 2D real raster of exactly \p nx by \p ny values.
/// \return values ordered j*nx + i.
std::vector<real_t> readRaster2D(const std::string& path, int nx, int ny);

/// Read a 2D integer raster (e.g. a soil-id map) of exactly \p nx by \p ny
/// values. Fractional values in the file are fatal.
std::vector<int> readIntRaster2D(const std::string& path, int nx, int ny);

/// Read a flat 3D field of exactly \p nx * \p ny * \p nz values ordered
/// (j*nx + i)*nz + k. Header rasters are not meaningful in 3D and are
/// rejected.
std::vector<real_t> readField3D(const std::string& path, int nx, int ny, int nz);

}  // namespace frehg::io

#endif  // FREHG_IO_GRIDDATAREADER_HPP
