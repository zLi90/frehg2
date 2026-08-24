/// \file GridDataReader.cpp
/// \brief Implementation of the gridded-input text readers.

#include "io/GridDataReader.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>

namespace frehg::io {

namespace {

/// Tokenized file content with any raster header decoded.
struct ParsedGridFile {
  bool hasHeader = false;
  long ncols = 0;
  long nrows = 0;
  bool hasNodata = false;
  real_t nodata = 0;
  std::vector<real_t> values;
};

bool isHeaderWord(const std::string& token) {
  return !token.empty() && (std::isalpha(static_cast<unsigned char>(token[0])) != 0 ||
                            token[0] == '_');
}

ParsedGridFile parseGridFile(const std::string& path) {
  std::ifstream in(path);
  if (!in) {
    log::fatal(log::msg() << "GridDataReader: cannot open '" << path << "'");
  }

  ParsedGridFile out;
  std::string line;
  bool inData = false;
  std::size_t lineNo = 0;
  while (std::getline(in, line)) {
    ++lineNo;
    const std::size_t hash = line.find('#');
    if (hash != std::string::npos) {
      line.erase(hash);
    }
    std::istringstream ls(line);
    std::string token;
    while (ls >> token) {
      if (!inData && isHeaderWord(token)) {
        std::string lowered = token;
        std::transform(lowered.begin(), lowered.end(), lowered.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        real_t headerValue = 0;
        if (!(ls >> headerValue)) {
          log::fatal(log::msg() << "GridDataReader: '" << path << "' line " << lineNo
                                << ": header key '" << token << "' has no value");
        }
        out.hasHeader = true;
        if (lowered == "ncols" || lowered == "n_cols") {
          out.ncols = static_cast<long>(headerValue);
        } else if (lowered == "nrows" || lowered == "n_rows") {
          out.nrows = static_cast<long>(headerValue);
        } else if (lowered == "nodata_value") {
          out.hasNodata = true;
          out.nodata = headerValue;
        } else if (lowered == "xllcorner" || lowered == "yllcorner" || lowered == "cellsize") {
          // Georeferencing is carried by the YAML domain section; these
          // header entries are accepted and ignored.
        } else {
          log::fatal(log::msg() << "GridDataReader: '" << path << "' line " << lineNo
                                << ": unknown header key '" << token << "'");
        }
        continue;
      }
      inData = true;
      try {
        std::size_t used = 0;
        const real_t value = std::stod(token, &used);
        if (used != token.size()) {
          log::fatal(log::msg() << "GridDataReader: '" << path << "' line " << lineNo
                                << ": bad number '" << token << "'");
        }
        out.values.push_back(value);
      } catch (const std::invalid_argument&) {
        log::fatal(log::msg() << "GridDataReader: '" << path << "' line " << lineNo
                              << ": bad number '" << token << "'");
      } catch (const std::out_of_range&) {
        log::fatal(log::msg() << "GridDataReader: '" << path << "' line " << lineNo
                              << ": number out of range '" << token << "'");
      }
    }
  }
  return out;
}

void applyNodata(ParsedGridFile& parsed) {
  if (!parsed.hasNodata) {
    return;
  }
  for (real_t& v : parsed.values) {
    if (v == parsed.nodata) {
      v = std::numeric_limits<real_t>::quiet_NaN();
    }
  }
}

}  // namespace

std::vector<real_t> readRaster2D(const std::string& path, int nx, int ny) {
  ParsedGridFile parsed = parseGridFile(path);
  const std::size_t expected = static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny);
  if (parsed.hasHeader && parsed.ncols > 0 && parsed.nrows > 0 &&
      (parsed.ncols != nx || parsed.nrows != ny)) {
    log::fatal(log::msg() << "GridDataReader: '" << path << "' header says " << parsed.ncols
                          << " x " << parsed.nrows << " but the domain is " << nx << " x " << ny);
  }
  if (parsed.values.size() != expected) {
    log::fatal(log::msg() << "GridDataReader: '" << path << "' has " << parsed.values.size()
                          << " values, expected " << expected << " (" << nx << " x " << ny
                          << ")");
  }
  applyNodata(parsed);
  return std::move(parsed.values);
}

std::vector<int> readIntRaster2D(const std::string& path, int nx, int ny) {
  const std::vector<real_t> raw = readRaster2D(path, nx, ny);
  std::vector<int> out(raw.size());
  for (std::size_t n = 0; n < raw.size(); ++n) {
    const real_t v = raw[n];
    if (std::isnan(v) || v != std::floor(v)) {
      log::fatal(log::msg() << "GridDataReader: '" << path << "' entry " << n << " (" << v
                            << ") is not an integer id");
    }
    out[n] = static_cast<int>(v);
  }
  return out;
}

std::vector<real_t> readField3D(const std::string& path, int nx, int ny, int nz) {
  ParsedGridFile parsed = parseGridFile(path);
  if (parsed.hasHeader) {
    log::fatal(log::msg() << "GridDataReader: '" << path
                          << "': raster headers are 2D-only; 3D fields are flat lists ordered "
                             "(j*nx + i)*nz + k");
  }
  const std::size_t expected =
      static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny) * static_cast<std::size_t>(nz);
  if (parsed.values.size() != expected) {
    log::fatal(log::msg() << "GridDataReader: '" << path << "' has " << parsed.values.size()
                          << " values, expected " << expected << " (" << nx << " x " << ny
                          << " x " << nz << ")");
  }
  return std::move(parsed.values);
}

}  // namespace frehg::io
