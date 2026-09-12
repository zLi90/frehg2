/// \file RunRecord.cpp
/// \brief Implementation of the persistent run record (v2 plan §2A).

#include "io/RunRecord.hpp"

#include "bc/BoundarySet.hpp"
#include "core/Logger.hpp"
#include "core/Timer.hpp"

#include <Kokkos_Core.hpp>
#include <yaml-cpp/yaml.h>

#include <unistd.h>

#include <array>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <iomanip>
#include <filesystem>
#include <fstream>
#include <sstream>

namespace frehg::io {

namespace {

// ---------------------------------------------------------------------------
// SHA-256 (FIPS 180-4), self-contained so the input-file hash needs no
// external dependency. Streaming form; used once per run on the input YAML.
// ---------------------------------------------------------------------------

struct Sha256 {
  std::array<std::uint32_t, 8> h{0x6a09e667u, 0xbb67ae85u, 0x3c6ef372u, 0xa54ff53au,
                                 0x510e527fu, 0x9b05688cu, 0x1f83d9abu, 0x5be0cd19u};
  std::array<unsigned char, 64> block{};
  std::size_t blockLen = 0;
  std::uint64_t totalBits = 0;

  static std::uint32_t rotr(std::uint32_t x, int n) { return (x >> n) | (x << (32 - n)); }

  void compress(const unsigned char* p) {
    static constexpr std::array<std::uint32_t, 64> k{
        0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u, 0x3956c25bu, 0x59f111f1u,
        0x923f82a4u, 0xab1c5ed5u, 0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u,
        0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u, 0xc19bf174u, 0xe49b69c1u, 0xefbe4786u,
        0x0fc19dc6u, 0x240ca1ccu, 0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau,
        0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u, 0xc6e00bf3u, 0xd5a79147u,
        0x06ca6351u, 0x14292967u, 0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u,
        0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u, 0xa2bfe8a1u, 0xa81a664bu,
        0xc24b8b70u, 0xc76c51a3u, 0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u,
        0x19a4c116u, 0x1e376c08u, 0x2748774cu, 0x34b0bcb5u, 0x391c0cb3u, 0x4ed8aa4au,
        0x5b9cca4fu, 0x682e6ff3u, 0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u,
        0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u};
    std::array<std::uint32_t, 64> w{};
    for (int i = 0; i < 16; ++i) {
      w[static_cast<std::size_t>(i)] =
          (static_cast<std::uint32_t>(p[4 * i]) << 24) |
          (static_cast<std::uint32_t>(p[4 * i + 1]) << 16) |
          (static_cast<std::uint32_t>(p[4 * i + 2]) << 8) |
          static_cast<std::uint32_t>(p[4 * i + 3]);
    }
    for (std::size_t i = 16; i < 64; ++i) {
      const std::uint32_t s0 = rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^ (w[i - 15] >> 3);
      const std::uint32_t s1 = rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^ (w[i - 2] >> 10);
      w[i] = w[i - 16] + s0 + w[i - 7] + s1;
    }
    std::uint32_t a = h[0], b = h[1], c = h[2], d = h[3];
    std::uint32_t e = h[4], f = h[5], g = h[6], hh = h[7];
    for (std::size_t i = 0; i < 64; ++i) {
      const std::uint32_t s1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25);
      const std::uint32_t ch = (e & f) ^ (~e & g);
      const std::uint32_t t1 = hh + s1 + ch + k[i] + w[i];
      const std::uint32_t s0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22);
      const std::uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
      const std::uint32_t t2 = s0 + maj;
      hh = g;
      g = f;
      f = e;
      e = d + t1;
      d = c;
      c = b;
      b = a;
      a = t1 + t2;
    }
    h[0] += a; h[1] += b; h[2] += c; h[3] += d;
    h[4] += e; h[5] += f; h[6] += g; h[7] += hh;
  }

  void update(const unsigned char* data, std::size_t len) {
    totalBits += static_cast<std::uint64_t>(len) * 8;
    while (len > 0) {
      const std::size_t take = std::min(len, block.size() - blockLen);
      std::memcpy(block.data() + blockLen, data, take);
      blockLen += take;
      data += take;
      len -= take;
      if (blockLen == block.size()) {
        compress(block.data());
        blockLen = 0;
      }
    }
  }

  std::string finish() {
    const unsigned char one = 0x80;
    const std::uint64_t bits = totalBits;
    update(&one, 1);
    const unsigned char zero = 0x00;
    while (blockLen != 56) {
      update(&zero, 1);
      totalBits -= 8;  // padding does not count
    }
    unsigned char lenBytes[8];
    for (int i = 0; i < 8; ++i) {
      lenBytes[i] = static_cast<unsigned char>(bits >> (56 - 8 * i));
    }
    totalBits = bits;
    update(lenBytes, 8);
    std::ostringstream hex;
    hex << std::hex << std::setfill('0');
    for (const std::uint32_t word : h) {
      hex << std::setw(8) << word;
    }
    return hex.str();
  }
};

std::string sha256File(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    return "unreadable";
  }
  Sha256 sha;
  std::array<char, 65536> buf{};
  while (in.read(buf.data(), static_cast<std::streamsize>(buf.size())) || in.gcount() > 0) {
    sha.update(reinterpret_cast<const unsigned char*>(buf.data()),
               static_cast<std::size_t>(in.gcount()));
  }
  return sha.finish();
}

std::string isoNow() {
  const std::time_t now = std::time(nullptr);
  std::tm utc{};
  gmtime_r(&now, &utc);
  char buf[32];
  std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%SZ", &utc);
  return std::string(buf);
}

std::string hostName() {
  char buf[256] = {0};
  if (gethostname(buf, sizeof(buf) - 1) != 0) {
    return "unknown";
  }
  return std::string(buf);
}

int ompThreads() {
  const char* env = std::getenv("OMP_NUM_THREADS");
  if (env != nullptr && env[0] != '\0') {
    return std::atoi(env);
  }
  return static_cast<int>(Kokkos::DefaultHostExecutionSpace().concurrency());
}

const char* targetName(BcTarget t) {
  switch (t) {
    case BcTarget::Surface: return "surface";
    case BcTarget::GroundwaterTop: return "groundwater_top";
    case BcTarget::GroundwaterBottom: return "groundwater_bottom";
    default: return "groundwater_side";
  }
}

const char* kindName(BcKind k) {
  switch (k) {
    case BcKind::Eta: return "eta";
    case BcKind::Discharge: return "discharge";
    case BcKind::Velocity: return "velocity";
    case BcKind::Outflow: return "outflow";
    case BcKind::Head: return "head";
    case BcKind::Flux: return "flux";
    default: return "scalar_value";
  }
}

}  // namespace

RunRecord::RunRecord(MPI_Comm comm, const FrehgConfig& config, const std::string& inputPath,
                     const BoundarySet& boundaries, int px, int py)
    : comm_(comm), began_(std::chrono::steady_clock::now()) {
  MPI_Comm_rank(comm_, &rank_);
  int size = 1;
  MPI_Comm_size(comm_, &size);

  const std::filesystem::path outDir =
      std::filesystem::path(config.output.filename).parent_path();
  path_ = (outDir / "run-record.yaml").string();

  if (rank_ != 0) {
    return;  // only rank 0 composes and writes
  }

  YAML::Node top;

  YAML::Node prov;
  prov["frehg_version"] = FREHG_VERSION;
  prov["git_sha"] = FREHG_GIT_SHA;
  prov["build_type"] = FREHG_BUILD_TYPE;
  prov["hostname"] = hostName();
  prov["mpi_ranks"] = size;
  YAML::Node decomp(YAML::NodeType::Sequence);
  decomp.SetStyle(YAML::EmitterStyle::Flow);
  decomp.push_back(px);
  decomp.push_back(py);
  prov["decomposition"] = decomp;
  prov["omp_threads"] = ompThreads();
  prov["kokkos_backend"] = Kokkos::DefaultExecutionSpace::name();
  prov["start_time"] = isoNow();
  prov["input_file"] = inputPath;
  prov["input_sha256"] = sha256File(inputPath);
  if (config.restart.enabled) {
    prov["restart_from"] = config.restart.file;
    prov["restart_time"] = config.restart.time;
  }
  prov["finished"] = false;
  top["provenance"] = prov;

  top["configuration"] = YAML::Load(resolvedConfigYaml(config));

  YAML::Node modules;
  modules["surface_water"] = config.modules.surfaceWater;
  modules["groundwater"] = config.modules.groundwater;
  modules["transport"] = config.modules.transport;
  modules["coupling_mode"] =
      (config.modules.surfaceWater && config.modules.groundwater)
          ? (config.coupling.mode == CouplingConfig::Mode::Subcycled ? "subcycled" : "sync")
          : "none";
  modules["density_coupling"] = config.groundwater.densityCoupling.enabled;
  top["modules"] = modules;

  YAML::Node bcs(YAML::NodeType::Sequence);
  for (std::size_t n = 0; n < boundaries.all().size(); ++n) {
    const BoundaryCondition& bc = boundaries.all()[n];
    const BoundaryConditionConfig& bcConfig = config.boundaryConditions[n];
    YAML::Node b;
    b["name"] = bc.name();
    b["target"] = targetName(bc.target());
    b["kind"] = kindName(bc.kind());
    YAML::Node value;
    switch (bc.valueForm()) {
      case BcValueConfig::Form::Constant:
        value["form"] = "constant";
        value["constant"] = bcConfig.value.constant;
        break;
      case BcValueConfig::Form::Series:
        value["form"] = "series";
        value["file"] = bcConfig.value.seriesFile;
        value["file_sha256"] = sha256File(config.resolvePath(bcConfig.value.seriesFile));
        break;
      case BcValueConfig::Form::Gravity:
        value["form"] = "gravity";
        break;
      case BcValueConfig::Form::Hydrostatic:
        value["form"] = "hydrostatic";
        value["eta"] = bcConfig.value.hydrostaticEta;
        break;
    }
    b["value"] = value;
    b["polygon_vertices"] = bcConfig.polygon.size();
    if (!bcConfig.polygon.empty()) {
      real_t xMin = bcConfig.polygon[0][0], xMax = xMin;
      real_t yMin = bcConfig.polygon[0][1], yMax = yMin;
      for (const std::array<real_t, 2>& pt : bcConfig.polygon) {
        xMin = std::min(xMin, pt[0]);
        xMax = std::max(xMax, pt[0]);
        yMin = std::min(yMin, pt[1]);
        yMax = std::max(yMax, pt[1]);
      }
      YAML::Node bbox(YAML::NodeType::Sequence);
      bbox.SetStyle(YAML::EmitterStyle::Flow);
      bbox.push_back(xMin);
      bbox.push_back(yMin);
      bbox.push_back(xMax);
      bbox.push_back(yMax);
      b["bounding_box"] = bbox;
    }
    b["global_cells"] = bc.globalCellCount();
    bcs.push_back(b);
  }
  top["boundary_conditions"] = bcs;

  YAML::Emitter emitter;
  emitter.SetDoublePrecision(17);
  emitter << top;
  staticText_ = emitter.c_str();
}

RunRecord::~RunRecord() = default;

void RunRecord::setSolver(const std::string& yamlText) { solverText_ = yamlText; }

void RunRecord::setClosure(const std::string& yamlText) { closureText_ = yamlText; }

void RunRecord::flush(bool finished) {
  Timer::Scoped timer("io/run_record");
  // Collective: every rank contributes its timer sections.
  const std::vector<TimerSection> sections = Timer::merged(comm_);
  if (rank_ != 0) {
    return;
  }

  YAML::Node top = YAML::Load(staticText_);
  top["provenance"]["finished"] = finished;
  top["provenance"]["end_time"] = isoNow();
  top["provenance"]["wall_seconds"] =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - began_).count();

  YAML::Node timers;
  for (const TimerSection& s : sections) {
    YAML::Node entry;
    entry["count"] = s.count;
    entry["min_s"] = s.minSeconds;
    entry["mean_s"] = s.meanSeconds;
    entry["max_s"] = s.maxSeconds;
    timers[s.path] = entry;
  }
  top["timers"] = timers;

  top["solver"] = solverText_.empty() ? YAML::Node(YAML::NodeType::Map)
                                      : YAML::Load(solverText_);
  top["closure"] = closureText_.empty() ? YAML::Node(YAML::NodeType::Map)
                                        : YAML::Load(closureText_);

  YAML::Emitter emitter;
  emitter.SetDoublePrecision(17);
  emitter << top;

  // Write-then-rename so a crash mid-write cannot truncate the previous
  // (still truthful) record.
  const std::string tmp = path_ + ".tmp";
  {
    std::ofstream out(tmp, std::ios::trunc);
    if (!out) {
      log::warn(log::msg() << "run record: cannot write '" << tmp << "'");
      return;
    }
    out << emitter.c_str() << "\n";
  }
  std::error_code ec;
  std::filesystem::rename(tmp, path_, ec);
  if (ec) {
    log::warn(log::msg() << "run record: rename to '" << path_ << "' failed: "
                         << ec.message());
  }
}

}  // namespace frehg::io
