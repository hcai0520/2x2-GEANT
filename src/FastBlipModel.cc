#include "FastBlipModel.hh"
#include "G4Exception.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {
[[noreturn]] void Fail(const std::string& message) {
    G4Exception("FastBlipModel", "FastBlip001", FatalException, message.c_str());
    throw std::runtime_error(message); // Also stop if an exception handler returns.
}
}

const FastBlipModel& FastBlipModel::Instance() {
    static const FastBlipModel model;
    return model;
}

FastBlipModel::FastBlipModel() {
    const char* mode = std::getenv("MCP_FAST_BLIPS");
    if (!mode || std::string(mode) == "0") return;
    if (std::string(mode) != "1") Fail("MCP_FAST_BLIPS must be 0 or 1.");
    enabled = true;

    // Saved blip_counts.ipynb: N1/N0 at 83.28 fiducial active cm.
    // Threshold and clustering are already included in these counts.
    rate = (105251. / 49894650.) / (83.28 * cm);
    const char* profile = std::getenv("MCP_BLIP_PROFILE");
    if (!profile || !*profile) {
        G4cout << "Fast blips: constant calibrated rate = " << rate * cm
               << " /cm (N1/N0 = 105251/49894650 at 83.28 cm)." << G4endl;
        return;
    }

    // CSV: cumulative active length in cm, P(at least one blip by that length).
    // Convert the first-arrival CDF to integrated intensity, not a position CDF.
    std::ifstream input(profile);
    if (!input) Fail("Cannot open MCP_BLIP_PROFILE: " + std::string(profile));
    std::string line;
    G4int lineNumber = 0;
    while (std::getline(input, line)) {
        ++lineNumber;
        const auto comment = line.find('#');
        if (comment != std::string::npos) line.erase(comment);
        if (line.find_first_not_of(" \t\r") == std::string::npos) continue;
        if (std::count(line.begin(), line.end(), ',') != 1)
            Fail("Profile requires two comma-separated numbers at line " + std::to_string(lineNumber));
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream row(line);
        G4double length, probability;
        std::string extra;
        if (!(row >> length >> probability) || (row >> extra)
            || !std::isfinite(length) || !std::isfinite(probability)
            || length < 0. || probability < 0. || probability >= 1.)
            Fail("Invalid profile row at line " + std::to_string(lineNumber));
        const G4double mean = -std::log1p(-probability);
        if (lengths.empty()) {
            if (length != 0. || probability != 0.) Fail("Profile must start at 0,0.");
        } else if (length * cm <= lengths.back() || mean < means.back()) {
            Fail("Profile lengths must increase and probabilities must not decrease.");
        }
        lengths.push_back(length * cm);
        means.push_back(mean);
    }
    if (lengths.size() < 2 || means.back() <= 0.)
        Fail("Profile needs at least two points and a positive final probability.");
    G4cout << "Fast blips: cumulative profile " << profile << "; support 0 to "
           << lengths.back() / cm << " active cm. No extrapolation." << G4endl;
}

void FastBlipModel::AddStep(const G4ThreeVector& pre, const G4ThreeVector& post,
                           std::vector<Segment>& path) const {
    // Caller has established that the pre-step volume is an active LAr pixel.
    // Clip its chord to the same fiducial bounds used by filter.py.
    const G4ThreeVector delta = post - pre;
    const G4double chord = delta.mag();
    if (chord <= 0.) return;
    std::vector<Segment> pieces;
    for (G4int slab = 0; slab < 2; ++slab) {
        const G4double low[] = {-63.931 * cm, -51.85 * cm,
                                (slab == 0 ? -54.32 : 12.68) * cm};
        const G4double high[] = {63.931 * cm, 51.85 * cm,
                                 (slab == 0 ? -12.68 : 54.32) * cm};
        G4double enter = 0., exit = 1.;
        for (G4int axis = 0; axis < 3 && enter < exit; ++axis) {
            if (delta[axis] == 0.) {
                if (pre[axis] < low[axis] || pre[axis] > high[axis]) exit = -1.;
            } else {
                G4double a = (low[axis] - pre[axis]) / delta[axis];
                G4double b = (high[axis] - pre[axis]) / delta[axis];
                if (a > b) std::swap(a, b);
                enter = std::max(enter, a);
                exit = std::min(exit, b);
            }
        }
        if (exit > enter)
            pieces.push_back({pre + enter * delta, pre + exit * delta,
                              (exit - enter) * chord});
    }
    // Preserve traversal order, including reverse-going tracks.
    std::sort(pieces.begin(), pieces.end(), [&pre](const Segment& a, const Segment& b) {
        return (a.start - pre).mag2() < (b.start - pre).mag2();
    });
    path.insert(path.end(), pieces.begin(), pieces.end());
}

G4double FastBlipModel::Mean(G4double length) const {
    if (lengths.empty()) return rate * length;
    // Only roundoff at the endpoint is tolerated. Never stretch a measured CDF.
    if (length > lengths.back() + 1.e-9 * cm)
        Fail("Track exceeds MCP_BLIP_PROFILE active-length coverage.");
    if (length >= lengths.back()) return means.back();
    const auto hi = std::upper_bound(lengths.begin(), lengths.end(), length) - lengths.begin();
    const auto lo = hi - 1;
    return means[lo] + (means[hi] - means[lo]) *
           (length - lengths[lo]) / (lengths[hi] - lengths[lo]);
}

G4double FastBlipModel::InverseMean(G4double mean) const {
    if (lengths.empty()) return mean / rate;
    // upper_bound skips zero-intensity plateaus.
    const auto it = std::upper_bound(means.begin(), means.end(), mean);
    if (it == means.end()) return lengths.back();
    const auto hi = it - means.begin();
    const auto lo = hi - 1;
    return lengths[lo] + (lengths[hi] - lengths[lo]) *
           (mean - means[lo]) / (means[hi] - means[lo]);
}

FastBlipModel::Pair FastBlipModel::Sample(const std::vector<Segment>& path) const {
    Pair result;
    for (const auto& segment : path) result.length += segment.length;
    if (result.length <= 0.) return result;
    result.mean = Mean(result.length);
    if (result.mean <= 0.) return result;
    result.weight = std::exp(-result.mean) * result.mean * result.mean / 2.;
    result.s1 = std::min(result.length, InverseMean(G4UniformRand() * result.mean));
    result.s2 = std::min(result.length, InverseMean(G4UniformRand() * result.mean));
    if (result.s1 > result.s2) std::swap(result.s1, result.s2);
    const auto position = [&path](G4double s) {
        for (const auto& segment : path) {
            if (s <= segment.length)
                return segment.start + (s / segment.length) * (segment.end - segment.start);
            s -= segment.length;
        }
        return path.back().end; // Accumulation roundoff only.
    };
    result.first = position(result.s1);
    result.second = position(result.s2);
    result.valid = true;
    return result;
}
