#pragma once

#include <array>
#include <cstdint>
#include <filesystem>
#include <string>

namespace astraxis {

// Marzok et al. 2022, Fig. 5: zenith OI 135.6 nm brightness in Rayleigh.
// Cell centers: west longitude 2.25..357.75, latitude -87.75..87.75 degrees.
struct GanymedeAurora {
    static constexpr int kWidth = 80;
    static constexpr int kHeight = 40;
    static constexpr double kTopKm = 70.0; // paper section 2.3, uniform column
    static constexpr double kTextureRange = 200.0;
    std::array<double, kWidth * kHeight> rayleigh{};

    bool load(const std::filesystem::path& file, std::string* error = nullptr);
    double brightness(double east_longitude, double latitude) const; // radians
    std::array<uint8_t, kWidth * kHeight * 4> texture_rgba() const;
};

} // namespace astraxis
