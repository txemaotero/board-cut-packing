#include <algorithm>
#include <array>
#include <atomic>
#include <execution>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <concepts>
#include <format>
#include <iostream>
#include <numeric>
#include <optional>
#include <print>
#include <string_view>
#include <utility>
#include <vector>
#include <ranges>
#include "./cuts.h"

namespace rg = std::ranges;
namespace rgv = std::ranges::views;
using uint = unsigned;

template<std::integral T>
bool isBetween(T number, T m, T M)
{
    return m <= number && number <= M;
}

struct CCOA
{
    size_t rectangleIdx;
    int x;
    int y;
    bool rotated;

    bool operator==(const CCOA& other) const = default;
    bool operator!=(const CCOA& other) const = default;
};

struct Rectangle
{
    uint width; // x dimension
    uint height; // y dimension
    int x{0};
    int y{0};

    int xmax() const
    {
        return x + width;
    }

    int ymax() const
    {
        return y + height;
    }

    void rotate()
    {
        std::swap(width, height);
    }

    void applyCCOA(const CCOA& ccoa)
    {
        x = ccoa.x;
        y = ccoa.y;
        if (ccoa.rotated)
            rotate();
    }

    uint area() const
    {
        return width * height;
    }

    bool overlaps(const Rectangle& other) const
    {
        return !(xmax() <= other.x || x >= other.xmax() || ymax() <= other.y || y >= other.ymax());
    }

    std::optional<float> distance(const Rectangle& other) const
    {
        return distanceSq(other).transform(std::sqrtf);
    }

    std::optional<uint> distanceSq(const Rectangle& other) const
    {
        if (overlaps(other))
            return {};

        int dx = std::max(x - other.xmax(), other.x - xmax());
        dx = std::max(dx, 0);
        int dy = std::max(y - other.ymax(), other.y - ymax());
        dy = std::max(dy, 0);
        return dx * dx + dy * dy;
    }

    std::array<uint, 4> distanceWithContainer(const Rectangle& container) const
    {
        // assert(container.contains(*this));
        return {static_cast<uint>(x - container.x),
                static_cast<uint>(container.xmax() - xmax()),
                static_cast<uint>(y - container.y),
                static_cast<uint>(container.ymax() - ymax())};
    }

    bool contains(const Rectangle& other) const
    {
        return (isBetween(other.x, x, xmax()) && isBetween(other.xmax(), x, xmax()) &&
                isBetween(other.y, y, ymax()) && isBetween(other.ymax(), y, ymax()));
    }
};

struct Configuration
{
    Rectangle container;
    std::vector<Rectangle> allRectangles;
    std::vector<size_t> packedRectangles;
    std::vector<size_t> unpackedRectangles;

    Configuration(const Rectangle& cont, const std::vector<Rectangle>& toPack):
        container{cont},
        allRectangles(toPack),
        unpackedRectangles(allRectangles.size())
    {
        std::iota(unpackedRectangles.begin(), unpackedRectangles.end(), 0);
    }

    bool isSuccessful() const
    {
        return unpackedRectangles.size() == 0;
    }

    float density() const
    {
        auto elementsAreas = packedRectangles | rgv::transform(
                                                    [this](const size_t id)
                                                    {
                                                        return allRectangles[id].area();
                                                    });
        return static_cast<float>(rg::fold_left(elementsAreas, 0, std::plus{})) / container.area();
    }

    uint noNormDensity() const
    {
        auto elementsAreas = packedRectangles | rgv::transform(
                                                    [this](const size_t id)
                                                    {
                                                        return allRectangles[id].area();
                                                    });
        return rg::fold_left(elementsAreas, 0, std::plus{});
    }

    /**
     * Returns the min distance between the rectangle and all the rectangles packed in the
     * configuration and container excluding the two edges that are shared
     */
    uint minDistance(const Rectangle& r) const
    {
        uint nZeros = 0;
        uint minDistSq = std::numeric_limits<uint>::max();
        // First with container
        for (const uint d: r.distanceWithContainer(container))
        {
            if (d == 0 && ++nZeros > 2)
            {
                return 0;
            }
            else
            {
                minDistSq = std::min(minDistSq, d * d);
            }
        }
        // Same with other rectangles
        for (const auto& packedRectsIdx: packedRectangles)
        {
            auto dOpt = r.distanceSq(allRectangles[packedRectsIdx]);
            if (!dOpt)
            {
                std::println("Trying to evaluate distance with overlapping rectangles");
                assert(false);
                return 0;
            }
            if (uint dSq = *dOpt; dSq == 0 && ++nZeros > 2)
            {
                return 0;
            }
            else
            {
                minDistSq = std::min(minDistSq, dSq);
            }
        }
        return std::sqrtf(minDistSq);
    }

    bool isInfeasible(const size_t placedRectIdx, const CCOA& ccoa) const
    {
        if (placedRectIdx == ccoa.rectangleIdx)
        {
            return true;
        }
        const Rectangle& placedRect = allRectangles[placedRectIdx];
        Rectangle candidate = allRectangles[ccoa.rectangleIdx];
        candidate.applyCCOA(ccoa);
        return candidate.overlaps(placedRect);
    }

    void write(const std::filesystem::path& filename) const
    {
        std::ofstream file(filename);

        if (!file.is_open())
        {
            std::cerr << "Failed to open file: " << filename << std::endl;
            return;
        }

        // Write the header
        file << "x;y;width;height\n";
        file << std::format("{};{};{};{}\n", container.x, container.y, container.width, container.height);

        // Write the data
        for (const auto& rectId: packedRectangles)
        {
            const auto& rect = allRectangles[rectId];
            file << std::format("{};{};{};{}\n", rect.x, rect.y, rect.width, rect.height);
        }

        file.close();
        if (!file)
        {
            std::cerr << "Failed to close file properly: " << filename << std::endl;
        }
    }
};

std::vector<CCOA> calculateInitialCCOAs(const Configuration& config)
{
    std::vector<CCOA> result;
    result.reserve(config.packedRectangles.size() * 4);
    auto processUnpacked = [&result, &config](const int rectIdx)
    {
        auto checkAndPush =
            [&result, &container = config.container, rectIdx](const Rectangle& r, bool rot)
        {
            if (container.contains(r))
                result.emplace_back(rectIdx, r.x, r.y, rot);
        };
        auto check4pos =
            [&checkAndPush, &container = config.container](Rectangle& r, bool rot)
        {
            r.x = container.x;
            r.y = container.y;
            checkAndPush(r, rot);
            r.y = container.ymax() - r.height;
            checkAndPush(r, rot);
            r.x = container.xmax() - r.width;
            checkAndPush(r, rot);
            r.y = container.y;
            checkAndPush(r, rot);
        };
        Rectangle r = config.allRectangles[rectIdx];
        check4pos(r, false);
        r.rotate();
        check4pos(r, true);
    };
    rg::for_each(config.unpackedRectangles, processUnpacked);
    return result;
}

float degree(const CCOA& ccoa, const Configuration& config)
{
    Rectangle r = config.allRectangles[ccoa.rectangleIdx];
    r.applyCCOA(ccoa);
    return 1. - config.minDistance(r) * 2. / (r.width + r.height);
}

struct Coordinate
{
    int x, y;
};

void addCandidateIfPossible(const Configuration& config,
                            const Rectangle& candidate,
                            std::vector<Coordinate>& candidateCoords)
{
    bool anyOverlap = rg::any_of(
        config.packedRectangles,
        [&candidate, &allRects = std::as_const(config.allRectangles)](const size_t packedIdx)
        {
            const Rectangle& packed = allRects[packedIdx];
            return packed.overlaps(candidate);
        });
    if (!anyOverlap)
    {
        candidateCoords.emplace_back(candidate.x, candidate.y);
    }
}

struct DimLimit
{
    uint min, max;

    bool contains(uint v) const
    {
        return min <= v && v <= max;
    }
};


enum class EFixedCornerPos
{
    bottomLeft,
    bottomRight,
    topRight,
    topLeft,
};

struct CandidateLimitFixed
{
    Coordinate c;
    EFixedCornerPos cPos;
    DimLimit widthLims;
    DimLimit heigthLims;

    Coordinate getPlacedRectanglePosition(const Rectangle& r) const
    {
        switch (cPos)
        {
        case EFixedCornerPos::bottomLeft:
            return c;
        case EFixedCornerPos::bottomRight:
            return {c.x - static_cast<int>(r.width), c.y};
        case EFixedCornerPos::topRight:
            return {c.x - static_cast<int>(r.width), c.y - static_cast<int>(r.height)};
        case EFixedCornerPos::topLeft:
            return {c.x, c.y - static_cast<int>(r.height)};
        }
        assert(false);
        return c;
    }

    bool isGoodCandidate(const Rectangle& r) const
    {
        return widthLims.contains(r.width) && heigthLims.contains(r.height);
    }
};

class PeripheralCandidateGenerator
{
public:
    PeripheralCandidateGenerator(const Rectangle& placed,
                                  const Rectangle& peripheral,
                                  const Rectangle& container)
    {
        // We are on top of placed. Do we see the peripheral?
        if (peripheral.ymax() > placed.ymax())
        {
            if (peripheral.xmax() < placed.xmax())
            {
                candLims.push_back({
                    .c = {peripheral.xmax(), placed.ymax()},
                    .cPos = EFixedCornerPos::bottomLeft,
                    .widthLims = {static_cast<uint>(std::max(0, placed.x - peripheral.xmax())),
                          static_cast<uint>(container.xmax() - peripheral.xmax())},
                    .heigthLims = {static_cast<uint>(std::max(0, peripheral.y - placed.ymax())),
                          static_cast<uint>(container.ymax() - placed.ymax())},
                });
            }
            if (peripheral.x > placed.x)
            {
                candLims.push_back({
                    .c = {peripheral.x, placed.ymax()},
                    .cPos = EFixedCornerPos::bottomRight,
                    .widthLims = {static_cast<uint>(std::max(0, peripheral.x - placed.xmax())),
                          static_cast<uint>(peripheral.x - container.x)},
                    .heigthLims = {static_cast<uint>(std::max(0, peripheral.y - placed.ymax())),
                          static_cast<uint>(container.ymax() - placed.ymax())},
                });
            }
        }
        // We are below the placed. Do we see the peripheral?
        if (peripheral.y < placed.y)
        {
            if (peripheral.xmax() < placed.xmax())
            {
                candLims.push_back({
                    .c = {peripheral.xmax(), placed.y},
                    .cPos = EFixedCornerPos::topLeft,
                    .widthLims = {static_cast<uint>(std::max(0, placed.x - peripheral.xmax())),
                          static_cast<uint>(container.xmax() - peripheral.xmax())},
                    .heigthLims = {static_cast<uint>(std::max(0, placed.y - peripheral.ymax())),
                          static_cast<uint>(placed.y - container.y)},
                });
            }
            if (peripheral.x > placed.x)
            {
                candLims.push_back({
                    .c = {peripheral.x, placed.y},
                    .cPos = EFixedCornerPos::topRight,
                    .widthLims = {static_cast<uint>(std::max(0, peripheral.x - placed.xmax())),
                          static_cast<uint>(peripheral.x - container.x)},
                    .heigthLims = {static_cast<uint>(std::max(0, placed.y - peripheral.ymax())),
                          static_cast<uint>(placed.y - container.y)},
                });
            }
        }

        // We are right to the placed
        if (peripheral.xmax() > placed.xmax())
        {
            if (peripheral.ymax() < placed.ymax())
            {
                candLims.push_back({
                    .c = {placed.xmax(), peripheral.ymax()},
                    .cPos = EFixedCornerPos::bottomLeft,
                    .widthLims = {static_cast<uint>(std::max(0, peripheral.x - placed.xmax())),
                          static_cast<uint>(container.xmax() - placed.xmax())},
                    .heigthLims = {static_cast<uint>(std::max(0, placed.y - peripheral.ymax())),
                          static_cast<uint>(container.ymax() - peripheral.ymax())},
                });
            }
            if (peripheral.y > placed.y)
            {
                candLims.push_back({
                    .c = {placed.xmax(), peripheral.y},
                    .cPos = EFixedCornerPos::topLeft,
                    .widthLims = {static_cast<uint>(std::max(0, peripheral.x - placed.xmax())),
                          static_cast<uint>(container.xmax() - placed.xmax())},
                    .heigthLims = {static_cast<uint>(std::max(0, peripheral.y - placed.ymax())),
                          static_cast<uint>(peripheral.y - container.y)},
                });
            }
        }

        // We are left to the placed
        if (peripheral.x < placed.x)
        {
            if (peripheral.ymax() < placed.ymax())
            {
                candLims.push_back({
                    .c = {placed.x, peripheral.ymax()},
                    .cPos = EFixedCornerPos::bottomRight,
                    .widthLims = {static_cast<uint>(std::max(0, placed.x - peripheral.xmax())),
                          static_cast<uint>(placed.x - container.x)},
                    .heigthLims = {static_cast<uint>(std::max(0, placed.y - peripheral.ymax())),
                          static_cast<uint>(container.ymax() - peripheral.ymax())},
                });
            }
            if (peripheral.y > placed.y)
            {
                candLims.push_back({
                    .c = {placed.x, peripheral.y},
                    .cPos = EFixedCornerPos::topRight,
                    .widthLims = {static_cast<uint>(std::max(0, placed.x - peripheral.xmax())),
                          static_cast<uint>(placed.x - container.x)},
                    .heigthLims = {static_cast<uint>(std::max(0, peripheral.y - placed.ymax())),
                          static_cast<uint>(peripheral.y - container.y)},
                });
            }
        }
    }

    std::vector<Coordinate> operator()(Rectangle r, const Configuration& config) const
    {
        std::vector<Coordinate> finalCandidates;
        for (const CandidateLimitFixed& cl : candLims)
        {
            if (cl.isGoodCandidate(r))
            {
                const auto [x, y] = cl.getPlacedRectanglePosition(r);
                r.x = x;
                r.y = y;
                addCandidateIfPossible(config, r, finalCandidates);
            }
        }
        return finalCandidates;
    }
private:
    std::vector<CandidateLimitFixed> candLims;
};

class ContainerCandidateGenerator
{
public:
    ContainerCandidateGenerator(const Rectangle& placed, const Rectangle& container)
    {
        const auto [dLeft, dRight, dBottom, dTop] = placed.distanceWithContainer(container);

        // We have potentially 8 candidate intervals
        // Above left
        CandidateLimitFixed al {
            .c = {container.x, placed.ymax()},
            .cPos = EFixedCornerPos::bottomLeft,
            .widthLims = {dLeft, container.width},
            .heigthLims = {0, dTop},
        };

        // Above right
        CandidateLimitFixed ar {
            .c = {container.xmax(), placed.ymax()},
            .cPos = EFixedCornerPos::bottomRight,
            .widthLims = {dRight, container.width},
            .heigthLims = {0, dTop},
        };

        // Bottom left
        CandidateLimitFixed bl {
            .c = {container.x, placed.y},
            .cPos = EFixedCornerPos::topLeft,
            .widthLims = {dLeft, container.width},
            .heigthLims = {0, dBottom},
        };

        // Bottom right
        CandidateLimitFixed br {
            .c = {container.xmax(), placed.y},
            .cPos = EFixedCornerPos::topRight,
            .widthLims = {dRight, container.width},
            .heigthLims = {0, dBottom},
        };

        // left top
        CandidateLimitFixed lt {
            .c = {placed.x, container.ymax()},
            .cPos = EFixedCornerPos::topRight,
            .widthLims = {0, dLeft},
            .heigthLims = {dTop, container.height},
        };

        // left bottom
        CandidateLimitFixed lb {
            .c = {placed.x, container.y},
            .cPos = EFixedCornerPos::bottomRight,
            .widthLims = {0, dLeft},
            .heigthLims = {dBottom, container.height},
        };

        // right top
        CandidateLimitFixed rt {
            .c = {placed.xmax(), container.ymax()},
            .cPos = EFixedCornerPos::topLeft,
            .widthLims = {0, dRight},
            .heigthLims = {dTop, container.height},
        };

        // right bottom
        CandidateLimitFixed rb {
            .c = {placed.xmax(), container.y},
            .cPos = EFixedCornerPos::bottomLeft,
            .widthLims = {0, dRight},
            .heigthLims = {dBottom, container.height},
        };

        candLims = {al, ar, bl, br, lt, lb, rt, rb};
    }

    std::vector<Coordinate> operator()(Rectangle r, const Configuration& config) const
    {
        std::vector<Coordinate> finalCandidates;
        for (const CandidateLimitFixed& cl : candLims)
        {
            if (cl.isGoodCandidate(r))
            {
                const auto [x, y] = cl.getPlacedRectanglePosition(r);
                r.x = x;
                r.y = y;
                addCandidateIfPossible(config, r, finalCandidates);
            }
        }
        return finalCandidates;
    }

private:
    std::array<CandidateLimitFixed, 8> candLims;
};

using CandidateGenerator = std::variant<PeripheralCandidateGenerator, ContainerCandidateGenerator>;

class PossibleNewPositionGenerator
{
public:
    PossibleNewPositionGenerator(const Configuration& configuration,
                                 const size_t lastPlacedRectIdx):
        config(configuration)
    {
        const Rectangle& lastPlaced = config.allRectangles[lastPlacedRectIdx];
        const Rectangle& container = config.container;
        auto packedRects = config.packedRectangles |
                           rgv::filter(
                               [lastPlacedRectIdx](auto i)
                               {
                                   return i != lastPlacedRectIdx;
                               }) |
                           rgv::transform(
                               [this](auto i) -> const auto&
                               {
                                   return config.allRectangles[i];
                               });
        for (const auto& rect: packedRects)
        {
            candGens.emplace_back(PeripheralCandidateGenerator(lastPlaced, rect, container));
        }
        candGens.emplace_back(ContainerCandidateGenerator(lastPlaced, container));
    }

    std::vector<Coordinate> operator()(const Rectangle& candidate) const
    {
        std::vector<Coordinate> result;
        for (const auto& cg : candGens)
        {
            const auto cands = std::visit(
                [this, &candidate](const auto& v)
                {
                    return v(candidate, config);
                },
                cg);
            result.insert(result.end(), cands.begin(), cands.end());
        }
        return result;
    }

private:
    const Configuration& config;
    std::vector<CandidateGenerator> candGens;
};

bool placedRectangleIsOk(const Configuration& config, const Rectangle& placed, size_t skip)
{
    uint nZeros = 0;
    for (const uint d: placed.distanceWithContainer(config.container))
    {
        if (d == 0 && ++nZeros == 2)
        {
            return true;
        }
    }
    // Same with other rectangles
    for (const auto& ri: config.packedRectangles)
    {
        if (ri == skip)
            continue;
        auto dOpt = placed.distanceSq(config.allRectangles[ri]);
        if (!dOpt)
            return false;
        if (*dOpt == 0 && ++nZeros == 2)
            return true;
    }
    return false;
}

bool placedRectangleIsOk(const Configuration& config, const size_t placedIDx)
{
    return placedRectangleIsOk(config, config.allRectangles[placedIDx], placedIDx);
}

bool ccoasAreOk(const Configuration& config, const std::vector<CCOA>& currentCCOAs)
{
    return rg::all_of(currentCCOAs,
                      [&config](const auto& ccoa)
                      {
                          Rectangle r = config.allRectangles[ccoa.rectangleIdx];
                          r.applyCCOA(ccoa);
                          bool isOK = placedRectangleIsOk(config, r, ccoa.rectangleIdx);
                          if (!isOK)
                              return false;
                          return isOK;
                      });
}

std::vector<CCOA> addNewPossibleCCOAs(Configuration& config,
                                      const size_t lastPlacedRectIdx,
                                      std::vector<CCOA>&& currentCCOAs)
{
    const PossibleNewPositionGenerator gen(config, lastPlacedRectIdx);
    auto processGen =
        [&gen, &currentCCOAs, &config = std::as_const(config)](const size_t candidateIdx)
    {
        Rectangle candidate = config.allRectangles[candidateIdx];
        for (auto& [x, y]: gen(candidate))
        {
            currentCCOAs.emplace_back(candidateIdx, x, y, false);
        }
        candidate.rotate();
        for (auto& [x, y]: gen(candidate))
        {
            currentCCOAs.emplace_back(candidateIdx, x, y, true);
        }
        // assert(ccoasAreOk(config, currentCCOAs));
    };
    rg::for_each(config.unpackedRectangles, processGen);
    // assert(ccoasAreOk(config, currentCCOAs));
    return currentCCOAs;
}

/**
 * Moves the rectangle rectIdx form unpacked to packed and returns the updated lists of CCOAs
 */
std::vector<CCOA> placeRectangle(Configuration& config,
                                 std::vector<CCOA>&& currentCCOAs,
                                 const size_t selectedCCOAIdx)
{
    CCOA ccoa = currentCCOAs[selectedCCOAIdx];
    currentCCOAs.erase(currentCCOAs.begin() + selectedCCOAIdx);
    const size_t selectedIdx = ccoa.rectangleIdx;
    Rectangle& r = config.allRectangles[selectedIdx];
    r.applyCCOA(ccoa);
    config.unpackedRectangles.erase(rg::find(config.unpackedRectangles, selectedIdx));
    config.packedRectangles.push_back(selectedIdx);

    // assert(placedRectangleIsOk(config, selectedIdx));

    auto [first, last] = rg::remove_if(currentCCOAs,
                                       [selectedIdx, &config = std::as_const(config)](const CCOA& c)
                                       {
                                           return config.isInfeasible(selectedIdx, c);
                                       });
    currentCCOAs.erase(first, last);

    return addNewPossibleCCOAs(config, selectedIdx, std::move(currentCCOAs));
}

Configuration A0(Configuration&& config, std::vector<CCOA>&& ccoas)
{
    auto distToOrgin = [x0 = config.container.x, y0 = config.container.y](int x, int y)
    {
        int dx = x - x0;
        int dy = y - y0;
        return dx * dx + dy * dy;
    };
    while (!ccoas.empty())
    {
        size_t maxCCOAIdx = 0;
        std::pair<int, int> maxCCOACoord{};
        float maxDegree = -1;
        for (size_t i = 0; const auto& ccoa: ccoas)
        {
            if (float deg = degree(ccoa, config);
                deg > maxDegree ||
                (deg == maxDegree &&
                 distToOrgin(ccoa.x, ccoa.y) < distToOrgin(maxCCOACoord.first, maxCCOACoord.second)))
            {
                maxDegree = deg;
                maxCCOAIdx = i;
                maxCCOACoord = {ccoa.x, ccoa.y};
            }
            ++i;
        }
        ccoas = placeRectangle(config, std::move(ccoas), maxCCOAIdx);
    }
    return config;
}

Configuration benefitA1(Configuration config, std::vector<CCOA> ccoas, const size_t selectedCCOAIdx)
{
    ccoas = placeRectangle(config, std::move(ccoas), selectedCCOAIdx);
    return A0(std::move(config), std::move(ccoas));
}

struct MaxBenefitInfo
{
    uint benefit = 0;
    size_t idx = 0;
    Coordinate c{};
    const Coordinate origin;

    MaxBenefitInfo(const Coordinate& orig) : origin(orig) {}

    uint distToOrgin(const int x, const int y)
    {
        int dx = x - origin.x;
        int dy = y - origin.y;
        return dx * dx + dy * dy;
    }

    void update(const uint& dens, const CCOA& ccoa, const size_t index)
    {
        if (dens > benefit ||
            (dens == benefit && distToOrgin(ccoa.x, ccoa.y) < distToOrgin(c.x, c.y)))
        {
            benefit = dens;
            idx = index;
            c = {ccoa.x, ccoa.y};
        }
    }
};

Configuration A1(const Rectangle& container, const std::vector<Rectangle>& toPack)
{
    Configuration initialConfig{container, toPack};
    auto ccoas = calculateInitialCCOAs(initialConfig);
    while (!ccoas.empty())
    {

        std::mutex mtx;
        Configuration solutionConfig = initialConfig;

        std::atomic_bool foundSol = false;
        std::vector<std::pair<size_t, CCOA>> idxsCcoa;
        idxsCcoa.reserve(ccoas.size());
        std::transform(ccoas.begin(),
                       ccoas.end(),
                   std::back_inserter(idxsCcoa),
                   [i = 0](const CCOA& ccoa) mutable -> std::pair<size_t, CCOA>
                       {
                           return {i++, ccoa};
                       });
        MaxBenefitInfo maxInfo({container.x, container.y});
        std::for_each(
            std::execution::par_unseq,
            idxsCcoa.begin(),
            idxsCcoa.end(),
            [&foundSol, &solutionConfig, &mtx, &initialConfig, &ccoas, &maxInfo](const auto& idxCcoa)
            {
                const auto& [index, ccoa] = idxCcoa;
                if (foundSol)
                    return;

                Configuration config = benefitA1(initialConfig, ccoas, index);
                if (config.isSuccessful())
                {
                    foundSol = true;
                    const std::lock_guard g(mtx);
                    solutionConfig = config;
                    return;
                }
                const uint dens = config.noNormDensity();
                maxInfo.update(dens, ccoa, index);
            });

        if (foundSol)
        {
            std::println("SUCCESS!!");
            return solutionConfig;
        }
        ccoas = placeRectangle(initialConfig, std::move(ccoas), maxInfo.idx);
    }
    std::println("Finished without packing all the rectangles");
    return initialConfig;
}

std::vector<std::vector<Rectangle>> placeAll(const Rectangle& container, std::vector<Rectangle> toPack)
{
    std::vector<std::vector<Rectangle>> result;
    while (!toPack.empty())
    {
        auto packed = A1(container, toPack);
        std::vector<Rectangle> partial;
        for (const size_t i : packed.packedRectangles)
        {
            partial.push_back(packed.allRectangles[i]);
        }
        result.push_back(std::move(partial));
        toPack.clear();
        for (const size_t i : packed.unpackedRectangles)
        {
            toPack.push_back(packed.allRectangles[i]);
        }
    }
    return result;
}

std::vector<std::vector<Rectangle>> placeAll(Rectangle container,
                                             std::vector<Rectangle> toPack,
                                             uint cutThick)
{
    cutThick += cutThick % 2;
    uint half = cutThick / 2;
    auto enlarge = [cutThick, half](Rectangle& r)
    {
        r.width += cutThick;
        r.height += cutThick;
        r.x -= half;
        r.y -= half;
    };
    auto reduce = [cutThick, half](Rectangle& r)
    {
        r.width -= cutThick;
        r.height -= cutThick;
        r.x += half;
        r.y += half;
    };
    enlarge(container);
    rg::for_each(toPack, enlarge);
    auto result = placeAll(container, std::move(toPack));
    rg::for_each(result,
                 [&reduce](auto& v)
                 {
                     rg::for_each(v, reduce);
                 });
    return result;
}
std::vector<std::vector<Rectangle>> placeAll(Rectangle container,
                                             std::vector<Rectangle> toPack,
                                             uint cutThick,
                                             uint containerBorder)
{
    // Reduce container
    uint half = containerBorder / 2;
    container.width -= containerBorder;
    container.height -= containerBorder;
    container.x += half;
    container.y += half;

    return placeAll(container, std::move(toPack), cutThick);
}

bool testCase1()
{
    Rectangle container{5, 5};
    std::vector<Rectangle> toPack{
        {2, 2},
        {2, 2},
        {2, 2},
        {2, 2},
        {2, 1},
        {2, 1},
        {2, 1},
        {2, 1},
        {2, 2},
    };
    auto result = A1(container, toPack);
    result.write("output.csv");
    return result.isSuccessful();
}

void writeBoards(const std::string& fnameTemplate, const Rectangle& container, const std::vector<std::vector<Rectangle>>& solutions)
{
    for (size_t i = 0; i < solutions.size(); ++i)
    {
        auto filename = std::format("{}_{}.csv", fnameTemplate, i);
        std::ofstream file(filename);

        if (!file.is_open())
        {
            std::cerr << "Failed to open file: " << filename << std::endl;
            return;
        }

        // Write the header
        file << "x;y;width;height\n";
        file << std::format("{};{};{};{}\n", container.x, container.y, container.width, container.height);

        // Write the data
        for (const auto& rect: solutions[i])
        {
            file << std::format("{};{};{};{}\n", rect.x, rect.y, rect.width, rect.height);
        }

        file.close();
        if (!file)
        {
            std::cerr << "Failed to close file properly: " << filename << std::endl;
        }
    }
}

void testPackCuts()
{
    const Rectangle board {2800, 2050};
    uint cutThick = 14;
    uint borderClear = 20;
    std::vector<Rectangle> allRects10;
    for (auto [w, h] : to10mmBoard)
    {
        allRects10.emplace_back(w, h);
    }
    auto boardsCut10 = placeAll(board, allRects10, cutThick, borderClear);
    std::println("We need {} 10mm thick boards", boardsCut10.size());
    writeBoards("outputs/board10mm", board, boardsCut10);

    std::vector<Rectangle> allRects16;
    for (auto [w, h] : to16mmBoard)
    {
        allRects16.emplace_back(w, h);
    }
    auto boardsCut16 = placeAll(board, allRects16, cutThick, borderClear);
    std::println("We need {} 16mm thick boards", boardsCut16.size());
    writeBoards("outputs/board16mm", board, boardsCut16);
}

void testLargeProfile()
{
    const Rectangle board {2800, 2050};
    uint cutThick = 14;
    uint borderClear = 20;

    std::vector<Rectangle> allRects16;
    for (int i = 0; i < 2; ++i)
    {
        for (auto [w, h]: to16mmBoard)
        {
            allRects16.emplace_back(w, h);
        }
    }
    auto boardsCut16 = placeAll(board, allRects16, cutThick, borderClear);
    std::println("We need {} 16mm thick boards", boardsCut16.size());
}

void testPackCatP1()
{
    const Rectangle board {20, 20};
    std::vector<Rectangle> allRectsCatP1;
    for (auto [w, h] : catP1)
    {
        allRectsCatP1.emplace_back(w, h);
    }
    auto boardsCut = placeAll(board, allRectsCatP1);
    std::println("We need {} thick boards", boardsCut.size());
    writeBoards("outputs/boardCatP1", board, boardsCut);
}

void testProblematicCase()
{
    const Rectangle board {2800, 2100};
    std::vector<Rectangle> allRects10;
    for (auto [w, h] : to10mmBoard)
    {
        allRects10.emplace_back(w, h);
        if (allRects10.size() == 10)
        {
            break;
        }
    }
    auto sol = A1(board, allRects10);
    sol.write("problematic.csv");
}

#include <iostream>

#define TEST_EQUAL(arg1, arg2) \
    do { \
        if ((arg1) != (arg2)) { \
            std::cerr << "Error: \n" \
                      << "    Expression 1: " << #arg1 << " = " << (arg1) << "\n" \
                      << "    Expression 2: " << #arg2 << " = " << (arg2) << "\n"; \
        } \
    } while (0)

#define TEST_TRUE(arg1) \
    do { \
        if (!(arg1)) { \
            std::cerr << "Error: \n" \
                      << "    Expression: " << #arg1 << " = false\n"; \
        } \
    } while (0)

#define TEST_FALSE(arg1) \
    do { \
        if (arg1) { \
            std::cerr << "Error: \n" \
                      << "    Expression: " << #arg1 << " = true\n"; \
        } \
    } while (0)

namespace test
{

    void testRectangle()
    {
        Rectangle r1 {10, 20, 0, 0};
        Rectangle r2 {10, 20, 5, 2};
        Rectangle r3 {10, 20, 5, 20};
        Rectangle cont {100, 200, -5, -5};

        TEST_TRUE(r1.overlaps(r2));
        TEST_FALSE(r1.overlaps(r3));
        TEST_TRUE(cont.contains(r1));
        auto d = r1.distanceSq(r2);
        TEST_FALSE(d);
        d = r1.distanceSq(r3);
        TEST_TRUE(d);
        TEST_EQUAL(*d, 0);
    }

    void testConfiguration()
    {
        Rectangle r1 {20, 10, 0, 0};
        Rectangle r2 {10, 20, 0, 0};
        Rectangle r3 {10, 20, 0, 0};
        Rectangle cont {100, 200, -10, -20};

        Configuration c {cont, {r1, r2, r3}};
        auto ccoas = calculateInitialCCOAs(c);
        TEST_EQUAL(ccoas.size(), 2*4*3);
        CCOA ccoa {0, -10, -20, false};
        TEST_TRUE(ccoas[0] == ccoa);

        ccoas = placeRectangle(c, std::move(ccoas), 0);
        TEST_EQUAL(ccoas.size(), (6 + 4) * 2);

        ccoa = {1, -10, -10, false};
        auto it = std::find(ccoas.begin(), ccoas.end(), ccoa);
        TEST_TRUE(it != ccoas.end());
        auto index = it - ccoas.begin();

        ccoas = placeRectangle(c, std::move(ccoas), index);
        TEST_EQUAL(ccoas.size(), 6 + 2 + 2 + 2);
    }
}

int main()
{
    // testPackCuts();
    testLargeProfile();
    // testPackCatP1();
    // test::testRectangle();
    // test::testConfiguration();
    return 0;
}
