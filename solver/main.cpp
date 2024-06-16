#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <fstream>
#include <functional>
#include <concepts>
#include <format>
#include <iostream>
#include <numeric>
#include <optional>
#include <print>
#include <utility>
#include <vector>
#include <ranges>

namespace rg = std::ranges;
namespace rgv = std::ranges::views;
using uint = unsigned;

template<std::integral T>
bool isBetween(T number, T m, T M)
{
    return m <= number && number <= M;
}

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
        assert(container.contains(*this));
        return {static_cast<uint>(x - container.x),
                static_cast<uint>(container.xmax() - xmax()),
                static_cast<uint>(y - container.y),
                static_cast<uint>(container.ymax() - ymax())};
    }

    bool contains(const Rectangle& other) const
    {
        return (isBetween(other.x, x, x + xmax()) && isBetween(other.xmax(), x, x + xmax()) &&
                isBetween(other.y, y, y + ymax()) && isBetween(other.ymax(), y, y + ymax()));
    }
};

struct CCOA
{
    int rectangleIdx;
    int x;
    int y;
    bool rotated;
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
                assert(false);
                std::println("Trying to evaluate distance with overlapping rectangles");
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
            return false;
        }
        const Rectangle& placedRect = allRectangles[placedRectIdx];
        Rectangle candidate = allRectangles[ccoa.rectangleIdx];
        candidate.x = ccoa.x;
        candidate.y = ccoa.y;
        if (ccoa.rotated)
            candidate.rotate();
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
            [&checkAndPush, &container = config.container, rectIdx](Rectangle& r, bool rot)
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
    if (ccoa.rotated)
        r.rotate();
    return 1. - config.minDistance(r) * 2. / (r.width + r.height);
}

void addCandidateIfPossible(const Configuration& config,
                            const Rectangle& candidate,
                            std::vector<std::pair<int, int>>& candidateCoords)
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

class PossibleNewPositionGeneratorAbove
{
public:
    const Configuration& config;

    PossibleNewPositionGeneratorAbove(const Configuration& configuration,
                                      const size_t lastPlacedRectIdx):
        config(configuration)
    {
        const Rectangle& lastPlaced = config.allRectangles[lastPlacedRectIdx];
        const Rectangle& container = config.container;

        yCoord = lastPlaced.ymax();
        if (yCoord == container.ymax())
        {
            return;
        }
        xToCheck.push_back(container.x);
        for (const size_t packedRectIdx: config.packedRectangles)
        {
            const auto& packedRect = config.allRectangles[packedRectIdx];
            int possibleX = packedRect.xmax();
            if (packedRect.ymax() > yCoord && possibleX < container.xmax())
            {
                xToCheck.push_back(possibleX);
            }
        }
    }

    std::vector<std::pair<int, int>> operator()(const Rectangle& candidate) const
    {
        if (yCoord + candidate.height > config.container.ymax())
            return {};
        std::vector<std::pair<int, int>> result;

        int lastPossibleX = config.container.xmax() - candidate.width;
        Rectangle candCopy = candidate;
        candCopy.y = yCoord;
        for (const int x: xToCheck)
        {
            if (x > lastPossibleX)
            {
                break;
            }
            candCopy.x = x;
            addCandidateIfPossible(config, candCopy, result);
        }
        // Need to check last one
        candCopy.x = lastPossibleX;
        addCandidateIfPossible(config, candCopy, result);
        return result;
    }

private:
    std::vector<int> xToCheck;
    int yCoord;
};

class PossibleNewPositionGeneratorBellow
{
public:
    const Configuration& config;

    PossibleNewPositionGeneratorBellow(const Configuration& configuration,
                                       const size_t lastPlacedRectIdx):
        config(configuration)
    {
        const Rectangle& lastPlaced = config.allRectangles[lastPlacedRectIdx];
        const Rectangle& container = config.container;

        yMaxCoord = lastPlaced.y;
        if (yMaxCoord == container.y)
        {
            return;
        }
        xToCheck.push_back(container.x);
        for (const size_t packedRectIdx: config.packedRectangles)
        {
            const auto& packedRect = config.allRectangles[packedRectIdx];
            int possibleX = packedRect.xmax();
            if (packedRect.y < yMaxCoord && possibleX < container.xmax())
            {
                xToCheck.push_back(possibleX);
            }
        }
    }

    std::vector<std::pair<int, int>> operator()(const Rectangle& candidate) const
    {
        int yMin = yMaxCoord - candidate.height;
        if (yMin < config.container.y)
            return {};
        std::vector<std::pair<int, int>> result;

        int lastPossibleX = config.container.xmax() - candidate.width;
        Rectangle candCopy = candidate;
        candCopy.y = yMin;
        for (const int x: xToCheck)
        {
            if (x > lastPossibleX)
            {
                break;
            }
            candCopy.x = x;
            addCandidateIfPossible(config, candCopy, result);
        }
        // Need to check last one
        candCopy.x = lastPossibleX;
        addCandidateIfPossible(config, candCopy, result);
        return result;
    }

private:
    std::vector<int> xToCheck;
    int yMaxCoord;
};

class PossibleNewPositionGeneratorRight
{
public:
    const Configuration& config;

    PossibleNewPositionGeneratorRight(const Configuration& configuration,
                                      const size_t lastPlacedRectIdx):
        config(configuration)
    {
        const Rectangle& lastPlaced = config.allRectangles[lastPlacedRectIdx];
        const Rectangle& container = config.container;

        xCoord = lastPlaced.xmax();
        if (xCoord == container.xmax())
        {
            return;
        }
        yToCheck.push_back(container.y);
        for (const size_t packedRectIdx: config.packedRectangles)
        {
            const auto& packedRect = config.allRectangles[packedRectIdx];
            int possibleY = packedRect.ymax();
            if (packedRect.xmax() > xCoord && possibleY < container.ymax())
            {
                yToCheck.push_back(possibleY);
            }
        }
    }

    std::vector<std::pair<int, int>> operator()(const Rectangle& candidate) const
    {
        if (xCoord + candidate.width > config.container.xmax())
            return {};
        std::vector<std::pair<int, int>> result;

        int lastPossibleY = config.container.ymax() - candidate.height;
        Rectangle candCopy = candidate;
        candCopy.x = xCoord;
        for (const int y: yToCheck)
        {
            if (y > lastPossibleY)
            {
                break;
            }
            candCopy.y = y;
            addCandidateIfPossible(config, candCopy, result);
        }
        // Need to check last one
        candCopy.y = lastPossibleY;
        addCandidateIfPossible(config, candCopy, result);
        return result;
    }

private:
    std::vector<int> yToCheck;
    int xCoord;
};

class PossibleNewPositionGeneratorLeft
{
public:
    const Configuration& config;

    PossibleNewPositionGeneratorLeft(const Configuration& configuration,
                                     const size_t lastPlacedRectIdx):
        config(configuration)
    {
        const Rectangle& lastPlaced = config.allRectangles[lastPlacedRectIdx];
        const Rectangle& container = config.container;

        xMaxCoord = lastPlaced.x;
        if (xMaxCoord == container.x)
        {
            return;
        }
        yToCheck.push_back(container.y);
        for (const size_t packedRectIdx: config.packedRectangles)
        {
            const auto& packedRect = config.allRectangles[packedRectIdx];
            int possibleY = packedRect.ymax();
            if (packedRect.x < xMaxCoord && possibleY < container.ymax())
            {
                yToCheck.push_back(possibleY);
            }
        }
    }

    std::vector<std::pair<int, int>> operator()(const Rectangle& candidate) const
    {
        int xMin = xMaxCoord - candidate.width;
        if (xMin < config.container.y)
            return {};
        std::vector<std::pair<int, int>> result;

        int lastPossibleY = config.container.ymax() - candidate.height;
        Rectangle candCopy = candidate;
        candCopy.x = xMin;
        for (const int y: yToCheck)
        {
            if (y > lastPossibleY)
            {
                break;
            }
            candCopy.y = y;
            addCandidateIfPossible(config, candCopy, result);
        }
        // Need to check last one
        candCopy.y = lastPossibleY;
        addCandidateIfPossible(config, candCopy, result);
        return result;
    }

private:
    std::vector<int> yToCheck;
    int xMaxCoord;
};

std::vector<CCOA> addNewPossibleCCOAs(Configuration& config,
                                      const size_t lastPlacedRectIdx,
                                      std::vector<CCOA>&& currentCCOAs)
{
    const PossibleNewPositionGeneratorAbove genAb(config, lastPlacedRectIdx);
    const PossibleNewPositionGeneratorRight genRi(config, lastPlacedRectIdx);
    const PossibleNewPositionGeneratorBellow genBe(config, lastPlacedRectIdx);
    const PossibleNewPositionGeneratorLeft genLe(config, lastPlacedRectIdx);
    auto processSide = [&currentCCOAs, &config = std::as_const(config)](const auto& sideGen,
                                                                        const size_t candidateIdx)
    {
        Rectangle candidate = config.allRectangles[candidateIdx];
        for (auto& [x, y]: sideGen(candidate))
        {
            currentCCOAs.emplace_back(candidateIdx, x, y, false);
        }
        candidate.rotate();
        for (auto& [x, y]: sideGen(candidate))
        {
            currentCCOAs.emplace_back(candidateIdx, x, y, true);
        }
    };
    auto processCandidate =
        [&genAb, &genRi, &genBe, &genLe, &processSide](const size_t candidateIdx)
    {
        processSide(genAb, candidateIdx);
        processSide(genRi, candidateIdx);
        processSide(genBe, candidateIdx);
        processSide(genLe, candidateIdx);
    };
    rg::for_each(config.unpackedRectangles, processCandidate);
    return currentCCOAs;
}

/**
 * Moves the rectangle rectIdx form unpacked to packed and returns the updated lists of CCOAs
 */
std::vector<CCOA> placeRectangle(Configuration& config,
                                 std::vector<CCOA>&& currentCCOAs,
                                 const std::vector<CCOA>::iterator selectedCCOA)
{
    CCOA ccoa = *selectedCCOA;
    currentCCOAs.erase(selectedCCOA);
    const size_t selectedIdx = ccoa.rectangleIdx;
    Rectangle& r = config.allRectangles[selectedIdx];
    r.x = ccoa.x;
    r.y = ccoa.y;
    if (ccoa.rotated)
        r.rotate();
    config.unpackedRectangles.erase(rg::find(config.unpackedRectangles, selectedIdx));
    config.packedRectangles.push_back(selectedIdx);

    auto [first, last] = rg::remove_if(currentCCOAs,
                                       [selectedIdx, &config = std::as_const(config)](const CCOA& c)
                                       {
                                           return c.rectangleIdx == selectedIdx ||
                                                  config.isInfeasible(selectedIdx, c);
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
        auto maxCCOAIter = ccoas.begin();
        float maxDegree = degree(*maxCCOAIter, config);
        for (auto it = ccoas.begin() + 1; it != ccoas.end(); ++it)
        {
            if (float deg = degree(*it, config);
                deg > maxDegree ||
                (deg == maxDegree &&
                 distToOrgin(it->x, it->y) < distToOrgin(maxCCOAIter->x, maxCCOAIter->y)))
            {
                maxDegree = deg;
                maxCCOAIter = it;
            }
        }
        ccoas = placeRectangle(config, std::move(ccoas), maxCCOAIter);
    }
    return config;
}

Configuration benefitA1(Configuration config, std::vector<CCOA> ccoas, const size_t selectedCCOAIdx)
{
    placeRectangle(config, std::move(ccoas), ccoas.begin() + selectedCCOAIdx);
    return A0(std::move(config), std::move(ccoas));
}

Configuration A1(const Rectangle& container, const std::vector<Rectangle>& toPack)
{
    Configuration initialConfig{container, toPack};
    auto ccoas = calculateInitialCCOAs(initialConfig);
    while (!ccoas.empty())
    {
        float maxBenefit = 0;
        auto maxBenefitIter = ccoas.begin();
        size_t index = 0;
        for (auto it = ccoas.begin(); it != ccoas.end(); ++it, ++index)
        {
            Configuration config = benefitA1(initialConfig, ccoas, index);
            if (config.isSuccessful())
            {
                std::println("SUCCESS!!");
                return config;
            }
            if (float dens = config.density(); dens > maxBenefit)
            {
                maxBenefit = dens;
                maxBenefitIter = it;
            }
        }
        ccoas = placeRectangle(initialConfig, std::move(ccoas), maxBenefitIter);
    }
    std::println("Finished without packing all the rectangles");
    return initialConfig;
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
    std::filesystem::path fname {"output.csv"};
    result.write(std::filesystem::absolute(fname));
    return result.isSuccessful();
}

int main()
{
    if (!testCase1())
    {
        std::println("Not found solution for testCase1");
    }
    return 0;
}
