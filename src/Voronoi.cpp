#include "Voronoi.h"
#include "PolygonUtils.h"
#include <random>
#include <iostream>

std::vector<Point> Voronoi::generateSites(int numSites, unsigned int seed, double minX, double maxX, double minY, double maxY) {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> distX(minX + 2, maxX - 2);
    std::uniform_real_distribution<double> distY(minY + 2, maxY - 2);

    std::vector<Point> sites;
    for (int i = 0; i < numSites; ++i) {
        float x = distX(gen);
        float y = distY(gen);
        sites.push_back({x, y});
    }
    return sites;
}

void Voronoi::computeVoronoiDiagram(std::vector<Point>& m_sites, std::vector<std::vector<Point>>& m_voronoiCells, 
                                         double minX, double maxX, double minY, double maxY) {
    m_voronoiCells.clear();
    m_voronoiCells.resize(m_sites.size());

    std::vector<Point> boundingPolygon = {
        { minX, minY },
        { maxX, minY },
        { maxX, maxY },
        { minX, maxY }
    };

    for (size_t i = 0; i < m_sites.size(); i++) {
        const Point& site = m_sites[i];
        std::vector<Point> cell = boundingPolygon;
        for (size_t j = 0; j < m_sites.size(); j++) {
            if (i == j) continue;
            const Point& otherSite = m_sites[j];
            Line bisector = GeometryUtils::perpendicularBisector(site, otherSite);
            double eval = GeometryUtils::evaluate(bisector, site);
            bool keepPositiveSide = eval > 0;
            auto clipped = PolygonUtils::clipPolygon(cell, bisector);
            cell = keepPositiveSide ? clipped.first : clipped.second;
            if (cell.empty()) {
                std::cout<<"\nCell is empty somefuckinghow";
                break;
            }
        }
        m_voronoiCells[i] = cell;
    }
}
