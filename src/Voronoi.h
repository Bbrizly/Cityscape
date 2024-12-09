#pragma once
#include <vector>
#include "GeometryUtils.h"

class Voronoi {
public:
    static std::vector<Point> generateSites(int numSites, unsigned int seed, double minX, double maxX, double minY, double maxY);
    static void computeVoronoiDiagram(std::vector<Point>& m_sites, std::vector<std::vector<Point>>& m_voronoiCells, 
                                      double minX, double maxX, double minY, double maxY);
};
