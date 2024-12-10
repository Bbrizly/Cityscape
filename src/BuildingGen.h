// BuildingGenerator.h
#pragma once
#include "PolygonUtils.h"
#include <algorithm>
#include <cmath>
#include <vector>
#include <array>
#include <random>
#include "Types.h" // For Building, District, Point, etc.
using namespace std;

class BuildingGen {
public:
    BuildingGen(float districtRadius, float minPolygonArea,
                      int IndustrialMinStories, int IndustrialMaxStories,
                      int CommercialMinStories, int CommercialMaxStories,
                      int ResidentialMinStories, int ResidentialMaxStories,
                      float baseStoryHeight, float baseAddition0, float baseAddition2,
                      float specialBuildingChance);

    Building generateBuildingDetails(const std::vector<Point>& polygon, const Point& districtCenter);

private:
    // Store parameters
    float m_districtRadius;
    float m_minPolygonArea;
    int m_IndustrialMinStories, m_IndustrialMaxStories;
    int m_CommercialMinStories, m_CommercialMaxStories;
    int m_ResidentialMinStories, m_ResidentialMaxStories;
    float m_baseStoryHeight;
    float m_baseAddition0;
    float m_baseAddition2;
    float m_specialBuildingChance;

    // Helper methods
    Building createIndustrialBuilding(const std::vector<Point>& polygon, double dist, float normalizedDist, std::mt19937& gen);
    Building createCommercialBuilding(const std::vector<Point>& polygon, double dist, float normalizedDist, std::mt19937& gen);
    Building createResidentialBuilding(const std::vector<Point>& polygon, double dist, float normalizedDist, std::mt19937& gen);
};
