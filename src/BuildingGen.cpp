#include "BuildingGen.h"

static vector<array<GLubyte,3>> m_industrialColors = {
    {50,50,50},   // Dark Gray
    {20,20,20},   // Black
    {100,100,100}, // Medium Gray
    {200,200,200},  // Light Gray
    {130,130,150}, // Darker steel blue
    {210,210,230}, // Very pale blue-gray
};

static vector<array<GLubyte,3>> m_commercialColors = {
    {240,240,240}, // White
    {180,180,180}, // Light Gray
    {100,100,100}, // Medium Gray
    {20,20,20},    // Black
    {100,100,100}, // Dark Gray
    {180,70,0}     // Dark Reddish Orange
};


static vector<array<GLubyte,3>> m_residentialColors = {
    {230,210,180}, // Beige
    {140,100,60},  // Dark Brown
    {210,180,140}, // Tan
    {165,42,42},   // Reddish Brown (Brownish Red)
    {139,69,19},   // Brown Brown (Saddle Brown)
    {85,107,47}    // Brownish Green (Olive Drab)
};

BuildingGen::BuildingGen(float districtRadius, float minPolygonArea,
                                     int IndustrialMinStories, int IndustrialMaxStories,
                                     int CommercialMinStories, int CommercialMaxStories,
                                     int ResidentialMinStories, int ResidentialMaxStories,
                                     float baseStoryHeight, float baseAddition0, float baseAddition2,
                                     float specialBuildingChance)
: m_districtRadius(districtRadius),
  m_minPolygonArea(minPolygonArea),
  m_IndustrialMinStories(IndustrialMinStories),
  m_IndustrialMaxStories(IndustrialMaxStories),
  m_CommercialMinStories(CommercialMinStories),
  m_CommercialMaxStories(CommercialMaxStories),
  m_ResidentialMinStories(ResidentialMinStories),
  m_ResidentialMaxStories(ResidentialMaxStories),
  m_baseStoryHeight(baseStoryHeight),
  m_baseAddition0(baseAddition0),
  m_baseAddition2(baseAddition2),
  m_specialBuildingChance(specialBuildingChance)
{
}

Building BuildingGen::generateBuildingDetails(const std::vector<Point>& polygon, const Point& districtCenter) {
    Building b;
    b.polygons = polygon;
    b.isSpecial = false;

    float area = PolygonUtils::calculatePolygonArea(polygon);
    if (area < m_minPolygonArea) {
        b.height = 0.0f;
        return b;
    }

    Point centroid = PolygonUtils::getCentroid(polygon);
    double dist = PolygonUtils::distanceBetweenPoints(centroid, districtCenter);
    float normalizedDist = std::clamp((float)(dist / m_districtRadius), 0.0f, 1.0f);

    float industrialThreshold = m_districtRadius / 3.0f;
    float commercialThreshold = (2.0f * m_districtRadius) / 3.0f;

    std::random_device rd;
    std::mt19937 gen(rd());

    if (dist < industrialThreshold) {
        return createIndustrialBuilding(polygon, dist, normalizedDist, gen);
    } else if (dist < commercialThreshold) {
        return createCommercialBuilding(polygon, dist, normalizedDist, gen);
    } else {
        return createResidentialBuilding(polygon, dist, normalizedDist, gen);
    }
}

Building BuildingGen::createIndustrialBuilding(const std::vector<Point>& polygon, double dist, float normalizedDist, std::mt19937& gen) {
    Building b;
    b.polygons = polygon;
    b.isSpecial = false;
    b.textureLayer = 2.0f; 

    std::uniform_int_distribution<int> storyDist(m_IndustrialMinStories, m_IndustrialMaxStories);
    int baseStories = storyDist(gen);

    // heightFactor example
    float industrialThreshold = m_districtRadius / 3.0f;
    float heightFactor = 2.0f - (float)(dist / industrialThreshold);

    float adjustedStories = baseStories * heightFactor;
    int finalStories = std::max((int)std::round(adjustedStories), 1);
    b.height = finalStories * m_baseStoryHeight + m_baseAddition2;

    // Color
    std::uniform_int_distribution<size_t> colorDist(0, m_industrialColors.size()-1);
    auto c = m_industrialColors[colorDist(gen)];
    b.r = c[0]; b.g = c[1]; b.b = c[2];

    // Special
    std::uniform_real_distribution<float> chanceDist(0.0f, 1.0f);
    if (chanceDist(gen) < m_specialBuildingChance) {
        b.isSpecial = true;
    }

    b.district = District::Industrial;
    return b;
}

Building BuildingGen::createCommercialBuilding(const std::vector<Point>& polygon, double dist, float normalizedDist, std::mt19937& gen) {
    Building b;
    b.polygons = polygon;
    b.isSpecial = false;

    std::uniform_int_distribution<int> storyDist(m_CommercialMinStories, m_CommercialMaxStories);
    int baseStories = storyDist(gen);

    float industrialThreshold = m_districtRadius / 3.0f;
    float commercialThreshold = (2.0f * m_districtRadius) / 3.0f;
    float range = commercialThreshold - industrialThreshold;
    float distFromInd = (float)(dist - industrialThreshold);
    float heightFactor = 1.8f - (0.8f * (distFromInd / range));

    float adjustedStories = baseStories * heightFactor;
    int finalStories = std::max((int)std::round(adjustedStories), 1);
    b.height = finalStories * m_baseStoryHeight;
    
    // Color
    std::uniform_int_distribution<size_t> colorDist(0, m_commercialColors.size()-1);
    auto c = m_commercialColors[colorDist(gen)];
    b.r = c[0]; b.g = c[1]; b.b = c[2];

    // Texture layer random
    std::uniform_real_distribution<float> textureChance(0.0f, 1.0f);
    if (textureChance(gen) < 0.5f) {
        b.textureLayer = 1.0f;
    } else {
        b.textureLayer = 0.0f;
        b.height += m_baseAddition0;
    }

    b.district = District::Commercial;
    return b;
}

Building BuildingGen::createResidentialBuilding(const std::vector<Point>& polygon, double dist, float normalizedDist, std::mt19937& gen) {
    Building b;
    b.polygons = polygon;
    b.isSpecial = false;

    std::uniform_int_distribution<int> storyDist(m_ResidentialMinStories, m_ResidentialMaxStories);
    int baseStories = storyDist(gen);

    float commercialThreshold = (2.0f * m_districtRadius) / 3.0f;
    float range = m_districtRadius - commercialThreshold;
    float distFromComm = (float)(dist - commercialThreshold);
    float heightFactor = 1.5f - (0.5f * (distFromComm / range));

    float adjustedStories = baseStories * heightFactor;
    int finalStories = std::max((int)std::round(adjustedStories), 1);
    b.height = finalStories * m_baseStoryHeight;

    // Color
    std::uniform_int_distribution<size_t> colorDist(0, m_residentialColors.size()-1);
    auto c = m_residentialColors[colorDist(gen)];
    b.r = c[0]; b.g = c[1]; b.b = c[2];

    // Texture layer random
    std::uniform_real_distribution<float> textureChance(0.0f, 1.0f);
    
    if (textureChance(gen) < 0.5f) {
        b.textureLayer = 1.0f;
    } else {
        b.textureLayer = 0.0f;
        b.height += m_baseAddition0;
    }

    b.district = District::Residential;
    return b;
}
