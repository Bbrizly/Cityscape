#include "CityGen.h"
#include "Building.h"
#include "Road.h"
#include "Cube.h"
#include <vector>
#include <random>
#include <cmath>
#include <iostream>
using namespace std;
using namespace glm;

vector<Building> buildings;
vector<Road> roads;

const int gridSize = 50;
const float spacing = 30.0f;
const float citySpacing = 1.25f;
const int sectorSize = 25;
const int skyScraperSector = 5;

const int numSites = 15;

struct Site {
    float x, z; // Using float for more precise positioning
};

vector<Site> generateVoronoiSites(int numSites, int gridSize, mt19937& gen)
{
    uniform_real_distribution<float> siteDist(-(gridSize * citySpacing), gridSize * citySpacing);
    vector<Site> sites;
    for (int i = 0; i < numSites; ++i) {
        sites.push_back({ siteDist(gen), siteDist(gen) });
    }
    return sites;
}

struct Cell {
    int siteIndex;
};

// Function to assign each grid point to the nearest Voronoi site
vector<vector<int>> assignGridToSites(int gridSize, const vector<Site>& sites, float spacing) {
    vector<vector<int>> grid(gridSize, vector<int>(gridSize, -1));
    float originX = - (gridSize / 2) * spacing;
    float originZ = - (gridSize / 2) * spacing;
    for (int x = 0; x < gridSize; ++x) {
        for (int z = 0; z < gridSize; ++z) {
            float worldX = originX + x * spacing;
            float worldZ = originZ + z * spacing;
            float minDist = numeric_limits<float>::max();
            int closestSite = -1;
            for (int i = 0; i < sites.size(); ++i) {
                float dx = worldX - sites[i].x;
                float dz = worldZ - sites[i].z;
                float dist = dx * dx + dz * dz;
                if (dist < minDist) {
                    minDist = dist;
                    closestSite = i;
                }
            }
            grid[x][z] = closestSite;
        }
    }
    return grid;
}

// Function to extract Voronoi edges based on grid assignments
vector<pair<vec3, vec3>> extractVoronoiEdges(const vector<vector<int>>& grid, int gridSize, float spacing) {
    vector<pair<vec3, vec3>> edges;
    float originX = - (gridSize / 2) * spacing;
    float originZ = - (gridSize / 2) * spacing;
    for (int x = 0; x < gridSize - 1; ++x) {
        for (int z = 0; z < gridSize - 1; ++z) {
            int current = grid[x][z];
            int right = grid[x + 1][z];
            int up = grid[x][z + 1];
            if (current != right) {
                // Vertical edge between (x, z) and (x+1, z)
                vec3 p1(originX + x * spacing, 0.0f, originZ + z * spacing);
                vec3 p2(originX + (x + 1) * spacing, 0.0f, originZ + z * spacing);
                edges.emplace_back(p1, p2);
            }
            if (current != up) {
                // Horizontal edge between (x, z) and (x, z+1)
                vec3 p1(originX + x * spacing, 0.0f, originZ + z * spacing);
                vec3 p2(originX + x * spacing, 0.0f, originZ + (z + 1) * spacing);
                edges.emplace_back(p1, p2);
            }
        }
    }
    return edges;
}

// Function to create roads from Voronoi edges
void createRoadsFromVoronoiEdges(const vector<pair<vec3, vec3>>& edges, vector<Road>& roads, GLuint program) {
    for (const auto& edge : edges) {
        // Depending on your Road implementation, you might need to create multiple Road segments
        // along the edge or represent the entire edge as a single Road.
        // Here, we'll create a Road at the midpoint for simplicity.

        vec3 midpoint = (edge.first + edge.second) / 2.0f;
        Road road(midpoint, 0);
        road.init(program);
        roads.push_back(road);
    }
}

// Function to place buildings around Voronoi sites within their respective cells
void placeBuildingsAroundSites(const vector<Site>& sites, float spacing, mt19937& gen, GLuint program, float scaleHeight) {
    uniform_real_distribution<float> offsetDist(-spacing / 2, spacing / 2);
    uniform_real_distribution<float> heightDist(1.0f, scaleHeight);
    
    for (const auto& site : sites) {
        // Define number of buildings per site
        int numBuildings = 10; // Adjust as needed based on site density
        
        for (int i = 0; i < numBuildings; ++i) {
            float offsetX = offsetDist(gen);
            float offsetZ = offsetDist(gen);
            vec3 buildingPosition = vec3(site.x + offsetX, 0.0f, site.z + offsetZ);
            
            // Determine building height
            float height = heightDist(gen);
            
            Building b(buildingPosition, height);
            b.init(program);
            buildings.push_back(b);
        }
    }
}

//SPLINE GRAPH TEST

// Function to compute Catmull-Rom spline points
std::vector<glm::vec3> computeCatmullRomSpline(const std::vector<Site>& sites, int pointsPerSegment)
{
    std::vector<glm::vec3> splinePoints;

    // Ensure there are enough points to form segments
    if (sites.size() < 4) {
        // Not enough points to form a spline
        return splinePoints;
    }

    for (size_t i = 0; i < sites.size() - 3; ++i) {

        cout<<"site "<<i<<" "<< sites[i].x<<", "<<sites[i].z<<endl;
        // Define the four control points for the current segment
        glm::vec3 p0(sites[i].x, 0.0f, sites[i].z);
        glm::vec3 p1(sites[i + 1].x, 0.0f, sites[i + 1].z);
        glm::vec3 p2(sites[i + 2].x, 0.0f, sites[i + 2].z);
        glm::vec3 p3(sites[i + 3].x, 0.0f, sites[i + 3].z);

        for (int j = 0; j < pointsPerSegment; ++j) {
            float t = j / static_cast<float>(pointsPerSegment);
            float t2 = t * t;
            float t3 = t2 * t;

            // Catmull-Rom spline formula
            glm::vec3 point = 0.5f * ((2.0f * p1) +
                                      (-p0 + p2) * t +
                                      (2.0f * p0 - 5.0f * p1 + 4.0f * p2 - p3) * t2 +
                                      (-p0 + 3.0f * p1 - 3.0f * p2 + p3) * t3);

            splinePoints.push_back(point);
        }
    }

    return splinePoints;
}

std::vector<glm::vec3> generateRoadMesh(const std::vector<glm::vec3>& splinePoints, float roadWidth)
{
    std::vector<glm::vec3> roadVertices;

    for (size_t i = 0; i < splinePoints.size() - 1; ++i) {
        glm::vec3 p0 = splinePoints[i];
        glm::vec3 p1 = splinePoints[i + 1];

        // Calculate the direction of the road segment
        glm::vec3 direction = glm::normalize(p1 - p0);

        // Calculate the perpendicular vector (assuming roads lie on the XZ plane)
        glm::vec3 normal = glm::vec3(-direction.z, 0.0f, direction.x);

        // Offset points to create the road width
        glm::vec3 left = p0 + normal * (roadWidth / 2.0f);
        glm::vec3 right = p0 - normal * (roadWidth / 2.0f);
        glm::vec3 leftNext = p1 + normal * (roadWidth / 2.0f);
        glm::vec3 rightNext = p1 - normal * (roadWidth / 2.0f);

        // Create two triangles for each road segment (forming a quad)
        // Triangle 1
        roadVertices.push_back(left);
        roadVertices.push_back(right);
        roadVertices.push_back(leftNext);

        // Triangle 2
        roadVertices.push_back(right);
        roadVertices.push_back(rightNext);
        roadVertices.push_back(leftNext);
    }

    return roadVertices;
}

// AAAAAAAAAAA SPLINE GRAPH TEST END

void CityGen::generate(GLuint program)
{
    m_program = program;

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> scaleDist(0.4f, 2.0f);
    uniform_real_distribution<> scaleHeight(1.0f, 2.0f);
    uniform_real_distribution<> offsetDist(-1.0f, 1.0f);

    vector<Site> sites = generateVoronoiSites(numSites, gridSize, gen);
    
    if (sites.size() >= 4) {
        sites.push_back(sites[0]);
        sites.push_back(sites[1]);
        sites.push_back(sites[2]);
    }

    /* AA */
    vector<vector<int>> grid = assignGridToSites(gridSize, sites, spacing);

    vector<pair<vec3, vec3>> voronoiEdges = extractVoronoiEdges(grid, gridSize, spacing);

    createRoadsFromVoronoiEdges(voronoiEdges, roads, m_program);
    
    /* AA */
    placeBuildingsAroundSites(sites, spacing, gen, m_program, 10.0f);
    
}

void CityGen::reGenerate()
{
    roads.clear();
    buildings.clear();
    generate(m_program);
}

void CityGen::render(const glm::mat4& proj, const glm::mat4& view, float m_timer)
{
    // glUseProgram(m_program);

    // glUniformMatrix4fv(glGetUniformLocation(m_program, "projection"), 1, GL_FALSE, glm::value_ptr(proj));
    // glUniformMatrix4fv(glGetUniformLocation(m_program, "view"), 1, GL_FALSE, glm::value_ptr(view));
    // glUniform1f(glGetUniformLocation(m_program, "u_time"), m_timer);

    // glBindVertexArray(buildings[0].getVAO());

    // glDrawArraysInstanced(GL_TRIANGLES, 0, buildings[0].getVertexCount(), instanceTransforms.size());

///
    // glBindVertexArray(m_vao);

    // glm::mat4 model = getTransformationMatrix();
    // GLint program;
    // glGetIntegerv(GL_CURRENT_PROGRAM,&program);

    // glUniformMatrix4fv(glGetUniformLocation(program, "projection"), 1, GL_FALSE, glm::value_ptr(proj));
    // glUniformMatrix4fv(glGetUniformLocation(program, "view"), 1, GL_FALSE, glm::value_ptr(view));
    // glUniformMatrix4fv(glGetUniformLocation(program, "world"), 1, GL_FALSE, glm::value_ptr(model));
    // glUniform1f(glGetUniformLocation(program, "u_time"), m_timer);
    
    // glDrawArrays(GL_TRIANGLES, 0, buildings.size() * 36);
///
    for (auto& road : roads)
    {
        road.draw(proj, view, m_timer);
    }
    for (auto& building : buildings)
    {
        building.draw(proj, view, m_timer);
    }
    // for (auto& cube : buildings)
    // {
    //     cube.draw(proj, view, m_timer);
    // }
}
