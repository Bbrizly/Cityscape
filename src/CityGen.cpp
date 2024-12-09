#include "CityGen.h"
#include "GeometryUtils.h"
#include "PolygonUtils.h"
#include "Debug.h"
#include "Voronoi.h"
#include "DrawRoad.h"
#include <functional>
#include <random>
#include <iostream>
#include <glm/gtc/matrix_transform.hpp>
#include <stb_image.h>

using namespace std;
using namespace glm;

void CityGen::computeChunks() {
    m_chunks = m_voronoiCells;
    unsigned int seed = static_cast<unsigned int>(time(nullptr));
    mt19937 gen(seed);
    std::uniform_int_distribution<int> chanceDist(1, 8);

    float minimumArea = 600.0f;
    float minimumEdgeLength = 50.0f;

    vector<vector<Point>> tempPolygonList;
    for(auto currentChunk : m_chunks) {
        if (PolygonUtils::calculatePolygonArea(currentChunk) <= minimumArea
         || currentChunk.size() <= 3)
        {
            tempPolygonList.push_back(currentChunk);
            continue;
        }

        pair<Point, Point> largestEdge = PolygonUtils::findLargestEdge(currentChunk);
        float edgeLength = (float)PolygonUtils::distanceBetweenPoints(largestEdge.first, largestEdge.second);

        if (edgeLength <= minimumEdgeLength) {
            tempPolygonList.push_back(currentChunk);
            continue;
        }

        uniform_int_distribution<int> splitDist(0, static_cast<int>(currentChunk.size() - 1));
        int index = splitDist(gen);

        Line cut = GeometryUtils::CreateLineFromPoints(currentChunk[index], currentChunk[(index+1)%currentChunk.size()]);
        cut = PolygonUtils::moveLineToCenter(cut, currentChunk);

        auto clipped = PolygonUtils::clipPolygon(currentChunk, cut);
        tempPolygonList.push_back(clipped.first);
        tempPolygonList.push_back(clipped.second);
    }
    m_chunks = tempPolygonList;
    m_voronoiCells = tempPolygonList;
}

void CityGen::sweepToBlocks()
{
    m_blocks.clear();
    
    m_blocks.resize(m_chunks.size());
    m_buildings.clear();
    vector<string> x;
    m_chunks = m_voronoiCells;

    for (size_t i = 0; i < m_chunks.size(); i++) {
        vector<Point> x = m_chunks[i];
        m_chunks[i] = PolygonUtils::scalePolygon(x, 0.8f);
    }
    
    for (size_t i = 0; i < m_chunks.size(); i++) {
        vector<Point> x = m_chunks[i];
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> dis(minMoveAmount, maxMoveAmount);

        pair<Point,Point> largestPointPair = PolygonUtils::findLargestEdge(x);
        double edgeLength = PolygonUtils::distanceBetweenPoints(largestPointPair.first,largestPointPair.second);

        Line largestEdge = GeometryUtils::CreateLineFromPoints(largestPointPair.first,largestPointPair.second);
        Point midpoint = GeometryUtils::findMidpoint(largestPointPair.first,largestPointPair.second);
        Point centerPoint = PolygonUtils::getCentroid(x);

        Point directionVector = GeometryUtils::getPerpendicularDirVector(largestPointPair,centerPoint);
        Point goToPoint = midpoint;
        largestEdge = GeometryUtils::moveToPoint(largestEdge,goToPoint);

        int iterationAmount = 0;
        double eval = GeometryUtils::evaluate(largestEdge, centerPoint);
        bool keepPositiveSide = eval > 0;
        m_strips.clear();

        while(true)
        {
            iterationAmount++;
            moveAmount = (float)dis(gen);
            goToPoint = GeometryUtils::movePointInDirection(goToPoint, directionVector, moveAmount);
            largestEdge = GeometryUtils::moveToPoint(largestEdge,goToPoint);
            auto clipped = PolygonUtils::clipPolygon(x, largestEdge);
            vector<Point> positive = clipped.first;
            vector<Point> negative = clipped.second;

            if((positive.empty() && negative.empty())) {break;}

            if((PolygonUtils::findSmallestEdgeAmount(positive, minEdge) && positive.size() <= 4) 
            ||(PolygonUtils::findSmallestEdgeAmount(negative, minEdge) && negative.size() <= 4))
            {
                m_strips.push_back(x);
                break;
            }
            else if(keepPositiveSide)
            {
                x = positive;
                if(!negative.empty()) {m_strips.push_back(negative);}
            }else
            {
                x = negative;
                if(!positive.empty()) {m_strips.push_back(positive);}
            }
        }

        float distance = (float)edgeLength * .7f;
        iterationAmount = 0;

        directionVector = GeometryUtils::perpendicularVector(directionVector.x, directionVector.y);
        Point oppVector = GeometryUtils::oppositeVector(directionVector.x, directionVector.y);

        goToPoint = centerPoint;
        goToPoint = GeometryUtils::movePointInDirection(goToPoint, oppVector, distance);
        largestEdge = GeometryUtils::makePerpendicularLine(largestEdge);
        largestEdge = GeometryUtils::moveToPoint(largestEdge, goToPoint);

        eval = GeometryUtils::evaluate(largestEdge, centerPoint);
        keepPositiveSide = eval > 0;

        vector<vector<Point>> tempStrips;
        while(!m_strips.empty())
        {
            iterationAmount++;
            moveAmount = (float)dis(gen);
            goToPoint = GeometryUtils::movePointInDirection(goToPoint, directionVector, moveAmount);
            largestEdge = GeometryUtils::moveToPoint(largestEdge, goToPoint);

            for (size_t j = 0; j < m_strips.size(); j++)
            {
                vector<Point> currPolygon = m_strips[j];
                auto clipped = PolygonUtils::clipPolygon(currPolygon, largestEdge);
                vector<Point> positive = clipped.first;
                vector<Point> negative = clipped.second;

                bool positiveMeetsCriteria = (PolygonUtils::findSmallestEdgeAmount(positive, minEdge2) && positive.size() <= 4);
                bool negativeMeetsCriteria = (PolygonUtils::findSmallestEdgeAmount(negative, minEdge2) && negative.size() <= 4);

                if((positiveMeetsCriteria) || (negativeMeetsCriteria))
                {
                    if(keepPositiveSide && !positive.empty())
                    {
                        currPolygon = positive;
                        tempStrips.push_back(currPolygon);
                    }
                    else if(!keepPositiveSide && !negative.empty())
                    {
                        currPolygon = negative;
                        tempStrips.push_back(currPolygon);
                    }
                    else
                    {
                        m_buildings.push_back(currPolygon);
                    }
                }
                else if(keepPositiveSide)
                {
                    if(!positive.empty()) {
                        currPolygon = positive;
                        tempStrips.push_back(currPolygon);
                    }
                    if(!negative.empty()) {
                        if(PolygonUtils::calculatePolygonArea(negative) > minPolygonArea)
                        {
                            m_buildings.push_back(negative);
                        }
                    }
                }
                else
                {
                    if(!negative.empty()) {
                        currPolygon = negative;
                        tempStrips.push_back(currPolygon);
                    }
                    if(!positive.empty()){
                        if(PolygonUtils::calculatePolygonArea(positive) > minPolygonArea)
                        {
                            m_buildings.push_back(positive);
                        }
                    }
                }
            }
            m_strips = tempStrips;
            tempStrips.clear();
        }
    }
}

glm::vec3 CityGen::calculateQuadNormal(const Point& p1, const Point& p2) {
    Point temp = GeometryUtils::getDirectionVector(p1,p2);
    return vec3(temp.y,0.0f,-temp.x);
}

void CityGen::pushVertexData(wolf::VertexBuffer*& vBuffer, wolf::VertexDeclaration*& vDecl, vector<Vertex>& vertices)
{
    vBuffer = wolf::BufferManager::CreateVertexBuffer(vertices.data(), vertices.size() * sizeof(Vertex));
    vDecl = new wolf::VertexDeclaration();
    vDecl->Begin();
    vDecl->AppendAttribute(wolf::AT_Position, 3, wolf::CT_Float);
    vDecl->AppendAttribute(wolf::AT_Color, 4, wolf::CT_UByte);
    vDecl->AppendAttribute(wolf::AT_TexCoord1, 2, wolf::CT_Float);
    vDecl->AppendAttribute(wolf::AT_TexCoord2, 1, wolf::CT_Float);
    vDecl->AppendAttribute(wolf::AT_Normal, 3, wolf::CT_Float);
    vDecl->SetVertexBuffer(vBuffer);
    vDecl->End();
}

void CityGen::buildVertexData() {
    m_vertices.clear();

    for (size_t i = 0; i < m_buildings.size(); i++) {
        vector<Point> x = m_buildings[i];
        m_buildings[i] = PolygonUtils::scalePolygon(x, 0.85f);
    }

    vector<GLubyte> colors = {
    255, 0, 0,
    0, 255, 0,
    0, 0, 255,
    255, 255, 0,
    255, 0, 255,
    0, 255, 255,
    191, 0, 191,
    255, 140, 0,
    0, 200, 0,
    200, 200, 200
    };

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> heightDist(30.0f, 100.0f); 
    std::uniform_int_distribution<int> chanceDist(1, 30);

    vec3 topNormal = vec3(0.0f,1.0f,0.0f);
    vector<Vertex> topVertices;
    for (size_t i = 0; i < m_buildings.size(); ++i) {
        topVertices.clear();
        const vector<Point>& cell = m_buildings[i];

        float wallHeight = (chanceDist(gen) == 1) ? 120.0f : heightDist(gen);
        GLubyte r = colors[(i * 3) % colors.size()];
        GLubyte g = colors[(i * 3 + 1) % colors.size()];
        GLubyte b = colors[(i * 3 + 2) % colors.size()];
        GLubyte a = 255;
        float buildingLayer = wall2;
        if(i % 2 == 0) {
            buildingLayer = wall;
        }

        float texWidth = 10.0f; 
        float texHeight = 10.0f;

        for (size_t j = 0; j < cell.size(); ++j) {
            Point p1 = cell[j];
            Point p2 = cell[(j + 1) % cell.size()];
            float width = (float)PolygonUtils::distanceBetweenPoints(p1,p2);
            float tileU = width / texWidth;
            float tileV = wallHeight / texHeight;
            glm::vec3 normal = calculateQuadNormal(p1, p2);

            Vertex v1 = { (GLfloat)p1.x, 0.0f, (GLfloat)p1.y, r, g, b, a,
            0.0f,0.0f,
            buildingLayer,
            normal.x, normal.y, normal.z};
            Vertex v2 = { (GLfloat)p2.x, 0.0f, (GLfloat)p2.y, r, g, b, a,
            tileU,0.0f,
            buildingLayer,
            normal.x, normal.y, normal.z};

            Vertex top1 = { (GLfloat)p1.x, wallHeight, (GLfloat)p1.y, r, g, b, a,
            0.0f,tileV,
            buildingLayer,
            normal.x, normal.y, normal.z};
            Vertex top2 = { (GLfloat)p2.x, wallHeight, (GLfloat)p2.y, r, g, b, a,
            tileU,tileV,
            buildingLayer,
            normal.x, normal.y, normal.z};

            m_vertices.push_back(v1);
            m_vertices.push_back(top1);
            m_vertices.push_back(v2);

            m_vertices.push_back(top1);
            m_vertices.push_back(top2);
            m_vertices.push_back(v2);

            topVertices.push_back(top1);
        }
        if(topVertices.empty()){continue;}
        
        auto roofVerts = PolygonUtils::fanTriangulatePolygon(cell, topNormal, wallHeight, 3.0f, r, g, b, a); //3 = sidewalk
        m_vertices.insert(m_vertices.end(), roofVerts.begin(), roofVerts.end());
    }

    m_lines.clear();
    for (size_t i = 0; i < m_voronoiCells.size(); ++i) {
        const vector<Point>& cell = m_voronoiCells[i];     
        GLubyte r = 255;
        GLubyte g = 255;
        GLubyte b = 255;
        GLubyte a = 255;
        for (size_t j = 0; j < cell.size(); ++j) {
            Point p1 = cell[j];
            Point p2 = cell[(j + 1) % cell.size()];
            // Add road decals
            DrawRoad::addRoadDecals(m_lines, p1, p2, sidewalk,
                                     PolygonUtils::distanceBetweenPoints,
                                     GeometryUtils::getDirectionVector,
                                     GeometryUtils::movePointInDirection);
        }
    }

    for (auto &cell : m_voronoiCells) {
        if (cell.size() < 3) continue;
        vector<Point> sidewalkPoly = PolygonUtils::scalePolygon(cell, 0.87f);
        GLubyte r=200,g=200,b=200,a=255;
        auto sidewalkVerts = PolygonUtils::fanTriangulatePolygon(sidewalkPoly, vec3((0.0f,1.0f,0.0f)), 0.01f, 3.0f, r, g, b, a);
        m_vertices.insert(m_vertices.end(), sidewalkVerts.begin(), sidewalkVerts.end());
    }

    m_numVertices = (int)m_vertices.size();
    pushVertexData(m_vertexBuffer,m_vertexDecl,m_vertices);
    m_numLines = (int)m_lines.size();
    pushVertexData(m_lineBuffer,m_lineDecl,m_lines);

    m_arrayTexture = wolf::TextureManager::CreateAutoArrayTexture({
        "data/building1.tga",
        "data/building2.tga",
        "data/building3.tga",
        "data/brickFloor.tga",
        "data/building1Windows.tga",
        "data/building2Windows.tga",
        "data/asphalt.tga",
    });

    if(!m_arrayTexture) {
        printf("Error: Array texture creation failed.\n");
        return;
    }
}

void CityGen::generate(GLuint program) {
    m_program = wolf::ProgramManager::CreateProgram("data/cube.vsh", "data/cube.fsh");
    unsigned int seed = (unsigned int)time(nullptr);
    m_sites = Voronoi::generateSites(numSites, seed, minX, maxX, minY, maxY);
    Voronoi::computeVoronoiDiagram(m_sites, m_voronoiCells, minX, maxX, minY, maxY);
    // computeChunks(); // Left commented as per original code
    sweepToBlocks();
    buildVertexData();
}

void CityGen::deGenerate() {
    if (m_program) {
        wolf::ProgramManager::DestroyProgram(m_program);
        m_program = nullptr;
    }
    if (m_buildingTexture) {
        wolf::TextureManager::DestroyTexture(m_buildingTexture);
        m_buildingTexture = nullptr;
    }
    if (m_arrayTexture) {
        wolf::TextureManager::DestroyTexture(m_arrayTexture);
        m_arrayTexture = nullptr;
    }
    if(m_lineBuffer || m_vertexBuffer) {
        wolf::BufferManager::DestroyBuffer(m_lineBuffer);
        wolf::BufferManager::DestroyBuffer(m_vertexBuffer);
    }

    m_vertices.clear();
    debugs.clear();
    debugs1.clear();
    debugs2.clear();
    m_sites.clear();
    m_voronoiCells.clear();
    m_buildings.clear();
}

void CityGen::reGenerate() {
    deGenerate();
    generate(0);
}

void CityGen::render(const glm::mat4& proj, const glm::mat4& view, float m_timer) {
    if (!m_program || !m_arrayTexture)
    {
        return;
    }
    m_program->Bind();
    m_program->SetUniform("projection", proj);
    m_program->SetUniform("view", view);
    m_program->SetUniform("world", glm::mat4(1.0f));
    m_program->SetUniform("u_texture", 0);
    m_arrayTexture->Bind(0);

    float hours = fmod(m_timer*2,24.0f);
    float normalizedTime = hours / 24.0f;
    float angle = (normalizedTime * 360.0f) - 180.0f;
    float rad = glm::radians(angle);
    rad += glm::radians(90.0f);
    glm::vec3 lightDir(cos(rad), sin(rad), 0.0f);
    float elevation = sin(rad);
    float dayFactor = glm::clamp(elevation,0.0f,1.0f);
    float ambientStrength = 0.2f + 0.2f * dayFactor; 
    glm::vec3 dynamicAmbient(ambientStrength);
    glm::vec3 dynamicLightColor = glm::mix(glm::vec3(0.5f,0.5f,0.7f), glm::vec3(1.0f,1.0f,0.9f), dayFactor);

    m_program->SetUniform("u_time", m_timer);
    m_program->SetUniform("u_lightDir", glm::normalize(lightDir));
    m_program->SetUniform("u_lightColor", dynamicLightColor);
    m_program->SetUniform("u_ambient", dynamicAmbient);

    glm::vec3 dayColor = glm::vec3(0.5f, 0.7f, 1.0f);
    glm::vec3 nightColor = glm::vec3(0.0f, 0.0f, 0.1f);
    glm::vec3 backgroundColor = glm::mix(nightColor, dayColor, dayFactor);
    glClearColor(backgroundColor.r, backgroundColor.g, backgroundColor.b, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    float desiredLitProbability = 0.3f;
    float u_windowLitFactor = (1.0f - dayFactor) * desiredLitProbability;
    u_windowLitFactor = glm::clamp(u_windowLitFactor, 0.0f, 1.0f);
    m_program->SetUniform("u_windowLitFactor", u_windowLitFactor);

    m_vertexDecl->Bind();
    glDrawArrays(GL_TRIANGLES, 0, m_numVertices);

    m_lineDecl->Bind();
    glDrawArrays(GL_LINES, 0, m_numLines);
    glBindVertexArray(0);
    glUseProgram(0);
}
