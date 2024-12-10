#include "CityGen.h"
using namespace std; using namespace glm;

Building CityGen::determineBuildingDetails(const vector<Point>& polygon)
{
    return m_buildingGen->generateBuildingDetails(polygon, districtCenter);
}

void CityGen::BuildingToVerticies(const Building& building, vector<Vertex>& m_vertices, float ground)
{
    if (building.polygons.size() < 3 || building.height <= 0.0f) { 
        return;
    }

    GLubyte r = building.r; GLubyte g = building.g; GLubyte b = building.b; GLubyte a = 255;

    float wallHeight = building.isSpecial ? building.height / 2 : building.height;
    float buildingLayer = building.textureLayer;
    vec3 topNormal = vec3(0.0f,1.0f,0.0f);
    float texWidth = 10.0f; 
    float texHeight = 10.0f;

    // Create walls
    for (size_t j = 0; j < building.polygons.size(); ++j) {
        Point p1 = building.polygons[j];
        Point p2 = building.polygons[(j + 1) % building.polygons.size()];
        float width = (float)PolygonUtils::distanceBetweenPoints(p1,p2);
        float tileU = width / texWidth;
        float tileV = wallHeight / texHeight;
        glm::vec3 normal = GeometryUtils::calculateQuadNormal(p1, p2);

        Vertex v1 = { (GLfloat)p1.x, ground, (GLfloat)p1.y, r, g, b, a,
        0.0f,0.0f,
        buildingLayer,
        normal.x, normal.y, normal.z};

        Vertex v2 = { (GLfloat)p2.x, ground, (GLfloat)p2.y, r, g, b, a,
        tileU,0.0f,
        buildingLayer,
        normal.x, normal.y, normal.z};

        Vertex top1 = { (GLfloat)p1.x, (ground + wallHeight), (GLfloat)p1.y, r, g, b, a,
        0.0f,tileV,
        buildingLayer,
        normal.x, normal.y, normal.z};

        Vertex top2 = { (GLfloat)p2.x, (ground + wallHeight), (GLfloat)p2.y, r, g, b, a,
        tileU,tileV,
        buildingLayer,
        normal.x, normal.y, normal.z};

        m_vertices.push_back(v1);
        m_vertices.push_back(top1);
        m_vertices.push_back(v2);

        m_vertices.push_back(top1);
        m_vertices.push_back(top2);
        m_vertices.push_back(v2);
    }

    // Create roof
    auto roofVerts = PolygonUtils::fanTriangulatePolygon(building.polygons, topNormal, (ground + wallHeight), roof, r, g, b, a);
    m_vertices.insert(m_vertices.end(), roofVerts.begin(), roofVerts.end());

    if (building.isSpecial) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<int> countDist(extraBuildingCountMin, extraBuildingCountMax);
        int extraCount = countDist(gen);

        float currentGround = ground + wallHeight;

        // Start from the building's original polygon
        vector<Point> topPolygon = building.polygons;

        for (int i = 0; i < extraCount; ++i) { //extraCount
        
            float area = PolygonUtils::calculatePolygonArea(topPolygon);
            if(area * extraBuildingScaleFactor < 5.0f) {
        cout << "Skipping building: Area too small after scaling.\n";return;}
            topPolygon = PolygonUtils::scalePolygon(topPolygon, extraBuildingScaleFactor);

            Building extraBuilding;
            extraBuilding.polygons = topPolygon;
            extraBuilding.height = building.height/2 * extraBuildingHeightMultiplier; // taller than original
            extraBuilding.textureLayer = building.textureLayer; 
            extraBuilding.isSpecial = false; 
            extraBuilding.district = building.district; // same district and style
            extraBuilding.r = r;
            extraBuilding.g = g;
            extraBuilding.b = b;

            // cout << "Building #" << i + 1 << " Details After Scaling:\n";
            // cout << "  Scaled Area: " << PolygonUtils::calculatePolygonArea(topPolygon) << "\n";
            // cout << "  Scaled Top Polygon Vertices:\n";
            // for (const auto& point : topPolygon) {
            //     cout << "    (" << point.x << ", " << point.y << ")\n";
            // }
            // cout << "  Height: " << extraBuilding.height << "\n";
            // cout << "  Texture Layer: " << extraBuilding.textureLayer << "\n";
            // cout << "  District: " << extraBuilding.district << "\n";
            // cout << "  Is Special: " << (extraBuilding.isSpecial ? "Yes" : "No") << "\n";

            
            // cout<<"A: "<<i<<endl;
            // Each extra building sits on top of the last
            BuildingToVerticies(extraBuilding, m_vertices, currentGround);
            // cout<<"B: "<<i<<endl;
            currentGround += extraBuilding.height; // Increase ground for the next building
        }
    }
}

void CityGen::sweepToBlocks()
{
    m_blocks.clear();
    
    m_blocks.resize(m_chunks.size());
    m_buildings.clear();
    
    m_chunks = m_voronoiCells;

    for (size_t i = 0; i < m_chunks.size(); i++) {
        vector<Point> x = m_chunks[i];
        m_chunks[i] = PolygonUtils::scalePolygon(x, 0.8f);
    }
    
    for (size_t i = 0; i < m_chunks.size(); i++) {
        //choose chunk district based on a random point within min and max
        //
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
                        Building b = determineBuildingDetails(currPolygon);
                        m_buildings.push_back(b);
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
                            Building b = determineBuildingDetails(negative);
                            m_buildings.push_back(b);
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
                            Building b = determineBuildingDetails(positive);
                            m_buildings.push_back(b);
                        }
                    }
                }
            }
            m_strips = tempStrips;
            tempStrips.clear();
        }
    }
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
        vector<Point> x = m_buildings[i].polygons;
        m_buildings[i].polygons = PolygonUtils::scalePolygon(x, 0.85f);
    }

    for(auto& b : m_buildings){
        BuildingToVerticies(b, m_vertices, 0.0f);
    }

    m_lines.clear(); //Roads
    for (size_t i = 0; i < m_voronoiCells.size(); ++i) {
        const vector<Point>& cell = m_voronoiCells[i];
        for (size_t j = 0; j < cell.size(); ++j) {
            Point p1 = cell[j];
            Point p2 = cell[(j + 1) % cell.size()]; 
            DrawRoad::addRoadDecals(m_lines, p1, p2, sidewalk,
                                     PolygonUtils::distanceBetweenPoints,
                                     GeometryUtils::getDirectionVector,
                                     GeometryUtils::movePointInDirection);
        }
    }

    for (auto &cell : m_voronoiCells) { //Sidewalk
        if (cell.size() < 3) continue;
        vector<Point> sidewalkPoly = PolygonUtils::scalePolygon(cell, 0.87f);
        auto sidewalkVerts = PolygonUtils::fanTriangulatePolygon(sidewalkPoly, vec3((0.0f,1.0f,0.0f)), 0.1f, sidewalk, 255,255,255,255);
        m_vertices.insert(m_vertices.end(), sidewalkVerts.begin(), sidewalkVerts.end());
    }

    {// GROUND FLOOR
        int div = 15;
        int padding = 10;
        float minUVX = minX / div;
        float maxUVX = maxX / div;
        float minUVY = minY / div;
        float maxUVY = maxY / div;
        Point bottomLeft = { minX-padding, minY-padding };
        Point bottomRight = { maxX+padding, minY-padding };
        Point topLeft = { minX-padding, maxY+padding };
        Point topRight = { maxX+padding, maxY+padding };
        float ground = -0.1f;
        float quadLayer = asphalt;
        glm::vec3 groundNormal = vec3(0.0f, 1.0f, 0.0f);
        GLubyte a,r,g,b = 255;
        int sex = 0;
        r=sex;g=sex;b=sex;a=255;
        
        Vertex bl = { (GLfloat)bottomLeft.x, ground, (GLfloat)bottomLeft.y, r, g, b, a,
            minUVX, minUVY, quadLayer, groundNormal.x, groundNormal.y, groundNormal.z };

        Vertex br = { (GLfloat)bottomRight.x, ground, (GLfloat)bottomRight.y, r, g, b, a,
            maxUVX, minUVY, quadLayer, groundNormal.x, groundNormal.y, groundNormal.z };

        Vertex tl = { (GLfloat)topLeft.x, ground, (GLfloat)topLeft.y, r, g, b, a,
            minUVX, maxUVY, quadLayer, groundNormal.x, groundNormal.y, groundNormal.z };

        Vertex tr = { (GLfloat)topRight.x, ground, (GLfloat)topRight.y, r, g, b, a,
            maxUVX, maxUVY, quadLayer, groundNormal.x, groundNormal.y, groundNormal.z };

        // Add quad vertices (two triangles)
        m_vertices.push_back(bl);
        m_vertices.push_back(tl);
        m_vertices.push_back(br);

        m_vertices.push_back(tl);
        m_vertices.push_back(tr);
        m_vertices.push_back(br);
    }


    m_numVertices = (int)m_vertices.size();
    pushVertexData(m_vertexBuffer,m_vertexDecl,m_vertices); 
    m_numLines = (int)m_lines.size();
    pushVertexData(m_lineBuffer,m_lineDecl,m_lines);

    m_arrayTexture = wolf::TextureManager::CreateAutoArrayTexture
    ({  "data/building1.tga",
        "data/building2.tga",
        "data/building3.tga",
        "data/building1Windows.tga",
        "data/building2Windows.tga",
        "data/building3Windows.tga",
        "data/ambigiousRoof.tga",
        "data/brickFloor.tga",
        "data/asphalt.tga" });

    if(!m_arrayTexture) {
        printf(" Array texture creation failed\n"); 
    }
}

void CityGen::generate(GLuint program) {
    m_program = wolf::ProgramManager::CreateProgram("data/mainShader.vsh", "data/mainShader.fsh");
    unsigned int seed = (unsigned int)time(nullptr);

    { //District point.
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> distX(minX, maxX);
        std::uniform_real_distribution<double> distY(minY, maxY);
        districtCenter = {distX(gen), distY(gen)};
        
        double width = maxX - minX;
        double height = maxY - minY;
        districtRadius = (float)((width + height)/2.0); 
    }

    {//BuildingGen init
        m_buildingGen = new BuildingGen(districtRadius, minPolygonArea,
        IndustrialMinStories, IndustrialMaxStories,
        CommercialMinStories, CommercialMaxStories,
        ResidentialMinStories, ResidentialMaxStories,
        baseStoryHeight, baseAddition0, baseAddition2, specialBuildingChance);

    }

    m_sites = Voronoi::generateSites(numSites, seed, minX, maxX, minY, maxY);
    Voronoi::computeVoronoiDiagram(m_sites, m_voronoiCells, minX, maxX, minY, maxY);
    // computeChunks();
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
    if (m_buildingGen) {
        delete m_buildingGen;
        m_buildingGen = nullptr;
    }
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

    #pragma region Day & NIght

    float hours = fmod(m_timer*2,cyceLength);
    float normalizedTime = hours / cyceLength;

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
    float u_windowLitFactor = ((1.0f - dayFactor) * desiredLitProbability) - 0.05f;
    u_windowLitFactor = glm::clamp(u_windowLitFactor, 0.0f, 1.0f);
    m_program->SetUniform("u_windowLitFactor", u_windowLitFactor);
    
    #pragma endregion

    m_vertexDecl->Bind();
    glDrawArrays(GL_TRIANGLES, 0, m_numVertices);

    m_lineDecl->Bind();
    glDrawArrays(GL_LINES, 0, m_numLines); //crosswalk lines
    glBindVertexArray(0);
    glUseProgram(0);
}
