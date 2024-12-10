#include "DrawRoad.h"

void DrawRoad::drawCrosswalk(std::vector<Vertex>& m_lines, Point base, Point direction, float lineLength, float spacing, float lineHeight, float sidewalk) {
    GLubyte r = 255, g = 255, b = 255, a = 255;
    Point perpendicular = GeometryUtils::perpendicularVector(direction.x, direction.y);
    for (float offset = -(spacing * 5); offset <= spacing * 5; offset += spacing) {
        Point lineStart = GeometryUtils::movePointInDirection(base, perpendicular, offset);
        Point lineEnd = GeometryUtils::movePointInDirection(lineStart, direction, lineLength);
        Vertex v1 = { static_cast<GLfloat>(lineStart.x), lineHeight, static_cast<GLfloat>(lineStart.y), r, g, b, a,
        0.0f,0.0f,
        sidewalk,
        0.0f,1.0f,0.0f };
        Vertex v2 = { static_cast<GLfloat>(lineEnd.x), lineHeight, static_cast<GLfloat>(lineEnd.y), r, g, b, a,
        1.0f,1.0f,
        sidewalk,
        0.0f,1.0f,0.0f };
        m_lines.push_back(v1);
        m_lines.push_back(v2);
    }
}

void DrawRoad::addRoadDecals(std::vector<Vertex>& m_lines, Point p1, Point p2, 
                              float sidewalk,
                              std::function<double(const Point&, const Point&)> distanceBetweenPoints,
                              std::function<Point(const Point&, const Point&)> getDirectionVector,
                              std::function<Point(const Point&, const Point&, double)> movePointInDirection)
{
    float lineHeight = 0.01f;
    float lineLength = 5.0f;
    float gapLength = 8.0f;
    float intersection = 15.0f;
    float crosswalkWidth = 8.0f;
    float crosswalkSpacing = 1.5f;
    float crosswalkLineLength = 5.0f;
    float minRoadLength = 30.0f;

    float totalLength = (float)distanceBetweenPoints(p1, p2);
    if ((totalLength * 0.65) <= (intersection * 2) || totalLength <= minRoadLength) {
        return;
    }

    Point p1Dir = getDirectionVector(p1, p2);
    Point p2Dir = getDirectionVector(p2, p1);

    Point startCrosswalk = movePointInDirection(p1, p1Dir, -crosswalkWidth);
    Point endCrosswalk = movePointInDirection(p2, p2Dir, -crosswalkWidth);

    p1 = movePointInDirection(p1, p1Dir, intersection);
    p2 = movePointInDirection(p2, p2Dir, intersection);

    // Add crosswalk at the start
    drawCrosswalk(m_lines, p1, p1Dir, crosswalkLineLength, crosswalkSpacing, lineHeight, sidewalk);
    // Add crosswalk at the end
    drawCrosswalk(m_lines, p2, p2Dir, crosswalkLineLength, crosswalkSpacing, lineHeight, sidewalk);

    totalLength = (float)distanceBetweenPoints(p1, p2);

    float segmentLength = lineLength + gapLength;
    int numSegments = static_cast<int>(totalLength / segmentLength);

    GLubyte r = 255, g = 255, b = 255, a = 255;
    Point direction = {
        (p2.x - p1.x) / totalLength,
        (p2.y - p1.y) / totalLength
    };

    for (int i = 0; i <= numSegments; ++i) {
        float startDist = i * segmentLength;
        float endDist = startDist + lineLength;
        if (startDist >= totalLength) break;
        if (endDist > totalLength) endDist = totalLength;
        Point interpStart = { p1.x + startDist * direction.x, p1.y + startDist * direction.y };
        Point interpEnd = { p1.x + endDist * direction.x, p1.y + endDist * direction.y };
        Vertex v1 = { static_cast<GLfloat>(interpStart.x), lineHeight, static_cast<GLfloat>(interpStart.y), r, g, b, a,
        0.0f,0.0f,
        sidewalk,
        0.0f,1.0f,0.0f};
        Vertex v2 = { static_cast<GLfloat>(interpEnd.x), lineHeight, static_cast<GLfloat>(interpEnd.y), r, g, b, a,
        1.0f,1.0f,
        sidewalk,
        0.0f,1.0f,0.0f };
        m_lines.push_back(v1);
        m_lines.push_back(v2);
    }
}
