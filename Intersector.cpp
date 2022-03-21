#include<Intersector.h>

void
Intersector::FindTrianglesInSubdomainScope()

void
Intersector::BuildKDTree()

void
Intersector::FindOcculudedNodes()

void
Intersector::FindIntersections()

void
Intersector::TagFirstLayerNodes();

void
Intersector::FloodFillStatus(); //by default, set sign = 1 (0 for occuluded) only do floodfill for sign = -1

//optional
void
Intersector::FindShortestDistanceForFirstLayer //don't forget to subtract half_distance

void
Intersector::FindShortestDistanceForOtherNodes //call level set reinitializer
