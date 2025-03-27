#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iostream>

struct Sphere {
    std::string name;
    double pos[3];
    double scale[3];
    double colour[3];
    double Ka, Kd, Ks, Kr;
    int specExp;
};

struct Light{
    std::string name;
    double pos[3];
    double intensity[3];
};

//scene data
double nearPlane;
double leftPlane, rightPlane, bottomPlane, topPlane;
int resX, resY;
std::vector<Sphere> spheres;
std::vector<Light> lights;
double background[3];
double ambient[3];
std::string outputFilename;

void parseInputFile(const std::string &filename);
bool finfNearestHit(const double origin[3], const double dir[3], double &thit, int &sphereIndex,
                     double hitPoint[3], double normal[3], bool isPrimaryRay);



void traceRay(const double origin[3], const double dir[3], int recursionDepth, bool isPrimaryRay, double colourOut[3]){
    if(recursionDepth > 3){
        colourOut[0]=colourOut[1]=colourOut[2] = 0.0;
        return;
    }
}