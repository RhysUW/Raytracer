#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include "ppm.cpp"
#include "invert.cpp"

// Data structures to store scene input

static const int MAX_SPHERES = 15;
static const int MAX_LIGHTS = 10;

struct Sphere
{
    std::string name;
    double pos[3];
    double scl[3];
    double color[3];
    double ka;
    double kd;
    double ks;
    double kr;
    int n;
    //store the sphere's forward transform and its inverse
    double M[4][4];
    double Minv[4][4];
    double MinvT[4][4];
};

struct Light
{
    std::string name;
    double pos[3];
    double intensity[3];
};

// Global scene parameters
static double g_nearPlane = 1.0;
static double g_left = -1.0, g_right = 1.0, g_bottom = -1.0, g_top = 1.0;
static int g_resX = 400, g_resY = 400;
static int g_numSpheres = 0, g_numLights = 0;
static Sphere g_spheres[MAX_SPHERES];
static Light g_lights[MAX_LIGHTS];
static double g_bgColor[3] = {0.0, 0.0, 0.0};      // background color
static double g_ambIntensity[3] = {0.2, 0.2, 0.2}; // ambient intensity
static char g_outputName[128] = {"output.ppm"};

// Max recursion depth for reflections
static const int MAX_BOUNCES = 3;

// Utility clamp double to [0,1]
static inline double clamp01(double x)
{
    if (x < 0.0)
        return 0.0;
    if (x > 1.0)
        return 1.0;
    return x;
}

// for each sphere we have:
//    - Translate to sphere's pos
//    - Scale by sphere's scl
static void buildSphereMatrix(Sphere &sph)
{

    printf("[buildSphereMatrix] Sphere %s\n", sph.name.c_str());
    fflush(stdout); 
    // Initialize to identity
    double M[4][4] = {
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    };
    printf("  After identity matrix.\n");
    fflush(stdout);

    // Scale
    M[0][0] = sph.scl[0];
    M[1][1] = sph.scl[1];
    M[2][2] = sph.scl[2];

    printf("  scale done: scl=%f,%f,%f\n", sph.scl[0], sph.scl[1], sph.scl[2]);
    fflush(stdout);

    // Then translate in a separate matrix T
    double T[4][4] = {
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    };
    T[0][3] = sph.pos[0];
    T[1][3] = sph.pos[1];
    T[2][3] = sph.pos[2];

    printf("  translate done: pos=%f,%f,%f\n", sph.pos[0], sph.pos[1], sph.pos[2]);
    fflush(stdout);


    // Combine final M = T * S
    double res[4][4];

    printf("  about to multiply T*S\n");
    fflush(stdout);

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            res[i][j] = 0.0;
            for (int k = 0; k < 4; k++)
            {
                res[i][j] += T[i][k] * M[k][j];
            }
        }
    }

    printf("  multiply done.\n");
    fflush(stdout);

    // Copy res into sph.M
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            sph.M[i][j] = res[i][j];
        }
    }
    
    printf("  M copied to sphere.\n");
    fflush(stdout);

    printf("  about to invert.\n");
    fflush(stdout);

    // Now invert that to get Minv
    invert_matrix(sph.M, sph.Minv);

    printf("  about to invert.\n");
    fflush(stdout);

    double tmp[4][4];
    // Transpose of Minv
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            tmp[i][j] = sph.Minv[j][i];
        }
    }

    printf("  MinvT done.\n");
    fflush(stdout);

    // store that in MinvT
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            sph.MinvT[i][j] = tmp[i][j];
        }
    }

     //Print out final M for debugging
     printf("Matrix M for sphere %s:\n", sph.name.c_str());
     for(int r=0; r<4; r++){
         printf("   ");
         for(int c=0; c<4; c++){
             printf("%f ", sph.M[r][c]);
         }
         printf("\n");
     }
     fflush(stdout);
}

// Ray-sphere intersection, using the sphere's inverse transform to transform
// the ray into object space, where the sphere is the "unit sphere" at origin.
// Return the t-value in world space. If no intersection, return -1.
static double intersectSphere(const Sphere &sph, const double eyeW[3], const double dirW[3])
{
    // Build hom coords for eye, direction
    double eye4[4] = {eyeW[0], eyeW[1], eyeW[2], 1};
    double dir4[4] = {dirW[0], dirW[1], dirW[2], 0};

    // Transform them using Minv
    double eyeL[4] = {0, 0, 0, 0};
    double dirL[4] = {0, 0, 0, 0};

    // localEye = Minv * eye4
    for (int i = 0; i < 4; i++)
    {
        eyeL[i] = sph.Minv[i][0] * eye4[0] + sph.Minv[i][1] * eye4[1] + sph.Minv[i][2] * eye4[2] + sph.Minv[i][3] * eye4[3];
    }
    // localDir = Minv * dir4
    for (int i = 0; i < 4; i++)
    {
        dirL[i] = sph.Minv[i][0] * dir4[0] + sph.Minv[i][1] * dir4[1] + sph.Minv[i][2] * dir4[2] + sph.Minv[i][3] * dir4[3];
    }

    // Now we have a local-space ray: eyeL + t*dirL
    // Intersection with unit sphere => solve: ||eyeL + t*dirL||^2 = 1
    double ox = eyeL[0], oy = eyeL[1], oz = eyeL[2];
    double dx = dirL[0], dy = dirL[1], dz = dirL[2];

    // Quadratic coefficients
    double A = dx * dx + dy * dy + dz * dz;
    double B = 2.0 * (ox * dx + oy * dy + oz * dz);
    double C = (ox * ox + oy * oy + oz * oz) - 1.0;

    double disc = B * B - 4.0 * A * C;

    if (disc < 0.0){
        return -1.0;
    }/*else{
        printf("Sphere %s: disc=%f (>=0) => possible intersection\n", 
           sph.name.c_str(), disc); // for debugging
    }*/

    double sqrtDisc = sqrt(disc);
    double t1 = (-B - sqrtDisc) / (2.0 * A);
    double t2 = (-B + sqrtDisc) / (2.0 * A);

    double tLocal = -1.0;
    if (t1 > t2 && t2 > 0.0){
        tLocal = t2;
    }else{
        tLocal = t1;
    }

    // no valid intersection in front of local eye
    if (tLocal < 0.0){
        return -1.0; 
    }
   

    double hitLocal[4] = {ox + tLocal * dx, oy + tLocal * dy, oz + tLocal * dz, 1.0};
    double hitWorld[4] = {0, 0, 0, 0};
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            hitWorld[i] += sph.M[i][j] * hitLocal[j];
        }
    }
    // The vector from eyeW to hitWorld:
    double vw[3] = {hitWorld[0] - eyeW[0],
                    hitWorld[1] - eyeW[1],
                    hitWorld[2] - eyeW[2]};
    // The direction dirW has length squared:
    double dirWLen2 = dirW[0] * dirW[0] + dirW[1] * dirW[1] + dirW[2] * dirW[2];
    // tWorld = (vw . dirW) / (dirW . dirW)
    double dotVD = vw[0] * dirW[0] + vw[1] * dirW[1] + vw[2] * dirW[2];
    double tWorld = dotVD / dirWLen2;

    if (tWorld < 0.0){
        return -1.0;
    }
    return tWorld;
}

// Compute the sphere's surface normal at a world-space point pW
static void getSphereNormal(const Sphere &sph, const double pW[3], double outNormalW[3])
{
    // Convert pW into local space
    double p4[4] = {pW[0], pW[1], pW[2], 1.0};
    double pLocal[4] = {0, 0, 0, 0};
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            pLocal[i] += sph.Minv[i][j] * p4[j];
        }
    }
    // local normal = pLocal (since it's a unit sphere at origin).
    // Then transform it by inverse-transpose for the sphere.
    // Then normalize in world space.
    double nLocal[4] = {pLocal[0], pLocal[1], pLocal[2], 0.0};

    double nWorld4[4] = {0, 0, 0, 0};
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            nWorld4[i] += sph.MinvT[i][j] * nLocal[j];
        }
    }
    // outNormalW = normalized(nWorld4)
    double nx = nWorld4[0], ny = nWorld4[1], nz = nWorld4[2];
    double len = sqrt(nx * nx + ny * ny + nz * nz);
    if (len < 1e-12)
    {
        // fallback
        outNormalW[0] = 0;
        outNormalW[1] = 1;
        outNormalW[2] = 0;
        return;
    }
    outNormalW[0] = nx / len;
    outNormalW[1] = ny / len;
    outNormalW[2] = nz / len;
}

// Scene trace function that returns color for a ray from eye at (0,0,0) with direction dirW
// bouncing up to bounceLevel times
static void traceRay(const double eyeW[3], const double dirW[3], int bounceLevel,
                     double outColor[3])
{
    outColor[0] = g_bgColor[0];
    outColor[1] = g_bgColor[1];
    outColor[2] = g_bgColor[2];

    double nearestT = 1e15;
    int hitSphereIdx = -1;

    // 1) Find nearest sphere intersection
    for (int i = 0; i < g_numSpheres; i++)
    {
        double t = intersectSphere(g_spheres[i], eyeW, dirW);
        if (t > 1e-6 && t < nearestT)
        {
            // debug: print t
            //printf("Sphere %s => t=%f\n", g_spheres[i].name.c_str(), t);
            nearestT = t;
            hitSphereIdx = i;
        }
    }
    if (hitSphereIdx < 0)
    {
        // no intersection => background
        return;
    }

    //We have an intersection
    Sphere &sph = g_spheres[hitSphereIdx];
    //Intersection point in world
    double hitPt[3] = {
        eyeW[0] + nearestT * dirW[0],
        eyeW[1] + nearestT * dirW[1],
        eyeW[2] + nearestT * dirW[2]};
    //Normal
    double normalW[3];
    getSphereNormal(sph, hitPt, normalW);

    //View vector is from hitPt toward eye => v = -dirW
    double viewVec[3] = {-dirW[0], -dirW[1], -dirW[2]};

    //Ambient contribution
    double outR = sph.ka * g_ambIntensity[0] * sph.color[0];
    double outG = sph.ka * g_ambIntensity[1] * sph.color[1];
    double outB = sph.ka * g_ambIntensity[2] * sph.color[2];

    //For each light => do shadow test, add diffuse & spec
    for (int li = 0; li < g_numLights; li++)
    {
        const Light &L = g_lights[li];
        //Light vector = L.pos - hitPt
        double Lvec[3] = {L.pos[0] - hitPt[0],
                          L.pos[1] - hitPt[1],
                          L.pos[2] - hitPt[2]};
        //distance to light
        double Ldist = sqrt(Lvec[0] * Lvec[0] + Lvec[1] * Lvec[1] + Lvec[2] * Lvec[2]);
        if (Ldist < 1e-12)
            continue;

        //Normalize
        double invLd = 1.0 / Ldist;
        Lvec[0] *= invLd;
        Lvec[1] *= invLd;
        Lvec[2] *= invLd;

        //cast a ray from hitPt + epsilon along Lvec
        double shadowEye[3] = {hitPt[0] + 1e-5 * Lvec[0],
                               hitPt[1] + 1e-5 * Lvec[1],
                               hitPt[2] + 1e-5 * Lvec[2]};
        //We want to see if we intersect ANY sphere at t>0 and t< Ldist
        bool inShadow = false;
        for (int si = 0; si < g_numSpheres; si++)
        {
            double tS = intersectSphere(g_spheres[si], shadowEye, Lvec);
            if (tS > 0.0 && tS < (Ldist - 1e-4))
            {
                inShadow = true;
                break;
            }
        }
        if (inShadow)
            continue; //no diffuse or spec for this light

        //Diffuse = kd * (N dot L) * color * LightIntensity
        double ndotl = normalW[0] * Lvec[0] + normalW[1] * Lvec[1] + normalW[2] * Lvec[2];
        if (ndotl < 0.0)
            ndotl = 0.0;
        outR += sph.kd * ndotl * sph.color[0] * L.intensity[0];
        outG += sph.kd * ndotl * sph.color[1] * L.intensity[1];
        outB += sph.kd * ndotl * sph.color[2] * L.intensity[2];

        //Specular = ks * (R dot V)^n * LightIntensity
        double minusL[3] = {-Lvec[0], -Lvec[1], -Lvec[2]};
        double ndotL = ndotl;
        double Rvec[3] = {minusL[0] + 2.0 * ndotL * normalW[0],
                          minusL[1] + 2.0 * ndotL * normalW[1],
                          minusL[2] + 2.0 * ndotL * normalW[2]};
        //Now R dot V
        double rdotv = Rvec[0] * viewVec[0] + Rvec[1] * viewVec[1] + Rvec[2] * viewVec[2];
        if (rdotv < 0.0)
            rdotv = 0.0;
        double specPow = pow(rdotv, (double)sph.n);
        outR += sph.ks * specPow * L.intensity[0];
        outG += sph.ks * specPow * L.intensity[1];
        outB += sph.ks * specPow * L.intensity[2];
    }

    //Reflection (Kr * colorFromReflectionRay)
    if (bounceLevel < MAX_BOUNCES && sph.kr > 0.0)
    {
        //reflect the view direction about normal
        double D[3] = {-viewVec[0], -viewVec[1], -viewVec[2]};
        double dDotN = D[0] * normalW[0] + D[1] * normalW[1] + D[2] * normalW[2];
        double R[3] = {
            D[0] - 2.0 * dDotN * normalW[0],
            D[1] - 2.0 * dDotN * normalW[1],
            D[2] - 2.0 * dDotN * normalW[2]};
        //cast reflection ray from hitPt + epsilon*R
        double reflEye[3] = {
            hitPt[0] + 1e-5 * R[0],
            hitPt[1] + 1e-5 * R[1],
            hitPt[2] + 1e-5 * R[2]};
        double reflColor[3] = {0, 0, 0};
        traceRay(reflEye, R, bounceLevel + 1, reflColor);
        outR += sph.kr * reflColor[0];
        outG += sph.kr * reflColor[1];
        outB += sph.kr * reflColor[2];
    }

    //printf("Final color for pixel => R=%.2f G=%.2f B=%.2f\n", outR, outG, outB);
    // clamp final
    outColor[0] = clamp01(outR);
    outColor[1] = clamp01(outG);
    outColor[2] = clamp01(outB);
}

//Parse the input file
static void parseFile(const char *filename)
{
    FILE *fp = fopen(filename, "r");
    if (!fp)
    {
        fprintf(stderr, "Could not open file '%s'\n", filename);
        exit(1);
    }
    char token[256];
    while (true)
    {
        if (1 != fscanf(fp, "%s", token))
        {
            //EOF or error
            break;
        }
        if (!strcmp(token, "NEAR"))
        {
            fscanf(fp, "%lf", &g_nearPlane);
        }
        else if (!strcmp(token, "LEFT"))
        {
            fscanf(fp, "%lf", &g_left);
        }
        else if (!strcmp(token, "RIGHT"))
        {
            fscanf(fp, "%lf", &g_right);
        }
        else if (!strcmp(token, "BOTTOM"))
        {
            fscanf(fp, "%lf", &g_bottom);
        }
        else if (!strcmp(token, "TOP"))
        {
            fscanf(fp, "%lf", &g_top);
        }
        else if (!strcmp(token, "RES"))
        {
            fscanf(fp, "%d %d", &g_resX, &g_resY);
        }
        else if (!strcmp(token, "SPHERE"))
        {
            //SPHERE <name> <px> <py> <pz> <sx> <sy> <sz> <r> <g> <b> <ka> <kd> <ks> <kr> <n>
            Sphere &s = g_spheres[g_numSpheres++];
            char namebuf[128];
            fscanf(fp, "%s", namebuf);
            s.name = namebuf;
            fscanf(fp, "%lf %lf %lf", &s.pos[0], &s.pos[1], &s.pos[2]);
            fscanf(fp, "%lf %lf %lf", &s.scl[0], &s.scl[1], &s.scl[2]);
            fscanf(fp, "%lf %lf %lf", &s.color[0], &s.color[1], &s.color[2]);
            fscanf(fp, "%lf", &s.ka);
            fscanf(fp, "%lf", &s.kd);
            fscanf(fp, "%lf", &s.ks);
            fscanf(fp, "%lf", &s.kr);
            fscanf(fp, "%d", &s.n);

            buildSphereMatrix(s);
        }
        else if (!strcmp(token, "LIGHT"))
        {
            //LIGHT <name> <px> <py> <pz> <ir> <ig> <ib>
            Light &L = g_lights[g_numLights++];
            char namebuf[128];
            fscanf(fp, "%s", namebuf);
            L.name = namebuf;
            fscanf(fp, "%lf %lf %lf", &L.pos[0], &L.pos[1], &L.pos[2]);
            fscanf(fp, "%lf %lf %lf", &L.intensity[0], &L.intensity[1], &L.intensity[2]);
        }
        else if (!strcmp(token, "BACK"))
        {
            fscanf(fp, "%lf %lf %lf", &g_bgColor[0], &g_bgColor[1], &g_bgColor[2]);
        }
        else if (!strcmp(token, "AMBIENT"))
        {
            fscanf(fp, "%lf %lf %lf", &g_ambIntensity[0], &g_ambIntensity[1], &g_ambIntensity[2]);
        }
        else if (!strcmp(token, "OUTPUT"))
        {
            fscanf(fp, "%s", g_outputName);
        }
        else
        {
            //unknown token, skip
            fprintf(stderr, "Skipping unknown token: %s\n", token);
        }
    }
    fclose(fp);
}

int main(int argc, char **argv)
{
    printf("beginning of main\n");
    if (argc < 2)
    {
        printf("Usage: %s inputFile.txt\n", argv[0]);
        return 0;
    }
    parseFile(argv[1]);
    printf("file parsed\n");

    // Allocate an RGB buffer for the final image
    unsigned char *pixels = new unsigned char[3 * g_resX * g_resY];

    double invNX = 1.0 / (double)g_resX;
    double invNY = 1.0 / (double)g_resY;

    int k = 0; // index into pixels
    for (int row = 0; row < g_resY; row++)
    {
        // pixel row goes from top=0 to bottom=g_resY-1
        //define "ty" from top to bottom:
        double sy = g_top - (row + 0.5) * ((g_top - g_bottom) / (double)g_resY);

        for (int col = 0; col < g_resX; col++)
        {
            double sx = g_left + (col + 0.5) * ((g_right - g_left) / (double)g_resX);

            // The ray from (0,0,0) to (sx, sy, -g_nearPlane)
            double dir[3] = {sx, sy, -g_nearPlane};
            double len = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
            if (len < 1e-12)
            {
                dir[0] = 0;
                dir[1] = 0;
                dir[2] = -1; 
            }
            else
            {
                dir[0] /= len;
                dir[1] /= len;
                dir[2] /= len;
            }
            double eye[3] = {0, 0, 0};
            double color[3];
            traceRay(eye, dir, 0, color);

            // store to pixels
            pixels[k] = (unsigned char)(255.0 * clamp01(color[0]));
            pixels[k + 1] = (unsigned char)(255.0 * clamp01(color[1]));
            pixels[k + 2] = (unsigned char)(255.0 * clamp01(color[2]));
            k += 3;
        }
    }

    // Save image
    save_imageP3(g_resX, g_resY, g_outputName, pixels);

    delete[] pixels;

    printf("Done. Saved image to '%s'\n", g_outputName);
    return 0;
}