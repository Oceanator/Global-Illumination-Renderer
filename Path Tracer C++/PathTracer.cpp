// g++ -O3 -fopenmp smallpt.cpp -o smallpt 
// Remove the "-fopenmp" flag to disable threading support

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

struct Vec
{
    double x, y, z; // position, also color (r,g,b)
	
    Vec(double x_=0, double y_=0, double z_=0){ x=x_; y=y_; z=z_; }
    Vec operator+(const Vec &b) const { return Vec(x+b.x,y+b.y,z+b.z); }
    Vec operator-(const Vec &b) const { return Vec(x-b.x,y-b.y,z-b.z); }
    Vec operator*(double b) const { return Vec(x*b,y*b,z*b); }
    Vec mult(const Vec &b) const { return Vec(x*b.x,y*b.y,z*b.z); }
    Vec& norm(){ return *this = *this * (1/sqrt(x*x+y*y+z*z)); } // normalize a vector
    double dot(const Vec &b) const { return x*b.x+y*b.y+z*b.z; } // dot product
    Vec operator%(Vec&b){return Vec(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);} // cross product
};

struct Ray // A ray is a parametric line with an origin o and a direction d
{ 
    Vec o, d; 
    Ray(Vec o_, Vec d_) : o(o_), d(d_) {} 
};

enum Refl_t { DIFF, SPEC, REFR };  // material types, used in radiance()

struct Sphere 
{
    double rad;       // radius
    Vec p, e, c;      // position, emission, color
    Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive)

    // constructor
    Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_):
        rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {}

    // returns distance, 0 if no hit
    double intersect(const Ray &r) const 
    {
        // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
        Vec op = p-r.o; 
        double t, eps = 1e-4;
        double b = op.dot(r.d);
        double det = b*b - op.dot(op) + rad * rad;
        if (det < 0) return 0; // ray misses sphere
        else det = sqrt(det); 
        return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0); // return smaller positive t
    }
};

// // Scene box
// Sphere spheres[] = {
// 	//Scene: radius, position, emission, color, material
//   	Sphere(1e5, Vec( 1e5+1,40.8,81.6), Vec(),Vec(.75,.25,.25),DIFF), //Left
//     Sphere(1e5, Vec(-1e5+99,40.8,81.6),Vec(),Vec(.25,.25,.75),DIFF), //Rght
//     Sphere(1e5, Vec(50,40.8, 1e5),     Vec(),Vec(.75,.75,.75),DIFF), //Back
//     Sphere(1e5, Vec(50,40.8,-1e5+170), Vec(),Vec(),           DIFF), //Frnt
//     Sphere(1e5, Vec(50, 1e5, 81.6),    Vec(),Vec(.75,.75,.75),DIFF), //Botm
//     Sphere(1e5, Vec(50,-1e5+81.6,81.6),Vec(),Vec(.75,.75,.75),DIFF), //Top
//     Sphere(16.5,Vec(27,16.5,47),       Vec(),Vec(1,1,1)*.999, SPEC), //Mirr
//     Sphere(16.5,Vec(73,16.5,78),       Vec(),Vec(1,1,1)*.999, REFR), //Glas
//     Sphere(600, Vec(50,681.6-.27,81.6),Vec(12,12,12),  Vec(), DIFF)  //Lite
// };

// Scene sky
Vec Cen(50,40.8,-860);
Sphere spheres[] = {
    Sphere(1600, Vec(1,0,2)*3000, Vec(1,.9,.8)*1.2e1*1.56*2,Vec(), DIFF), // sun
    Sphere(1560, Vec(1,0,2)*3500,Vec(1,.5,.05)*4.8e1*1.56*2, Vec(),  DIFF), // horizon sun2
    Sphere(10000,Cen+Vec(0,0,-200), Vec(0.0627, 0.188, 0.569)*6e-2*8, Vec(.7,.7,1)*.25,  DIFF), // sky
    Sphere(10000,Cen+Vec(0,0,-200), Vec(0.00063842, 0.02001478, 0.28923243)*6e-2*8, Vec(.7,.7,1)*.25,  DIFF), // sky

    Sphere(100000, Vec(50, -100000, 0),  Vec(),Vec(.3,.3,.3),DIFF), // grnd
    Sphere(110000, Vec(50, -110048.5, 0),  Vec(.9,.5,.05)*4,Vec(),DIFF),// horizon brightener
    Sphere(4e4, Vec(50, -4e4-30, -3000),  Vec(),Vec(.2,.2,.2),DIFF),// mountains
    Sphere(3.99e4, Vec(50, -3.99e4+20.045, -3000),  Vec(),Vec(.7,.7,.7),DIFF),// mountains snow

    Sphere(26.5,Vec(22,26.5,42),   Vec(),Vec(1,1,1)*.596, SPEC), // white Mirr
    Sphere(13,Vec(75,13,82),   Vec(),Vec(.96,.96,.96)*.96, REFR),// Glas
    Sphere(22,Vec(87,22,24),   Vec(),Vec(.6,.6,.6)*.696, REFR)    // Glas2
};

// convert colors (output of radiance) to displayable range
inline double clamp(double x){ return x<0 ? 0 : x>1 ? 1 : x; }
inline int toInt(double x){ return int(pow(clamp(x),1/2.2)*255+.5); } // gamma correction of 2.2

// check each sphere, one at a time. Keep the closest intersection
inline bool intersect(const Ray &r, double &t, int &id)
{
    double n = sizeof(spheres)/sizeof(Sphere);
    double d; 
    double inf = t = 1e20;
    for(int i = int(n); i--;)
    {
        if((d = spheres[i].intersect(r))&&d<t)
        {
            t = d; 
            id = i;
        }
    }
    return t < inf;
}

// compute the radiance estimate along ray
// return value - Vec the radiance estimate
// r - the ray we are casting
// depth - the ray depth  
// Xi - random number seed
// E - whether to include emissive color
Vec radiance(const Ray &r, int depth, unsigned short *Xi){
    double t;                               // distance to intersection
    int id = 0;                               // id of intersected object
    if (!intersect(r, t, id)) return Vec(); // if miss, return black
    const Sphere &obj = spheres[id];        // the hit object

    // surface properties
    Vec x = r.o+r.d*t; // ray intersection point
    Vec n = (x-obj.p).norm(); // sphere normal
    Vec nl = n.dot(r.d) < 0 ? n : n*-1; // property oriented surface normal - for refractions
    Vec f = obj.c; // object color (BRDF modulator)

    // use maximum reflectivity amount for Russian Roulette
    double p = f.x>f.y && f.x>f.z ? f.x : f.y>f.z ? f.y : f.z; // use the maximum component (r,g,b) of the surface color.
    if (++depth>5) if ((rand()/(double)RAND_MAX)<p) f=f*(1/p); else return obj.e;  // don't do Russian Roulette until after depth 5

    // ideal diffuse reflection
    if (obj.refl == DIFF)
    {                  
        double r1 = 2*M_PI*(rand()/(double)RAND_MAX); // angle around
        double r2 = (rand()/(double)RAND_MAX), r2s = sqrt(r2); // distance from center
        
        Vec w = nl; // w = normal
        Vec u = ((fabs(w.x)>.1 ? Vec(0,1) : Vec(1))%w).norm(); // u is perpendicular to w
        Vec v = w % u; // v is perpendicular to u and w
        // sampling unit hemisphere
        Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm(); // d is random reflection ray
        
        return obj.e + f.mult(radiance(Ray(x,d),depth,Xi));
    } else if (obj.refl == SPEC){
        return obj.e + f.mult(radiance(Ray(x,r.d-n*2*n.dot(r.d)),depth,Xi)); // Ideal SPECULAR reflection
    }

    Ray reflRay(x, r.d-n*2*n.dot(r.d));     // Ideal dielectric REFRACTION
    bool into = n.dot(nl)>0;                // Ray from outside going in?
    double nc = 1, nt = 1.5, nnt = into ? nc/nt : nt/nc, ddn = r.d.dot(nl), cos2t;

    if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0)    // Total internal reflection
        return obj.e + f.mult(radiance(reflRay,depth,Xi));

    // otherwise choose reflection or  refraction
    Vec tdir = (r.d*nnt - n*((into ? 1 : -1)*(ddn*nnt + sqrt(cos2t)))).norm();
    double a = nt-nc, b = nt+nc, R0 = a*a/(b*b), c = 1-(into ? -ddn : tdir.dot(n));
    double Re = R0+(1-R0)*c*c*c*c*c,Tr=1-Re, P = .25+.5*Re, RP = Re/P, TP = Tr/(1-P);
    return obj.e + f.mult(depth>2 ? ((rand()/(double)RAND_MAX) < P ?   // Russian roulette
            radiance(reflRay,depth,Xi)*RP:radiance(Ray(x,tdir),depth,Xi)*TP) :
            radiance(reflRay,depth,Xi)*Re+radiance(Ray(x,tdir),depth,Xi)*Tr);
}

int main(int argc, char *argv[])
{
    // set up image
    int w = 1024, h = 768;
    int samps = argc == 2 ? atoi(argv[1])/4 : 1; // # samples

    // set up camera
    Ray cam(Vec(50,52,295.6), Vec(0,-0.042612,-1).norm()); // camera pos, dir
    // x direction increment (uses implicit 0 for y, z)
    // Assume upright camera 
    // 0.5135 defines FOV
    Vec cx = Vec(w*.5135/h);  // horizontal camera vector direction
    // cross product gets vector perpendicular to both cx and gaze direction
    Vec cy = (cx%cam.d).norm()*.5135; // vertical camera direction

    Vec r; // used for colors of samples
    Vec *c = new Vec[w*h]; // the image

    #pragma omp parallel for schedule(dynamic, 1) private(r) // OpenMP - states each loop should be run in its own thread

    // creat image by looping over all image pixels
    for (int y=0; y<h; y++) // Loop over image rows
    {                       
        fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps*4,100.*y/(h-1)); // print progress
        for (unsigned short x=0, Xi[3]={0,0,y*y*y}; x<w; x++)   // Loop cols
        {
            // for each pixel do 2X2 subsamples, and samps samples per subsample. The subpixel color will be averaged
            // calculate array index for (x, y): i=(h-y-1)*w+x
            for (int sy=0, i=(h-y-1)*w+x; sy<2; sy++)     // 2x2 subpixel rows
            {
                for (int sx=0; sx<2; sx++, r=Vec()) 	  // 2x2 subpixel cols
                {        
                    for (int s=0; s<samps; s++) 
                    {
                        // r1 and r2 are random values of a tent filter
                        // Determine location of sample within pixel
                        double r1=2*(rand()/(double)RAND_MAX), dx=r1<1 ? sqrt(r1)-1: 1-sqrt(2-r1);
                        double r2=2*(rand()/(double)RAND_MAX), dy=r2<1 ? sqrt(r2)-1: 1-sqrt(2-r2);
                        // compute ray direction
                        Vec d = cx*( ( (sx+.5 + dx)/2 + x)/w - .5) +
                                cy*( ( (sy+.5 + dy)/2 + y)/h - .5) + cam.d;
                        r = r + radiance(Ray(cam.o+d*140,d.norm()),0,Xi)*(1./samps); // estimate radiance
                    } // Camera rays are pushed ^^^^^ forward to start in interior
                    // add the gamma-corrected subpixel color estimate to the Pixel color c[i]
                    c[i] = c[i] + Vec(clamp(r.x),clamp(r.y),clamp(r.z))*.25;
                }
            }
        }
    }

    FILE *f = fopen("image.ppm", "w"); // Write image to PPM file.
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    for (int i=0; i<w*h; i++)
    {
        fprintf(f,"%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
    }
}


