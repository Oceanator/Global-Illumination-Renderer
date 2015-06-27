using System;
using System.IO;

namespace Path_Tracer_CS
{
    class Program
    {
        struct Vec
        {
            public double x, y, z;  // position,also color (r,g,b)

            public Vec(double x_, double y_, double z_) { x = x_; y = y_; z = z_; }
            public static Vec operator +(Vec a, Vec b) { return new Vec(a.x + b.x, a.y + b.y, a.z + b.z); }
            public static Vec operator -(Vec a, Vec b) { return new Vec(a.x - b.x, a.y - b.y, a.z - b.z); }
            public static Vec operator *(Vec a, double b) { return new Vec(a.x * b, a.y * b, a.z * b); }
            public Vec mult(Vec b) { return new Vec(x * b.x, y * b.y, z * b.z); }
            public Vec norm() { return this = this * (1 / Math.Sqrt(x * x + y * y + z * z)); }
            public double dot(Vec b) { return x * b.x + y * b.y + z * b.z; }
            public static Vec operator %(Vec a, Vec b) { return new Vec(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x); }
        }


        struct Ray // A ray is a parametric line with an origin o and a direction d
        { 
            public Vec o, d; 
            public Ray(Vec o_, Vec d_) { o = o_; d = d_; } 
        }

        enum Refl_t { DIFF, SPEC, REFR }; // material types, used in radiance()

        class Sphere
        {
            public double rad;      // radius
            public Vec p, e, c;     // position,emission,color
            public Refl_t refl;     // reflection type (DIFFuse,SPECular,REFRactive)
            
            public Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_)
            {
                rad = rad_; p = p_; e = e_; c = c_; refl = refl_;
            }

            // returns distance, 0 if no hit
            public double intersect(Ray r) 
            { 
                Vec op = p - r.o; // Solve t^2*d.d+2*t*(o-p).d+(o-p).(o-p)-R^2=0
                double t, eps = 1e-4, b = op.dot(r.d), det = b * b - op.dot(op) + rad * rad;

                if (det < 0) // ray  misses sphere
                    return 0; 
                else 
                    det = Math.Sqrt(det);

                return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
            }
        };

        class smallpt
        {
            static Random random = new Random();
            static double erand48() { return random.NextDouble(); }
            //Scene: radius,position,emission,color,material
            static Sphere[] spheres ={
                new Sphere(1e5,new Vec( 1e5+1,40.8,81.6), new Vec(),new Vec(.75,.25,.25),Refl_t.DIFF),//Left
                new Sphere(1e5,new Vec(-1e5+99,40.8,81.6),new Vec(),new Vec(.25,.25,.75),Refl_t.DIFF),//Rght
                new Sphere(1e5,new Vec(50,40.8,1e5),      new Vec(),new Vec(.75,.75,.75),Refl_t.DIFF),//Back
                new Sphere(1e5,new Vec(50,40.8,-1e5+170), new Vec(),new Vec(),           Refl_t.DIFF),//Frnt
                new Sphere(1e5,new Vec(50,1e5,81.6),      new Vec(),new Vec(.75,.75,.75),Refl_t.DIFF),//Botm
                new Sphere(1e5,new Vec(50,-1e5+81.6,81.6),new Vec(),new Vec(.75,.75,.75),Refl_t.DIFF),//Top
                new Sphere(16.5,new Vec(27,16.5,47),      new Vec(),new Vec(1,1,1)*.999, Refl_t.SPEC),//Mirr
                new Sphere(16.5,new Vec(73,16.5,78),      new Vec(),new Vec(1,1,1)*.999, Refl_t.REFR),//Glas
                new Sphere(600,new Vec(50,681.6-.27,81.6),new Vec(12,12,12), new Vec(),  Refl_t.DIFF) //Lite
            };

            // convert colors (output of radiance) to displayable range
            static double clamp(double x) { return x < 0 ? 0 : x > 1 ? 1 : x; }
            static int toInt(double x) { return (int)(Math.Pow(clamp(x), 1 / 2.2) * 255 + .5); } // gamma correction of 2.2

            // check each sphere, one at a time. Keep the closest intersection
            static bool intersect(Ray r, ref double t, ref int id)
            {
                double d, inf = t = 1e20;
                for (int i = spheres.Length - 1; i >= 0; i--)
                    if ((d = spheres[i].intersect(r)) != 0 && d < t) 
                    { 
                        t = d; 
                        id = i; 
                    }
                return t < inf;
            }

            // compute the radiance estimate along ray
            // return value - Vec the radiance estimate
            // r - the ray we are casting
            // depth - the ray depth  
            // Xi - random number seed
            // E - whether to include emissive color
            static Vec radiance(Ray r, int depth)
            {
                double t = 0;                            // distance to intersection
                int id = 0;                              // id of intersected object
                if (!intersect(r, ref t, ref id)) return new Vec(); // if miss,return black
                Sphere obj = spheres[id];                // the hit object

                // surface properties
                Vec x = r.o + r.d * t; // ray intersection point
                Vec n = (x - obj.p).norm(); // sphere normal
                Vec nl = n.dot(r.d) < 0 ? n : n * -1; // property oriented surface normal - for refractions
                Vec f = obj.c; // object color (BRDF modulator)

                // use maximum reflectivity amount for Russian Roulette
                double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z; // use the maximum component (r,g,b) of the surface color.
                if (++depth > 5) if (erand48() < p) f = f * (1 / p); else return obj.e; // Russian Roulette
                if (depth > 100) return obj.e;

                // Ideal diffuse reflection
                if (obj.refl == Refl_t.DIFF)
                {
                    double r1 = 2 * Math.PI * erand48(), r2 = erand48(); // angle around
                    double r2s = Math.Sqrt(r2); // distance from center

                    Vec w = nl; // w = normals
                    Vec u = ((Math.Abs(w.x) > .1 ? new Vec(0, 1, 0) : new Vec(1, 0, 0)) % w).norm();
                    Vec v = w % u; // v is perpendicular to u and w
                    // sampling unit hemisphere
                    Vec d = (u * Math.Cos(r1) * r2s + v * Math.Sin(r1) * r2s + w * Math.Sqrt(1 - r2)).norm(); // d is random reflection ray
                    return obj.e + f.mult(radiance(new Ray(x, d), depth));
                }
                else if (obj.refl == Refl_t.SPEC)           // Ideal SPECULAR reflection
                {
                    return obj.e + f.mult(radiance(new Ray(x, r.d - n * 2 * n.dot(r.d)), depth));
                }

                Ray reflRay = new Ray(x, r.d - n * 2 * n.dot(r.d));    // IdealdielectricREFRACTION
                bool into = n.dot(nl) > 0;                             // Ray from outside going in?
                double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.d.dot(nl), cos2t;

                if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0)    // Total internal reflection
                    return obj.e + f.mult(radiance(reflRay, depth));

                // otherwise choose reflection or  refraction
                Vec tdir = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + Math.Sqrt(cos2t)))).norm();
                double a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : tdir.dot(n));
                double Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re, P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);
                return obj.e + f.mult(depth > 2 ? (erand48() < P ?  // Russian roulette
                  radiance(reflRay, depth) * RP : radiance(new Ray(x, tdir), depth) * TP) :
                  radiance(reflRay, depth) * Re + radiance(new Ray(x, tdir), depth) * Tr);
            }

            public static void Main(string[] args)
            {
                DateTime start = DateTime.Now;

                // set up image
                int w = 256, h = 256;
                int samps = args.Length == 2 ? int.Parse(args[1]) / 4 : 25;// # samples

                // set up camera
                Ray cam = new Ray(new Vec(50, 52, 295.6), new Vec(0, -0.042612, -1).norm()); //cam pos,dir

                // x direction increment (uses implicit 0 for y, z)
                // Assume upright camera 
                // 0.5135 defines FOV
                Vec cx = new Vec(w * .5135 / h, 0, 0); // horizontal camera vector direction
                // cross product gets vector perpendicular to both cx and gaze direction
                Vec cy = (cx % cam.d).norm() * .5135; // vertical camera direction

                Vec r; // used for colors of samples
                Vec[] c = new Vec[w * h]; // the image

                // creat image by looping over all image pixels
                for (int y = 0; y < h; y++)
                {                        
                    Console.Write("\rRendering ({0}spp) {1:F2}%", samps * 4, 100.0 * y / (h - 1)); // print progress
                    for (int x = 0; x < w; x++)   // Loop cols
                    {
                        // for each pixel do 2X2 subsamples, and samps samples per subsample. The subpixel color will be averaged
                        // calculate array index for (x, y): i=(h-y-1)*w+x
                        for (int sy = 0, i = (h - y - 1) * w + x; sy < 2; sy++)     // 2x2 subpixel rows
                        {
                            for (int sx = 0; sx < 2; sx++)                          // 2x2 subpixel cols
                            {
                                r = new Vec();
                                for (int s = 0; s < samps; s++)
                                {
                                    // r1 and r2 are random values of a tent filter
                                    // Determine location of sample within pixel
                                    double r1 = 2 * erand48(), dx = r1 < 1 ? Math.Sqrt(r1) - 1 : 1 - Math.Sqrt(2 - r1);
                                    double r2 = 2 * erand48(), dy = r2 < 1 ? Math.Sqrt(r2) - 1 : 1 - Math.Sqrt(2 - r2);
                                    // compute ray direction
                                    Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
                                          cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;
                                    r = r + radiance(new Ray(cam.o + d * 140, d.norm()), 0) * (1.0 / samps); // estimate radiance
                                } // Camera rays are pushed ^^^^^ forward to start in interior
                                // add the gamma-corrected subpixel color estimate to the Pixel color c[i]
                                c[i] = c[i] + new Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * .25;
                            }
                        }
                    }
                }

                Console.WriteLine("\n{0} sec", (DateTime.Now - start).TotalSeconds);
                using (StreamWriter sw = new StreamWriter("image.ppm"))
                {
                    sw.Write("P3\r\n{0} {1}\r\n{2}\r\n", w, h, 255);
                    for (int i = 0; i < w * h; i++)
                        sw.Write("{0} {1} {2}\r\n", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
                    sw.Close();
                }
            }
        }
    }
}
