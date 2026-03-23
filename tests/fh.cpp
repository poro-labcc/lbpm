
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <chrono>
#include <iostream>
#include "common/Array.h"

static inline float intersection(int q, int vk, float fq, float fvk)
{
        float qq = (float)q  * (float)q;
        float vv = (float)vk * (float)vk;
        return ((fq + qq) - (fvk + vv)) / (2.0 * ((float)q - (float)vk));
}

void edt_1d(const int *f, int *g, int n)
{
    int *v = (int *) malloc((size_t)n * sizeof(int));
    float *z = (float *) malloc((size_t)(n + 1) * sizeof(float));

    int k = 0;
    v[0] = 0;
    z[0] = -INFINITY;
    z[1] =  INFINITY;

    for (int q = 1; q < n; q++) {
        float s;

        while (1) {
            int vk = v[k];
            s = intersection(q,vk,f[q],f[vk]);

            if (s <= z[k]) {
                k--;
                if (k < 0) {
                    k = 0;
                    break;
                }
            } else {
                break;
            }
        }

        if (k == 0) {
            int vk = v[k];
            s = intersection(q,vk,f[q],f[vk]);

            if (s <= z[k]) {
                v[0] = q;
                z[0] = -INFINITY;
                z[1] =  INFINITY;
                continue;
            }
        }

        k++;
        v[k] = q;
        z[k] = s;
        z[k + 1] = INFINITY;
    }

    k = 0;
    for (int x = 0; x < n; x++) {
        while (z[k + 1] < (float)x) {
            k++;
        }

        int vk = v[k];
        int dx = x - vk;
        g[x] = dx * dx + f[vk];
    }

    free(v);
    free(z);
}

void process_line(int* edt2, const size_t first, const size_t stride, int n)
{

            IntArray f(n);
            IntArray g(n);

            size_t pos = first;
            for (int i = 0; i < n; i++) {
                f(i) = edt2[pos];
                pos += stride;
            }

            edt_1d(f.data(), g.data(), n);

            pos = first;
            for (int i = 0; i < n; i++) {
                edt2[pos] = g(i);
                pos += stride;
            }
}

void edt_3d( unsigned char target, IntArray &image, IntArray &distance2)
{

    const int nx = image.size(0);
    const int ny = image.size(1);
    const int nz = image.size(2);

    size_t nvox =  image.length();
    int BIG = static_cast<float>(nx*nx + ny*ny + nz*nz) + 1;

    int* img = image.data();
    int* edt2 = distance2.data();

    for (int i = 0; i < (int) nvox; i++) {
        edt2[i] = (img[i] == target) ? 0.0 : BIG;
    }

    size_t maxdim = std::max( std::max(nz, ny), nx );

    for (int z = 0; z < nz; z++) {
        for (int y = 0; y < ny; y++) {
            process_line( edt2 , nx * (y + z * ny) , 1, nx );
        }
    }

    for (int z = 0; z < nz; z++) {
        for (int x = 0; x < nx; x++) {
            process_line( edt2 , z * nx * ny + x ,  nx, ny);
        }
    }

    for (int y = 0; y < ny; y++) {
        for (int x = 0; x < nx; x++) {
            process_line( edt2 , y * nx + x , nx * ny, nz);
        }
    }
}

#include <iostream>

using namespace std;

void generate_random_binary(IntArray& A,
                            double prob_one)
{
     size_t n = A.length();
     int* data = A.data();

      for (size_t i = 0; i < n; i++) {
          double r = (double)rand() / (double)RAND_MAX;
          data[i] = (r < prob_one) ? 1 : 0;
      }
}


int main()
{
        int samples = 1;
        double tFH = 0;

        int nx = 300;
        int ny = 300;
        int nz = 300;

        IntArray  img (nx,ny,nz);
        IntArray dist2( img.size() );

        for (int k =0; k < samples; k++)
        {
            generate_random_binary(img, 0.8);
            auto t0 = std::chrono::steady_clock::now();
            edt_3d( 0, img, dist2);
            auto t1 = std::chrono::steady_clock::now();
            tFH += std::chrono::duration<double>(t1 - t0).count();
        }

        std::cout << "Nx x Ny x Nz: " << nx << " x " << ny << " x " << nz <<  " pixels\n";
        std::cout << "Tempo EDT FH = " << tFH/samples << " s\n";

        FILE* f = fopen("distance.raw","w");
        fwrite( dist2.data()  , sizeof(int), dist2.length() , f);
        fclose(f);

        f = fopen("geo.raw","w");
        fwrite(img.data() , sizeof(int), img.length() , f);
        fclose(f);

        return 0;
}
