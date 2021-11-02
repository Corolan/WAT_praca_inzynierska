#include <cstdio>
#include <cstdlib>
#include <cstring>

#define SETIERROR( err, errmsg ) do { \
    exit (1); \
} while (0)

#define BAD_DECODE 1

typedef float sah_complex[2];

#define MEM_ALIGN 128

void *malloc_a(size_t size, size_t alignment) {
  void *palignedMem;
  if (posix_memalign(&palignedMem,alignment,size)) palignedMem=NULL;
  return palignedMem;
}

int decode(unsigned char* bin, int nbytes, FILE* f) {
  unsigned char buf[256], *p, c0, c1, c2;
  int i, n, m, nleft = nbytes, offset=0;
  int nbadlines = 0;
  while (1) {
    if (nleft <= 0) break;
    p = (unsigned char*)fgets((char*)buf, 256, f);
    if (!p) {
      SETIERROR(BAD_DECODE, "file ended too soon");
    }
    n = (int)strlen((char*)buf)-1;
    if (n%4) {
      nbadlines++;
      if (nbadlines == 100) {
        SETIERROR(BAD_DECODE, "too many bad lines - rejecting file");
      }
      n = 64;
      fprintf(stderr, "encoded line has bad length: %s\n", buf);
    }
    m = (n/4)*3;
    if (m > nleft+2) {
      n = ((nleft+2)/3)*4;
    }
    p = buf;
    for (i=0; i<n/4; i++) {
      p[0] -= 0x20;
      p[1] -= 0x20;
      p[2] -= 0x20;
      p[3] -= 0x20;
      c0 = (p[0]&0x3f) | p[1]<<6; // 6 + 2
      c1 = p[1]>>2 | p[2]<<4; // 4 + 4
      c2 = p[2]>>4 | p[3]<<2; // 2 + 6
      p += 4;
      bin[offset++] = c0;
      bin[offset++] = c1;
      bin[offset++] = c2;
      nleft -= 3;
    }
  }
  return 0;
}

void bits_to_floats(unsigned char* raw, sah_complex* data, int nsamples, int bits_per_sample) {
  int i, j, k=0;
  unsigned char c;
  const float v2[2]={-1.0,1.0};
  const float v4[4]={-3.3358750,-1.0,1.0,3.3358750};
  switch (bits_per_sample) {
    case 2:
        for (i=0; i<nsamples/4; i++) {
            j = (i&1) ? i-1 : i+1;
            c = raw[j];
            for (j=0; j<4; j++) {
                data[k][0] = v2[(c>>1)&1];
                data[k][1] = v2[(c&1)];
                k++;
                c >>= 2;
            }
        }
        break;
    case 4:
        for (i=0; i<nsamples/2; i++) {
            j = (i&1) ? i-1 : i+1;
            c = raw[j];
            for (j=0; j<2; j++) {
                data[k][0] = v4[(c>>2)&3];
                data[k][1] = v4[c&3];
                k++;
                c >>= 4;
            }
        }
        break;
    case 8:
        for (i=0; i<nsamples; i++) {
            j = (i&1) ? i-1 : i+1;
            c = raw[j];
            data[k][0] = static_cast<float>((c>>4)&15)-7.5f;
            data[k][1] = static_cast<float>(c&15)-7.5f;
            k++;
        }
        break;
    case 16:
        for (i=0; i<nsamples*2; i+=2) {
            signed char sc = raw[i+1];
            data[k][0] = static_cast<float>(sc);
            sc = raw[i];
            data[k][1] = static_cast<float>(sc);
            k++;
        }
        break;
    default: 
        fprintf(stderr,"Unsupported bit depth (%d)\n",bits_per_sample);
        throw BAD_DECODE;
  }
}

int main (int argc, char **argv) {

    unsigned int nsamples = 1048576;

    //size_t buffsize = nsamples / 4 + 1 + 2; 
    size_t buffsize = nsamples / 2  + 2;
  
    unsigned char *buffer = (unsigned char*) malloc (buffsize);
    sah_complex *complex_data = (sah_complex *)malloc_a (nsamples * sizeof (sah_complex), MEM_ALIGN);

    FILE *f = fopen ("work_unit.data", "r");
    decode (buffer, buffsize - 4, f);
    //decode (buffer, buffsize - 2, f);

    //bits_to_floats (buffer, complex_data, nsamples, 2);//here was a mistake. I hardcoded 2 as "bits per sample"
    bits_to_floats (buffer, complex_data, nsamples, 4); //but with this, i get "Segmentation fault"

    sah_complex *c = complex_data;
    for (unsigned int i = 0; i < nsamples; ++i, ++c) {
        printf ("%f %f\n", (*c)[0], (*c)[1]); 
    }

    return 0;
}
