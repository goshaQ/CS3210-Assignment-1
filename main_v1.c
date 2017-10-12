#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>
#include <math.h>

typedef struct {
    unsigned char *data;
    float *dst;
    int width;
    int height;
    float sigma;
    int blur_size;
    int n;
} gaussian_blur_arg;

typedef struct {
    unsigned char *src;
    float *dst;
    float *kernel;
    int ksize;
    int width;
    int height;
    int start;
    int stop;
} blur_helper_arg;

int read_BMP(char* filename, unsigned char *info, unsigned char **dataR, unsigned char **dataG, unsigned char **dataB, int *size, int *width, int *height, int *offset, int *row_padded)
{
    int i = 0, j, k, read_bytes, h, w, o, p;
    unsigned char *data;

    FILE* f = fopen(filename, "rb");

    if(f == NULL)
    {
        printf ("Invalid filename: %s\n", filename);
        return -1;
    }


    read_bytes = fread(info, sizeof(unsigned char), 54, f); // read the 54-byte header
    if (read_bytes != 54)
    {
        printf ("Error at read: %d instead of 54 bytes", read_bytes);
        return -1;
    }


    // extract image data from header
    *width = *(int*)&info[18];
    *height = *(int*)&info[22];
    *size = *(int*)&info[2];
    *offset = *(int*)&info[10];
    *row_padded = (*width*3 + 3) & (~3);


    //printf ("Filename: %s, Width: %d, Row_padded: %d, Height: %d, Size:  %d, Offset: %d\n", filename, *width, *row_padded, *height, *size, *offset);
    w = *width;
    p = *row_padded;
    h = *height;
    o = *offset;

    data = (unsigned char*) malloc (p * h);
    *dataR = (unsigned char*) malloc (w * h);
    *dataG = (unsigned char*) malloc (w * h);
    *dataB = (unsigned char*) malloc (w * h);

    fseek(f, sizeof(unsigned char) * o, SEEK_SET);
    read_bytes = fread(data, sizeof(unsigned char), p * h, f);
    if (read_bytes != p * h)
    {
        printf ("Error at read: %d\n", read_bytes);
        free (data);
        return -1;
    }
    for (k = 0; k < h; k++)
    {
        i = k * p;
        for (j = 0; j < w; j++)
        {
            (*dataB)[k*w + j] = data[i];
            (*dataG)[k*w + j] = data[i + 1];
            (*dataR)[k*w + j] = data[i + 2];

            //printf ("BGR %d %d i= %d: %d %d %d\n", k, j, i, data[i], data[i+1], data[i+2]);
            i+= 3;
        }
    }

    free (data);
    fclose(f);
    return 0;
}

int write_BMP(char* filename, float *dataB, float *dataG, float *dataR, unsigned char *header, int offset, int width,  int row_padded, int height)
{
    int write_bytes = 0, i, pad_size;
    FILE* f = fopen(filename, "wb");
    unsigned char null_byte = 0, valR, valB, valG;

    write_bytes = fwrite (header, sizeof(unsigned char), offset, f);
    if (write_bytes < offset)
    {
        printf( "Error at writing the header\n");
        return -1;
    }


    for (i = 0; i< width*height; i++)
    {
        if ( dataB[i] > 256.0f || dataR[i] > 256.0f || dataG[i] > 256.0f ){
            printf( "Error: invalid value %f %f %f", dataB[i], dataG[i], dataR[i]);
            return -1;
        }

        valB = dataB[i];
        valG = dataG[i];
        valR = dataR[i];
        write_bytes = fwrite(&valB, sizeof( unsigned char ), 1, f );
        if (write_bytes != 1)
        {
            printf ("Error at write: i = %d %d\n", i, valB);
            return -1;
        }
        write_bytes = fwrite(&valG, sizeof( unsigned char ), 1, f );
        if (write_bytes != 1)
        {
            printf ("Error at write: i = %d %d\n", i, valG);
            return -1;
        }
        write_bytes = fwrite(&valR, sizeof( unsigned char ), 1, f );
        if (write_bytes != 1)
        {
            printf ("Error at write: i = %d %d\n", i, valR);
            return -1;
        }

        if ((i + 1) % width == 0 ) {
            pad_size = row_padded - width *3 ;
            while( pad_size-- > 0 ) {
                fwrite(&null_byte, sizeof( unsigned char), 1, f );
            }
        }
    }

    fclose (f);
    return 0;
}



float convolve(const float *kernel, const float *buffer, const int ksize) {
    float sum = 0.0f;
    int i;
    for(i=0; i<ksize; i++)
    {
        sum += kernel[i]*buffer[i];
    }
    return sum;
}

void *row_blur_helper(void *ptr)
{
    blur_helper_arg *arg = (blur_helper_arg *) ptr;
    unsigned char *src = arg->src;
    float *dst = arg->dst, *kernel = arg->kernel;
    int width = arg->width, ksize = arg->ksize;
    int start = arg->start, stop = arg->stop;

    int x, y, i;
    int halfksize = ksize / 2;
    float *buffer = (float*)malloc(ksize * sizeof(float));

    for (y = start; y < stop; y++)
    {
        for (x = 0; x < halfksize  ; x++)
        {
            buffer[x] = (float) src[y*width];
        }

        for (x = halfksize; x < ksize-1; x++)
        {
            buffer[x] = (float) src[y*width + x-halfksize];
        }

        for (x = 0; x < width; x++)
        {

            i = (x+ksize-1)%ksize;

            if(x < width - halfksize)
            {
                buffer[i] = (float) src[ y * width + x + halfksize ];
            }
            else
            {
                buffer[i] = (float) src[ y * width + width - 1 ];
            }

            dst[y*width + x] = convolve(kernel, buffer, ksize);
        }
    }

    free(buffer);
}

void *column_blur_helper(void *ptr)
{
    blur_helper_arg *arg = (blur_helper_arg *) ptr;
    float *dst = arg->dst, *kernel = arg->kernel;
    int width = arg->width, height = arg->height, ksize = arg->ksize;
    int start = arg->start, stop = arg->stop;

    int x, y, i;
    int halfksize = ksize / 2;
    float *buffer = (float*)malloc(ksize * sizeof(float));

    for (x = start; x < stop; x++)
    {
        for (y = 0; y < halfksize  ; y++)
        {
            buffer[y] = dst[0*width + x];
        }
        for (y = halfksize; y < ksize-1; y++)
        {
            buffer[y] = dst[(y-halfksize)*width + x];
        }

        for (y = 0; y < height; y++)
        {
            i = (y+ksize-1)%ksize;
            if(y < height - halfksize)
            {
                buffer[i] = dst[(y+halfksize)*width + x];
            }
            else
            {
                buffer[i] = dst[(height-1)*width + x];
            }

            dst[y*width + x] = convolve(kernel, buffer, ksize);
        }
    }

    free(buffer);
}

void distribute_array(short n_processors, int array_size, int *l_bound, int *u_bound)
{
    int i, interval;

    for (i = 0; i < array_size % n_processors; i++)
    {
        interval = (int) ceilf(array_size / n_processors);

        l_bound[i] = (i) ? u_bound[i - 1] : 0;
        u_bound[i] = (i + 1) * interval;
    }

    for (i = array_size % n_processors; i < n_processors; i++)
    {
        interval = (int) floor(array_size / n_processors);

        l_bound[i] = u_bound[i - 1];
        u_bound[i] = (i + 1) * interval;
    }
}

void *gaussian_blur(void *ptr)
{
    gaussian_blur_arg *arg = (gaussian_blur_arg *) ptr;

    unsigned char *src = arg->data;
    float *dst = arg->dst;
    int width = arg->width;
    int height = arg->height;
    float sigma = arg->sigma;
    int ksize = arg->blur_size;
    short n_processors = (short) arg->n;

    int x, y, i;

    int halfksize = ksize / 2;
    float sum = 0.f, t;
    float *kernel;

    // create Gaussian kernel
    kernel = (float*)malloc(ksize * sizeof(float));

    if (!kernel)
    {
        printf ("Error in memory allocation!\n");
        return NULL;
    }

    // if sigma too small, just copy src to dst
    if (ksize <= 1)
    {
        for (y = 0; y < height; y++)
            for (x = 0; x < width; x++)
                dst[y*width + x] = src[y*width + x];
        return NULL;
    }


    //compute the Gaussian kernel values
    for (i = 0; i < ksize; i++)
    {
        x = i - halfksize;
        t = expf(- x * x/ (2.0f * sigma * sigma)) / (sqrt(2.0f * M_PI) * sigma);
        kernel[i] = t;
        sum += t;
    }
    for (i = 0; i < ksize; i++)
    {
        kernel[i] /= sum;
        //printf ("Kernel [%d] = %f\n", i, kernel[i]);
    }

    int l_bound[n_processors], u_bound[n_processors], ret_code = 0;

    pthread_t threads[n_processors];
    blur_helper_arg args[n_processors];

    distribute_array(n_processors, height, l_bound, u_bound);

    // blur each row
    for (i = 0; i < n_processors; i++)
    {
        args[i].src = src;
        args[i].dst = dst;
        args[i].kernel = kernel;
        args[i].width = width;
        args[i].height = height;
        args[i].ksize = ksize;
        args[i].start = l_bound[i];
        args[i].stop = u_bound[i];

        ret_code = pthread_create(&threads[i], NULL, row_blur_helper, (void *) &args[i]);
        if (ret_code < 0)
        {
            free(kernel);

            exit(EXIT_FAILURE);
        }
    }

    for (i = 0; i < n_processors; i++)
    {
        pthread_join(threads[i], NULL);
    }

    distribute_array(n_processors, width, l_bound, u_bound);

    // blur each column
    for (i = 0; i < n_processors; i++)
    {
        args[i].start = l_bound[i];
        args[i].stop = u_bound[i];

        ret_code = pthread_create(&threads[i], NULL, column_blur_helper, (void *) &args[i]);
        if (ret_code < 0)
        {
            free(kernel);

            exit(EXIT_FAILURE);
        }
    }

    for (i = 0; i < n_processors; i++)
    {
        pthread_join(threads[i], NULL);
    }

    // clean up
    free(kernel);
}

void safe_exit_failure(gaussian_blur_arg *argB, gaussian_blur_arg *argR, gaussian_blur_arg *argG)
{
    free (argB->dst);
    free (argR->dst);
    free (argG->dst);
    free (argB->data);
    free (argR->data);
    free (argG->data);

    exit(EXIT_FAILURE);
}

int main(int argc, char ** argv)
{
    unsigned char info[54];
    int blur_size, ret_code = 0, size, width, height, offset, row_padded, n;
    char *in_filename, *out_filename;
    float sigma;

    pthread_t threadB, threadR, threadG;

    gaussian_blur_arg argB, argR, argG;
    argB.data = argG.data = argR.data = NULL;

    if (argc != 6)
    {
        printf ("Usage: %s <filename.bmp> <sigma> <blur_size> <output_filename.bmp>  <num_threads>", argv[0]);
        return -1;
    }
    in_filename = argv[1];
    out_filename = argv[4];
    n = atoi(argv[3]);
    blur_size = atoi (argv[3]);
    sigma = atof (argv[2]);
    ret_code = read_BMP(in_filename, info, &argR.data, &argG.data, &argB.data, &size, &width, &height, &offset, &row_padded);
    if (ret_code < 0)
    {
        safe_exit_failure(&argB, &argR, &argG);
    }

    argB.dst = (float*)malloc (width*height* sizeof(float));
    argR.dst = (float*)malloc (width*height* sizeof(float));
    argG.dst = (float*)malloc (width*height* sizeof(float));

    argB.width = argR.width = argG.width = width;
    argB.height = argR.height = argG.height = height;
    argB.sigma = argR.sigma = argG.sigma = sigma;
    argB.blur_size = argR.blur_size = argG.blur_size = blur_size;
    argB.n = argR.n = argG.n = n;

    ret_code = pthread_create(&threadB, NULL, gaussian_blur, (void *) &argB);
    if(ret_code)
    {
        fprintf(stderr,"Error while creating threadB; return code: %d\n", ret_code);

        safe_exit_failure(&argB, &argR, &argG);
    }

    ret_code = pthread_create(&threadR, NULL, gaussian_blur, (void *) &argR);
    if(ret_code)
    {
        fprintf(stderr,"Error while creating threadR; return code: %d\n", ret_code);

        safe_exit_failure(&argB, &argR, &argG);
    }

    ret_code = pthread_create(&threadG, NULL, gaussian_blur, (void *) &argG);
    if(ret_code)
    {
        fprintf(stderr,"Error while creating threadG; return code: %d\n", ret_code);

        safe_exit_failure(&argB, &argR, &argG);
    }

    pthread_join(threadB, NULL);
    pthread_join(threadR, NULL);
    pthread_join(threadG, NULL);


    ret_code = write_BMP (out_filename, argB.dst, argG.dst, argR.dst, info, offset, width, row_padded, height);

    free (argB.dst);
    free (argR.dst);
    free (argG.dst);
    free (argB.data);
    free (argR.data);
    free (argG.data);

    return ret_code;
}