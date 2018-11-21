#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <fcntl.h>
#include <malloc.h>
#include <string.h>
#include <unistd.h>

#define pi 3.142857
#define C 8

const int m = 8, n = 8;

// Function to perform Discrete Cosine Transform
double** dctTransform(int **matrix)
{
    int i, j, k, l;
  
    // dct will store the discrete cosine transform
    double **dct = (double **)malloc(m * sizeof(double *));
    for (i=0; i<8; i++)
         dct[i] = (double *)malloc(n * sizeof(double));
 
    double ci, cj, dct1, sum;
 
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
 
            // ci and cj depends on frequency as well as
            // number of row and columns of specified matrix
            if (i == 0)
                ci = 1 / sqrt(m);
            else
                ci = sqrt(2) / sqrt(m);
            if (j == 0)
                cj = 1 / sqrt(n);
            else
                cj = sqrt(2) / sqrt(n);
 
            // sum will temporarily store the sum of 
            // cosine signals
            sum = 0;
            for (k = 0; k < m; k++) {
                for (l = 0; l < n; l++) {
                    dct1 = matrix[k][l] * 
                           cos((2 * k + 1) * i * pi / (2 * m)) * 
                           cos((2 * l + 1) * j * pi / (2 * n));
                    sum = sum + dct1;
                }
            }
            dct[i][j] = ci * cj * sum;
        }
    }
    
    return dct;
}

// Function to perform Inverse Discrete Cosine Transform
double** inverseDct(double** in)
{
  int i, j, u, v;
  double s;

  double **idct = (double **)malloc(8 * sizeof(double *));
    for (i=0; i<8; i++)
         idct[i] = (double *)malloc(8 * sizeof(double));

  for (i = 0; i < 8; i++)
    for (j = 0; j < 8; j++)
    {
      s = 0;

      for (u = 0; u < 8; u++)
        for (v = 0; v < 8; v++)
          s += in[u][v] * cos((2 * i + 1) * u * M_PI / 16) *
                          cos((2 * j + 1) * v * M_PI / 16) *
               ((u == 0) ? 1 / sqrt(2) : 1.) *
               ((v == 0) ? 1 / sqrt(2) : 1.);

      idct[i][j] = s / 4;
    }
    return idct;
}

// Function made for debugging purpose
void print_matrix(int **matrix ,int x, int y)
{
  for(int i=0;i<x;i++)
  {
    for(int j=0;j<y;j++)
    {
      printf("%d\t" ,matrix[i][j]);
    }
    printf("\n");
  }
  printf("\n\n\n");
}

// Function made for debugging purpose
void print_matrix_float(double **matrix ,int x, int y)
{
  for(int i=0;i<x;i++)
  {
    for(int j=0;j<y;j++)
    {
      printf("%f\t" ,matrix[i][j]);
    }
    printf("\n");
  }
  printf("\n\n\n");
}

// Function to store image in ZigZag Order
int* getZigZagMatrix( int **arr, int n, int m)
{
    int row = 0, col = 0;
    
    // Boolean variable that will true if we
    // need to increment 'row' value otherwise
    // false- if increment 'col' value
    int row_inc = 0;
    int array_size=m*n;
    int* w;
    w = (int *)malloc(array_size * sizeof(int));
    int array_counter=0;
 
    // Print matrix of lower half zig-zag pattern
    int mn = (m<n)?m:n;
    for (int len = 1; len <= mn; ++len) {
        for (int i = 0; i < len; ++i) {
            w[array_counter]=arr[row][col];
            array_counter++;
 
            if (i + 1 == len)
                break;
            // If row_increment value is true
            // increment row and decrement col
            // else decrement row and increment
            // col
            if (row_inc)
                ++row, --col;
            else
                --row, ++col;
        }
 
        if (len == mn)
            break;
 
        // Update row or col vlaue according
        // to the last increment
        if (row_inc)
            ++row, row_inc = 0;
        else
            ++col, row_inc = 1;
    }
 
    // Update the indexes of row and col variable
    if (row == 0) {
        if (col == m - 1)
            ++row;
        else
            ++col;
        row_inc = 1;
    }
    else {
        if (row == n - 1)
            ++col;
        else
            ++row;
        row_inc = 0;
    }
 
    // Print the next half zig-zag pattern
    int MAX = ((m>n)?m:n) - 1;
    for (int len, diag = MAX; diag > 0; --diag) {
 
        if (diag > mn)
            len = mn;
        else
            len = diag;
 
        for (int i = 0; i < len; ++i) {
            w[array_counter]=arr[row][col];
            array_counter++;
 
            if (i + 1 == len)
                break;
 
            // Update row or col value according
            // to the last increment
            if (row_inc)
                ++row, --col;
            else
                ++col, --row;
        }
 
        // Update the indexes of row and col variable
        if (row == 0 || col == m - 1) {
            if (col == m - 1)
                ++row;
            else
                ++col;
 
            row_inc = 1;
        }
 
        else if (col == 0 || row == n - 1) {
            if (row == n - 1)
                ++col;
            else
                ++row;
 
            row_inc = 0;
        }
    }

    return w;
}

// Function to get ZigZag Order indices(pattern)
int* getZigZagOrder(int n, int m)
{
    int row = 0, col = 0;

    // Boolean variable that will true if we
    // need to increment 'row' value otherwise
    // false- if increment 'col' value
    int row_inc = 0;
    int array_size=m*n;
    int* order;
    order = (int *)malloc(array_size * sizeof(int));
    int array_counter=0;
 
    // Print matrix of lower half zig-zag pattern
    int mn = (m<n)?m:n;
    for (int len = 1; len <= mn; ++len) {
        for (int i = 0; i < len; ++i) {
            order[array_counter]=row*m+col;
            array_counter++;
 
            if (i + 1 == len)
                break;
            // If row_increment value is true
            // increment row and decrement col
            // else decrement row and increment
            // col
            if (row_inc)
                ++row, --col;
            else
                --row, ++col;
        }
 
        if (len == mn)
            break;
 
        // Update row or col vlaue according
        // to the last increment
        if (row_inc)
            ++row, row_inc = 0;
        else
            ++col, row_inc = 1;
    }
 
    // Update the indexes of row and col variable
    if (row == 0) {
        if (col == m - 1)
            ++row;
        else
            ++col;
        row_inc = 1;
    }
    else {
        if (row == n - 1)
            ++col;
        else
            ++row;
        row_inc = 0;
    }
 
    // Print the next half zig-zag pattern
    int MAX = ((m>n)?m:n) - 1;
    for (int len, diag = MAX; diag > 0; --diag) {
 
        if (diag > mn)
            len = mn;
        else
            len = diag;
 
        for (int i = 0; i < len; ++i) {
            order[array_counter]=row*m+col;
            array_counter++;
 
            if (i + 1 == len)
                break;
 
            // Update row or col vlaue according
            // to the last increment
            if (row_inc)
                ++row, --col;
            else
                ++col, --row;
        }
 
        // Update the indexes of row and col variable
        if (row == 0 || col == m - 1) {
            if (col == m - 1)
                ++row;
            else
                ++col;
 
            row_inc = 1;
        }
 
        else if (col == 0 || row == n - 1) {
            if (row == n - 1)
                ++col;
            else
                ++row;
 
            row_inc = 0;
        }
    }

    return order;
}

// Main function
void main()
{
    int fp1,fp2,i,j,k,l;
    FILE* fp;
    double a;
    char input_image[100],output_image[100];
    unsigned char image[1000],image1[1000];
    int npix,nscan,pix,scan,i1,i2;

    // Reading input image
    printf("Give input image name \n");
    scanf("%s",input_image);
    printf("Give size of image \n");
    scanf("%d%d",&npix,&nscan);
    
    // opening input image
    fp1= open(input_image,O_RDONLY);
    if (fp1<0)
    {
        printf("Error in opening %s image \n",input_image);
        exit(1);
    }

    int **img = (int **)malloc(nscan * sizeof(int *));
    for (i=0; i<nscan; i++)
         img[i] = (int *)malloc(npix * sizeof(int));
    
    double **dct = (double **)malloc(8 * sizeof(double *));
    for (i=0; i<8; i++)
         dct[i] = (double *)malloc(8 * sizeof(double));

    int **quantized_value = (int **)malloc(8 * sizeof(int *));
    for (i=0; i<8; i++)
         quantized_value[i] = (int *)malloc(8 * sizeof(int));
   
    // Quantization Matrix 
    int q[8][8] = { { 16, 11, 10, 16, 24, 40, 51, 61 },
        						{ 12, 12, 14, 19, 26, 58, 60, 55 },
        						{ 14, 13, 16, 24, 40, 57, 69, 56 },
        						{ 14, 17, 22, 29, 51, 87, 80, 62},
        						{ 18, 22, 37, 56, 68, 109, 103, 77 },
        						{ 24, 35, 55, 64, 81, 104, 113, 92},
        						{ 49, 64, 78, 87, 103, 121, 120, 101 },
        						{ 72, 92, 95, 98, 112, 100, 103, 99 } };

    for (scan=0;scan<nscan;scan++)
    {
        read(fp1,&image[0],npix*sizeof(unsigned char));
        for (pix=0;pix<npix;pix++)
        {
            a=(double)image[pix];
            img[scan][pix] = a-128;
        }
    }

    // file to store compressed image (nonzero coefficients only)
    fp = fopen ("compressed_image", "w+");

    double **d = (double **)malloc(8 * sizeof(double *));
    for (i=0; i<8; i++)
         d[i] = (double *)malloc(8 * sizeof(double));

    int **im = (int **)malloc(8 * sizeof(int *));
    for (i=0; i<8; i++)
         im[i] = (int *)malloc(8 * sizeof(int));

    printf("Compressing and reconstructing.....\n");

    for(scan=0;scan<nscan;scan+=8)
    {
        for(pix=0;pix<npix;pix+=8)
        {
            for(i=scan;i<scan+8;i++)
            {
              for(j=pix;j<pix+8;j++)
              {
                im[i-scan][j-pix]=img[i][j];
              }
            }

            // finding DCT of 8*8 images
            d=dctTransform(im);
            for(i=scan;i<scan+8;i++)
            {
              for(j=pix;j<pix+8;j++)
              {
                dct[i-scan][j-pix]=d[i-scan][j-pix];
                quantized_value[i-scan][j-pix]=dct[i-scan][j-pix]/q[i-scan][j-pix];
              }
            }

            int* w;
            int write_elements=m*n;
            w = (int *)malloc((m*n) * sizeof(int));
            w = getZigZagMatrix(quantized_value,8,8);

            for(i=m*n-1;i>=0;i--)
            {
              if(w[i]==0) write_elements--;
              else break;
            }
            
            fprintf(fp, "%d ",write_elements);
            for(i=0;i<write_elements;i++)
            {
              fprintf(fp, "%d ",w[i]);
            }          
        }
    }
    fclose(fp);
    close(fp1);

    //-------------JPEG Image Reconstruction------------------//

    int* order;
    int write_elements=m*n;
    order = (int *)malloc((m*n) * sizeof(int));
    order = getZigZagOrder(8,8);

    unsigned char img_rc[npix][nscan];

    fp=fopen("compressed_image","r");
    if(fp<0)
    {
      printf("unable to open file\n");
      exit(1);
    }

    int buffer[64];
    memset(buffer,0,64*sizeof(int));
    int nzTotal;
    
    double **block = (double **)malloc(8 * sizeof(double *));
    for (i=0; i<8; i++)
         block[i] = (double *)malloc(8 * sizeof(double));

    double **idct = (double **)malloc(8 * sizeof(double *));
    for (i=0; i<8; i++)
         idct[i] = (double *)malloc(8 * sizeof(double));
    
    // Logic for reconstructiong image from compressed image
    for(scan=0;scan<nscan;scan+=8)
    {
        for(pix=0;pix<npix;pix+=8)
        {
          memset(buffer,0,64*sizeof(int));
          fscanf(fp,"%d",&nzTotal);
          for(i=0;i<nzTotal;i++)
          {
            fscanf(fp,"%d",&buffer[i]); 
          }

          for (i = 0; i < 64; ++i)
          {
            block[order[i]/8][order[i]%8]=buffer[i];
            block[order[i]/8][order[i]%8]*=q[order[i]/8][order[i]%8];
          }

          idct=inverseDct(block);
          for(i=scan;i<scan+8;i++)
          {
            for(j=pix;j<pix+8;j++)
            {
              img_rc[i][j]=(unsigned char)((int)(idct[i-scan][j-pix]+128));
            }
          }
        }
    }

    // reconstructing output image
    sprintf(output_image,"%s_jpeg_reconstructed",input_image);
    fp2= creat(output_image,0667);
    if (fp2<0)
    {
        printf("Error in opening %s image \n",output_image);
        exit(1);
    }

    write(fp2,&img_rc,(npix)*(nscan)*sizeof(unsigned char));
    printf("check Reconstructed image: %s\n",output_image);
    close(fp2);
}