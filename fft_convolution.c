/**
Author: Kimiya Saadat
This program will do convolution in frequency domain.

*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define SWAP(a,b)  tempr=(a);(a)=(b);(b)=tempr

// Function to pad zeros to the input array to make its length M
void pad_zeros_to(double *arr, int current_length, int M) {
    int padding = M - current_length;
    for (int i = 0; i < padding; ++i) {
        arr[current_length + i] = 0.0;
    }
}


// struct to hold all data up until the end of subchunk1
typedef struct {
    char chunk_ID[4];           // "RIFF"
    int chunk_Size;             // Size of the entire file in bytes minus 8 bytes
    char format[4];            // "WAVE"
    char subchunk1_ID[4];       // "fmt "
    int subchunk1_Size;         // PCM = 16
    short audioFormat;         // PCM = 1
    short numChannels;         // Mono = 1, Stereo = 2, etc.
    int sampleRate;            // 8000, 44100, etc.
    int byteRate;              // == SampleRate * NumChannels * BitsPerSample/8
    short blockAlign;          // == NumChannels * BitsPerSample/8
    short bitsPerSample;       // 8 bits = 8, 16 bits = 16, etc.
} WavHeader;

void writeWAVHeader(int numChannels, int numSamples, int outputRate, int bitsPerSample, FILE *outputFile){

    //calculations for header fields

    //subchunk sizes
    int BytesPerSample = bitsPerSample/8;
    int dataChunkSize = numChannels * numSamples * BytesPerSample; 
    int formSize = 36 + dataChunkSize; //36 is sum of all field sizes before data
    
    //number of bytes for one sample after accounting all channels (frame size)
    short frameSize = numChannels * BytesPerSample;
    
    //bytes per second  = sample rate * total bytes per sample
    int byteRate = (int) outputRate * frameSize;

    //write header to file
    //use fputs for big endian fields, use respective LSB methods for little endian

    //RIFF chunk descriptor
    fputs("RIFF", outputFile);
    fwriteIntLSB(formSize, outputFile);
    fputs("WAVE", outputFile);
    

    fputs("fmt ", outputFile);
    fwriteIntLSB(16, outputFile); 
    fwriteShortLSB(1, outputFile);
    fwriteShortLSB((short)numChannels, outputFile);
    fwriteIntLSB(outputRate, outputFile);
    fwriteIntLSB(byteRate, outputFile);
    fwriteShortLSB(frameSize, outputFile);
    fwriteShortLSB(bitsPerSample, outputFile);

    fputs("data", outputFile);
    fwriteIntLSB(dataChunkSize, outputFile);

}




//  The four1 FFT from Numerical Recipes in C,
//  p. 507 - 508.
//  Note:  changed float data types to double.
//  nn must be a power of 2, and use +1 for
//  isign for an FFT, and -1 for the Inverse FFT.
//  The data is complex, so the array size must be
//  nn*2. This code assumes the array starts
//  at index 1, not 0, so subtract 1 when
//  calling the routine (see main() below).
void four1(double data[], int nn, int isign)
{
    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;

    n = nn << 1;
    j = 1;

    for (i = 1; i < n; i += 2) {
	if (j > i) {
	    SWAP(data[j], data[i]);
	    SWAP(data[j+1], data[i+1]);
	}
	m = nn;
	while (m >= 2 && j > m) {
	    j -= m;
	    m >>= 1;
	}
	j += m;
    }

    mmax = 2;
    while (n > mmax) {
	istep = mmax << 1;
	theta = isign * (6.28318530717959 / mmax);
	wtemp = sin(0.5 * theta);
	wpr = -2.0 * wtemp * wtemp;
	wpi = sin(theta);
	wr = 1.0;
	wi = 0.0;
	for (m = 1; m < mmax; m += 2) {
	    for (i = m; i <= n; i += istep) {
		j = i + mmax;
		tempr = wr * data[j] - wi * data[j+1];
		tempi = wr * data[j+1] + wi * data[j];
		data[j] = data[i] - tempr;
		data[j+1] = data[i+1] - tempi;
		data[i] += tempr;
		data[i+1] += tempi;
	    }
	    wr = (wtemp = wr) * wpr - wi * wpi + wr;
	    wi = wi * wpr + wtemp * wpi + wi;
	}
	mmax = istep;
    }
    
}

// Function to reverse an array
void reverse_array(double *arr, int length) {
    int i, j;
    double temp;

    for (i = 0, j = length - 1; i < j; ++i, --j) {
        temp = arr[i];
        arr[i] = arr[j];
        arr[j] = temp;
    }
}

// Function to find the next power of 2 greater than or equal to n
int next_power_of_2(int n) {
    return pow(2, (int)(log2(n - 1) + 1));
}


void convolution(double *x, int K, double *h, double *y) {


    // Perform the DFT
    for (int k = 0, nn = 0; k < K; k++, nn += 2)
    {
	    y[nn] = ((x[nn] * h[nn]) - (x[nn+1] * h [nn+1]));
	    y[nn+1] = ((x[nn] * h[nn+1]) + (x[nn+1] * h[nn]));
	}
    
    }


int main (int argc, char *argv[]){

    char *inputFile, *IRfile, *outputFile;

    if (argc != 4){
        fprintf(stderr, "Usage : %s \t3 arguments required: (input file) (IR file) (output file)\n", argv[0]);
    }

    inputFile = argv[1];
    IRfile = argv[2];
    outputFile = argv[3];

    readWavFile(inputFile, IRfile, outputFile);

    
    return 0;
}

int readWavFile(char* sampleTone, char* impulseTone, char*outputTone ){
    FILE *inputFile = fopen(sampleTone, "rb");
    FILE *IRfile = fopen(impulseTone, "rb");
    FILE *outputFile = fopen(outputTone, "wb");

    if (inputFile == NULL){
        printf("Couldn't load input file, Error!!\n");
        perror("");
        exit(-1);
    }
    if (IRfile == NULL){
        printf("Couldn't load IR file, Error!!\n");
        perror("");
        exit(-1);
    }
    if (outputFile == NULL){
        printf("Couldn't load ouput file, Error!!\n");
        perror("");
        fclose(inputFile);
        exit(-1);
    }

    WavHeader inputHeader, IRheader;
    // read header subchunk 1
    fread(&inputHeader, sizeof(WavHeader), 1, inputFile);
    fread(&IRheader, sizeof(WavHeader), 1, IRfile);



}
int convolveTone(float x[], int N, float y[], int M, float h[], int P){

    // Example usage
    int inputLength = sizeof(x);
    int IRlength = sizeof(y);
    int M = inputLength + IRlength - 1; // Length of x
    int K = next_power_of_2(2*M);
    double *inputSignal = (double *)calloc(K*2, sizeof(double));
    double *IRsignal = (double *)calloc(K , sizeof(double));

    for(int i=0; i< 2*M; i++)
    {
        x[i] = 0.0;
    }
    // sample array pass instead of x
    x[0] = 1.0;
    x[2] = 2.0;
    x[4] = 3.0;
    x[6] = 4.0;
    x[8] = 5.0;
    //double x[] = {1.0, 0, 2.0, 0, 3.0, 0, 4.0, 0, 5.0, 0};

    int N = 3; // Length of h
    double *h = (double *)calloc(K*2, sizeof(double));
    //double h[] = {0.5, 0, 0.5, 0};

    for(int i=0; i< 2*N; i++)
    {
        h[i] = 0.0;
    }
    // impulse array
    h[0] = 0.5;
    h[1] = 0.0;
    h[2] = 0.5;
    h[3] = 0.0;
    h[4] = 1.0;
    h[5] = 0.0;
    // Calculate the output size
   
    printf("k: %d \n", K);
    pad_zeros_to(x, 2*M, K*2);
    pad_zeros_to(h, 2*N, K*2);


    four1(x-1, K, 1);
    four1(h-1, K, 1);

    double y[K*2];
    // Perform convolution in frequency domain
    convolution(x, K, h, y);
    //  for (int i = 0; i < 2 * K; ++i) {
    //     printf("%f ", y[i]);
    // }
        

printf("\n ");
 //int u = K;
    four1(y-1, K, -1);
    
  

    for(int k = 0, i=0; k<M+N-1; k++, i += 2)
    {
        y[i] /= (double)K;
        y[i+1] /= (double)K;
    }
    printf("\n Result: \n");
    for (int i = 0; i < 2*K; i=i+2) {
        printf("%f ", y[i]);
    }

    //free(y);

    return 0;
}
// extracting real part into imaginary array

