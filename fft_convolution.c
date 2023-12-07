/**
Author: Kimiya Saadat
This program will do convolution in frequency domain.

*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define SWAP(a,b)  tempr=(a);(a)=(b);(b)=tempr

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
    char subChunk2ID[4];
    int subChunk2Size;
} WavHeader;

// function declaration
void pad_zeros_to(double *arr, int current_length, int M);
void writeWAVHeader(int numChannels, int numSamples, int outputRate, int bitsPerSample, FILE *outputFile);
size_t fwriteShortLSB(short int data, FILE *stream);
size_t fwriteIntLSB(int data, FILE *stream);
void four1(double data[], int nn, int isign);
void reverse_array(double *arr, int length);
int next_power_of_2(int n);
void convolution(double *x, int K, double *h, double *y) ;
void readWav(double **data, int *length, WavHeader *header, FILE *file);
int convolveTone(char* sampleTone, char* impulseTone, char* outputTone);

// Function to pad zeros to the input array to make its length M
void pad_zeros_to(double *arr, int current_length, int M) {
    int padding = M - current_length;
    for (int i = 0; i < padding; ++i) {
        arr[current_length + i] = 0.0;
    }
}




void writeWAVHeader(int numChannels, int numSamples, int outputRate, int bitsPerSample, FILE *outputFile){

    /*  Calculate the total number of bytes for the data chunk  */
    int BYTES_PER_SAMPLE = bitsPerSample/8;
    int dataChunkSize = numChannels * numSamples * BYTES_PER_SAMPLE;
	
    int formSize = 36 + dataChunkSize;
    short int frameSize = numChannels * BYTES_PER_SAMPLE;
    int bytesPerSecond = (int)ceil(outputRate * frameSize);

  
    fputs("RIFF", outputFile);
    fwriteIntLSB(formSize, outputFile);
    fputs("WAVE", outputFile);

    fputs("fmt ", outputFile);  
    fwriteIntLSB(16, outputFile);
    fwriteShortLSB(1, outputFile);
    fwriteShortLSB((short)numChannels, outputFile);
    fwriteIntLSB(outputRate, outputFile);
    fwriteIntLSB(bytesPerSecond, outputFile);
    fwriteShortLSB(frameSize, outputFile);
    fwriteShortLSB(bitsPerSample, outputFile);

    fputs("data", outputFile);
    fwriteIntLSB(dataChunkSize, outputFile);
}

size_t fwriteIntLSB(int data, FILE *stream)
{
    unsigned char array[4];

    array[3] = (unsigned char)((data >> 24) & 0xFF);
    array[2] = (unsigned char)((data >> 16) & 0xFF);
    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 4, stream);
}

size_t fwriteShortLSB(short int data, FILE *stream)
{
    unsigned char array[2];

    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 2, stream);
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
        exit(-1);
    }

    inputFile = argv[1];
    IRfile = argv[2];
    outputFile = argv[3];

    convolveTone(inputFile, IRfile, outputFile);

    
    return 0;
}


void readWav(double **data, int *length, WavHeader *header, FILE *file) {
    // Read the header

    fread(header, sizeof(WavHeader), 1, file);

    // Ensure the WAV file is of the expected format
    if (header->audioFormat != 1 || header->numChannels != 1 || header->bitsPerSample != 16) {
        fprintf(stderr, "Unsupported WAV format\n");
        fclose(file);
        exit(-1);
    }

    // Calculate the number of samples
    *length = header->subChunk2Size / (header->numChannels * (header->bitsPerSample / 8));

    // Allocate memory for the audio data
    *data = (double *)malloc(*length * sizeof(double));
    if (*data == NULL) {
        perror("Error allocating memory for audio data");
        fclose(file);
        exit(-1);
    }

    // Read the audio data and convert to floating point
    for (int i = 0; i < *length; ++i) {
        short shortSample;  // random initialization
        size_t readFile = fread(&shortSample, sizeof(short), 1, file);
        
        if (readFile< 1) {
            perror("Error reading audio data");
            free(*data);
            fclose(file);
            exit(-1);
        }
        (*data)[i] = (double)shortSample / 32767.0;
    }
    // for (int i = 0; i < *length; i++) {
    //     // Normalize 16-bit PCM data to [-1.0, 1.0]
    //     (*data)[i] = 32767.0;
    // }

    fclose(file);
}



int convolveTone(char* sampleTone, char* impulseTone, char* outputTone) {
    FILE *inputFile = fopen(sampleTone, "rb");
    FILE *IRfile = fopen(impulseTone, "rb");
    //FILE *outputFile = fopen(outputTone, "wb");

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
    

    
    double *inputSignal, *IRsignal;
    int inputLength, IRlength;

    WavHeader inputHeader, IRheader;
    // fread(&inputHeader, sizeof(WavHeader), 1, inputFile);
    // fread(&IRheader, sizeof(WavHeader), 1, IRfile);

    printf("Reading input and IR files...\n");
    readWav(&inputSignal, &inputLength, &inputHeader, inputFile);
    printf("input length: %d\n", inputLength);
    printf("IR length: %d\n", IRlength);
    readWav(&IRsignal, &IRlength, &IRheader, IRfile);
    printf("IR length: %d\n", IRlength);

    // size needed for fft
    int M = inputLength + IRlength - 1; 
    int K = next_power_of_2(M);
    // int K = next_power_of_2(2*M);
    double *x = (double *)calloc(K*2, sizeof(double));
    double *h = (double *)calloc(K*2, sizeof(double));

    for(int i=0; i< inputLength; i++)
    {
        x[2*i] = inputSignal[i]; // real part
        x[2*i + 1] = 0.0; // initially imaginary is 0
    }
    for(int i=0; i< IRlength; i++)
    {
        h[2*i] = IRsignal[i]; // real part
        h[2*i + 1] = 0.0; // initially imaginary is 0
    }

    // pad zeros
    pad_zeros_to(x, 2*inputLength, K*2);
    pad_zeros_to(h, 2*IRlength, K*2);

    // FFT
    four1(x-1, K, 1);
    four1(h-1, K, 1);

    // Convolution
    double *y = (double *)calloc(K*2, sizeof(double));
    convolution(x, K, h, y);

    // Inverse FFT
    four1(y-1, K, -1);

    for(int i=0; i < K*2; i += 2)
    {
        y[i] /= (double)K; // only real is reqd since imaginary is 0
    }

    for (int i = 0; i < 2*K; i=i+2) {
        printf("%f ", y[i]);
    }

    int numChannels = inputHeader.numChannels;
    int bitsPerSample = inputHeader.bitsPerSample;
    int outputRate = inputHeader.sampleRate;

    FILE *outputFile = fopen(outputTone, "wb");
    
    if (outputFile == NULL){
        printf("Couldn't load ouput file, Error!!\n");
        perror("");
        fclose(inputFile);
        free(inputSignal);
        free(IRsignal);
        exit(-1);
    }

    writeWAVHeader(numChannels, M, outputRate, bitsPerSample , outputFile);

    for (int i = 0; i < K; i++) {
        short int sample = (short int)(y[2 * i] * 32767.0); // Convert double to 16-bit sample
        fwriteShortLSB(sample, outputFile);
    }

    // Clean up
    free(inputSignal);
    free(IRsignal);
    free(x);
    free(h);
    free(y);

    fclose(inputFile);
    fclose(IRfile);
    fclose(outputFile);
    

    return 0;
}
        
    

//  now extract real part into imaginary array, imaginary part will be 0. calculate properly

