

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
    // char subChunk2ID[4];
    // int subChunk2Size;
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


int main(int argc, char *argv[]) {
    
    char *sampleTone, *impulseTone, *outputTone;
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <input_wave_file> <impulse_response_wave_file> <output_wave_file>\n", argv[0]);
        return -1;
    }

    sampleTone = argv[1];
    impulseTone = argv[2];
    outputTone = argv[3];

    FILE *inputFile = fopen(sampleTone, "rb");
    FILE *IRfile = fopen(impulseTone, "rb");
    printf("inputFile: %s\n", sampleTone);
    printf("IRfile: %s\n", impulseTone);
    //FILE *outputFile = fopen(outputTone, "wb");

    if (inputFile == NULL){
        printf("Couldn't load input file, Error!!\n");
        perror("");
        return -1;
    }
    if (IRfile == NULL){
        printf("Couldn't load IR file, Error!!\n");
        perror("");
        return -1;
    }
    

    WavHeader inputHeader, IRheader;
    
    // Reading WAV header
    
    if (fread(&inputHeader, sizeof(WavHeader), 1, inputFile) < 1 || fread(&IRheader, sizeof(WavHeader), 1, IRfile) < 1) {
        perror("Error reading WAV header");
        fclose(inputFile);
        fclose(IRfile);
        return -1;
    }
    if (inputHeader.subchunk1_Size != 16){
        // eliminate null bytes
        int remainder = inputHeader.subchunk1_Size - 16;
        char randomVar[remainder];
        fread(randomVar, remainder  , 1, inputFile);
    }

    if (IRheader.subchunk1_Size != 16){
        // eliminate null bytes
        int remainder = IRheader.subchunk1_Size - 16;
        char randomVar[remainder];
        fread(randomVar, remainder  , 1, IRfile);
    }

    char subchunk2_id_sample[4];
    char subchunk2_id_impulse[4];
    int subchunk2_size_sample; // an integer is 4 bytes
    int subchunk2_size_impulse; // an integer is 4 bytes
    fread(&subchunk2_id_sample, sizeof(subchunk2_id_sample), 1, inputFile);
    fread(&subchunk2_size_sample, sizeof(subchunk2_size_sample), 1, inputFile);
    fread(&subchunk2_id_impulse, sizeof(subchunk2_id_impulse), 1, IRfile);
    fread(&subchunk2_size_impulse, sizeof(subchunk2_size_impulse), 1, IRfile);

    

    int inputLength = subchunk2_size_sample / (inputHeader.bitsPerSample / 8);
    int IRlength = subchunk2_size_impulse / (IRheader.bitsPerSample / 8);

    // After reading the WAV header
    printf("Chunk Size: %d, Format: %s\n", inputHeader.chunk_Size, inputHeader.format);
    printf("Subchunk1 Size: %d, Audio Format: %d\n", inputHeader.subchunk1_Size, inputHeader.audioFormat);
    printf("Byte Rate: %d, Block Align: %d\n", inputHeader.byteRate, inputHeader.blockAlign);
    printf("Sub Chunk 2 Size: %d\n", subchunk2_size_sample);

    // After reading the WAV header
    printf("IR Chunk Size: %d, Format: %s\n", IRheader.chunk_Size, IRheader.format);
    printf("IR Subchunk1 Size: %d, Audio Format: %d\n", IRheader.subchunk1_Size, IRheader.audioFormat);
    printf("IR Byte Rate: %d, Block Align: %d\n", IRheader.byteRate, IRheader.blockAlign);
    printf("IR Sub Chunk 2 Size: %d\n", subchunk2_size_impulse);



    // Add a check for negative lengths
    if (inputLength < 0 || IRlength < 0) {
        fprintf(stderr, "Error: Negative length found in WAV headers.\n");
        fclose(inputFile);
        fclose(IRfile);
        return -1;
    }

    printf("inputLength: %d\n", inputLength);
    printf("IRlength: %d\n", IRlength);

    // M is length of output signal
    int M = inputLength + IRlength - 1; 
    int K = next_power_of_2(M);   // K is length of input signal

    double *inputSignal = (double *)malloc(2 * K * sizeof(double));
    double *IRsignal = (double *)malloc(2 * K * sizeof(double));

    

    printf("K: %d\n", K);
    printf("M: %d\n", M);
    

    if (!inputSignal || !IRsignal) {
        perror("Error allocating memory for audio data");
        fclose(inputFile);
        fclose(IRfile);
        if (inputSignal) 
            free(inputSignal);
        if (IRsignal) 
            free(IRsignal);
        return -1;
    }


    for (int i = 0; i < inputLength; i++) {
        inputSignal[2 * i] = inputSignal[i]; // real part
        inputSignal[2 * i + 1] = 0.0; // imaginary part
    }
    for (int i = 0; i < IRlength; i++) {
        IRsignal[2 * i] = IRsignal[i]; // real part
        IRsignal[2 * i + 1] = 0.0; // imaginary part
    }

    printf("inputSignal[0]: %f\n", inputSignal[0]);
    printf("IRsignal[0]: %f\n", IRsignal[0]);
    // pad zeros
    short tempSample;
    for (int i = 0; i < inputLength; i++) {
        if (fread(&tempSample, sizeof(short), 1, inputFile) == 1) {
            inputSignal[2 * i] = (double)tempSample/32767.0;
            inputSignal[2 * i + 1] = 0.0; // imaginary part is zero
        } else {
            fprintf(stderr, "Error reading input file\n");
            fclose(inputFile);
            fclose(IRfile);
            free(inputSignal);
            free(IRsignal);
            return -1;
        }
    }
    for (int i = 0; i < IRlength; i++) {
        if (fread(&tempSample, sizeof(short), 1, IRfile) == 1) {
            IRsignal[2 * i] = (double)tempSample/32767.0;
            IRsignal[2 * i + 1] = 0.0; // imaginary part is zero
        } else {
            fprintf(stderr, "Error reading IR file\n");
            fclose(inputFile);
            fclose(IRfile);
            free(inputSignal);
            free(IRsignal);
            return -1;
        }
    }

    pad_zeros_to(inputSignal, 2*inputLength, 2*K);
    pad_zeros_to(IRsignal, 2*IRlength, 2*K);
    
    
    

    

    printf("\nFFT...\n");
    // FFT
    four1(inputSignal-1, K, 1);
    four1(IRsignal-1, K, 1);

    // Convolution
    double *outputSignal = (double *)malloc(2*K * sizeof(double));
    for (int i = 0; i < 2 * M; i++) {
        outputSignal[i] = 0.0;
    }


    printf("\nConvolving...\n");
    convolution(inputSignal, M, IRsignal, outputSignal);
 
    
    printf("inputHeader.numChannels: %d\n", inputHeader.numChannels);
    // Inverse FFT
    four1(outputSignal-1, K, -1);
  

    for(int i=0; i < M; i++)
    {
        outputSignal[2*i] /= (double)K; 
        outputSignal[2*i+1] /= (double)K;
    }


    int numChannels = inputHeader.numChannels;
    int bitsPerSample = inputHeader.bitsPerSample;
    int outputRate = inputHeader.sampleRate;

    FILE *outputFile = fopen(outputTone, "wb");
    
    if (outputFile == NULL){
        printf("Couldn't load ouput file, Error!!\n");
        perror("");
        fclose(inputFile);
        fclose(IRfile);
        free(inputSignal);
        free(IRsignal);
        free(outputSignal);
        exit(-1);
    }

    writeWAVHeader(numChannels, M, outputRate, bitsPerSample , outputFile);

    for (int i = 0; i < M; i++) {
        short int sample = (short int)(outputSignal[2 * i] * 32767.0); // Convert double to 16-bit sample
        fwriteShortLSB(sample, outputFile);
    }
    
    printf("\nOutput written to %s\n", outputTone);

    // Clean up
    free(inputSignal);
    free(IRsignal);
    free(outputSignal);

    fclose(inputFile);
    fclose(IRfile);
    fclose(outputFile);
    

    return 0;
}
        


   


