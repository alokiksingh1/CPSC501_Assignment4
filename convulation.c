#include <stdlib.h>
#include <stdio.h>
#include <math.h>


/*  CONSTANTS  ***************************************************************/
#define PI                3.14159265358979

/*  Standard sample rate in Hz  */
#define SAMPLE_RATE       44100.0

/*  Number of channels  */
#define MONOPHONIC        1

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

// function definition
void convolve(float x[], int N, float h[], int M, float y[], int P);
// void readWaveFileHeader(WavHeader *header, FILE *inputFile);
float bytesToFloat(char firstByte, char secondByte);
void readTone(char *sampleTone, char *impulseTone, char *outputTone);
size_t fwriteShortLSB(short int data, FILE *stream);
size_t fwriteIntLSB(int data, FILE *stream);
void writeWaveFileHeader(int channels, int numberSamples, double outputRate, int bitsPerSample,  FILE *outputFile);

// Function to convert two bytes to one float in the range -1 to 1
// This is used as .wav files store data in short format (typically 16 bits, can also be extracted from the bits_per_sample header)
// assumes the data is read in bytes
float bytesToFloat(char firstByte, char secondByte) {
    // Convert two bytes to one short (little endian)
    short s = (secondByte << 8) | firstByte;
    // Convert to range from -1 to (just below) 1
    return s / 32768.0;
}




// Writing the header in WAVE format to the output file.
void writeWaveFileHeader(int channels, int numberSamples, double outputRate, int bitsPerSample,  FILE *outputFile)
{
    /*  Calculate the total number of bytes for the data chunk  */
    int BYTES_PER_SAMPLE = bitsPerSample/8;
    int dataChunkSize = channels * numberSamples * BYTES_PER_SAMPLE;
	
    /*  Calculate the total number of bytes for the form size  */
    int formSize = 36 + dataChunkSize;
	
    /*  Calculate the total number of bytes per frame  */
    short int frameSize = channels * BYTES_PER_SAMPLE;
	
    /*  Calculate the byte rate  */
    int bytesPerSecond = (int)ceil(outputRate * frameSize);

    /*  Write header to file  */
    /*  Form container identifier  */
    fputs("RIFF", outputFile);
      
    /*  Form size  */
    fwriteIntLSB(formSize, outputFile);
      
    /*  Form container type  */
    fputs("WAVE", outputFile);

    /*  Format chunk identifier (Note: space after 't' needed)  */
    fputs("fmt ", outputFile);
      
    /*  Format chunk size (fixed at 16 bytes)  */
    fwriteIntLSB(16, outputFile);

    /*  Compression code:  1 = PCM  */
    fwriteShortLSB(1, outputFile);

    /*  Number of channels  */
    fwriteShortLSB((short)channels, outputFile);

    /*  Output Sample Rate  */
    fwriteIntLSB((int)outputRate, outputFile);

    /*  Bytes per second  */
    fwriteIntLSB(bytesPerSecond, outputFile);

    /*  Block alignment (frame size)  */
    fwriteShortLSB(frameSize, outputFile);

    /*  Bits per sample  */
    fwriteShortLSB(bitsPerSample, outputFile);

    /*  Sound Data chunk identifier  */
    fputs("data", outputFile);

    /*  Chunk size  */
    fwriteIntLSB(dataChunkSize, outputFile);
      

    // /*  Bits per sample  */
    // unsigned char array[5];
    // array[0] = (unsigned char)(bitsPerSample);
    // array[1] = (unsigned char)(dataChunkSize);
    // array[2] = (unsigned char)(bytesPerSecond);
    // array[3] = (unsigned char)((int)outputRate);
    // array[4] = (unsigned char)((short)channels);
    // fwrite(array, 1, sizeof(unsigned char), outputFile);
    fputs("data", outputFile);

}



int main(int argc, char *argv[]) {
    
    char *inputFile, *IRfile, outputFile;

    if (argc != 4){
        fprintf(stderr, "Usage : %s \t3 arguments required: (input file) (IR file) (output file)\n", argv[0]);
    }

    inputFile = argv[1];
    IRfile = argv[2];
    outputFile = argv[3];

    readTone(inputFile, IRfile, outputFile);
    
    return 0;
}


// function to convolve x[n](Input) with h[n](IR) to produce y[n](Output)
void convolve(float x[], int N, float y[], int M, float h[], int P){
    int n, m;

    /* Clear Output Buffer y[] */
    for (n = 0; n < P; n++)
        y[n] = 0.0;
    
    // outer loop:  for each input vale x[n] to be processed in turn
    for (n = 0; n < N; n++){
        // inner loop: for all x[n] to be processed with impulse sample h[n]
        for (m = 0; m < M; m++){
            // here * means multiply
            y[n + m] += x[n] * h[m];
            
        }
    }
}


void readTone(char *sampleTone, char *impulseTone, char *outputTone){
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

    
   
    // fscanf(inputFile, "%s", inputFile); // fscanf is used to read the input file
    // fscanf(IRfile, "%s", IRfile); // fscanf to read IR file
    WavHeader inputHeader, IRheader;
    // read header subchunk 1
    fread(&inputHeader, sizeof(WavHeader), 1, inputFile);
    fread(&IRheader, sizeof(WavHeader), 1, IRfile);

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

    int num_samples = subchunk2_size_sample / (inputHeader.bitsPerSample / 8); // number of data points in the sample
    int num_impulse = subchunk2_size_impulse / (IRheader.bitsPerSample / 8); // number of data points in the impulse
    // frequency = 1 / (duration / numberOfSamples)
    // duration = numberOfSamples / frequency
    // numberOfSamples = duration * SAMPLE_RATE
    

    


    /*
    Let's print out some header values
    */
    printf("audio format : %s\n", inputHeader.format);
    printf("number of channels : %s\n", inputHeader.numChannels);
    printf("size of subchunk1 : %s\n", inputHeader.subchunk1_Size);
    printf("number of samples : %s\n", num_samples);
    printf("\n***************************************************\n");
    printf("IR audio format : %s\n", IRheader.format);
    printf("IR number of channels : %s\n", IRheader.numChannels);
    printf("IR size of subchunk1 : %s\n", IRheader.subchunk1_Size);
     printf("IR number of samples : %s\n", num_impulse);

    int outputSize = num_samples + num_impulse - 1;
    float *x = (float*)malloc(num_samples * sizeof(float));
    float *y = (float*)malloc(num_impulse * sizeof(float));
    float *h = (float*)malloc(outputSize * sizeof(float));
    // *h = calloc(outputSize, sizeof(float));
    if (!x || !y || !h) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        if (x) free(x);
        if (y) free(y);
        if (h) free(h);
        // Close any opened files here as well
        fclose(inputFile);
        fclose(IRfile);
        fclose(outputFile);
        exit(-1);
    }


    for (int i = 0; i < num_samples; ++i) {
        char buffer[2];
        fread(buffer, sizeof(char), 2, inputFile);
        x[i] = bytesToFloat(buffer[0], buffer[1]);
    }

    for (int i = 0; i < num_impulse; ++i) {
        char buffer[2];
        fread(buffer, sizeof(char), 2, IRfile);
        y[i] = bytesToFloat(buffer[0], buffer[1]);
    }

    // convolve(x, subchunk2_size_sample, y, subchunk2_size_impulse, h, outputSize);

    convolve(x, num_samples, h, num_impulse, y, outputSize);

    int numChannels =inputHeader.numChannels;
    int outputRate = inputHeader.sampleRate;
    int bitsPerSample = inputHeader.bitsPerSample;
    int NumSamples = subchunk2_size_sample/(bitsPerSample/8);
    
    //writeWaveFileHeader(numChannels, NumSamples, outputRate , bitsPerSample, outputFile);

    int IR_numChannels =IRheader.numChannels;
    int IR_outputRate = IRheader.sampleRate;
    int IR_bitsPerSample = IRheader.bitsPerSample;
    int IR_NumSamples = subchunk2_size_impulse/(IR_bitsPerSample/8);
    
    
    //fwrite(h, sizeof(float), outputSize, outputFile);
    fwrite(h, outputSize, 1, outputFile);

    writeWaveFileHeader(num_samples, num_impulse, IR_outputRate , IR_bitsPerSample, outputFile);

    for (int i = 0; i < outputSize; ++i) {
        // Convert h[i] from float to short
        short sample = (short)(h[i] * 32768.0);
        fwrite(&sample, sizeof(short), 1, outputFile);
    }

    fclose(inputFile);
    fclose(outputFile);
    fclose(IRfile);

    free(h);
    free(x);
    free(y);

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