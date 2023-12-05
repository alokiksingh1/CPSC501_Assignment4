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
void readWaveFileHeader(WavHeader *header, FILE *inputFile);
void readTone(char *sampleTone, char *impulseTone, char *outputTone);
float bytesToFloat(char firstByte, char secondByte);
size_t fwriteShortLSB(short int data, FILE *stream);
size_t fwriteIntLSB(int data, FILE *stream);








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
void writeWaveFileHeader(int channels, int numberSamples,
                         double outputRate, int bitsPerSample, FILE *outputFile)
{
    /*  Calculate the total number of bytes for the data chunk  */
    int Bytes_per_Sample = bitsPerSample/8;
    int dataChunkSize = channels * numberSamples * Bytes_per_Sample;
	
    /*  Calculate the total number of bytes for the form size  */
    int formSize = 36 + dataChunkSize;
	
    /*  Calculate the total number of bytes per frame  */
    short int frameSize = channels * Bytes_per_Sample;
	
    /*  Calculate the byte rate  */
    int bytesPerSecond = (int)ceil(outputRate * frameSize);

    /*  Write header to file  */
    /*  Form container identifier  */
    fputs("RIFF", outputFile);
      
    /*  Form size  */
    //fwriteIntLSB(formSize, outputFile);
      
    /*  Form container type  */
    fputs("WAVE", outputFile);

    /*  Format chunk identifier (Note: space after 't' needed)  */
    fputs("fmt ", outputFile);
      
    /*  Format chunk size (chunk_Size)  */
    
    //fwriteIntLSB(16, outputFile);

    /*  Compression code:  1 = PCM  */
    //fwriteShortLSB(1, outputFile);

    /*  Number of channels  */
    //fwriteShortLSB((short)channels, outputFile);

    /*  Output Sample Rate  */
    //fwriteIntLSB((int)outputRate, outputFile);

    /*  Bytes per second  */
    //fwriteIntLSB(bytesPerSecond, outputFile);

    /*  Block alignment (frame size)  */
    //fwriteShortLSB(frameSize, outputFile);

    /*  Bits per sample  */
    unsigned char array[5];
    array[0] = (unsigned char)(bitsPerSample);
    array[1] = (unsigned char)(dataChunkSize);
    array[2] = (unsigned char)(bytesPerSecond);
    array[3] = (unsigned char)((int)outputRate);
    array[4] = (unsigned char)((short)channels);
    fwrite(array, 1, sizeof(unsigned char), outputFile);

    /*  Sound Data chunk identifier  */
    fputs("data", outputFile);

    /*  Chunk size  */
    //fwriteIntLSB(dataChunkSize, outputFile);
}



int main(int argc, char *argv[]) {
    
    char *inputFile, *IRfile, outputFile;

    if (argc != 4){
        printf(stderr, "Usage : %s \t3 arguments required: (input file) (IR file) (output file)\n", argv[0]);
    }

    if (argc == 4){
        inputFile = argv[1];
        IRfile = argv[2];
        outputFile = argv[3];

    }

    
    readTone(inputFile, IRfile, outputFile);
    
    
    
    fscanf(inputFile, "%s", inputFile); // fscanf is used to read the input file
    // now let's calculate frequency and duration of the input .wav file
    // we need to use the formula: frequency = 1 / (duration / numberOfSamples)
    // we need to use the formula: duration = numberOfSamples / frequency
    // we need to use the formula: numberOfSamples = duration * SAMPLE_RATE

    // subchunk2_Size chunk size
    int subchunk2_Size;
    fread(&subchunk2_Size, 4, 1, inputFile);
    int chunkSize = 36 + subchunk2_Size;
    
    // find duration
    double duration = chunkSize / SAMPLE_RATE;

    // find number of samples
    int numberOfSamples = duration * SAMPLE_RATE;
    
    double Frequency = 1 / (duration / numberOfSamples);

    // **********************************************************************

    // now similarly, let's calculate for IR file as well

    fscanf(IRfile, "%s", IRfile); // fscanf is used to read the input file
    
    // subchunk2_Size chunk size
    int IR_subchunk2_Size;
    fread(&IR_subchunk2_Size, 4, 1, IRfile);
    int IR_chunkSize = 36 + IR_subchunk2_Size;
    
    // find duration
    double IR_duration = IR_chunkSize / SAMPLE_RATE;

    // find number of samples
    int IR_numberOfSamples = IR_duration * SAMPLE_RATE;

    double IR_Frequency = 1 / (IR_duration / IR_numberOfSamples);


    WavHeader inputHeader, IRheader;
    // read header subchunk 1
    fread(&inputHeader, sizeof(inputHeader), 1, inputFile);
    fread(&IRheader, sizeof(IRheader), 1, IRfile);


    int outputSize = subchunk2_Size * IR_subchunk2_Size;
    float *x , *y, *h = malloc(outputSize);
    // *h = calloc(outputSize, sizeof(float));

    convolve(x, subchunk2_Size, y, IR_subchunk2_Size, h, outputSize);

    


    int numChannels =inputHeader.numChannels;
    int outputRate = inputHeader.sampleRate;
    int bitsPerSample = inputHeader.bitsPerSample;
    int NumSamples = subchunk2_Size/(bitsPerSample/8);
    
    writeWaveFileHeader(numChannels, NumSamples, outputRate , bitsPerSample, outputFile);

    int IR_numChannels =IRheader.numChannels;
    int IR_outputRate = IRheader.sampleRate;
    int IR_bitsPerSample = IRheader.bitsPerSample;
    int IR_NumSamples = IR_subchunk2_Size/(IR_bitsPerSample/8);
    
    // writeWaveFileHeader(IR_numChannels, IR_NumSamples, IR_outputRate , IR_bitsPerSample, outputFile);
    fwrite(h, sizeof(float), outputSize, outputFile);

    fclose(inputFile);
    fclose(outputFile);
    fclose(IRfile);

    free(h);
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

}