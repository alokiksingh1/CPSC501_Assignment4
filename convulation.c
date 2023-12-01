#include <stdlib.h>
#include <stdio.h>


/*  CONSTANTS  ***************************************************************/
#define PI                3.14159265358979

/*  Test tone frequency in Hz  */
#define FREQUENCY         440.0

/*  Standard sample rate in Hz  */
#define SAMPLE_RATE       44100.0

/*  Number of channels  */
#define MONOPHONIC        1

// function definition
void convolve(float x[], int N, float h[], int M, float y[], int P);
void readWaveFileHeader(WavHeader *header, FILE *inputFile);
void readTone(char *sampleTone, char *impulseTone);
float bytesToFloat(char firstByte, char secondByte); 

// struct to hold all data up until the end of subchunk2
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
    char subchunk2_ID[4];       // "data"
    int subchunk2_Size;         // == NumSamples * NumChannels * BitsPerSample/8
} WavHeader;



// read wavefileheader
// function for reading headers from input files

// void readWaveFileHeader(WavHeader *header, FILE *inputFile){
//     fread(header -> chunk_ID, sizeof(header -> chunk_Size), 1, inputFile);
//     fread(&header->chunk_Size, sizeof(header->chunk_Size), 1, inputFile);
//     fread(header->format, sizeof(header->format), 1, inputFile);
//     fread(header->subchunk1_ID, sizeof(header->subchunk1_ID), 1, inputFile);
//     fread(&header->subchunk1_Size, sizeof(header->subchunk1_Size), 1, inputFile);
//     fread(&header->audioFormat, sizeof(header->audioFormat), 1, inputFile);
//     fread(&header->numChannels, sizeof(header->numChannels), 1, inputFile);
//     fread(&header->sampleRate, sizeof(header->sampleRate), 1, inputFile);
//     fread(&header->byteRate, sizeof(header->byteRate), 1, inputFile);
//     fread(&header->blockAlign, sizeof(header->blockAlign), 1, inputFile);
//     fread(&header->bitsPerSample, sizeof(header->bitsPerSample), 1, inputFile);
//     fread(header->subchunk2_ID, sizeof(header->subchunk2_ID), 1, inputFile);
//     fread(&header->subchunk2_Size, sizeof(header->subchunk2_Size), 1, inputFile);

// }

void rreadWaveFileHeader(WavHeader *header, FILE *inputFile){
    fread(header, sizeof(WavHeader), 1, inputFile);
    // fread(header, sizeof(header), 1, inputFile);
}


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
    fwriteIntLSB(formSize, outputFile);
      
    /*  Form container type  */
    fputs("WAVE", outputFile);

    /*  Format chunk identifier (Note: space after 't' needed)  */
    fputs("fmt ", outputFile);
      
    /*  Format chunk size (chunk_Size)  */
    
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
}



int main(int argc, char *argv[]) {
    

    if (argc != 4){
        printf("Usage : %s \t3 arguments required: (input .wav file) (IR file) (output file)\n", argv[0]);
    }

    FILE *inputFile = fopen(argv[1], "rb");
    FILE *IRfile = fopen(argv[2], "rb");
    FILE *outputFile = fopen(argv[3], "wb");
    

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
    readWaveFileHeader(&inputHeader, inputFile);
    readWaveFileHeader(&IRheader, IRfile);

    
    int outputSize = IRheader.subchunk2_Size * inputHeader.chunk_Size;
    float *x , *y, *h = calloc(outputSize, sizeof(float));

    convolve(x, )

    


    int numChannels =inputHeader.numChannels;
    int NumSamples = inputHeader.subchunk2_Size/(inputHeader.bitsPerSample/8);
    int outputRate = inputHeader.sampleRate;
    int bitsPerSample = inputHeader.bitsPerSample;
    writeWaveFileHeader(numChannels, NumSamples, outputRate , bitsPerSample, outputFile);

    fclose(inputFile);
    fclose(outputFile);
    fclose(IRfile);

    clear(inputHeader);
    clear(outputFile);

    free();
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


/*      purpose:        Writes a 4-byte integer to the file stream, starting
*                       with the least significant byte (i.e. writes the int
*                       in little-endian form).  This routine will work on both
*                       big-endian and little-endian architectures.
*/

size_t fwriteIntLSB(int data, FILE *stream)
{
    unsigned char array[4];

    array[3] = (unsigned char)((data >> 24) & 0xFF);
    array[2] = (unsigned char)((data >> 16) & 0xFF);
    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 4, stream);
}




/*       purpose:       Writes a 2-byte integer to the file stream, starting
*                       with the least significant byte (i.e. writes the int
*                       in little-endian form).  This routine will work on both
*                       big-endian and little-endian architectures.
*/

size_t fwriteShortLSB(short int data, FILE *stream)
{
    unsigned char array[2];

    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 2, stream);
}
