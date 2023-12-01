/******************************************************************************
*
*     Program:       testtone
*     
*     Description:   This program generates a two-second 440 Hz test tone and
*                    writes it to a .wav file.  The sound file has 16-bit
*                    samples at a sample rate of 44.1 kHz, and is monophonic. 
*
*     Author:        Leonard Manzara
*
*     Date:          November 21, 2009
*
******************************************************************************/


/*  HEADER FILES  ************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/*  CONSTANTS  ***************************************************************/
#define PI                3.14159265358979

/*  Test tone frequency in Hz  */
#define FREQUENCY         440.0

/*  Test tone duration in seconds  */
#define DURATION          2.0				

/*  Standard sample rate in Hz  */
#define SAMPLE_RATE       44100.0

/*  Standard sample size in bits  */
#define BITS_PER_SAMPLE   16

/*  Standard sample size in bytes  */		
#define BYTES_PER_SAMPLE  (BITS_PER_SAMPLE/8)

/*  Number of channels  */
#define MONOPHONIC        1
#define STEREOPHONIC      2


/*  FUNCTION PROTOTYPES  *****************************************************/
void createTestTone(double frequency, double duration,
                    int numberOfChannels, char *filename);
void writeWaveFileHeader(int channels, int numberSamples,
                         double outputRate, FILE *outputFile);
size_t fwriteIntLSB(int data, FILE *stream);
size_t fwriteShortLSB(short int data, FILE *stream);



/******************************************************************************
*
*       function:       main
*
*       purpose:        Creates the test tone and writes it to the
*                       specified .wav file.
*			
*       arguments:      argv[1]:  the filename for the output .wav file
*                       
******************************************************************************/

int main (int argc, char *argv[])
{
    char *outputFilename = NULL;
	
    /*  Process the command line arguments  */
    if (argc == 2) {
        /*  Set a pointer to the output filename  */
        outputFilename = argv[1];
    }
    else {
        /*  The user did not supply the correct number of command-line
            arguments.  Print out a usage message and abort the program.  */
        fprintf(stderr, "Usage:  %s outputFilename\n", argv[0]);
        exit(-1);
    }
	
    /*  Create the sine wave test tone, using the specified
        frequency, duration, and number of channels, writing to
        a .wav file with the specified output filename  */
    createTestTone(FREQUENCY, DURATION, MONOPHONIC, outputFilename);
	
    /*  End of program  */
    return 0;
}



/******************************************************************************
*
*       function:       createTestTone
*
*       purpose:        Calculates and writes out a sine test tone to file
*
*       arguments:      frequency:  frequency of the test tone in Hz
*                       duration:  length of the test tone in seconds
*                       numberOfChannels:  number of audio channels
*                       filename:  name of the file to create
*                       
*       internal
*       functions:      writeWaveFileHeader, fwriteShortLSB
*
*       library
*       functions:      ceil, pow, fopen, fprintf, sin, rint, fclose
*
******************************************************************************/

void createTestTone(double frequency, double duration,
                    int numberOfChannels, char *filename)
{
    int i;
		
    /*  Calculate the number of sound samples to create,
        rounding upwards if necessary  */
    int numberOfSamples = (int)ceil(duration * SAMPLE_RATE);
	
    /*  Calculate the maximum value of a sample  */
    int maximumValue = (int)pow(2.0, (double)BITS_PER_SAMPLE - 1) - 1;
	
    /*  Open a binary output file stream for writing */
    FILE *outputFileStream = fopen(filename, "wb");
    if (outputFileStream == NULL) {
        fprintf(stderr, "File %s cannot be opened for writing\n", filename);
        return;
    }
	
    /*  Write the WAVE file header  */
    writeWaveFileHeader(numberOfChannels, numberOfSamples,
                        SAMPLE_RATE, outputFileStream);
						
    /*  Create the sine tone and write it to file  */
    /*  Since the frequency is fixed, the angular frequency
        and increment can be precalculated  */
    double angularFrequency = 2.0 * PI * frequency;
    double increment = angularFrequency / SAMPLE_RATE;
    for (i = 0; i < numberOfSamples; i++) {
        /*  Calculate the sine wave in the range -1.0 to + 1.0  */
        double value = sin(i * increment);
		
        /*  Convert the value to a 16-bit integer, with the
            range -maximumValue to + maximumValue.  The calculated
            value is rounded to the nearest integer  */
        short int sampleValue = rint(value * maximumValue);
		
        /*  Write out the sample as a 16-bit (short) integer
            in little-endian format  */
        fwriteShortLSB(sampleValue, outputFileStream);
		
        /*  If stereo output, duplicate the sample in the right channel  */
        if (numberOfChannels == STEREOPHONIC)
            fwriteShortLSB(sampleValue, outputFileStream);
    }
						
    /*  Close the output file stream  */
    fclose(outputFileStream);
}



/******************************************************************************
*
*       function:       writeWaveFileHeader
*
*       purpose:        Writes the header in WAVE format to the output file.
*
*       arguments:      channels:  the number of sound output channels
*                       numberSamples:  the number of sound samples
*                       outputRate:  the sample rate
*                       outputFile:  the output file stream to write to
*                       
*       internal
*       functions:      fwriteIntLSB, fwriteShortLSB
*
*       library
*       functions:      ceil, fputs
*
******************************************************************************/

void writeWaveFileHeader(int channels, int numberSamples,
                         double outputRate, FILE *outputFile)
{
    /*  Calculate the total number of bytes for the data chunk  */
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
    fwriteShortLSB(BITS_PER_SAMPLE, outputFile);

    /*  Sound Data chunk identifier  */
    fputs("data", outputFile);

    /*  Chunk size  */
    fwriteIntLSB(dataChunkSize, outputFile);
}



/******************************************************************************
*
*       function:       fwriteIntLSB
*
*       purpose:        Writes a 4-byte integer to the file stream, starting
*                       with the least significant byte (i.e. writes the int
*                       in little-endian form).  This routine will work on both
*                       big-endian and little-endian architectures.
*
*       internal
*       functions:      none
*
*       library
*       functions:      fwrite
*
******************************************************************************/

size_t fwriteIntLSB(int data, FILE *stream)
{
    unsigned char array[4];

    array[3] = (unsigned char)((data >> 24) & 0xFF);
    array[2] = (unsigned char)((data >> 16) & 0xFF);
    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 4, stream);
}



/******************************************************************************
*
*       function:       fwriteShortLSB
*
*       purpose:        Writes a 2-byte integer to the file stream, starting
*                       with the least significant byte (i.e. writes the int
*                       in little-endian form).  This routine will work on both
*                       big-endian and little-endian architectures.
*
*       internal
*       functions:      none
*
*       library
*       functions:      fwrite
*
******************************************************************************/

size_t fwriteShortLSB(short int data, FILE *stream)
{
    unsigned char array[2];

    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 2, stream);
}
