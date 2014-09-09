/* Copyright (c) 2013 Y. William Yu. All rights reserved. */
/* threshold.cpp
 * Changes quality values above [quality score] in a file to [quality score]
*/

#if !defined(NDEBUG)
#define BOOST_MULTI_INDEX_ENABLE_INVARIANT_CHECKING
#define BOOST_MULTI_INDEX_ENABLE_SAFE_MODE
#endif

#define BOOST_MULT_INDEX_DISABLE_SERIALIZATION

//#define VERBOSE

#include "global.h"

#include <iterator>
#include <fstream>


int main( int argc, char *argv[])
{
	if (argc < 3) {
		std::cerr << "Smoothes all quality scores above a threshold." << std::endl;
		std::cerr << "Usage: " << argv[0] << " [quality score] [input_file]" << std::endl;
		std::cerr << "\tWill output a modified file named '[input_file].reduced' with all" << std::endl
				  << "\tascii values above [quality score] in column 11 replaced with" << std::endl
				  << "\t[quality score]." << std::endl;
		exit(-1);
	}
	if (strlen(argv[1])>1) {
		std::cerr << "[quality score] must be a single ascii character." << std::endl;
		exit(-50);
	}

	
	std::FILE * inFile;
	std::FILE * outFile;
	char outfile_name[4096];
	unsigned int name_length;
	int col=1;

	char quality_threshold;
	quality_threshold = argv[1][0];

	int c; // character input
	for (int argwalker = 2; argwalker < argc; ++argwalker) {
		name_length = strlen(argv[argwalker]);
		if (name_length>4000) {
			std::cerr << "Filepaths cannot be longer than 4000 characters." << std::endl;
			exit(-2);
		}
		strcpy(outfile_name,argv[argwalker]);
		outfile_name[name_length++]='.';
		outfile_name[name_length++]='r';
		outfile_name[name_length++]='e';
		outfile_name[name_length++]='d';
		outfile_name[name_length++]='u';
		outfile_name[name_length++]='c';
		outfile_name[name_length++]='e';
		outfile_name[name_length++]='d';
		outfile_name[name_length++]='\0';
		std::cout << argv[argwalker] << " --> " << outfile_name << std::endl;
		inFile=fopen(argv[argwalker],"r");
		if (inFile==NULL) {
			perror ("Error opening input file.");
			exit(-3);
		}
		outFile=fopen(outfile_name,"w");
		if (outFile==NULL) {
			perror ("Error opening output file.");
			exit(-4);
		}

		while ( (c=getc(inFile)) !=EOF)
		{
			if (col == 11) {
				//printf(".");
				switch (c) {
					case '\t':
						++col;
						putc(c,outFile);
						break;
					case '\n':
						col = 1;
						putc(c,outFile);
						break;
					default:
						if (c>quality_threshold) {
							putc(quality_threshold,outFile);
						}
						else
							putc(c,outFile);
						break;
				}

			} else {
				switch (c) {
					case '\n':
						col=0;
						// Deliberate fall-through
					case '\t':
						++col;
						// Deliberate fall-through
					default:
						break;
				}
				putc(c,outFile);
			}
		}
		fclose(inFile);
		fclose(outFile);
	}
	return 0;
};

