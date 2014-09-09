/* Copyright (c) 2013 Y. William Yu. All rights reserved. */
/* sparsify.cpp
 * Sparsifies the quality vector for a read based on k-mer distance from dictionary
*/

#if !defined(NDEBUG)
#define BOOST_MULTI_INDEX_ENABLE_INVARIANT_CHECKING
#define BOOST_MULTI_INDEX_ENABLE_SAFE_MODE
#endif

#define BOOST_MULT_INDEX_DISABLE_SERIALIZATION

//#define VERBOSE

#include "global.h"

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/member.hpp>
#include <iterator>
#include <fstream>

#include <ctime>

using boost::multi_index_container;
using namespace boost::multi_index;

struct read_entry
{
public:
	readseq encoded;
	unsigned long m1, m2;
	unsigned long m3, m4;
	read_entry() {};
	read_entry(
		const readseq encoded_):
		encoded(encoded_)
	{
//		m1 = encoded & 0x0000FFFFFFFFFFFF;
//		m2 = encoded & 0xFFFF0000FFFFFFFF;
//		m3 = encoded & 0xFFFFFFFF0000FFFF;
//		m4 = encoded & 0xFFFFFFFFFFFF0000;

		m1 = encoded & 0x3F3F3F3F3F3F3F3F;
		m2 = encoded & 0xCFCFCFCFCFCFCFCF;
		m3 = encoded & 0xF3F3F3F3F3F3F3F3;
		m4 = encoded & 0xFCFCFCFCFCFCFCFC;
	}

	friend std::ostream& operator<<(std::ostream& os,const read_entry& e)
	{
		os << decode_read(e.encoded);
		os << std::endl;
		return os;
	}
};

//struct id{};
struct m1{};
struct m2{};
struct m3{};
struct m4{};

typedef multi_index_container<
	read_entry,
	indexed_by<
		//hashed_non_unique<tag<encoded>,member<read_entry,std::bitset,&read_entry::encoded::data> >
		//hashed_unique<tag<id>,member<read_entry,uint64_t,&read_entry::encoded> >,
		hashed_non_unique<tag<m1>,member<read_entry,uint64_t,&read_entry::m1> >,
		hashed_non_unique<tag<m2>,member<read_entry,uint64_t,&read_entry::m2> >,
		hashed_non_unique<tag<m3>,member<read_entry,uint64_t,&read_entry::m3> >,
		hashed_non_unique<tag<m4>,member<read_entry,uint64_t,&read_entry::m4> >//,
	>
> read_entry_database;

template<typename Tag,typename MultiIndexContainer>
void print_out_by(
	const MultiIndexContainer& s
)
{
	const typename boost::multi_index::index<MultiIndexContainer,Tag>::type& i=get<Tag>(s);
	typedef typename MultiIndexContainer::value_type value_type;

	/* dump the elements of the index to cout */
	std::copy(i.begin(),i.end(),std::ostream_iterator<value_type>(std::cout));
}

// Function returns a "readseq" with "A/00" in all positions except those which
// are SNPs (i.e. the database has an entry exactly one substitution away from
// tested), which are coded as "G/10" instead.
// 
// In the event that there are no hits in the database, return -1=1111...111
readseq check_red(const read_entry_database &red, const std::unordered_set<readseq> &red_whole, const readseq tested)
{
	std::set<readseq> tmp;
	readseq ans = 0;
	int curr;
	const readseq constant_one = 1;
	const readseq constant_max = -1;

	read_entry work = {tested};
	//assert(work.encoded == tested);
	//std::cerr  << decode_read(tested) << std::endl;

	int numhits = red.get<m1>().count(work.m1) + red.get<m2>().count(work.m2) + red.get<m3>().count(work.m3) + red.get<m4>().count(work.m4);
	if (numhits == 0) {
		return constant_max;
	} else if (numhits < 96) {
		typedef read_entry_database::index<m1>::type read_by_m1;
		std::pair <read_by_m1::iterator, read_by_m1::iterator> range_m1 = red.get<m1>().equal_range(work.m1);
		while (range_m1.first != range_m1.second) {
			tmp.insert((*range_m1.first).encoded);
			range_m1.first++;
		}

		typedef read_entry_database::index<m2>::type read_by_m2;
		std::pair <read_by_m2::iterator, read_by_m2::iterator> range_m2 = red.get<m2>().equal_range(work.m2);
		while (range_m2.first != range_m2.second) {
			tmp.insert((*range_m2.first).encoded);
			range_m2.first++;
		}
	
		typedef read_entry_database::index<m3>::type read_by_m3;
		std::pair <read_by_m3::iterator, read_by_m3::iterator> range_m3 = red.get<m3>().equal_range(work.m3);
		while (range_m3.first != range_m3.second) {
			tmp.insert((*range_m3.first).encoded);
			range_m3.first++;
		}

		typedef read_entry_database::index<m4>::type read_by_m4;
		std::pair <read_by_m4::iterator, read_by_m4::iterator> range_m4 = red.get<m4>().equal_range(work.m4);
		while (range_m4.first != range_m4.second) {
			tmp.insert((*range_m4.first).encoded);
			range_m4.first++;
		}
	
		bool one_subst = false;
		bool exact_match = false;
		for (auto itr = tmp.begin(); itr != tmp.end(); ++itr) {
			if ((curr = subst_find((*itr),(work.encoded))) == -1) {
				exact_match = true;
			}
			else if (curr >= 0) {
				ans |= (constant_one<<(2*curr));
				one_subst = true;
			}
		}
		if (!(one_subst || exact_match))
			return constant_max;
		return ans;
	} else {
		const readseq orig = work.encoded;
		readseq x, y;
		unsigned int curr_value;
		for (int i=0; i<32; ++i) {
			x = orig & ~(3UL << (2*i));
			curr_value = (orig >> (2*i) & 3);
			for (unsigned long j=0; j<4; ++j) {
				if (curr_value == j) {
					// silently do nothing
				} else {
					y = x | (j << (2*i));
					if (red_whole.count(y)==1) {
						ans |= (constant_one << (2*(31-i)));
					}
				}
			}
		}
		if ((ans != 0) || (red_whole.count(orig)==1))
			return ans;
		else
			return constant_max;
	}

}

// Constructs the read_entry_database hash table
void construct_red(read_entry_database &red, std::unordered_set<readseq> &red_whole, std::string dict_filename) {
	std::deque<readseq> mer_list;
	char line[4096];
	std::string line2;
	{
		FILE * dictionary;
		dictionary = fopen (dict_filename.c_str(), "r");
		setvbuf (dictionary, NULL, _IOFBF, BUFSIZ);

		readseq *array;
		unsigned long arraysize = 1024;
		array = (readseq*) malloc(arraysize*sizeof(readseq));
		unsigned long arraycounter = 0;
		std::string encoded;
		readseq reversed;

		if (dictionary == NULL) perror ("Error opening file");
		else
		{
		
			printf("Reading database file: %016d",0);
			while ( fgets (line, 4096, dictionary) != NULL ) {
				encoded=line;
				encoded.resize(32);
				mer_list = encode_read(encoded.c_str());
				array[arraycounter++] = mer_list[0];
				if (arraycounter == arraysize) {
					array = (readseq*) realloc(array,arraysize*2*sizeof(readseq));
					arraysize *= 2;
				}
				if (!(arraycounter%4096)) {
					printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%016ld", arraycounter);
					fflush(stdout);
				}
			}
			printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%016ld", arraycounter);
			printf("\n");

			red.rehash(2*arraycounter);
			red_whole.reserve(2*arraycounter);
			printf("Generating keys table: %016d",0);
			unsigned long i;
			for ( i = 0; i<arraycounter; ++i) {
				red.insert(read_entry(array[i]));
				red_whole.insert(array[i]);
				reversed = rev_compl(array[i]);
				red.insert(reversed);
				red_whole.insert(reversed);
				if (!(i%4096)) {
					printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%016ld", i);
					fflush(stdout);
				}
			}
			printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%016ld", i);
			printf("\n");
			free(array);

		}
		fclose(dictionary);
	}

}

int main( int argc, char *argv[])
{
	int64checker();
	if (argc < 3) {
		std::cerr << "Discards non-SNP quality values for known reads." << std::endl;
		std::cerr << "Usage: " << argv[0] << " dictionary_file input_file(s)" << std::endl;
		std::cerr << "\tInput is assumed to be SAM file without headers. Because this program" << std::endl
				  << "\tonly modifies column 11 if column 10 contains a valid read, however," << std::endl
				  << "\tit may or may not work on SAM files with headers. This behaviour" << std::endl
				  << "\tshould *not* be depended on." << std::endl << std::endl
				  << "\tOutput will be modified SAM files named [input_file].filtered," <<std::endl
				  << "\tidentical to the original except in column 11." << std::endl << std::endl
				  << "\tThe quality of nearly all bases corresponding to a 32-mer listed in the 1st" << std::endl
				  << "\tcolumn of [dictionary_file] will be set to '~'. Correspondance shall be" << std::endl
				  << "\tdefined as a Hamming distance of less than or equal to 1 in the 32-mer." << std::endl
				  << "\tNote that a 32-mer can correspond to multiple 32-mers in the dictionary." << std::endl
				  << "\tThe exceptions to the setting of quality values above will be all bases" << std::endl
				  << "\tthat are different to any of the corresponding 32-mers in the dictionary," << std::endl
				  << "\twhich will retain their original quality value data." << std::endl;
		exit(-1);
	}


	std::string dict_filename = argv[1];
	read_entry_database red;
	std::unordered_set<readseq> red_whole;
	construct_red(red, red_whole, dict_filename);

	std::FILE * inFile;
	std::FILE * outFile;
	char outfile_name[4096];
	char infile_name[4096];
	unsigned int name_length;
	int col=1;
	char str[4096];
	char timestr[4096];

	//int c; // character input
	int j; // string counter
	std::vector<readseq> mer_list;
	for (int argwalker = 2; argwalker < argc; ++argwalker) {
		strcpy(infile_name,argv[argwalker]);
		name_length = strlen(infile_name);
		if (name_length>4000) {
			std::cerr << "Filepaths cannot be longer than 4000 characters." << std::endl;
			exit(-2);
		}
		strcpy(outfile_name,infile_name);
		outfile_name[name_length++]='.';
		outfile_name[name_length++]='f';
		outfile_name[name_length++]='i';
		outfile_name[name_length++]='l';
		outfile_name[name_length++]='t';
		outfile_name[name_length++]='e';
		outfile_name[name_length++]='r';
		outfile_name[name_length++]='e';
		outfile_name[name_length++]='d';
		outfile_name[name_length++]='\0';

		std::cout << infile_name << " ==> " << outfile_name << ":\t";
		// Timer
		std::time_t result = std::time(NULL);
		strcpy(timestr,std::asctime(std::localtime(&result)));
		timestr[strlen(timestr)-1]='\0';
		std::cerr << timestr << " --> " ;

		inFile=fopen(infile_name,"r");
		if (inFile==NULL) {
			perror ("Error opening input file.");
			exit(-3);
		}
		setvbuf (inFile, NULL, _IOFBF, BUFSIZ);
		outFile=fopen(outfile_name,"w");
		if (outFile==NULL) {
			perror ("Error opening output file.");
			exit(-4);
		}
		setvbuf (outFile, NULL, _IOFBF, BUFSIZ );


		const readseq constant_max = -1;
		readseq a;
		readseq b;
		//int counter = 0;
		int lasty;
		bool alter_this[4096] = { 0 };
		bool pass_through_line;
		char line[4096];
		char *c;
        while ( fgets (line, 4096, inFile) != NULL ) {
			c = &line[0];
			pass_through_line = false;
			j = 0;
			col = 1;
			for (int x=0; x<4096; ++x)
				alter_this[x]=false;
			while ( *c != '\0' )
			{
				if (pass_through_line == true) {
					break;
				} else if (col == 10) {
					switch (*c) {
						case 'A': case 'a':
						case 'C': case 'c':
						case 'G': case 'g':
						case 'T': case 't':
							str[j++]=(char) *c;
							break;
						case '\t': case '\0':
							str[j]='\0';
							mer_list = encode_read_vector(str);
							lasty=-100;
							for (unsigned int y=0; y<mer_list.size(); ++y) {
								a = mer_list[y];
								if (((alter_this[y])||(y-lasty<16))&&(!(y==(mer_list.size()-1)))) {
									// Silently do nothing
								} else {
									lasty=y;
									b = check_red(red,red_whole,a);
									if (b==constant_max) {
										// Silently do nothing	
									} else {
										for (int x=0; x<32; ++x) {
											alter_this[x+y] = alter_this[x+y] || !(b&1);
											b>>=2;
										}
									}
									
								/*	b = check_red(red,red_whole,rev_compl(a));
									if (b==constant_max) {
										// Silently do nothing	
									} else {
										for (int x=0; x<32; ++x) {
											alter_this[31-x+y] = alter_this[31-x+y] || !(b&1);
											b>>=2;
										}
									} */
								}
							}
							j=0;
							++col;
							break;
						case '\n':
							// This should NOT happen
							std::cerr << "Malformed read. Not enough columns." << std::endl;
							// Deliberate fall-through here.
						default:
							// If something other than A,T,C,G (perhaps N) appears, just skip line.
							pass_through_line = true;
							break;
					}
				} else if (col == 11) {
					switch (*c) {
						case '\t':
							pass_through_line = true;
							break;
						case '#':
							break;
						default:
							if ((alter_this[j++]))
								*c = 126;
							break;
					}

				} else {
					if (*c == '\t')
						++col;
				}
				c++;
			} 
			fputs(line, outFile);
		}

		// Timer
		std::time_t result2 = std::time(NULL);
		std::cerr << std::asctime(std::localtime(&result2));

		fclose(inFile);
		fclose(outFile);
	}
	return 0;
};

