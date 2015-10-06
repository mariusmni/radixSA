/*
 * runtimes.cpp
 *
 *  Created on: Oct 3, 2011
 *      Author: marius
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "radix.h"

//#define BENCHMARK


uchar* readData(const char * const filename, uint& n) {
	struct stat fileInfo;
	FILE *file;
	if (stat(filename, &fileInfo)) {
		printf("Unable to get stat of file %s \n", filename);
		return NULL;
	}
	n = fileInfo.st_size;
	uchar *result = new uchar[n];
	if (!(file = fopen(filename, "r"))) {
		printf("Unable to open file %s \n", filename);
		return NULL;
	}
	rewind(file);
	if (n > fread(result, sizeof(uchar), n, file)) {
		printf("Error reading file %s \n", filename);
		fclose(file);
		delete[] result;
		return NULL;
	}
	fclose(file);
	return result;
}

void printIntData(uint *data, uint n, char *outputFile) {
	FILE *f = fopen(outputFile, "w");
	for (uint i = 0; i < n; ++i) {
		fprintf(f, "%d ", data[i]);
	}
	fprintf(f, "\n");
	fclose(f);
}

int main(int argc, char *argv[]) {
	int i = 1;
	char *c = (char *) &i;
#ifdef LITTLE_ENDIAN_FLAG
	bool fail = (*c == 0);
#else
	bool fail = (*c == 1);
#endif
	if (fail) {
		printf(
				"Failed to detect endianness! \nPlease manually undefine LITTLE_ENDIAN_FLAG in radix.h and recompile\n");
		exit(0);
	}

#ifdef BENCHMARK
	compareWithBPR(argc, argv);
	return 0;
#endif

	if (argc < 3) {
		printf("Expected 2 arguments: input_file output_file\n");
		return 0;
	}
	char *inputFile = argv[1];
	char *outputFile = argv[2];
	uint n = 0;
	uchar *ustr = readData(inputFile, n);
	if (ustr == NULL)
		return 0;

	if (ustr[n - 1] == 10) { // remove line feed
		ustr[--n] = 0;
	}

	cout << "Read " << n << " bytes from file " << inputFile << endl;
	cout << "Computing Suffix Array..." << endl;
	clock_t startTime = clock();
	unsigned int *radixSA = Radix(ustr, n).build();
	clock_t endTime = clock();
	float mseconds = clock_diff_to_msec(endTime - startTime);
	printf("RadixSA took [%.2fs]", mseconds/1000.0);
	cout << endl << "Writing output to file " << outputFile << "..." << endl;
	printIntData(radixSA, n, outputFile);
	delete[] radixSA;

	return 0;
}

#ifdef BENCHMARK
extern "C" {
#include "kbs_Error.h"
#include "kbs_Math.h"
#include "kbs_String.h"
#include "kbs_SuffixArray.h"
//#include "kbs_SuffixArrayAnnotated.h"
#include "kbs_SuffixArrayChecker.h"
#include "kbs_SuffixArrayConstDStepAndPre.h"
#include "kbs_Types.h"
}
#endif
#ifdef BENCHMARK

uint *findSA(Kbs_Ustring *ustr) {
	Radix r(ustr->str, ustr->strLength);
	return r.build();
}

void findSAAndDelete(Kbs_Ustring *ustr) {
	uint *sa = findSA(ustr);
	delete[] sa;
}

Kbs_SuffixArray *bprFindSA(Kbs_Ustring *ustr) {
	Kbs_Ulong q = 3;
	kbs_get_AlphabetForUstring(ustr);
	if (ustr == NULL) {
		KBS_ERROR( KBS_ERROR_NULLPOINTER);
		exit(1);
	}
	//		if (argc == 2) {
	if (ustr->alphabet->alphaSize <= 9) {
		q = 7;
	}
	if (9 < ustr->alphabet->alphaSize && ustr->alphabet->alphaSize <= 13) {
		q = 6;
	}
	if (13 < ustr->alphabet->alphaSize && ustr->alphabet->alphaSize <= 21) {
		q = 5;
	}
	if (13 < ustr->alphabet->alphaSize && ustr->alphabet->alphaSize <= 21) {
		q = 5;
	}
	if (21 < ustr->alphabet->alphaSize && ustr->alphabet->alphaSize <= 46) {
		q = 4;
	}
	if (46 < ustr->alphabet->alphaSize) {
		q = 3;
	}
	// implementation using direct pointers as bucket pointers
	return kbs_buildDstepUsePrePlusCopyFreqOrder_SuffixArray(ustr, q);
}

void bprFindSAAndDelete(Kbs_Ustring *ustr) {
	Kbs_SuffixArray *sa;
	sa = bprFindSA(ustr);
	kbs_delete_SA_IncludingString(sa);
}

void compareAnswers(unsigned char* data, unsigned int* mySA,
		unsigned long* bprSA, int bases) {
	int i;
	int j = 0;
	unsigned int m = bases;
	for (i = 0; i < bases; ++i) {
		if (mySA[i] != bprSA[i]) {
			if (mySA[i] < m)
				m = mySA[i], j = i;
		}
	}
	i = j;
	if (m < bases) {
		printf("NO GOOD! %d\n", i);
		printf("Suffix at %d = %d ...\n", mySA[i], data[mySA[i]]);
		printf("Suffix at %d = %d ...\n", bprSA[i], data[bprSA[i]]);
		for (int j = 0; j <= 40; ++j) {
			printf("%d ", data[mySA[i] + j]);
		}
		printf("\n");
		for (int j = 0; j <= 40; ++j) {
			printf("%d ", data[bprSA[i] + j]);
		}
		printf("\n");
	} else {
		printf("OK\n");
	}

}

int compareWithBPR(int argc, char *argv[]) {
	srand(123);
	Kbs_Ustring* ustr = NULL;

	bool timeOnlyMode = false;

	double runTime[argc][2];
	for (int i = 1; i < argc; ++i) {
		ustr = kbs_getUstring_FromFile(argv[i]);
//		if (ustr->str[ustr->strLength - 1] == 10) { // remove line feed
//			ustr->str[--ustr->strLength] = 0;
//		}
		if (!timeOnlyMode) {
			cout << "===============================================" << endl;
			cout << "string length " << ustr->strLength << " of file " << endl
					<< argv[i] << endl;
		}

		double seconds;
		unsigned int *radixSA;
		Kbs_SuffixArray *sa;

		{
			if (timeOnlyMode) {
				//				seconds = time_it(findSAAndDelete, ustr, runs);
				clock_t startTime = clock();
				findSAAndDelete(ustr);
				clock_t endTime = clock();
				seconds = clock_diff_to_msec(endTime - startTime);
			} else {
				clock_t startTime = clock();
				radixSA = findSA(ustr);
				clock_t endTime = clock();
				seconds = clock_diff_to_msec(endTime - startTime);
			}

			runTime[i][0] = seconds;
			printf("Radix \t %s \t %.3fs ", runTime[i][0] / 1000.0, argv[i]);
			cout << endl;
		}

//		if (false)
		{
			if (timeOnlyMode) {
				//				seconds = time_it(bprFindSAAndDelete, ustr, runs);
				clock_t startTime = clock();
				bprFindSAAndDelete(ustr);
				clock_t endTime = clock();
				seconds = clock_diff_to_msec(endTime - startTime);
			} else {
				clock_t startTime = clock();
				sa = bprFindSA(ustr);
				clock_t endTime = clock();
				seconds = clock_diff_to_msec(endTime - startTime);
			}
			runTime[i][1] = seconds;
			printf("BPR \t %s \t %.3fs", argv[i], runTime[i][1] / 1000.0);
			cout << endl;
		}

		if (!timeOnlyMode) {
			int n = ustr->strLength;
			compareAnswers(ustr->str, radixSA, sa->posArray, n);
			kbs_delete_SA_IncludingString(sa);
			delete[] radixSA;
			cout
					<< "===== Dataset ============= \t Radix ==== \t BPR ==============="
					<< endl;
			for (int j = 1; j <= i; ++j) {
				printf("%-25s \t %.2f \t %.2f\n", argv[j], runTime[j][0],
						runTime[j][1]);
			}
		}
		cout << endl;

		//exit(0);
	}
	return 0;
}
#endif
