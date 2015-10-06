/*
 * utils.h
 *
 *  Created on: Oct 17, 2011
 *      Author: marius
 */

#ifndef MY_RADIX_UTILS_H_
#define MY_RADIX_UTILS_H_
#include <time.h>
#include <iostream>
#include <cstdio>
#include <ctime>
using namespace std;

typedef unsigned int uint;
typedef unsigned char uchar;
typedef unsigned short int ushort;
//typedef unsigned long long int ulong;
const int uintBits = 8 * sizeof(uint);

#ifdef WORD64
typedef unsigned long long word;
#define DIV_BY_WORDSIZE(p) ((p) >> 6)
#define MOD_WORDSIZE(p) ((p) & 63)
#else
typedef uint word;
#define DIV_BY_WORDSIZE(p) ((p) >> 5)
#define MOD_WORDSIZE(p) ((p) & 31)
#endif

#ifdef WORD64BUCK
typedef unsigned long long bword;
#else
typedef uint bword;
#endif

const int bitsPerWord = 8 * sizeof(word);
const int bitsPerBWord = 8 * sizeof(bword);
static const int maxChar = 256;

#ifdef LITTLE_ENDIAN_FLAG
#define getBufferStart(buffer, len) ((buffer) + (len) - 1)
#define getMSDOffset(increment, remainingKey) ((remainingKey) - 1)
#define getLSDOffset(increment, remainingKey) ((increment) - (remainingKey))
#else
#define getBufferStart(buffer, len) (buffer)
#define getMSDOffset(increment, remainingKey) ((increment) - (remainingKey))
#define getLSDOffset(increment, remainingKey) ((remainingKey) - 1)
#endif

word mask[bitsPerWord][bitsPerWord]; // mask of i 1's shifted to left j positions (e.g. mask[3][2] = 0..011100)
word maskLeft[bitsPerWord]; // mask of i 1's as most signif bits (e.g. mask[3] = 11100..0)
word maskRight[bitsPerWord]; // mask of i 1's  (e.g. mask[3] = 0...00111 = 7)

void init() {
	for (int i = 0; i < bitsPerWord; ++i) {
		word ones = (((word) 1) << i) - 1;
		maskRight[i] = ones;
		maskLeft[i] = ones << (bitsPerWord - i);
		for (int j = 0; j < bitsPerWord; ++j, ones <<= 1)
			mask[i][j] = ones;
	}
}

#define GET_LAST(w, b) ((w) & maskRight[b])
#define GET_FIRST(w, b) (((w) & maskLeft[b]) >> ((bitsPerWord) - (b)))
#define GET_MIDDLE(w, b, o) (((w) & mask[b][o]) >> (o))

inline int max(int a, int b) {
	return a > b ? a : b;
}
inline int min(int a, int b) {
	return a < b ? a : b;
}

/**
 * Take the last1 bits from w1 and the first2 bits from w2 and form a new word
 */
inline word stitch(word w1, word w2, int last1, int first2) {
	return (GET_LAST(w1, last1) << first2) | GET_FIRST(w2, first2);
}

int bitsFor(int n) {
	int b = 1;
	if (n >= 0x10000) {
		b += 16;
		n >>= 16;
	}
	if (n >= 0x100) {
		b += 8;
		n >>= 8;
	}
	if (n >= 0x10) {
		b += 4;
		n >>= 4;
	}
	if (n >= 4) {
		b += 2;
		n >>= 2;
	}
	if (n >= 2)
		b += 1;
	return b;
}

template<class T>
inline void swapItems(T *a, T *b) {
	T tmp = *a;
	*a = *b;
	*b = tmp;
}

template<class T>
void prefixSum(T *s, T *d, int length) {
	*d++ = *s++;
	for (; --length; ++d, ++s)
		*d = *s + *(d - 1);
}

int renameAlphabet(const uchar *in, uint length, uint *code, uint& codeLength) {
	const int maxAlpha = 256;
	memset(code, 0, maxAlpha * sizeof(uint));

	for (; length--;)
		code[*in++] = 1;

//	for (int i = 0; i < 256; ++i)
//		if (code[i])
//			cout <<"char " << (char)i << endl;

	code[0] -= 1;
	prefixSum(code, code, maxAlpha);

	int alphaSize = code[maxAlpha - 1] + 1;
	codeLength = bitsFor(alphaSize-1);

	return alphaSize;
}

clock_t _startTime = clock();
clock_t _endTime;

void resetTime() {
	_startTime = clock();
}

float getTime() {
	_endTime = clock();
	return (float) (_endTime - _startTime) / 1000; //CLOCKS_PER_SEC;
}

#ifdef DEBUG
void printTime(char *msg) {
	float seconds = getTime();
	cout << msg << " " << seconds << endl;
	//	resetTime();
}
#else
#define printTime(msg) {}
#endif

uint max(uint *a, uint n) {
	uint m = *a;
	for (; --n;)
		if (*++a > m)
			m = *a;
	return m;
}

template<class Word>
void insertSort(uint* sa, uint length, Word* key) {
	if (length == 2) {
		if (key[1] < key[0]) {
			swapItems(key, key + 1);
			swapItems(sa, sa + 1);
		}
		return;
	}

	for (uint i = 1; i < length; ++i) {
		Word kp = key[i];
		uint t = sa[i];
		uint j;
		for (j = i; j > 0 && kp < key[j - 1]; --j) {
			key[j] = key[j - 1];
			sa[j] = sa[j - 1];
		}
		key[j] = kp;
		sa[j] = t;
	}
}

#define PRINT_LIMIT 500
void printCharData(const uchar *array, int nKey) {
	int e = (nKey > PRINT_LIMIT) ? PRINT_LIMIT : nKey;
	for (int i = 0; i < e; ++i) {
		uchar c = array[i];
		if (c == '\n')
			c = 'N';
		if (c == '\r')
			c = 'n';
		if (c == '\t')
			c = 'T';
		if (c == 0)
			printf("0");
		else
			printf("%c", c);
	}
	if (e < nKey)
		printf("... %d more", nKey - e);
}

template<typename T>
void FreeAll(T & t) {
	T tmp;
	t.swap(tmp);
}

double clock_diff_to_msec(long clock_diff) {
	return double(clock_diff) / CLOCKS_PER_SEC * 1000;
}

template<class Proc, class Arg>
double time_it(Proc proc, Arg a, int N) // returns time in microseconds
{
	std::clock_t const start = std::clock();
	for (int i = 0; i < N; ++i)
		proc(a);
	std::clock_t const end = std::clock();
	if (clock_diff_to_msec(end - start) < 200)
		return time_it(proc, a, N * 5);
	return clock_diff_to_msec(end - start) / N;
}

template<class T>
void reverseArray(T *a, int n) {
	for (int i = 0, j = n - 1; i < j; ++i, --j) {
		T t = a[i];
		a[i] = a[j];
		a[j] = t;
	}
}
#endif
