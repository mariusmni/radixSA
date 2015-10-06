/*
 * RadixLSDCache.h
 *
 *  Created on: Jul 19, 2012
 *      Author: marius
 */

#ifndef RADIXLSDCACHE_H_
#define RADIXLSDCACHE_H_
#include "Sorter.h"
#include "MergeSorter.h"
#include <vector>

template<class Word>
class RadixLSDCache: public Sorter<Word> {
	static const int charBuckets = 1 << (8 * sizeof(uchar));
	static const int shortBuckets = 1 << (8 * sizeof(ushort));
	static const uint mergeSortThreshold = 512;
	static const uint doubleCharSortThreshold = shortBuckets / 2;

	uint bucketSizeShort[sizeof(Word) / sizeof(ushort)][shortBuckets];
	uint bucketSizeChar[sizeof(Word) / sizeof(uchar)][charBuckets];
	uint bucketStart[shortBuckets];
	vector<Word> tKey;
	vector<uint> tSa;
	MergeSorter<Word> ss;

public:
	RadixLSDCache() {

	}

	~RadixLSDCache() {
	}

	void sort(uint length, uint *data, Word *wKey) {
		//	return;
		if (length < mergeSortThreshold)
			ss.sort(length, data, wKey);
		else if (length > doubleCharSortThreshold)
			sortP<ushort, shortBuckets> (length, data, wKey, bucketSizeShort);
		else
			sortP<uchar, charBuckets> (length, data, wKey, bucketSizeChar);
	}

	template<class T, int N>
	void sortP(uint length, uint *data, Word *wKey, uint bucketSize[][N]) {
		tSa.reserve(length);
		tKey.reserve(length);
		uint *sa = data;
		uint *t_sa = &tSa[0];
		Word *key = wKey;
		Word *t_key = &tKey[0];

		const int passes = sizeof(Word) / sizeof(T);

		for (int r = passes; r; --r)
			memset(bucketSize[r - 1], 0, N * sizeof(uint));

		uint i = 0;
		for (T *k = (T*) key; i < length; ++i, k += passes) {
			for (int r = passes; r; --r) {
				int o = getLSDOffset(passes, r);
				++bucketSize[r - 1][k[o]];
			}
		}

		for (uint r = passes; r; --r) {
			uint *bSize = bucketSize[r - 1];

			for (i = 0; bSize[i] == 0; ++i)
				;
			if (bSize[i] == length)
				continue;

			bucketStart[0] = 0;
			for (i = 1; i < N; ++i)
				bucketStart[i] = bucketStart[i - 1] + bSize[i - 1];

			int o = getLSDOffset(passes, r);
			T *k = ((T*) key) + o;
			for (i = 0; i < length; ++i, k += passes) {
				uint b = bucketStart[*k]++;

				t_sa[b] = sa[i];
				t_key[b] = key[i];
			}

			uint * t = t_sa;
			t_sa = sa;
			sa = t;

			Word * tk = t_key;
			t_key = key;
			key = tk;

		}
		if (t_sa == data) {
			memcpy(t_sa, sa, length * sizeof(uint));
			memcpy(t_key, key, length * sizeof(Word));
		}
	}

	int maxCapacity() {
		return tSa.capacity();
	}
	void freeBuffers() {
		FreeAll(tSa);
		FreeAll(tKey);
	}

};

#endif /* RADIXLSDCACHE_H_ */
