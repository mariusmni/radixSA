/*
 * radix.h
 *
 *  Created on: Oct 21, 2011
 *      Author: marius
 */

#ifndef RADIX_H_
#define RADIX_H_

#define NDEBUG
#include <cassert>
#include <vector>

#ifndef NDEBUG
#define DEBUG
#endif

#include <math.h>
#define LITTLE_ENDIAN_FLAG
#define WORD64
#define WORD64BUCK

#ifdef WORD64BUCK
#define BWORD64FLAG 1
#else
#define BWORD64FLAG 0
#endif

#include "utils.h"
#include "RadixLSDCache.h"

class Radix {
private:
	const uchar *originalInput;
	const uint length;

	uint *sa;

	vector<uchar> bucketFlagVec;
	uchar *bucketFlag;

	vector<uint> bucketVec;
	uint *bucket;

	vector<uint> uintBuff1;
	vector<uint> uintBuff2;

	word *input;
	vector<bword> bucketBuffVec;

	RadixLSDCache<bword> sorter;

	uint bitsPerChar;
	uint charCode[maxChar];
	word shiftedCode[maxChar];
	uint expectedSorted;
	const uint bucketPiggyBackBits;
	const uint bucketPiggyBackLimit;
	const uint bucketPiggyBackMask;

	void prepareInput(word *input) {
		int remainBits = bitsPerWord;
		word l = 0;
		int nWords = 0;
		for (uint i = 0; i < length; ++i) {
			uint code = charCode[originalInput[i]];
			remainBits -= bitsPerChar;
			if (remainBits >= 0) {
				l = (l << bitsPerChar) | code;
			} else {
				int bitsInNextWord = -remainBits;
				int bitsAtEnd = bitsPerChar - bitsInNextWord;
				l = (l << bitsAtEnd) | (code >> bitsInNextWord);
				input[nWords++] = l;
				l = GET_LAST(code, bitsInNextWord);
				remainBits = bitsPerWord - bitsInNextWord;
			}
		}
		if (remainBits < bitsPerWord)
			input[nWords++] = l << remainBits;
		input[nWords] = 0;
	}

	void shiftCode(int bitsAtATime) {
		for (int i = 0; i < maxChar; ++i)
			shiftedCode[i] = ((word) charCode[i]) << bitsAtATime;
	}

	void firstPass(int bitsPerFirstPass, int n, uint *bStart, uint *bSize) {
		shiftCode(bitsPerFirstPass - bitsPerChar);

		memset(bSize, 0, n * sizeof(uint));

		{
			uint w = 0;
			for (int i = length - 1; i >= 0; --i) {
				int c = originalInput[i];
				w >>= bitsPerChar;
				w |= shiftedCode[c];
				bSize[w]++;
			}
		}

		bStart[0] = 0;
		prefixSum(bSize, bStart + 1, n - 1);

		{
			uint w = 0;
			for (int i = length - 1; i >= 0; --i) {
				int c = originalInput[i];
				w >>= bitsPerChar;
				w |= shiftedCode[c];
				sa[bStart[w]++] = i;
			}

		}

		for (int i = 0; i < n; ++i)
			bStart[i] -= bSize[i];

		int bitVectorSize = length + 1;
		createBucketFlags(bitVectorSize);
		bucketFlag[length] = 1;

		uint zeros = fixEndingZeros();
		bStart[0] += zeros;
		bSize[0] -= zeros;
	}

	word getWordAtBit(int p) {
		int w = DIV_BY_WORDSIZE(p);
		int b = MOD_WORDSIZE(p);
		return b ? stitch(input[w], input[w + 1], bitsPerWord - b, b)
				: input[w];
	}

	void copyWords(uint* lo, int bLen, word* buffer, int bitsToSkip) {
		for (; bLen--;) {
			int bit = (*lo++) * bitsPerChar + bitsToSkip;
			*buffer++ = getWordAtBit(bit);
		}
	}

	void updateBucketStart(word *key, int length, int bStart) {
		bucketFlag[bStart] = 1;
		for (int i = 1; i < length; ++i)
			if (key[i] != key[i - 1])
				bucketFlag[bStart + i] = 1;
	}

	void updateBucketStart(bword *key, int length, int bStart, uchar flag) {
		bucketFlag[bStart] = flag;
		for (int i = 1; i < length; ++i)
			if (key[i] != key[i - 1])
				bucketFlag[bStart + i] = flag;
	}

	void createBucketFlags(int bitVectorSize) {
		bucketFlagVec.reserve(bitVectorSize);
		bucketFlag = &bucketFlagVec[0];
		memset(bucketFlag, 0, bitVectorSize * sizeof(uchar));
	}

	int secondPass(uint bitsToSkip, uint nBuckets, uint *bStart, uint *bSize) {
		RadixLSDCache<word> sorter;
		uint maxSize = max(bSize, nBuckets);
		vector<word> bufferVec(maxSize);
		word *buffer = &bufferVec[0];
		for (uint i = 0; i < nBuckets; ++i) {
			uint start = bStart[i];
			bucketFlag[start] = 1;

			uint bLen = bSize[i];
			if (bLen > 1) {
				uint * lo = sa + start;
				copyWords(lo, bLen, buffer, bitsToSkip);
				sorter.sort(bLen, lo, buffer);
				updateBucketStart(buffer, bLen, start);
			}
		}

		return bitsPerWord;
	}

	bword getBucketWord(uint s1, uint D) {
#ifdef WORD64BUCK
		return (((bword) bucket[s1]) << uintBits) | bucket[s1 + D];
#else
		return bucket[s1];
#endif
	}

	void copyBucketNumbers(uint *sa, uint D, uint D2, bword *buf, int bLen) {
		for (; bLen--;)
			*buf++ = getBucketWord((*sa++) + D, D2);
	}

	void assignSubBucketNumbersAndLen(uint bStart, uint bLen) {
		uint currentBStart = bStart;
		for (uint e = bStart + bLen, i = bStart + 1; i <= e; ++i) {
			if (bucketFlag[i]) {
				uint prevBucketLen = i - currentBStart;
				uint buck = ((currentBStart + 1) << bucketPiggyBackBits);
				if (prevBucketLen < bucketPiggyBackLimit)
					buck |= prevBucketLen;

				for (; currentBStart < i;)
					bucket[sa[currentBStart++]] = buck;
			}
		}
	}

	void detectPeriods(uint *lo, uint bLen, uint D, uint &periodLength,
			int &maxNbPeriods) {
		periodLength = 0;
		maxNbPeriods = 0;
		int currentPeriodSize = 1;
		for (uint i = 1; i < bLen; ++i) {
			uint diff = lo[i - 1] - lo[i];
			if (diff <= D) { // period
				if (periodLength && diff != periodLength) {
					periodLength = 0;
					maxNbPeriods = 0;
					return; // no good
				}
				periodLength = diff;
				currentPeriodSize++;
			} else {
				if (currentPeriodSize > maxNbPeriods)
					maxNbPeriods = currentPeriodSize;
				currentPeriodSize = 1;
			}
		}
		if (currentPeriodSize > maxNbPeriods)
			maxNbPeriods = currentPeriodSize;
	}

	bool detectAndTreatPeriods(uint bStart, uint bLen, uint D, uint D2,
			uchar bucketStartFlag) {
		uint period;
		int maxPeriodLen;
		detectPeriods(sa + bStart, bLen, D, period, maxPeriodLen);

		if (maxPeriodLen * period > 2 * D)
			return treatPeriodLeaders(bStart, bLen, D, D2, bucketStartFlag,
					period);
		else
			return false;
	}

	bool treatPeriodLeaders(uint bStart, uint bLen, uint D, uint D2,
			uchar bucketStartFlag, uint periodLength) {

		bucketBuffVec.reserve(bLen);
		bword *bufferP = &bucketBuffVec[0];

		uintBuff1.reserve(bLen);
		uint *suffixBuff3P = &uintBuff1[0];

		uint *lo = sa + bStart;
		int leaders[2] = { 0, bLen - 1 };
		for (uint i = 0; i < bLen;) {
			int l = lo[i];
			bword leaderEnd = getBucketWord(l + D - periodLength, D2);
			bword rest = getBucketWord(l + D, D2);

			if (leaderEnd == rest)
				return false; // can happen with the copy from next optimization

			if (leaderEnd < rest) {
				bufferP[leaders[1]] = rest;
				suffixBuff3P[leaders[1]] = i;
				--leaders[1];
			} else {
				bufferP[leaders[0]] = rest;
				suffixBuff3P[leaders[0]] = i;
				++leaders[0];
			}

			for (++i; i < bLen && lo[i - 1] - lo[i] == periodLength; ++i)
				;
		}

		uintBuff2.reserve(bLen);
		uint *dest = &uintBuff2[0];
		int k = 0;
		if (leaders[0] > 1)
			sorter.sort(leaders[0], suffixBuff3P, bufferP);

		if (leaders[0])
			k += populate<false> (leaders[0], periodLength, bStart, bLen,
					bucketStartFlag, suffixBuff3P, dest, bufferP);

		int buffStart = leaders[1] + 1;
		leaders[1] = bLen - buffStart;
		if (leaders[1]) {
			bword *bf = bufferP + buffStart;
			uint *sb = suffixBuff3P + buffStart;
			reverseArray(bf, leaders[1]);
			reverseArray(sb, leaders[1]);
			if (leaders[1] > 1)
				sorter.sort(leaders[1], sb, bf);
			k += populate<true> (leaders[1], periodLength, bStart, bLen,
					bucketStartFlag, sb, dest, bf);
		}

		memcpy(sa + bStart, dest, bLen * sizeof(uint));
		assignSubBucketNumbersAndLen(bStart, bLen);
		return true;
	}

	template<bool isTypeS>
	int populate(int leaders, uint periodLength, uint bStart, uint bLen,
			uint bucketStartFlag, uint *suffixBuff, uint *dest, bword *buffer) {

		uint *lo = sa + bStart;
		int end = isTypeS ? bLen - leaders : 0;
		while (leaders > 0) {
			int start = end;
			int s = 0;
			for (int j = 0; j < leaders; ++j) {
				if (j == 0 || buffer[j] != buffer[j - 1]) {
					bucketFlag[bStart + end] = bucketStartFlag;
				}

				uint i = suffixBuff[j];
				uint suf = lo[i];
				dest[end] = suf;
				++end;

				if (i + 1 < bLen && suf - lo[i + 1] == periodLength) {
					suffixBuff[s] = i + 1;
					buffer[s] = buffer[j];
					s++;
				}
			}
			leaders = s;

			if (isTypeS)
				end = start - leaders;
		}

		return isTypeS ? bLen - end : end;
	}

	void handleBucketRepeatsOnly(uint bStart, uint bLen, uint D, uint D2,
			uchar bucketStartFlag) {
		detectAndTreatPeriods(bStart, bLen, D, D2, bucketStartFlag);
	}

	void handleBucketNoRepeats(uint bStart, uint bLen, uint D, uint D2,
			uchar bucketStartFlag) {
		bucketBuffVec.reserve(bLen);
		bword *buffer = &bucketBuffVec[0];
		uint *lo = sa + bStart;
		copyBucketNumbers(lo, D, D2, buffer, bLen);
		sorter.sort(bLen, lo, buffer);
		updateBucketStart(buffer, bLen, bStart, bucketStartFlag);
		assignSubBucketNumbersAndLen(bStart, bLen);
	}

	void handleBucketWithRepeats(uint bStart, uint bLen, uint D, uint D2,
			uchar bucketStartFlag) {
		if (!detectAndTreatPeriods(bStart, bLen, D, D2, bucketStartFlag))
			handleBucketNoRepeats(bStart, bLen, D, D2, bucketStartFlag);
	}

	bool adjacentSuffixes(uint p1, uint p2, uint n) {
		uint *a = sa + p1;
		uint *b = sa + p2;
		for (; n--;)
			if (*a++ != *b++ - 1)
				return false;
		return true;
	}

	void copyAdjacentSuffs(uint d, uint s, uint n) {
		uint *a = sa + d;
		uint *b = sa + s;
		for (; n--;)
			*a++ = *b++ - 1;
	}

	void copyFlags(uint d, uint s, uint n, uint f) {
		uchar *a = bucketFlag + d;
		uchar *b = bucketFlag + s;
		for (; n--; ++a)
			if (*b++)
				*a = f;
	}

	uint depthForFlag(uint flag) {
		return ((1 + BWORD64FLAG) * (flag - 1) + 1) * expectedSorted;
	}

	void finalTouches() {
		assignSubBucketNumbersAndLen(0, length);

		const int allowedRepeats = 128;
//		int totalAccess = 2 * length;
		do {
			bool done = true;
			uint prevBLen = 0;
			uint prevBStart = 0;
			uint nextBStart = 0;
			uint bStart = 0;
			uint bLen = 0;
			bool copyFromNext = false;
			uint copyDepth = expectedSorted;

			for (int i = length - 1; i >= 0; --i) {

				bool prevShouldCopy = false;
				{
					uint buck = bucket[i];
					prevBLen = buck & bucketPiggyBackMask;
					prevBStart = (buck >> bucketPiggyBackBits) - 1;
					if (prevBLen != 1) {
						if (prevBStart != bStart) {
							if (!prevBLen)
								prevBLen = getBucketLength(prevBStart);
							prevShouldCopy = prevBLen == bLen
									&& adjacentSuffixes(prevBStart, bStart,
											bLen);
						}
					}
				}

				if (bLen > 1) {
//					totalAccess += bLen;
					int flag = bucketFlag[bStart];
					if (flag <= allowedRepeats) {
						if (copyFromNext) {
							copyDepth += 1;
							bool treatPeriods = copyDepth >= sa[bStart]
									- sa[bStart + 1];

							copyAdjacentSuffs(bStart, nextBStart, bLen);
							copyFlags(bStart, nextBStart, bLen, flag + 1);
							assignSubBucketNumbersAndLen(bStart, bLen);

							if (treatPeriods) {
								uint subLen = 1;
								for (int i = bStart + bLen - 1; i >= bStart; --i) {
									if (bucketFlag[i]) {
										if (subLen > 1) {
											handleBucketRepeatsOnly(i, subLen,
													copyDepth, expectedSorted,
													flag + 1);
										}
										subLen = 1;
									} else
										++subLen;
								}
							}
						} else {
							uint D = depthForFlag(flag);
							handleBucketWithRepeats(bStart, bLen, D,
									expectedSorted, flag + 1);
							copyDepth = depthForFlag(flag + 1);
						}
					} else
						done = false;
				}

				nextBStart = bStart;
				if (bStart == prevBStart) {
					uint buck = bucket[i];
					prevBLen = buck & bucketPiggyBackMask;
					prevBStart = (buck >> bucketPiggyBackBits) - 1;
					if (!prevBLen)
						prevBLen = getBucketLength(prevBStart);
				}
				bLen = prevBLen;
				bStart = prevBStart;
				copyFromNext = prevShouldCopy;
			}

			if (done)
				break;

			for (uint i = 0; i < length; ++i) {
				uchar f = bucketFlag[i];
				if (f > 1)
					bucketFlag[i] = 1;
			}
			expectedSorted *= (2 + BWORD64FLAG);

		} while (1);
//		cout << "Avg access per suffix " << (totalAccess * 1.0 / length)
//				<< endl;
	}

	int inputBasedSort(uint bitsPerFirstPass) {
		vector<word> inputVec(2 + (length * bitsPerChar / bitsPerWord));
		input = &inputVec[0];
		prepareInput(input);

		const int n = 1 << bitsPerFirstPass;
		vector<uint> bSizeVec(n);
		vector<uint> bStartVec(n);
		firstPass(bitsPerFirstPass, n, &bStartVec[0], &bSizeVec[0]);

		uint bits =
				secondPass(bitsPerFirstPass, n, &bStartVec[0], &bSizeVec[0]);
		return bits;
	}

public:
	Radix(uchar *input, uint n) :
		originalInput(input), length(n), bucketPiggyBackBits(uintBits
				- bitsFor(length)), bucketPiggyBackLimit(1
				<< bucketPiggyBackBits), bucketPiggyBackMask(
				bucketPiggyBackLimit - 1) {
	}

	uint* build() {
		init();

//		uint alphaSize =
				renameAlphabet(originalInput, length, charCode,
				bitsPerChar);

		const int bitsPerFirstPass = (16 / bitsPerChar) * bitsPerChar;

		sa = new uint[length];

		int minBitsSortedInSecondPass = inputBasedSort(bitsPerFirstPass);
		expectedSorted = (bitsPerFirstPass + minBitsSortedInSecondPass)
				/ bitsPerChar;

//		printTime("First two passes");
//		printBucketStats();

		uint sentinels = 100; // add sentinels to avoid if's in getBucketWord
		bucketVec.reserve(length + sentinels);
		bucket = &bucketVec[0];
		memset(bucket + length, 0, sentinels * sizeof(uint));

		finalTouches();

		return sa;
	}

private:
	/* To avoid using a special character for the terminal $, we
	 * solve any suffixes of the form 00...0$ first
	 */
	uint fixEndingZeros() {
		uint p = 0;
		for (uint i = length - 1; i && charCode[originalInput[i]] == 0; --i)
			bucketFlag[p++] = 1;
		return p;
	}

	uint getBucketLength(int bStart) {
		uint j = bStart + bucketPiggyBackLimit;
		while (!bucketFlag[j])
			++j;
		return j - bStart;
	}

#ifdef DEBUG
	void printBucketStats() {
		int n = 32;
		int bb[n];
		for (int i = 0; i < n; ++i)
		bb[i] = 0;

		int b = 0;
		int mb = 0;
		int singletons = 0;
		for (uint i = 0; i <= length; ++i) {
			if (bucketFlag[i]) {
				if (b > mb) {
					mb = b;
				}
				//				/*
				for (int k = 0x10000, j = 0; k > 0; k >>= 1, ++j)
				if (b > k) {
					bb[j]++;
				}
				//				 */
				if (b == 1) {
					singletons++;
				}

				b = 1;
			} else {
				++b;
			}
		}
		printf("LARGEST bucket %d\n", mb);
		//		/*
		for (int i = 0x10000, j = 0; i > 0; i >>= 1, ++j)
		if (bb[j] > 0)
		printf("buckets larger than %d : %d\n", i, bb[j]);
		//		 */
		printf("singletons %d\n", singletons);
	}

	uint largestBuck() {
		uint b = 0;
		uint mb = 0;
		for (uint i = 0; i <= length; ++i)
		if (bucketFlag[i]) {
			if (b > mb)
			mb = b;
			b = 1;
		} else
		++b;
		return mb;
	}

	void printSuffix(int s, int len) {
		printf(" s=%d len=%d [", s, len);
		printCharData(originalInput + s, len);
		printf("]\n");
	}

	void printSuffixTranslated(int s, int len) {
		printf("Translated s=%d len=%d: ", s, len);
		for (int i = 0; i < len; ++i)
		printf("%3d ", charCode[originalInput[s + i]]);
		printf("\n");
	}

	void printLargestBucket(uint expectedSorted) {
		uint lb = largestBuck();
		printf("Largest %d\n", lb);
		uint bLen = 0;
		int bStart = 0;
		if (expectedSorted > 1000) {
			expectedSorted = 1000;
		}
		for (uint i = 0; i <= length; ++i) {
			if (bucketFlag[i]) {
				if (bLen > 1 && lb == bLen) {
					uint *lo = sa + bStart;
					uint s = *lo;
					printf("Bucket of size %d starts with", bLen);
					printSuffix(s, expectedSorted);
					break;
				}
				bLen = 1;
				bStart = i;
			} else {
				++bLen;
			}
		}
	}

	bool checkIncreasingSuffixes(int bStart, int bLen, uint expectedSorted) {
		for (int i = bStart + 1; i < bStart + bLen; ++i) {
			int j;
			for (j = 0; j < expectedSorted; ++j) {
				if (sa[i] + j >= length || sa[i - 1] + j >= length) {
					break;
				} else {
					if (originalInput[sa[i] + j]
							!= originalInput[sa[i - 1] + j])
					break;
				}
			}
			if (j < expectedSorted && originalInput[sa[i] + j]
					< originalInput[sa[i - 1] + j]
					//				&& originalInput[mySA[i]] != originalInput[mySA[i] + 1]
			) {
				printf(
						"OHO suffix %d before suffix %d but they differ at %d\n",
						sa[i - 1], sa[i], j);
				printf("OHO suffix %d: ", sa[i - 1]);
				printSuffix(sa[i - 1], expectedSorted);
				printSuffixTranslated(sa[i - 1], expectedSorted);
				printf("OHO suffix %d: ", sa[i]);
				printSuffix(sa[i], expectedSorted);
				printSuffixTranslated(sa[i], expectedSorted);

				return false;
				break;
			}
		}
		//
		return true;
	}

	void checkSortedBucket(uint bStart, uint bLen, uint expectedSorted) {
		uint *lo = sa + bStart;
		for (int i = 1; i < bLen; ++i)
		if (memcmp(originalInput + lo[i - 1], originalInput + lo[i],
						expectedSorted) != 0) {
			printf("These should not be in the same bucket\n");
			printSuffix(lo[i - 1], expectedSorted);
			printSuffix(lo[i], expectedSorted);
		}

	}

#endif

	/*	/////////////////////

	 template<bool forwardDirection>
	 int Radix::populateWithDepth(BinObj& b, uint leaders, uint periodLength,
	 uint *leaderIndex, word *leaderWords, uint *leaderDepth,
	 uint *suffixDest) {

	 int periodDepth = b.D;
	 const uint *lo = sa + b.start;
	 int end = forwardDirection ? 0 : b.len - leaders;
	 while (leaders > 0) {
	 int start = end;
	 end = populateOneRound(end, periodDepth, b, leaders, leaderIndex,
	 leaderWords, leaderDepth, suffixDest);

	 int s = 0;
	 for (uint j = 0; j < leaders; ++j) {
	 uint i = leaderIndex[j];
	 if (i + 1 < b.len && lo[i] - lo[i + 1] == periodLength) {
	 leaderIndex[s] = i + 1;
	 leaderWords[s] = leaderWords[j];
	 leaderDepth[s] = leaderDepth[j];
	 s++;
	 }
	 }
	 leaders = s;

	 periodDepth += periodLength;

	 if (!forwardDirection)
	 end = start - leaders;
	 }

	 return forwardDirection ? end : b.len - end;
	 }

	 answer handleBucketRepeatsOnly2(BinObj& b, WordProvider& wp,
	 uchar bucketStartFlag) {

	 answer a = detectAndTreatPeriods(b, bucketStartFlag, wp);
	 if (a != NotHandledPeriods)
	 assignBucketNumbers(b.start, b.len);
	 return a;
	 }

	 answer handleBucketWithRepeats2(BinObj& b, WordProvider& wp,
	 uchar bucketStartFlag) {

	 answer a = detectAndTreatPeriods2(b, bucketStartFlag, wp);
	 if (a == NotHandledPeriods) {

	 wordBuffer.reserve(b.len);
	 wp.copyWords(b, &wordBuffer[0]);

	 uint *lo = sa + b.start;
	 sorter.sort(b.len, lo, &wordBuffer[0]);
	 updateBucketStart(&wordBuffer[0], b.len, bucketFlag, b.start,
	 bucketStartFlag);
	 }
	 int depthFlag = bucketFlagForDepth(b.D);
	 ForNonSingleSubBucket(b.start, b.len, {
	 updateDepth(_subStart, b.D, depthFlag);
	 //				checkSortedBucket(_subStart, _subLen,
	 //						bucketDepthForFlag(getBucketDepth(_subStart)));
	 });
	 assignBucketNumbers(b.start, b.len);

	 return a;
	 }
	 */
	/*

	 uint populateOneRound(uint end, uint bStart, uint bLen,
	 uchar bucketStartFlag, uint leaders, uint *leaderIndex,
	 word *leaderWords, uint *suffixDest) {
	 const uint *lo = sa + bStart;
	 //		uint currentBucket;
	 for (uint j = 0; j < leaders; ++j) {
	 if (j == 0 || leaderWords[j] != leaderWords[j - 1]) {
	 bucketFlag[bStart + end] = bucketStartFlag;
	 //				currentBucket = bStart + end + 1;
	 }

	 uint i = leaderIndex[j];
	 uint s = lo[i];
	 suffixDest[end] = s;
	 //			bucket[s] = currentBucket;
	 ++end;
	 }
	 return end;
	 }

	 template<bool forwardDirection>
	 int populateWithDepth(uint leaders, uint periodLength, uint bStart,
	 uint bLen, uchar bucketStartFlag, uint *leaderIndex,
	 word *leaderWords, uint *suffixDest) {

	 const uint *lo = sa + bStart;
	 int end = forwardDirection ? 0 : bLen - leaders;
	 while (leaders > 0) {
	 int start = end;
	 end = populateOneRound(end, bStart, bLen, bucketStartFlag, leaders,
	 leaderIndex, leaderWords, suffixDest);

	 int s = 0;
	 for (uint j = 0; j < leaders; ++j) {
	 uint i = leaderIndex[j];
	 if (i + 1 < bLen && lo[i] - lo[i + 1] == periodLength) {
	 leaderIndex[s] = i + 1;
	 leaderWords[s] = leaderWords[j];
	 s++;
	 }
	 }
	 leaders = s;

	 if (!forwardDirection)
	 end = start - leaders;
	 }

	 return forwardDirection ? end : bLen - end;
	 }
	 answer treatPeriodLeaders3(uint bStart, uint bLen, uint D, uint D2,
	 uchar bucketStartFlag, uint period) {

	 vector<word>& leaderWordsVec = bufferVec;
	 leaderWordsVec.clear();
	 vector<uint>& leaderIndex = uintBuff1;
	 leaderIndex.clear();

	 const uint *lo = sa + bStart;
	 const uint periodThreshold = D;
	 uint prevSuff = length + D + 1;
	 uint periodFound = 0;
	 uint periodLength = 0;
	 for (uint i = 0; i < bLen; ++i) {
	 int l = lo[i];
	 bool shouldAdd = true;
	 uint potentialPeriod = prevSuff - l;
	 if (potentialPeriod <= periodThreshold) {
	 if (!periodFound) {
	 periodFound = i;
	 periodLength = potentialPeriod;
	 } else
	 shouldAdd = (potentialPeriod != periodLength);
	 }

	 if (shouldAdd) {
	 word next = getBucketWord(l + D, D2);
	 leaderWordsVec.push_back(next);
	 leaderIndex.push_back(i);
	 }
	 prevSuff = l;
	 }

	 uint nLeaders = leaderWordsVec.size();

	 word *leaderWords = &leaderWordsVec[0];
	 sorter.sort(nLeaders, &leaderIndex[0], leaderWords);


	 uintBuff2.reserve(bLen);
	 uint *suffixDest = &uintBuff2[0];

	 if (nLeaders < bLen) {
	 uint i = 0;
	 for (; leaderIndex[i] != periodFound; ++i)
	 assert(i < nLeaders);
	 if ((i > 0 && leaderWords[i - 1] == leaderWords[i]) || (i + 1
	 < bLen && leaderWords[i] == leaderWords[i + 1])) {
	 //				cout << "Horror blen " << bLen << " leaders " << nLeaders
	 //						<< endl;
	 return NotHandledPeriods;
	 }

	 uint populated = populateWithDepth<true> (i, periodLength, bStart,
	 bLen, bucketStartFlag, &leaderIndex[0], leaderWords,
	 suffixDest);

	 populated += populateWithDepth<false> (nLeaders - (i + 1),
	 periodLength, bStart, bLen, bucketStartFlag, &leaderIndex[i
	 + 1], leaderWords + i + 1, suffixDest);

	 assert(populated == bLen);
	 } else {
	 assert(nLeaders == bLen);
	 uint populated = populateOneRound(0, bStart, bLen, bucketStartFlag,
	 nLeaders, &leaderIndex[0], leaderWords, suffixDest);
	 assert(populated == bLen);
	 }
	 memcpy(sa + bStart, suffixDest, bLen * sizeof(uint));
	 assignBucketNumbers(bStart, bLen);
	 return HandledPeriodsFullFit;
	 }
	 */
	/*
	 uint getBucketLen(uint bStart) {
	 uint code = bucketLen[bStart + 1];
	 if (!(code & 128))
	 return 1;

	 uint bLen = (code ^ 128);
	 bStart++;
	 for (uint bits = 7; (code = bucketLen[++bStart]) & 128; bits += 7)
	 bLen = ((code ^ 128) << bits) | bLen;
	 return bLen;
	 }

	 void cleanupLength(uint bStart) {
	 while (bucketLen[++bStart] & 128)
	 bucketLen[bStart] = 0;
	 }

	 void setBucketLen(uint bStart, uint bLen) {
	 if (bLen > 1) {
	 uint oldBLen = bLen;
	 uint oldBStart = bStart;
	 while (bLen) {
	 bucketLen[++bStart] = (uchar) ((bLen & 127) | 128);
	 bLen >>= 7;
	 }
	 if (getBucketLen(oldBStart) != oldBLen) {
	 cout << "Hallou" << endl;
	 uint a = getBucketLen(oldBStart);
	 setBucketLen(oldBStart, oldBLen);
	 uint b = getBucketLen(oldBStart);
	 }
	 }
	 }*/
	/*
	 void assignBucketNumbers(uint bStart, uint bLen) {
	 int currentBStart = bStart;
	 bucket[sa[bStart]] = currentBStart + 1;
	 for (uint e = bStart + bLen; ++bStart < e;) {
	 if (bucketFlag[bStart])
	 currentBStart = bStart;
	 bucket[sa[bStart]] = currentBStart + 1;
	 }
	 }*/

};
#endif /* RADIX_H_ */
