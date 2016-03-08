# RadixSA - Suffix array construction algorithm

As described in the paper:

Rajasekaran, Sanguthevar, and Marius Nicolae. "An elegant algorithm for the construction of suffix arrays." Journal of Discrete Algorithms 27 (2014): 21-28.

http://www.sciencedirect.com/science/article/pii/S1570866714000173

**NOTE:** This algorithm uses 32 bit integers for the suffix array, therefore the maximum length of the string should fit within 32 bits. **We recommend everyone to get [RadixSA64](https://github.com/mariusmni/radixSA64) instead**, which uses configurable integer type: 32 bits for shorter strings, 64 bits for longer strings.


## To compile:

Make sure you have ```g++``` and ```make``` installed. 
Open a command line and type: 

```
make
```

## To run:

```
./radixSA inputFile outputFile
```

Please report any bugs to mariusmni@gmail.com




