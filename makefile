all:
	make -C Debug/
	mv Debug/radixSA .

clean:
	make -C Debug/ clean
	rm -f radixSA Debug/radixSA
