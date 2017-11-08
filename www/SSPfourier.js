// ****************************************************************
// Fourier processing procedures to be used in the SSP package
// i.t. young
// Wednesday, 21 September 2016


// ****************************************************************
// From Mathematica: f2D1A = Abs[f2D1]; MagOnlyImage = Chop[InverseFourier[f2D1A]]
// checked Sunday, 18 December 2016

function magOnlySignal(data, nn, ndim)
	{
		var theData = data.slice();
		var temp = fft(theData, nn, ndim, FORWARD);
		var magOnly = abs(temp);
		temp = fft(magOnly, nn, ndim, BACKWARD);
		return chop(temp);
	};


function testMagOnly(testChoice)
	{
		if(testChoice)
		{
		console.log("\nStart testMagOnly");
		
		// prepare unit impulse signal
		var myLength = 6;
		var data = createArray(2*myLength);
		for(var i = 0; i < 2*myLength; i += 2)
			{
			data[i] = 0; data[i+1] = 0;
			};
		data[0] = 1;
		
		var nn = [myLength];
		var ndim = dimensions(nn);
/*
		isign = FORWARD;
	
		var spectrum = fft(data, nn, ndim, isign);
		console.log("Input to & output from fft: ",data,spectrum);
		
		isign = BACKWARD;
		var inverse = fft(spectrum, nn, ndim, isign);
		console.log("output from inverse transform: ",spectrum,inverse);
*/
		var tlist1 = magOnlySignal(data,nn,ndim);
		console.log("mag-only input & output A: ",data,round(tlist1,rounding));

		var places = 3;
		var rotData = rotateRight(data,2*places);
		var tlist2 = magOnlySignal(rotData,nn,ndim);
		console.log("mag-only input & output B: ",rotData,round(tlist2,rounding));

		console.log("Stop testMagOnly");
		console.log("");
		};
	};


// ****************************************************************
// For 1D & 2D signals
// From Mathematica, PhaseOnlyImage: f2D1P = f2D1/Abs[f2D1] -->
// PhaseOnlyImage = Chop[InverseFourier[f2D1P]]
// checked Sunday, 18 December 2016

function phaseOnlySignal(data, nn, ndim)
	{
		var n = data.length;
		var phaseOnly = createArray(n);
		var theData = data.slice();
		var temp = fft(theData, nn, ndim, FORWARD);
		var mag = abs(temp);
		for(var i = 0; i < n; i += 2)
			{
			if(mag[i] === 0)
				{
				phaseOnly[i] = 0;
				phaseOnly[i+1] = 0;
				}
			else
				{
				phaseOnly[i] = temp[i]/mag[i];
				phaseOnly[i+1] = temp[i+1]/mag[i];
				}
			};
		temp = fft(phaseOnly, nn, ndim, BACKWARD);
		return chop(temp);
	};


function testPhaseOnly(testChoice)
	{
		if(testChoice)
		{
		console.log("\nStart testPhaseOnly");
	
		// prepare input test signal
		var myLength = 9;
		var data = createArray(2*myLength);
		for(var i = 0; i < 2*myLength; i += 2)
			{
			data[i] = 0; data[i+1] = 0;
			};
		data[0] = -2;
//		var real = [-2,0,0,0,0,0];
//		data[i] = 4 - Math.abs(4 - i/2);
//		var real = nlist1;
		
		var nn = [myLength];
		var ndim = dimensions(nn);

		console.log("Input A: ",data);
		var tlist1 = phaseOnlySignal(data,nn,ndim);
		console.log("Phase-only input & output A: ",data,round(tlist1,rounding));

		var places = 3;
		var rotData = rotateRight(data,2*places);
		var tlist2 = phaseOnlySignal(rotData,nn,ndim);
		console.log("Phase-only input & output B: ",rotData,round(tlist2,rounding));

		console.log("Stop testPhaseOnly");
		console.log("");
		};
	};


// ****************************************************************
// Copied and slightly modified from fft.js
// define() has been removed as has the complicated return statement

/* 
	* Free FFT and convolution (JavaScript)
 * 
 * Copyright (c) 2014 Project Nayuki
 * http://www.nayuki.io/page/free-small-fft-in-multiple-languages
 * 
 * (MIT License)
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 * the Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 * - The above copyright notice and this permission notice shall be included in
 *   all copies or substantial portions of the Software.
 * - The Software is provided "as is", without warranty of any kind, express or
 *   implied, including but not limited to the warranties of merchantability,
 *   fitness for a particular purpose and noninfringement. In no event shall the
 *   authors or copyright holders be liable for any claim, damages or other
 *   liability, whether in an action of contract, tort or otherwise, arising from,
 *   out of or in connection with the Software or the use or other dealings in the
 *   Software.
 */
// checked Sunday, 18 December 2016

	var cosTable = [];
	var sinTable = [];
	var lutInitialised = undefined;

	"use strict";

	function initLut(n) {
		cosTable = new Array(n / 2);
		sinTable = new Array(n / 2);
		for (var i = 0; i < n / 2; i++) {
			cosTable[i] = Math.cos(2 * Math.PI * i / n);
			sinTable[i] = Math.sin(2 * Math.PI * i / n);
		}
		lutInitialised = n;
	}


	/* 
	 * Computes the discrete Fourier transform (DFT) of the given complex vector,
	 * storing the result back into the vector.
	 * The vector can have any length. This is a wrapper function.
	 */
	function transform(real, imag) {
		if (real.length != imag.length)
			throw "Mismatched lengths";

		var n = real.length;
		if (n == 0)
			return;
		else if ((n & (n - 1)) == 0)  // Is power of 2
			transformRadix2(real, imag);
		else  // More complicated algorithm for arbitrary sizes
			transformBluestein(real, imag);
	}


	/* 
	 * Computes the inverse discrete Fourier transform (IDFT) of the given complex
	 * vector, storing the result back into the vector.
	 * The vector can have any length. This is a wrapper function. This
	 * transform does not perform scaling, so the inverse is not a true inverse.
	 */

	function inverseTransform(real, imag)
		{
		var ntot = real.length;
		transform(imag, real); 
		}

	/* 
	 * Computes the discrete Fourier transform (DFT) of the given complex vector,
	 * storing the result back into the vector.
	 * The vector's length must be a power of 2. Uses the Cooley-Tukey
	 * decimation-in-time radix-2 algorithm.
	 */
	function transformRadix2(real, imag) {
		// Initialization
		if (real.length != imag.length)
			throw "Mismatched lengths";
		var n = real.length;
		if (n == 1)  // Trivial transform
			return;
		var levels = -1;
		for (var i = 0; i < 32; i++) {
			if (1 << i == n)
				levels = i;  // Equal to log2(n)
		}
		if (n !== lutInitialised) { // has the LUT been initialised with the right samplesize?
			if (levels == -1)
				throw "Length is not a power of 2";
			initLut(n);
		}

		// Bit-reversed addressing permutation
		for (var i = 0; i < n; i++) {
			var j = reverseBits(i, levels);
			if (j > i) {
				var temp = real[i];
				real[i] = real[j];
				real[j] = temp;
				temp = imag[i];
				imag[i] = imag[j];
				imag[j] = temp;
			}
		}

		// Cooley-Tukey decimation-in-time radix-2 FFT
		for (var size = 2; size <= n; size *= 2) {
			var halfsize = size / 2;
			var tablestep = n / size;
			for (var i = 0; i < n; i += size) {
				for (var j = i, k = 0; j < i + halfsize; j++, k += tablestep) {
					var tpre = real[j + halfsize] * cosTable[k] + imag[j + halfsize] * sinTable[k];
					var tpim = -real[j + halfsize] * sinTable[k] + imag[j + halfsize] * cosTable[k];
					real[j + halfsize] = real[j] - tpre;
					imag[j + halfsize] = imag[j] - tpim;
					real[j] += tpre;
					imag[j] += tpim;	
				}
			}
		}

		// Returns the integer whose value is the reverse of the lowest 'bits' bits of the integer 'x'.
		function reverseBits(x, bits) {
			var y = 0;
			for (var i = 0; i < bits; i++) {
				y = (y << 1) | (x & 1);
				x >>>= 1;
			}
			return y;
		}
	}


	/* 
	 * Computes the discrete Fourier transform (DFT) of the given complex vector,
	 * storing the result back into the vector.
	 * The vector can have any length. This requires the convolution function,
	 * which in turn requires the radix-2 FFT function.
	 * Uses Bluestein's chirp z-transform algorithm.
	 */
	function transformBluestein(real, imag) {
		// Find a power-of-2 convolution length m such that m >= n * 2 + 1
		if (real.length != imag.length)
			throw "Mismatched lengths";
		var n = real.length;
		var m = 1;
		while (m < n * 2 + 1)
			m *= 2;

		// Trignometric tables
		var cosTable = new Array(n);
		var sinTable = new Array(n);
		for (var i = 0; i < n; i++) {
			var j = i * i % (n * 2);  // This is more accurate than j = i * i
			cosTable[i] = Math.cos(Math.PI * j / n);
			sinTable[i] = Math.sin(Math.PI * j / n);
		}

		// Temporary vectors and preprocessing
		var areal = new Array(m);
		var aimag = new Array(m);
		for (var i = 0; i < n; i++) {
			areal[i] = real[i] * cosTable[i] + imag[i] * sinTable[i];
			aimag[i] = -real[i] * sinTable[i] + imag[i] * cosTable[i];
		}
		for (var i = n; i < m; i++)
			areal[i] = aimag[i] = 0;
		var breal = new Array(m);
		var bimag = new Array(m);
		breal[0] = cosTable[0];
		bimag[0] = sinTable[0];
		for (var i = 1; i < n; i++) {
			breal[i] = breal[m - i] = cosTable[i];
			bimag[i] = bimag[m - i] = sinTable[i];
		}
		for (var i = n; i <= m - n; i++)
			breal[i] = bimag[i] = 0;

		// Convolution
		var creal = new Array(m);
		var cimag = new Array(m);
		convolveComplex(areal, aimag, breal, bimag, creal, cimag);

		// Postprocessing
		for (var i = 0; i < n; i++) {
			real[i] = creal[i] * cosTable[i] + cimag[i] * sinTable[i];
			imag[i] = -creal[i] * sinTable[i] + cimag[i] * cosTable[i];
		}
	}


	/* 
	 * Computes the circular convolution of the given real vectors. Each vector's length must be the same.
	 */
	function convolveReal(x, y, out) {
		if (x.length != y.length || x.length != out.length)
			throw "Mismatched lengths";
		var zeros = new Array(x.length);
		for (var i = 0; i < zeros.length; i++)
			zeros[i] = 0;
		convolveComplex(x, zeros, y, zeros.slice(0), out, zeros.slice(0));
	}


	/* 
	 * Computes the circular convolution of the given complex vectors.
	 * Each vector's length must be the same.
	 */
	function convolveComplex(xreal, ximag, yreal, yimag, outreal, outimag) {
		if (xreal.length != ximag.length || xreal.length != yreal.length || yreal.length != yimag.length || xreal.length != outreal.length || outreal.length != outimag.length)
			throw "Mismatched lengths";

		var n = xreal.length;
		xreal = xreal.slice(0);
		ximag = ximag.slice(0);
		yreal = yreal.slice(0);
		yimag = yimag.slice(0);

		transform(xreal, ximag);
		transform(yreal, yimag);
		for (var i = 0; i < n; i++) {
			var temp = xreal[i] * yreal[i] - ximag[i] * yimag[i];
			ximag[i] = ximag[i] * yreal[i] + xreal[i] * yimag[i];
			xreal[i] = temp;
		}
		inverseTransform(xreal, ximag);

		for (var i = 0; i < n; i++)
		{  // Scaling (because this FFT implementation omits it)
			outreal[i] = xreal[i] / n;
			outimag[i] = ximag[i] / n;
		}

	}


// ****************************************************************
// Following is based upon Numerical Recipe's's procedure fourn.c (2nd Edition).
// checked Sunday, 18 December 2016

function fourND(data, nn, ndim, isign)
// Replaces data by its ndim-dimensional discrete Fourier transform
// if isign is input as 1 (FORWARD). nn[1..ndim] is an integer array containing
// the lengths of each dimension (number of complex values), which MUST all be
// powers of 2.

// Data is a real array of length twice the product of these lengths, in which the data
// are stored as in a multidimensional complex array: real and imaginary parts of each
// element are in consecutive locations, and the rightmost index of the array increases
// most rapidly as one proceeds along data.

// For a two-dimensional array, this is equivalent to storing the array by rows.
// If isign is input as −1 (BACKWARD), data is replaced by its inverse transform times the
// product of the lengths of all dimensions.
	{
	// Compute total number of complex values
	var ntot = 1;
	var idim = 1;
	for (idim = 0; idim < ndim; idim++)
		ntot *= nn[idim];
	
	var nprev = 1;
	var n = 1;
	var nrem = 1;
	var i1 = 1; var i2 = 1; var i3 = 1;
	var ip1 = 1; var ip2 = 1; var ip3 = 1;
	var i2rev = 1; var i3rev = 1;
	var ibit = 1;
	var ifp1 = 1; var ifp2 = 1;
	var theta = 1;
	var wr = 1; var wi = 1;
	var wpr = 1; var wpi = 1;
	var wtemp = 1;
	var tempr = 1; var tempi = 1;
	var k1 = 1; var k2 = 1;
	
	// Main loop over the dimensions.
	for (idim = ndim; idim >= 1; idim--)
		{
		n = nn[idim - 1];
		nrem = ntot / (n * nprev);
		ip1 = nprev << 1;
		ip2 = ip1 * n;
		ip3 = ip2 * nrem;
		i2rev = 1;
		for (i2 = 1; i2 <= ip2; i2 += ip1)
			// This is the bit-reversal section of the routine.
			{
			if (i2 < i2rev)
				{
				for (i1 = i2; i1 <= i2 + ip1 - 2; i1 += 2)
					{
					for (i3 = i1; i3 <= ip3; i3 += ip2)
						{
						i3rev = i2rev + i3 - i2;
						// inline swap routine
						temp = data[i3 - 1];
						data[i3 - 1] = data[i3rev - 1];
						data[i3rev - 1] = temp;
						// inline swap routine
						temp = data[i3];
						data[i3] = data[i3rev];
						data[i3rev] = temp;
						}
					}
				}
			ibit = ip2 >>> 1;
			while (ibit >= ip1 && i2rev > ibit)
				{
				i2rev -= ibit;
				ibit >>>= 1;
				}
			i2rev += ibit;
			}
			
		// Here begins the Danielson - Lanczos section of the routine
		ifp1 = ip1;
		while (ifp1 < ip2)
			{
			ifp2 = ifp1 << 1;
			// Initialize for the trig. recurrence.
			theta = isign * 2 * Math.PI / (ifp2 / ip1); 
			wtemp = Math.sin(0.5 * theta);
			wpr = -2.0 * wtemp * wtemp;
			wpi = Math.sin(theta);
			wr = 1.0;
			wi = 0.0;
			for (i3 = 1; i3 <= ifp1; i3 += ip1)
				{
				for (i1 = i3; i1 <= i3 + ip1 - 2; i1 += 2)
					{
					for (i2 = i1; i2 <= ip3; i2 += ifp2)
						{
						// Danielson - Lanczos formula:
						k1 = i2;
						k2 = k1 + ifp1;
						tempr = wr * data[k2-1] - wi * data[k2];
						tempi = wr * data[k2] + wi * data[k2-1];
						data[k2-1] = data[k1-1] - tempr;
						data[k2] = data[k1] - tempi;
						data[k1-1] += tempr;
						data[k1] += tempi;
						}
					}
				// Trigonometric recurrence.
				wtemp = wr;
				wr = wr * wpr - wi * wpi + wr;
				wi = wi * wpr + wtemp * wpi + wi;
				}
			ifp1 = ifp2;
			}
		nprev *= n;
		};
	// following is NOT in Numerical Recipes but it should be, right?
	if (isign === BACKWARD)
		for(var i = 0; i < data.length; i++) data[i] /= ntot;
	return data};

function testFourND(testChoice)
	{
		if(testChoice)
		{
		console.log("\nStart testFourND");
		
		var myLength = 16; // number of complex samples
		var real = createArray(myLength);
		var imag = createArray(myLength);
		var data1D = createArray(2*myLength);
		
		for(var i = 0; i< myLength; i++)
			{
			real[i] = 0; imag[i] = 0;
			};
		real[0] = 1;
		
		var rcopy = real.slice();
		var icopy = imag.slice();
		var myData = [rcopy,icopy];
	
		var spectrum = myFourier(myData);
		console.log("myFourier 1D: ",myData, spectrum, round(spectrum[0],rounding),
			round(spectrum[1],rounding));
		var inverse = myInverseFourier(spectrum);
		console.log("myInverseFourier 1D: ",round(inverse[0],rounding),
			round(inverse[1],rounding));

		var nn = [myLength];
		var data = merge(myData, nn);
		var ndim = 1;
		var isign = FORWARD;
		
		var result = fourND(data,nn,ndim,isign);
		console.log("testFourND 1D forward: myData, result",myData,
			round(result,rounding));
		var isign = BACKWARD;
		var invresult = fourND(data,nn,ndim,isign);
		console.log("testFourND 1D backward: myData, invresult",myData,
			round(invresult,rounding));
			
		for(var i = 0; i< myLength; i++)
			{
			real[i] = Math.cos(i*2*Math.PI / myLength);
			imag[i] = 0;
			};

		var rcopy = real.slice();
		var icopy = imag.slice();
		var myData = [rcopy,icopy];
	
		var spectrum = myFourier(myData);
		console.log("myFourier 1D cos(): ",myData,
			round(spectrum[0],rounding),round(spectrum[1],rounding));
		var inverse = myInverseFourier(spectrum);
		console.log("myInverseFourier 1D cos(): ",round(inverse[0],rounding),
			round(inverse[1],rounding));

		nn = [myLength];
		data = merge(myData, nn);
		ndim = 1;
		isign = FORWARD;
		
		result = fourND(data,nn,ndim,isign);
		console.log("testFourND 1D forward cos(): myData, data",
			myData,round(result,rounding));
			
		var isign = BACKWARD;
		invresult = fourND(result,nn,ndim,isign);
		splitData = split(round(invresult,rounding), [1,myLength])
		console.log("testFourND 1D backward cos(): myData, invresult",
			myData,round(invresult,rounding),splitData);

		var dims = arrayDimension(image1);
		nn = [dims[0]/2, dims[1]];  // correction for complex numbers
		ndim = nn.length;
		var data2D = merge(image1,nn);
		isign = FORWARD;
		
		var spectrum2D = fourND(data2D,nn,ndim,isign);
		var splitSpect2D = split(spectrum2D, nn)
		console.log("testFourND 2D forward impulse: origData, splitSpect2D",
			image1, splitSpect2D);
			
		var isign = BACKWARD;
		var invresult2D = fourND(spectrum2D,nn,ndim,isign);
		var splitInv2D = split(invresult2D, nn)
		console.log("testFourND 2D backward impulse: origData, splitInv2D",
			splitSpect2D,splitInv2D);

		var dataL2D = dims[0] * dims[1] / 2;
		for(var i = 0; i< dataL2D; i++)
			{
			data2D[2*i] = Math.cos(i*4*2*Math.PI / dataL2D);
			data2D[2*i+1] = 0;
			};
		isign = FORWARD;
		var dcopy = data2D.slice();
		
		var spectrum2D = chop(fourND(data2D,nn,ndim,isign));
		var splitSpect2D = split(spectrum2D, nn)
		console.log("testFourND 2D forward impulse: origData, splitSpect2D",
			chop(split(dcopy,nn)), splitSpect2D);


		console.log("Stop testFourND");
		console.log("");
		};
	};


function timeFourND(testChoice)
	{
		if(testChoice)
		{
		console.log("\nStart timeFourND");
		var minPower = 10;
		var nSamples = [];
		for (var a = 0; a < 6; a++)
			nSamples[a] = Math.pow(2,2*a+minPower);
		var time1 = createArray(nSamples.length);
		var time2 = createArray(nSamples.length);
		var iterations = 10;
	
		for(var i = 0; i < nSamples.length; i++)
			{
			time1[i] = 0;
			time2[i] = 0;
			};

// First check the old Fourier routine myFourier(*)

		for(var i = 0; i < nSamples.length; i++)
			{
			var samples = nSamples[i];
			var rnoise = Array(samples);
			var inoise = Array(samples);
			var data1D =createArray(2*samples);
			
			for(var j = 0; j < samples; j++)
				{
				rnoise[j] = Math.random();
				inoise[j] = Math.random();
				};
				
			var noise = [rnoise,inoise];  // complex noise input
			var spectrum = [];
			var start = new Date().getTime();
			
			for(var k = 0; k < iterations; k++) spectrum = myFourier(noise);
			
			var finish = new Date().getTime(); // in ms
			time1[i] = finish - start;
			};
		
		for(var i = 0; i < nSamples.length; i++)
			{
			time1[i] /= iterations;  // time per iteration in ms
			};
			
		console.log("1D - First check the old Fourier routine myFourier(*): ");
		console.log("myFourier Samples & Time [ms]: total 1D data, time/fft, ",
		nSamples,round(time1,rounding));

// Now check the new Fourier routine fourND(*) with shuffling
		var cases = 6;
		var caseSample = createArray(cases);
		var totalLength = createArray(cases);
		var time3 = createArray(cases);
		for(var ii = 0; ii < cases; ii++)
			{
			caseSample[ii] = Math.pow(2,ii+5);
			var samples = caseSample[ii];
			// square array of complex numbers
			var rows = samples;
			var cols = samples;
			totalLength[ii] = rows*cols;
			var noiseLinear = createArray(2*rows*cols);
			
			for (var m = 0; m < noiseLinear.length; m++)
				{
				noiseLinear[m] = randomInteger(0,9);
				};
			
			var noise = createArray(rows,2*cols);
			var nn = [rows,cols];
			var ndim = nn.length;
			noise = split(noiseLinear,nn);
			var isign = FORWARD;
			var spectrum = [];
			var data3D = [];
			
			var start = new Date().getTime();
			for(var k = 0; k < iterations; k++)
				{
				data3D = merge(noise,nn);
				spectrum = fourND(data3D,nn,ndim,isign);
				};
			var finish = new Date().getTime(); // in ms
			time3[ii] = finish - start;
			};
		
		for(var n = 0; n < cases; n++)
			{
			time3[n] /= iterations;  // time per iteration in ms
			};
			
			console.log("2D - Now check the new Fourier routine fourND(*) with merge: ");
			console.log("fourND Samples & Time [ms]: total 2D data, time/fft, rows & cols ",totalLength,round(time3,rounding),caseSample);

		console.log("Stop timeFourND");
		console.log("");
		};

	};

// ****************************************************************
// Following is the general fft routine for 1D and 2D data.
// Both 1D and 2D versions are working Wednesday, 7 December 2016

// Remember 1D can be ANY length but 2D MUST be powers of two lengths in both dimensions.
// A 2D 512 x 512 fft takes about 75 ms (determined by the performance command below).
// checked Sunday, 18 December 2016


// The calling sequence is:

function fft(tempData, nn, ndim, isign)
// Computes n-dim-dimensional discrete Fourier transform starting from "data"
// Does NOT overwrite input "data". (Assumption: Memory is cheap and plentiful.)
// If isign = 1 (FORWARD), this is a transform from time (or space) to frequency. 
// If isign = –1 (BACKWARD), this is a transform from frequency to time (or space).
// The parameter ndim is the number of dimensions , e.g. ndim = 1 or ndim = 2.
// nn[1,2, ..., ndim] is an integer array containing the lengths of each dimension,
// For example:
//     if ndim = 1, then nn = [197] means 197 samples.
//     if ndim = 2, then nn = [128,256] means 128 rows by 256 columns leading
// to 128*256 samples and each sample is assumed to be a COMPLEX value.
// Note that a 1D signal can have ANY length but a 2D signal MUST have rows and
// columns that are each a power of two.

// The array "data" is a REAL array of a length that is twice the product of the
// total signal length. In the examples above, the 1D array has a total length
// of 2*197 and the 2D array has a total length of 2*128*256.
// The real and imaginary parts of each element are in consecutive locations.
// As one indexes along the array one is moving along a row of the data.

// This fft routine automatically determines which Fourier subroutine should be used.
// For 1D data this is the routine "transform" and for 2D data this is "fourND".
	
	{
	// Compute total number of complex values
	var ntot = 1;
	var idim = 1;
	for (idim = 0; idim < ndim; idim++) ntot *= nn[idim];

	if (ndim === 1)
		{
		var splitData = split(tempData,nn);
		if (isign === FORWARD)
			{
			transform(splitData[0], splitData[1]);
			var mergeData = merge(splitData,nn);
			}
		else
			{
			inverseTransform(splitData[0],splitData[1]);
			var mergeData = merge(splitData,nn);
			// Scaling (because this FFT implementation omits it)
			var n = (mergeData.length)/2;
			for (var i = 0; i < mergeData.length; i++)
				mergeData[i] /= n;
			};
		return mergeData;
		}
	else if (ndim === 2)
		{
		var data = tempData.slice();
		var nprev = 1;
		var n = 1;
		var nrem = 1;
		var i1 = 1; var i2 = 1; var i3 = 1;
		var ip1 = 1; var ip2 = 1; var ip3 = 1;
		var i2rev = 1; var i3rev = 1;
		var ibit = 1;
		var ifp1 = 1; var ifp2 = 1;
		var theta = 1;
		var wr = 1; var wi = 1;
		var wpr = 1; var wpi = 1;
		var wtemp = 1;
		var tempr = 1; var tempi = 1;
		var k1 = 1; var k2 = 1;
	
		// Main loop over the dimensions.
		for (idim = ndim; idim >= 1; idim--)
			{
			n = nn[idim - 1];
			nrem = ntot / (n * nprev);
			ip1 = nprev << 1;
			ip2 = ip1 * n;
			ip3 = ip2 * nrem;
			i2rev = 1;
			for (i2 = 1; i2 <= ip2; i2 += ip1)
				// This is the bit-reversal section of the routine.
				{
				if (i2 < i2rev)
					{
					for (i1 = i2; i1 <= i2 + ip1 - 2; i1 += 2)
						{
						for (i3 = i1; i3 <= ip3; i3 += ip2)
							{
							i3rev = i2rev + i3 - i2;
							// inline swap routine
							temp = data[i3 - 1];
							data[i3 - 1] = data[i3rev - 1];
							data[i3rev - 1] = temp;
							// inline swap routine
							temp = data[i3];
							data[i3] = data[i3rev];
							data[i3rev] = temp;
							}
						}
					}
				ibit = ip2 >>> 1;
				while (ibit >= ip1 && i2rev > ibit)
					{
					i2rev -= ibit;
					ibit >>>= 1;
					}
				i2rev += ibit;
				}
			
			// Here begins the Danielson - Lanczos section of the routine
			ifp1 = ip1;
			while (ifp1 < ip2)
				{
				ifp2 = ifp1 << 1;
				// Initialize for the trig. recurrence.
				theta = isign * 2 * Math.PI / (ifp2 / ip1); 
				wtemp = Math.sin(0.5 * theta);
				wpr = -2.0 * wtemp * wtemp;
				wpi = Math.sin(theta);
				wr = 1.0;
				wi = 0.0;
				for (i3 = 1; i3 <= ifp1; i3 += ip1)
					{
					for (i1 = i3; i1 <= i3 + ip1 - 2; i1 += 2)
						{
						for (i2 = i1; i2 <= ip3; i2 += ifp2)
							{
							// Danielson - Lanczos formula:
							k1 = i2;
							k2 = k1 + ifp1;
							tempr = wr * data[k2-1] - wi * data[k2];
							tempi = wr * data[k2] + wi * data[k2-1];
							data[k2-1] = data[k1-1] - tempr;
							data[k2] = data[k1] - tempi;
							data[k1-1] += tempr;
							data[k1] += tempi;
							}
						}
					// Trigonometric recurrence.
					wtemp = wr;
					wr = wr * wpr - wi * wpi + wr;
					wi = wi * wpr + wtemp * wpi + wi;
					}
				ifp1 = ifp2;
				}
			nprev *= n;
			};
		// following is NOT in Numerical Recipes but it should be, right?
		if (isign == BACKWARD)
			for(var i = 0; i < data.length; i++) data[i] /= ntot;
		return data;
		}
	else
		throw "Wrong number of dimensions!";
	};

function testFFT1D(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testFFT1D");
	
		var ndim = 1; // number of dimensions
		var nn = createArray(ndim);
		nn[0] = 16; // number of complex samples
		var myLength = nn[0];
		var testData1D = createArray(2*myLength);
	
		// prepare unit impulse signal
		for(var i = 0; i < 2*myLength; i += 2)
			{
			testData1D[i] = 0; testData1D[i+1] = 0;
			};
		testData1D[0] = 1;
	
		// forward transform
		var isign = FORWARD;
		var spectrum = fft(testData1D, nn, ndim, isign);
		console.log("fft 1D impulse: ",testData1D, round(spectrum,rounding));

		isign = BACKWARD;
		var inverse = fft(spectrum, nn, ndim, isign);
		console.log("inverse fft 1D: ",round(spectrum,rounding), round(inverse,rounding));

		// prepare sinusoidal signal
		var harmonic = 1;
		for(var i = 0; i< 2*myLength; i += 2)
			{
			testData1D[i] = Math.cos(i*harmonic*Math.PI / myLength);
			testData1D[i+1] = 0;
			};

		// forward transform
		isign = FORWARD;
		spectrum = fft(testData1D, nn, ndim, isign);
		console.log("fft 1D cos(): ",round(testData1D,rounding), round(spectrum,rounding));
		
		isign = BACKWARD;
		inverse = fft(spectrum, nn, ndim, isign);
		console.log("inverse fft 1D: ",round(spectrum,rounding), round(inverse,rounding));

		console.log("Stop testFFT1D");
		console.log("");
		};
	};


function testLigteringen(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testLigteringen");
	
		var rdata = [1, 2, 3, 4, 5, 4, 3, 2, 1];  // the real data
		var dLength = rdata.length;  // number of samples, any positive number
		var ndim = dimensions(rdata); // number of dimensions
	
		var idata = Array(dLength);  // the imaginary data
		idata.fill(0);  // fill imaginary data with zeroes
		var nn = [];
		nn[0] = dLength; // number of complex samples
		var data = merge([rdata, idata], nn);  // merge real & imag into universal format
		console.log("fft audio: ndim, nn,rdata, idata, data",ndim, nn,rdata, idata, data);

		// forward transform
		var isign = FORWARD;  // defined in SSPconstants.js, FORWARD = 1
		var spectrum = fft(data, nn, ndim, isign);
		console.log("fft audio: ",data, spectrum);
		
		// absolute value of spectrum to display
		var magnitude = abs(spectrum);
		console.log("fft audio: spectrum, magnitude",spectrum,magnitude);
		
		// normalized absolute value of spectrum to display
		var normed = abs(ampNormalize(spectrum));
		console.log("fft audio: spectrum, normed",spectrum,normed);
		
		// phase of spectrum to display
		var phase = arg(spectrum);
		console.log("fft audio: spectrum, phase",spectrum,phase);

		isign = BACKWARD;  // defined in SSPconstants.js, BACKWARD = -1
		var inverse = fft(spectrum, nn, ndim, isign);
		console.log("inverse audio: ",spectrum, inverse);

		console.log("Stop testLigteringen");
		console.log("");
		};
	};


function testFFT2D(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testFFT2D");
	
		// prepare unit impulse 2D signal (image)
		var dims = arrayDimension(image1);
		nn = [dims[0]/2, dims[1]];  // correction for complex numbers
		ndim = nn.length;
		var rawImage = flatten(image1);
		var myLength = (rawImage.length)/2;
		var data2D = rawImage.slice();
		
		// forward transform
		var isign = FORWARD;
		var spectrum = data2D.slice();
		// the 2D fft overwrites the data it receives
		var spectrum = fft(spectrum, nn, ndim, isign);
//		console.log("fft 2D impulse: ",data2D, round(spectrum,rounding));

		isign = BACKWARD;
		var inverse = spectrum.slice();
		// the 2D inverse fft overwrites the data it receives
		var inverse = fft(inverse, nn, ndim, isign);
		console.log("inverse fft 2D: ",round(spectrum,rounding), round(inverse,rounding));

		// prepare sinusoidal signal
		var harmonic = 3;
		for(var i = 0; i< myLength; i++)
			{
			data2D[2*i] = Math.cos(i*harmonic*Math.PI / (myLength/dims[1]));
			data2D[2*i+1] = 0;
			};

		// forward transform
		isign = FORWARD;
		var spectrum = data2D.slice();
		// the 2D fft overwrites the data it receives
		spectrum = fft(spectrum, nn, ndim, isign);
		console.log("fft 2D cos(): ",round(data2D,rounding), round(spectrum,rounding));
		
		isign = BACKWARD;
		var inverse = spectrum.slice();
		// the 2D inverse fft overwrites the data it receives
		var inverse = fft(inverse, nn, ndim, isign);
		console.log("inverse fft 2D: ",round(spectrum,rounding), round(inverse,rounding));

		// prepare unit impulse 2D full-size image
		var rows = 512
		var cols = 512;
		nn = [rows, cols];  // correction for complex numbers
		ndim = nn.length;
		myLength = rows*cols;
		rawImage = createArray(2*myLength);
		harmonic = 3;
		for(var i = 0; i< myLength; i++)
			{
			rawImage[2*i] = Math.cos(i*harmonic*Math.PI / (myLength/cols));
			rawImage[2*i+1] = 0;
			};
		rawImage[0] += 1;

		data2D = rawImage.slice();
		
		// forward transform
		var start = performance.now();
		isign = FORWARD;
		spectrum = data2D.slice();
		// the 2D fft overwrites the data it receives
		spectrum = fft(spectrum, nn, ndim, isign);
		var stopfft = performance.now();

		console.log("fft 2D impulse + cos(): ",round(data2D,rounding), round(spectrum,rounding));

		var middle = performance.now();
		isign = BACKWARD;
		inverse = spectrum.slice();
		// the 2D inverse fft overwrites the data it receives
		inverse = fft(inverse, nn, ndim, isign);
		var end = performance.now();
		console.log("inverse fft 2D: ",round(spectrum,rounding), round(inverse,rounding));
		console.log("fft 2D timing – dimensions, fft, ifft [ms]: ",nn,Math.round(stopfft-start),Math.round(end-middle));

		console.log("Stop testFFT2D");
		console.log("");
		};
	};


// ****************************************************************
// ****************************************************************
// For 1D signals
// data is assumed to be complex with the structure data = [real_array,imag_array]
// transform & inverseTransform mutate the arrays sent to the subroutine
// Caveat emptor

/* No longer used replacde by fft()

function myFourier(data)
	{
	var real = data[0].slice();
	var imag = data[1].slice();
	transform(real, imag);
	return [real,imag];
	};

function myInverseFourier(data)
	{
	var real = data[0].slice();
	var imag = data[1].slice();
	inverseTransform(real, imag);
	var n = real.length;
	
	// install the missing scaling
	for(var i = 0; i < n; i++)
		{
			real[i] /= n;
			imag[i] /= n;
		};
	return [real,imag];
	};


function testFourier(testChoice)
	{
		if(testChoice)
		{
		console.log("\nStart");
		var myLength = 16;
		var real = Array(myLength);
		var imag = createArray(real.length);

//		for (var i = 0; i < myLength; i++)	real[i] = randomInteger(0,10);
//		var imag = createArray(real.length);
//		for(i = 0; i < imag.length; i++)	imag[i] = randomInteger(0,2);

		for(i = 0; i < imag.length; i++)	imag[i] = 0;
		var real = imag.slice();
		real[0] = 1;

		for(var i = 0; i< myLength; i++)
			real[i] = Math.sin(i*2*Math.PI / myLength);
		var imag = createArray(real.length);
		for(i = 0; i < imag.length; i++)	imag[i] = 0;

		rcopy = real.slice();
		icopy = imag.slice();
		var data = [rcopy,icopy];
	
		console.log("Input to transform: ",rcopy,icopy);
	
		var spectrum = myFourier(data);
		console.log("output from transform: ",spectrum,round(spectrum[0],rounding),
			round(spectrum[1],rounding));
		
		var inverse = myInverseFourier(spectrum);
		console.log("output from inverse transform: ",round(inverse[0],rounding),
			round(inverse[1],rounding));
		
		console.log("Stop");
		console.log("");
		};
	};


function timeFourier(testChoice)
	{
		if(testChoice)
		{
		console.log("\nStart");
		var nSamples = Array(42-8+1);
		var time = Array(nSamples.length);
		var iterations = 10;
	
		for(var i = 0; i < nSamples.length; i++)
			{
			nSamples[i] = Math.floor(Math.sqrt(Math.pow(2,i+1)));
			time[i] = 0;
			};
	
		for(var i = 0; i < nSamples.length; i++)
			{
			var samples = nSamples[i];
			var rnoise = Array(samples);
			var inoise = Array(samples);
			for(var j = 0; j < samples; j++)
				{
				rnoise[j] = Math.random();
				inoise[j] = Math.random();
				};
			var noise = [rnoise,inoise];  // complex noise input
			var spectrum = [];
			var start = new Date().getTime();
			for(var k = 0; k < iterations; k++) spectrum = myFourier(noise);
			var finish = new Date().getTime(); // in ms
			time[i] = finish - start;
			};
		
		for(var i = 0; i < nSamples.length; i++)
			{
			time[i] /= iterations;  // time per iteration in ms
			};
		console.log("myFourier Samples & Time [ms]: ",nSamples,time);

		console.log("Stop");
		console.log("");
		};
	};

*/

