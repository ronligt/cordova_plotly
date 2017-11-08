// ****************************************************************
// Probability and statistical procedures to be used in the SSP package
//
// All of the routines for generating random numbers according to a
// specific, desired distribution provide a single value of the random
// variate. That is, to generate a list of 10,000 numbers you must call
// the given routine 10,000 times.
//
// Statistics routines, however, assume that an array of data is presented
// that consists of (preferably) more than one sample. This array must be
// in the "universal" format described in the header of "SSPconstants.js"
// so that it can be used on 1D and 2D complex data.
//
// i.t. young
// Tuesday, 5 July 2016


// ****************************************************************
function randomInteger(lower,upper)
// Uniform integer distribution on x for lower ≤ x ≤ upper
// Checked Wednesday, 11 January 2017

	{
	var ulimit = upper - lower + 1;
	return (Math.floor(ulimit*Math.random() + lower));
	};


function testRandomInteger(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testRandomInteger");
		
		var myLength = 10000;
		var numbers = Array(myLength);
		var lower = 1;
		var upper = 9;
		for (var i = 0; i < myLength; i++)	numbers[i] = randomInteger(lower,upper);
		console.log("randomInteger - list, mean, average, sigma: ",
			numbers,(upper+lower)/2,mean(numbers),stdev(numbers));
			
		console.log("Stop testRandomInteger");
		};
	};


// ****************************************************************
function randomReal(lower,upper)
// Uniform distribution on x for lower ≤ x < upper
// Checked Wednesday, 11 January 2017

	{
	var ulimit = upper - lower;
	return (ulimit*Math.random() + lower);
	};


function testRandomReal(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testRandomReal");
		
		var myLength = 10000;
		var uniform = Array(myLength);
		var lower = 0;
		var upper = 10;
		for (var i = 0; i < myLength; i++)
			uniform[i] = randomReal(lower,upper);
		console.log("randomReal - list, mean, average, th_var, meas_var: ",
			uniform,(upper-lower)/2,mean(uniform),
			(upper-lower)*(upper-lower)/12,variance(uniform));
			
		console.log("Stop testRandomReal");
		console.log("");
		};
	};


// ****************************************************************
// Based upon remark in Wikipedia about Triangular distribution
// <https://en.wikipedia.org/wiki/Triangular_distribution>.
// See "Generating Triangular-distributed random variates"
// Checked Wednesday, 11 January 2017

function randomTriangular(lower,upper,mode)
// Triangular distribution on x for lower ≤ x ≤ upper with peak at mode
	{
	var f = (mode - lower)/(upper - lower);
	var u = Math.random();
	var result = 0;
	
	if (0 < u && u < f)
		{
		result = lower + Math.sqrt(u*(upper - lower)*(mode - lower));
		};
	if (f <= u && u < 1)
		{result = upper - Math.sqrt((1-u)*(upper - lower)*(upper - mode));
		};
	return result;
	};


function testRandomTriangular(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testRandomTriangular");
		
		var myLength = 100000;
		var triangular = Array(myLength);
		var lower = 0;
		var upper = 10;
		var mode = (lower + upper) / 9;
		for (var i = 0; i < myLength; i++)
			triangular[i] = randomTriangular(lower,upper,mode);
		console.log("randomTriangular - list, mean, average, th_var, meas_var: ",
			triangular,(mode+lower+upper)/3,
			mean(triangular),
			((mode*mode)-(mode*upper)+(upper*upper)-
			(mode*lower)-(lower*upper)+(lower*lower))/18,
			variance(triangular));
			
		console.log("Stop testRandomTriangular");
		console.log("");
		};
	};


// ****************************************************************
// Copied from:
// <http://www.ollysco.de/2012/04/gaussian-normal-functions-in-javascript.html>

// Returns a Gaussian Random Number around a normal distribution defined by the mean
// and standard deviation parameters.

// Uses the algorithm used in Java's random class, which in turn comes from
// Donald Knuth's implementation of the Box–Muller transform.

// @param {Number} [mean = 0.0] The mean value, default 0.0
// @param {Number} [standardDeviation = 1.0] The standard deviation, default 1.0
// @return {Number} A random number

// Checked Wednesday, 11 January 2017

function randomGaussian(mean, standardDeviation)
	{	
	if (randomGaussian.nextGaussian !== undefined)
		{
		var nextGaussian = randomGaussian.nextGaussian;
		delete randomGaussian.nextGaussian;
		return (nextGaussian * standardDeviation) + mean;
		}
	else
		{
		var v1, v2, s, multiplier;
		do
			{
			v1 = 2 * Math.random() - 1; // between -1 and 1
			v2 = 2 * Math.random() - 1; // between -1 and 1
			s = v1 * v1 + v2 * v2;
			} while (s >= 1 || s == 0);
		multiplier = Math.sqrt(-2 * Math.log(s) / s);
		randomGaussian.nextGaussian = v2 * multiplier;
		return (v1 * multiplier * standardDeviation) + mean;
		};
	};


function testRandomGaussian(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testRandomGaussian");
		
		var myLength = 10000;
		var normal = Array(myLength);
		var average = 5;
		var sigma = Math.SQRT2;
		for (var i = 0; i < myLength; i++)
			normal[i] = randomGaussian(average,sigma);
		console.log("randomGaussian - list, mean, average, th_sigma, meas_sigma: ",
			normal,average,mean(normal),sigma,stdev(normal));
			
		console.log("Stop testRandomGaussian");
		console.log("");
		};
	};


// ****************************************************************
// Based upon remark in Wikipedia about Rayleigh distribution
// <https://en.wikipedia.org/wiki/Rayleigh_distribution>.
// "[If] random complex numbers whose real and imaginary components are
// independently and identically distributed Gaussian with equal variance
// and zero mean, [then] the absolute value of the complex number is
// Rayleigh-distributed"
// Checked Wednesday, 11 January 2017

function randomRayleigh(standardDeviation)
	{	
	var v1, v2, s, multiplier;
	do
		{
		v1 = 2 * Math.random() - 1; // between -1 and 1
		v2 = 2 * Math.random() - 1; // between -1 and 1
		s = v1 * v1 + v2 * v2;
		} while (s >= 1 || s == 0);
	multiplier = Math.sqrt(-2 * Math.log(s) / s);
	v1 *= multiplier;
	v2 *= multiplier;
	return (standardDeviation*Math.hypot(v1,v2));
	};


function testRandomRayleigh(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testRandomRayleigh");
		
		var myLength = 100000;
		var rayleigh = Array(myLength);
		var lambda = 8;
		for (var i = 0; i < myLength; i++)	rayleigh[i] = randomRayleigh(lambda);
		console.log("randomRayleigh - list, lambda, th_mean, average, th_var, meas_var: ",
			rayleigh,lambda,Math.sqrt(Math.PI/2)*lambda,mean(rayleigh),
			(2-(Math.PI/2))*lambda*lambda,variance(rayleigh));
			
		console.log("Stop testRandomRayleigh");
		console.log("");
		};
	};


// ****************************************************************
// Based upon Knuth, Semi-Numerical Algorithsms, Vol. 2, p. 117, (1969)
// Returns a Poisson-distributed random integer whose mean is given.
// Checked Wednesday, 11 January 2017

function randomPoisson(mean)
	{ 
	var p = Math.exp(-mean);
	var q = 1.0;
	var n = 0;

	while (q >= p)
		{
		n++;
		q *= Math.random();
		};
	return (n - 1);
	};


function testRandomPoisson(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testRandomPoisson");
		
		var myLength = 10000;
		var poisson = Array(myLength);
		var lambda = Math.SQRT2;
		for (var i = 0; i < myLength; i++)	poisson[i] = randomPoisson(lambda);
		console.log("randomPoisson - list, lambda, th_avg, average, th_var, meas_var: ",
			poisson,lambda,lambda,mean(poisson),lambda,variance(poisson));
			
		console.log("Stop testRandomPoisson");
		console.log("");
		};
	};


// ****************************************************************
// Based upon Knuth, Semi-Numerical Algorithsms, Vol. 2, p. 114, (1969)
// Returns an Exponential-distributed random real whose mean is given.
// Checked Wednesday, 11 January 2017

function randomExponential(lambda)
	{ 
	var p = Math.random();

	if (p < tooSmall) p = tooSmall;
	
	return (-Math.log(p)*lambda);
	};


function testRandomExponential(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testRandomExponential");
		
		var myLength = 10000;
		var exponential = Array(myLength);
		var lambda = 10;
		for (var i = 0; i < myLength; i++)
			exponential[i] = randomExponential(lambda);
		console.log("randomExponential - list, MinOfArray, MaxOfArray, lambda, th_avg,  average, lambda^2, th_var, meas_var: ",exponential,getMinOfArray(exponential), getMaxOfArray(exponential),lambda,lambda,mean(exponential),lambda*lambda,variance(exponential));
		
		console.log("Stop testRandomExponential");
		console.log("");
		};
	};


// ****************************************************************
// Method #1 Knuth, "Seminumerical Algorithms", Vol2 Sec 3.4.1 D,
// the "logarithm method"

// Method #2 See
// <https://en.wikipedia.org/wiki/Laplace_distribution#Generating_random_variables_according_to_the_Laplace_distribution>
// "The difference between two independent identically distributed
// exponential random variables is governed by a Laplace distribution"
//
// The following is a slightly modified Method #1.

// Checked Wednesday, 11 January 2017

function randomLaplace(mean,b)
	{ 
	var u = randomReal(-1/2,+1/2);
	if (u === 0) u = Math.sign(u)*tooSmall;  // non illigitimum carborundum
	var v = mean - b * Math.sign(u) * Math.log(1-(2*Math.abs(u)));
	
	return v;
	};


function testRandomLaplace(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testRandomLaplace");
		
		var myLength = 100000;
		var laplace = Array(myLength);
		var beta = 10;
		var average = 73;
		for (var i = 0; i < myLength; i++)	laplace[i] = randomLaplace(average,beta);
		console.log("randomLaplace - list, lambda, beta, th_avg, mean, th_var, meas_var: ",laplace,average,beta,mean(laplace),2*beta*beta,variance(laplace));
			
		console.log("Stop testRandomLaplace");
		console.log("");
		};
	};


// ****************************************************************
// Based upon flipping a coin with prob(Heads) = p and
// prob(Tails) = 1-p. The former returns 1 and the latter 0.
// Checked Wednesday, 11 January 2017

function randomBinary(p)
	{ 
	var value = Math.random();
	
	return ((value < p) ? 1 : 0);
	};


function testRandomBinary(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testRandomBinary");
		
		var myLength = 100000;
		var binary = Array(myLength);
		var p = Math.SQRT1_2;
		for (var i = 0; i < myLength; i++)	binary[i] = randomBinary(p);
		console.log("randomBinary - list, p, mean, th_var, meas_var: ",
			binary,p,mean(binary),p*(1-p),variance(binary));
			
		console.log("Stop testRandomBinary");
		console.log("");
		};
	};


// ****************************************************************
// All of the routines for generating random numbers according to a
// specific, desired distribution provide a single value of the random
// variate. That is, to generate a list of 10,000 numbers you must call
// the given routine 10,000 times.
//
// Statistics routines, however, assume that an array of data is presented
// that consists of (preferably) more than one sample. This array must be
// in the "universal" format described in the header of "SSPconstants.js"
// so that it can be used on 1D and 2D complex data.
//
// ****************************************************************
// Checked Wednesday, 11 January 2017

function mean(data)
	{
	var totalr = 0;
	var totali = 0;
	var dLength = data.length;
	for (var i = 0; i < dLength; i += 2)
		{
		totalr += data[i];
		totali += data[i+1];
		};
		return [2*totalr/dLength, 2*totali/dLength];
	};

function testMean(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testMean");
		var tdata = nlist1;
		console.log("mean - list, samples, average:",tdata,tdata.length,mean(tdata));
	
		var myLength = 10000
		var tlist1 = Array(2*myLength); /* temporary list */
		for (var i = 0; i < 2*myLength; i += 2)
			{
			tlist1[i] = randomInteger(0,1);
			tlist1[i+1] = 2*Math.PI*randomInteger(0,1);
			};
		console.log("mean - samples, th_avg, average:",myLength,1/2,Math.PI,mean(tlist1));
	
		var myLength = 10000
		var tlist2 = Array(2*myLength); /* temporary list */
		for (var i = 0; i < 2*myLength; i += 2)
			{
			tlist2[i] = randomReal(2,23);
			tlist2[i+1] = 0;
			};

		console.log("mean - samples, th_avg, average:",myLength,(2+23)/2,0,mean(tlist2));
		
		console.log("Stop testMean");
		console.log("");
		};
	};


// ****************************************************************
// Checked Thursday, 12 January 2017

function variance(data)
	{
	var z = Array(data.length);
	var dLength = (data.length)/2;
	var tmean = mean(data);
	var meanr = tmean[0];
	var meani = tmean[1];
	// computed the old fashioned way
	for (var i = 0; i < 2*dLength; i += 2)
		{
		z[i] = data[i] - meanr;
		z[i+1] = data[i+1] - meani;
		};
	var sqz = dot(z,z);
	return [sqz[0]/(dLength-1), 0];
	};

function stdev(data)
	{
	var temp = variance(data);
	var std = Math.sqrt(temp[0]);
	return [std, 0];
	};

function testVarStd(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testVarStd");
		
		var tdata = [1,0,1,0,1,0];
		console.log("variation - data, samples, th_avg, average, th_var, variance, th_sigma, sigma:",
			tdata,tdata.length/2,1,mean(tdata),0,variance(tdata)[0],Math.sqrt(0),
			stdev(tdata)[0]);
	
		tdata = [1,0,2,0,3,0,4,0];
		console.log("variation - data, samples, th_avg, average, th_var, variance, th_sigma, sigma:",
			tdata,tdata.length/2,2.5,mean(tdata),5/3,variance(tdata)[0],Math.sqrt(5/3),
			stdev(tdata)[0]);
	
		tdata = nlist1;
		console.log("variation - data, samples, average, variance, sigma:",
			tdata,tdata.length/2,mean(tdata),variance(tdata)[0],stdev(tdata)[0]);
	
		var myLength = 100000
		var tlist = Array(myLength); /* temporary list */
		var mu = -1;
		var sigma = 7;
		for (var i=0; i<myLength; i++)
			{
			tlist[i] = randomGaussian(mu,sigma);
			};

		console.log("variation Gauss - samples, mean, average, variance, th_sigma, meas_sigma:",
			myLength/2,mu,mean(tlist),variance(tlist)[0],
			sigma*Math.SQRT2,stdev(tlist)[0]);
	
		var myLength = 100000
		var tlist = Array(myLength); /* temporary list */
		var lambda = 17;
		for (var i=0; i<myLength; i++)
			{
			tlist[i] = randomPoisson(lambda);
			};

		console.log("variation Possson - samples, mean, average, th_var, variance, meas_sigma:",
			myLength/2,lambda,mean(tlist),2*lambda,variance(tlist)[0],
			stdev(tlist)[0]);

		console.log("Stop testVarStd");
		console.log("");
		};
	};


// ****************************************************************
// Mathematica function Quantile
// Returns the value at a given percentile in a sorted numeric array.
// "Linear interpolation between closest ranks" method. Take from:
// <https://gist.github.com/IceCreamYou/6ffa1b18c4c8f6aeaad2>
//
// We assume the data are REAL as one cannot sort complex numbers.
// Thus we only look at the real part of the array data.
//
// If the array is not already sorted, it will be sorted. If it is already
// sorted then no harm will be done.
// Checked Monday, 16 January 2017

function quantile(arr, p)
	{
		if (typeof p !== 'number') throw new TypeError('p must be a number');
		if (arr.length === 0) return [0,0];
		var sdata = sort(arr);
		if (p <= 0) return [sdata[0], 0];
		if (p >= 1) return [sdata[sdata.length] - 2, 0];
		if (p === 1/2) return median(sdata);

		var index = (sdata.length/2) * p;
		var lower = Math.floor(index);
		var upper = lower + 1;
		var weight = index % 1;
		if (2*upper >= sdata.length) return sdata[2*lower];
		return [sdata[2*lower] * (1 - weight) + sdata[2*upper] * weight, 0];
	};

// Returns the percentile of the given value in a sorted numeric array.
function percentRank(arr, v)
	{
		if (typeof v !== 'number') throw new TypeError('v must be a number');
		for (var i = 0, l = arr.length; i < l; i++)
			{
				if (v <= arr[i])
					{
						while (i < l && v === arr[i]) i++;
							if (i === 0) return 0;
							if (v !== arr[i-1])
								{
									i += (v - arr[i-1]) / (arr[i] - arr[i-1]);
								}
					return i / l;
					};
			};
		return 1;
	};

function testQuantile(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testQuantile");
		
		var tdata = nlist4.slice();  // no mutation
		var iter = tdata.length;
		var p = 1/3;
		console.log("sort 1: ",tdata,p,quantile(tdata,p));

		tdata = nlist5.slice();  // no mutation
		iter = tdata.length;
		p = 1/2;
		console.log("sort 2: ",tdata,p,quantile(tdata,p));
	
		var myLength = 100000
		tdata = Array(myLength); /* temporary list */
		var lambda = 9;
		for (var i=0; i<myLength; i++)
			{
			tdata[i] = randomExponential(lambda);
			};
			
		p = 0.4999;
		console.log("quantile - samples, p, lambda, th_median, meas_median: ",
			tdata.length,p,lambda,lambda*Math.LN2,quantile(tdata,p));

		p = 1/2;
		console.log("quantile - samples, p, lambda, th_median, meas_median: ",
			tdata.length,p,lambda,lambda*Math.LN2,quantile(tdata,p));

		console.log("Stop");
		console.log("");
		};
	};


// ****************************************************************
// Copied from <https://gist.github.com/caseyjustus/1166258>
// See description from quantile (above) for more details
//
// Checked Monday, 16 January 2017

function median(data)
	{
	if (data.length === 0) return [0,0];
	var values = sort(data);
	var dLength = data.length/2

	if(oddQ(dLength)) // test if length is even or odd
		return [values[dLength-1], 0];  // odd
	else
		return [(values[dLength-2] + values[dLength]) / 2.0, 0];  // even
}


function testMedian(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testMedian");
		
		var tdata = [1,0,2,0,3,0];
		console.log("median - list, mean, median: ",tdata,median(tdata));
	
		var tdata = [1,0,2,0,3,0,4,0,5,0,6,0,7,0];
		console.log("median - list, mean, median: ",tdata,median(tdata));
	
		var tdata = [1,0,2,0,3,0,4,0];
		console.log("median - list, mean, median: ",tdata,median(tdata));
	
		var tdata = nlist4.slice();
		console.log("median - list, mean, median: ",tdata,median(tdata));
	
		tdata = nlist5.slice();
		console.log("median - list, mean, median: ",tdata,median(tdata));
	
		var myLength = 100000
		var tdata = Array(myLength); /* temporary list */
		var lambda = 9;
		for (var i=0; i<myLength; i++)
			{
			tdata[i] = randomRayleigh(lambda);
			};
		console.log("median Rayleigh - samples, lambda, th_median, meas_median: ",
			myLength,lambda,lambda*Math.sqrt(Math.log(4)),median(tdata));
	
		console.log("Stop");
		console.log("");
		};
	};


// ****************************************************************
// Based upon Mathematica's HistogramList function. Given an array of
// REAL numbers, produce a histogram list of the occurrences of the values.
// The bin width for the histograms is determined by "width". If width
// is TRUE, then the width of the bin, used for the construction of the
// histogram, will be 1. An image with 256 grey levels, for example, will
// have 256 bins each of width 1 irrespective of the number of pixels.
// (The actual Mathematica routine is more general but enough is enough.)
// If the width is FALSE, then the width will be computed automatically.
// The range of the data will be defined as 0.9*minval to 1.1*maxval.
// If the number of data samples (e.g. pixels) is N, the number of bins
// will be approximately sqrt(N).
//
// Checked Friday, 27 January 2017

function histogramList(data,width)
	{
	var myData = data.slice();
	var totalN = myData.length/2;  // complex data format
	var range = 0;
	var bins = 0;
	var binWidth = 0;
	var value = 0;
	var frac = 0;
	var index = 0;
	var dLength = 0;
	var minval = Infinity;
	var maxval = -Infinity;
	
	// The following is for bin width = 1.
	if (width)
		{
		// to avoid stack issues the min and max are computed here
		// on already flattened data
		for (var i = 0; i < 2*totalN; i += 2)
			{
			value = myData[i];
			if (value < minval) minval = value;
			if (value > maxval) maxval = value;
			};
		
		range = Math.ceil(maxval) - Math.floor(minval) + 1;
		var dbins = range;
		if (Math.ceil(maxval) === maxval)
			var hbins = range;
		else
			var hbins = range - 1;
		binWidth = 1;
		
		// Create list of bin delimiters
		var delimitersList = Array(2*dbins);  // complex data format
		dLength = delimitersList.length;
		for (var i = 0; i < dLength; i += 2)
			{
			delimitersList[i] = binWidth*i/2 + Math.floor(minval);
			delimitersList[i+1] = 0;
			};
		
		// Initialize histogram array
		var myHistogramList = Array(2*hbins);  // complex data format
		myHistogramList.fill(0);
		
		// Compute histogram
		for (var i = 0; i < 2*totalN; i += 2)
			{
			value = Math.floor(myData[i] - minval);
			index = 2*value;
			myHistogramList[index] += 1;
			};
		}
	// The following is for automatically computed bin width.
	else
		{
		var sorted = sort(myData);
		minval = sorted[0];
		maxval = sorted[sorted.length-2];
		range = (maxval - minval);
		bins = Math.floor(Math.sqrt(totalN));
		binWidth = range / bins;
		
		// Create list of bin delimiters
		var delimitersList = Array(2*(bins+1));
		dLength = delimitersList.length;
		for (var i = 0; i < dLength; i += 2)
			{
			delimitersList[i] = binWidth*(i/2) + minval;
			delimitersList[i+1] = 0;
			}
		
		// Initialize histogram array
		var myHistogramList = Array(2*bins);
		myHistogramList.fill(0);
		
		for (var i = 0; i < 2*totalN; i += 2)
			{
			frac = (myData[i] - minval)/range;  // 0 ≤ frac ≤ 1
			index = Math.floor(bins*frac);
			if (index === bins) index -= 1;
			myHistogramList[2*index] += 1;
			};

		};
	return [delimitersList, myHistogramList];
	}

function testHistogramList(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testHistogramList");
		
		var data = [2,0,3,0,4,0,5,0,7,0];
		var endList = histogramList(data,true);
		console.log("histogramList TRUE, data, delimiters, histogram, total N: ", sort(data), endList[0], endList[1], total(endList[1]));
		
		data = [15,0,16,0,18,0,17,0,1,0,2,0,3,0,4,0,5,0];
		endList = histogramList(data,true);
		console.log("histogramList TRUE, data, delimiters, histogram, total N: ", sort(data), endList[0], endList[1], total(endList[1]));
		
		endList = histogramList(data,false);
		console.log("histogramList FALSE, data, delimiters, histogram, total N: ", sort(data), round(endList[0],rounding), endList[1], total(endList[1]));
		
		var data = [27.1976, 0, 16.2564, 0, 23.5566, 0, 9.60064, 0, 12.2811, 0,
		 16.3955, 0, 13.6815, 0, 0.408964, 0, 4.46461, 0, 12.0228, 0, 7.18286, 0,
		 5.88196, 0, 10.0249, 0, 3.83276, 0, 10.0125, 0, 1.3949, 0, 11.4326, 0,
		 2.11597, 0, 4.50059, 0, 8.2486, 0, 3.79871, 0, 15.123, 0, 9.02716, 0,
		 10.5584, 0, 17.3638, 0, 5.19686, 0, 11.0076, 0, 8.92813, 0, 13.9632, 0, 
		 11.466, 0, 14.1743, 0, 1.06765, 0, 9.55211, 0, 19.4067, 0, 16.6116, 0,
		 20.7735, 0, 2.08838, 0, 18.9789, 0, 14.9563, 0, 4.00338, 0, 21.7476, 0,
		 6.14579, 0, 14.4046, 0, 8.75238, 0, -3.07091, 0, -2.07081, 0, 11.1903, 0,
		 21.5998, 0, 12.6512, 0, 4.90017, 0];
		endList = histogramList(data,false);
		console.log("histogramList FALSE, sort(data), delimiters, histogram, total N: ", sort(data), round(endList[0],rounding), endList[1], total(endList[1]));
		
		endList = histogramList(data,true);
		console.log("histogramList TRUE, sort(data), delimiters, histogram, total N: ", sort(data), round(endList[0],rounding), endList[1], total(endList[1]));
		
		var myLength = 1000;
		var gauss = Array(2*myLength);
		var mu = 0;
		var sigma = 1;
		for (var i = 0; i < 2*myLength; i += 2)
			{
			gauss[i] = randomGaussian(mu,sigma);
			gauss[i+1] = 0;
			};
		endList = histogramList(gauss,true);
		console.log("histogramList TRUE, sort(gauss), delimiters, histogram, total N: ", round(sort(gauss),rounding),endList[0], endList[1], total(endList[1]));

		endList = histogramList(gauss,false);
		console.log("histogramList FALSE, sort(gauss), delimiters, histogram, total N: ", 
round(sort(gauss),rounding),round(endList[0],rounding),endList[1], total(endList[1]));

		console.log("Stop");
		console.log("");
		};
	};


// ****************************************************************
// Based upon Mathematica's PearsonChiSquareTest function. There is a 
// hypothesis that data collected are from a certain probability
// distribution / density function. We examine the normalized mean square
// difference between the actual data histogram and the histogram from the
// hypothesized distribution with the same number of samples N.
// Based upon procedure found at <http://simplestatistics.org>

// The χ2 statistic uses the χ2 probability distribution which is an
// incomplete gamma function. See Numerical Recipes 6.2.18 which points
// to the procedure "gammp", "gammq", and others
//
// Checked Friday, 27 January 2017

function gammp(a, x) 
	{  // Note: no checking of variables
	var result = [];
	if (x < a+1)  // use series approximation
		{
		result = gser(a, x);
		return result[0];
		}
	else  // use continued fraction representation
		{
		result = gcf(a, x);
		return 1-result[0];  // and take complement
		};
	};
	
function gammq(a, x) 
	{  // Note: no checking of variables
	var result = [];
	if (x < a+1)  // use series approximation
		{
		result = gser(a, x);
		return 1-result[0];
		}
	else  // use continued fraction representation
		{
		result = gcf(a, x);
		return result[0];  // and take complement
		};
	};
	
function gser(a,x)
	{  // Note: no checking of variables
	var gamser = 0;
	var gln = gammln(a);
	var ap = a;
	var sum = 1/a;
	var del = sum;
	for (var n = 1; n <= iTMAX; n++)
		{
		++ap;
		del *= (x/ap);
		sum += del;
		if (Math.abs(del) < Math.abs(sum)*ePS)
			{
			gamser = sum*Math.exp(-x + (a*Math.log(x)) - gln);
			return [gamser, gln];
			};
		};
	throw("gser: a too large, iTMAX too small");
	};

function gcf(a,x)
	{  // Note: no checking of variables
	var gln = gammln(a);
	var b = x+1-a;
	var c = 1/fPMIN;
	var d = 1/b;
	var h = d;
	var an = 0;
	var del = 0;
	for (var i = 1; i <= iTMAX; i++)
		{
		an = -i*(i-a);
		b += 2;
		d = an*d + b;
		if (Math.abs(d) < fPMIN) d = fPMIN;
		c = b+(an/c);
		if (Math.abs(c) < fPMIN) c = fPMIN;
		d = 1/d;
		del = d*c;
		h *= del;
		if (Math.abs(del-1) < ePS) break;
		};
	if (i > iTMAX) throw("gcf: a too large, iTMAX too small");
	var gammcf = Math.exp(-x + (a*Math.log(x)) - gln)*h;
	return [gammcf, gln];
	
	};

function gammln(xx)
	{  // Note: no checking of variables
	var cof = [76.18009172947146, -86.50532032941677, 24.01409824083091,
		-1.231739572450155, 0.1208650973866179e-2,-0.5395239384953e-5];
	var x = xx;
	var y = x;
	var tmp = x+5.5;
	tmp = tmp - (x+0.5)*Math.log(tmp);
	var ser = 1.000000000190015;
	for (var j = 0; j < cof.length; j++)
		{
		++y;
		ser += (cof[j]/y);
		};
	return -tmp+Math.log(2.506628274631005*ser/x);
	};
	
// Here begin the various density functions currently supported in this test
function gaussDensity(mu,sigma,totalN,delimiters)
	{
	var x = 0;
	var z = 0;

	var xmin = delimiters[0];  // not -Infinity
	var xmax = delimiters[delimiters.length - 2];  // not +Infinity

	var binWidth = delimiters[4] - delimiters[2];

	// the expected histogram cells to be returned
	var bins = delimiters.length/2 - 1;
	var probE = 0;
	var probP = 0;
	var cells = [];
	var normalize = 0;
	var aconst = 1/Math.sqrt(2*Math.PI);
	for (var i = 0; i < 2*bins; i += 2)
		{

		// in the middle of each bin
		x = delimiters[i] + binWidth/2;
		z = (x - mu)/sigma;

		// probability of being in this bin
		// probE = Math.exp(-(z*z)/2);
		probP = cdfGauss(delimiters[i+2], mu, sigma) - cdfGauss(delimiters[i], mu, sigma);

		normalize += probP;
		cells[i] = probP;
		cells[i+1] = 0;
		};
		
	// correct for the use of a limited number of bins of finite width
	var bconst = totalN/normalize;
	for (var i = 0; i < 2*bins; i += 2)  cells[i] *= bconst;

	return cells;
};

function laplaceDensity(mid,b,totalN,delimiters)
	{
	var x = 0;
	var z = 0;
	var prob = 0;
	
	var xmin = delimiters[0];  // not -Infinity
	var xmax = delimiters[delimiters.length - 2];  // not +Infinity
	var binWidth = delimiters[4] - delimiters[2];
	var bins = delimiters.length/2 - 1;
	// the expected histogram cells to be returned
	var cells = [];
	
	var normalize = 0;
	for (var i = 0; i < 2*bins; i += 2)
		{
		// probability of each bin
		prob = cdfLaplace(delimiters[i+2],mid,b) - cdfLaplace(delimiters[i],mid,b);
		normalize += prob;
		cells[i] = prob;
		cells[i+1] = 0;
		};
		
	// correct for the use of a limited number of bins
	var bconst = totalN/normalize;
	for (var i = 0; i < 2*bins; i += 2)  cells[i] *= bconst;
	
	return cells;
};

function uniformDensity(left,right,totalN,delimiters)
	{
	// the expected histogram cells to be returned
	var cells = Array(delimiters.length-2);
	var bins = cells.length/2;
	var prob = 1/bins;
	var expect = totalN*prob;
	for (var i = 0; i < 2*bins; i += 2)
		{
		cells[i] = expect;
		cells[i+1] = 0;
		};
	return cells;
};

function chiSquaredGoodnessOfFit(data, distributionType)
	{
	'use strict';
	// Estimate mean & standard deviation from the sample data
	var histoResults = histogramList(data, false);
	var delimiters = histoResults[0].slice();  // delimiters on histogram bins
	var observedN = histoResults[1].slice();  // binned histogram
	var observed = Array(observedN.length/2);
	var expected = Array(observed.length);  // same length as observed
	var totalN = data.length/2;

	var avg = mean(data)[0];  // ML estimate of mean
	var std = stdev(data)[0];  // ML estimate of sigma
	var middle = median(data)[0];  // ML estimate of middle
	var ttotal = 0;
	for (var i = 0; i < data.length; i += 2)
		ttotal += Math.abs(data[i] - middle);
	var inputWidth = ttotal / totalN;  // ML estimate of width
	var value = 0;
	var lower = Infinity;
	var upper = -Infinity;
	for (var i = 0; i < 2*totalN; i += 2)
		{
		value = data[i];
		if (value < lower) lower = value;
		if (value > upper) upper = value;
		};

	var nParams = 0;
	var hypothesizedN = [];
	// Generate the hypothesized distribution.
	switch(distributionType)
		{
		case "Gauss":
			{
			hypothesizedN =
				gaussDensity(avg,std,totalN,delimiters);
			nParams = 2;  // mean & sigma
			break;
			}
	
	// See <https://en.wikipedia.org/wiki/Laplace_distribution#Parameter_estimation>
		case "Laplace":
			{
			hypothesizedN =
				laplaceDensity(middle,inputWidth,totalN,delimiters);
			nParams = 2;  // middle & width
			break;
			}
	
		case "Uniform":
			{
			hypothesizedN =
				uniformDensity(lower,upper,totalN,delimiters);
			nParams = 2;  // lower & upper
			break;
			}
	
		default: throw("Non-existent probability choice");
		};

	// Working backward and forwards through the hypothesized (expected)
	// frequencies, collapse bins if less than three observations are expected
	// for a bin. This transformation is applied to the observed frequencies
	// as well.
	for(var i = 0; i < observed.length; i++)
		{
		observed[i] = observedN[2*i];
		expected[i] = hypothesizedN[2*i];
		};

	// see W.G. Cochran, "Some Methods for Strengthening the Common Chi-Square
	// Tests", Biometrics,10, 417-450, 1954.
	var cutoff = 1;
	for (var k = expected.length - 1; k > 0; k--)
		{
		if (expected[k] < cutoff)
			{
			expected[k-1] += expected[k];
			expected.splice(k, 1);

			observed[k-1] += observed[k];
			observed.splice(k, 1);
			};
		};
		
	for (var k = 0; k < expected.length - 1; k++)
		{
		if (expected[k] < cutoff)
			{
			expected[k+1] += expected[k];
			expected.splice(k, 1);

			observed[k+1] += observed[k];
			observed.splice(k, 1);
			};
		};

	// Degrees of freedom = (number of class intervals -
	// the number of independent parameters for the distribution - 1).
	var n = observed.length;
	var degreesOfFreedom = n-nParams-1;
	var chiSquared = 0;
	var o = 0;
	var e = 0;
	var q = 0;
	for (var i = 0; i < n; i++)
		{
		o = observed[i];
		e = expected[i];
		q = Math.abs(o-e);
		chiSquared += (q*q)/e;
		};

	// We want probChiSquared close to 1.0 if the null hypothesis
	// (good fit) is not to be rejected.
var probChiSquared = gammq(degreesOfFreedom/2,chiSquared/2);
	return [degreesOfFreedom, chiSquared, probChiSquared]
	};

function testChiSquared(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testChiSquared");
		
		// choose a number of degrees-of-freedom f and the χ2 value g from a table,
		// the probability associated with that is given by gammq()
		var x = Math.PI;
		console.log("x, gammln(x), true value: ", x, gammln(x), 0.827695);
		console.log("");


		// choose a number of degrees-of-freedom f and the χ2 value g from a table,
		// the probability associated with that is given by gammq()
		var dof = 40;
		var x = 26.51;
		console.log("dof, x, gammq(dof/2,x/2), true value: ",dof,x,gammq(dof/2,x/2),0.949989);
		console.log("");


		// choose a number of degrees-of-freedom f and the χ2 value g from a table,
		// the probability associated with that is given by gammq()
		var dof = 1;
		var x = 0.5;
		console.log("dof, x, gammq(dof/2,x/2), true value: ",dof,x,gammq(dof/2,x/2),0.4795);


		// these data (N = 50) are from a Gaussian random generator with
		// mu = 10 and sigma = 7
		var testdata = [27.1976, 0, 16.2564, 0, 23.5566, 0, 9.60064, 0, 12.2811, 0,
		 16.3955, 0, 13.6815, 0, 0.408964, 0, 4.46461, 0, 12.0228, 0, 7.18286, 0,
		 5.88196, 0, 10.0249, 0, 3.83276, 0, 10.0125, 0, 1.3949, 0, 11.4326, 0,
		 2.11597, 0, 4.50059, 0, 8.2486, 0, 3.79871, 0, 15.123, 0, 9.02716, 0,
		 10.5584, 0, 17.3638, 0, 5.19686, 0, 11.0076, 0, 8.92813, 0, 13.9632, 0, 
		 11.466, 0, 14.1743, 0, 1.06765, 0, 9.55211, 0, 19.4067, 0, 16.6116, 0,
		 20.7735, 0, 2.08838, 0, 18.9789, 0, 14.9563, 0, 4.00338, 0, 21.7476, 0,
		 6.14579, 0, 14.4046, 0, 8.75238, 0, -3.07091, 0, -2.07081, 0, 11.1903, 0,
		 21.5998, 0, 12.6512, 0, 4.90017, 0];
		var myLength = testdata.length/2;

		var result = chiSquaredGoodnessOfFit(testdata, "Gauss");
		console.log("\nGaussian goodness-of-fit test");
		console.log("cannot REJECT Ho at 5% level if P(dof, χ2) > 0.05");
		console.log("χ2 Gauss 1 - total N, mean, sigma, [dof, χ2, P(dof, χ2]): ", myLength,round(mean(testdata)[0],rounding),round(stdev(testdata)[0],rounding),round(result,rounding));


		myLength = 1000;
		var testdataG1 = Array(2*myLength);
		var mu = 0;
		var sigma = 1;
		for (var i = 0; i < 2*myLength; i += 2)
			{
			testdataG1[i] = randomGaussian(mu,sigma);
			testdataG1[i+1] = 0;
			}

		var result = chiSquaredGoodnessOfFit(testdataG1, "Gauss");
		console.log("\nGaussian goodness-of-fit test");
		console.log("cannot REJECT Ho at 5% level if P(dof, χ2) > 0.05");
		console.log("χ2 Gauss 2 - total N, mean, sigma, [dof, χ2, P(dof, χ2]): ", myLength,round(mean(testdataG1)[0],rounding),round(stdev(testdataG1)[0],rounding),round(result,rounding));


		myLength = 1000;
		var testdataU1 = Array(2*myLength);
		var lower = 1;
		var upper = 6;
		for (var i = 0; i < 2*myLength; i += 2)
			{
			testdataU1[i] = randomReal(lower,upper);
			testdataU1[i+1] = 0;
			}

		var result = chiSquaredGoodnessOfFit(testdataU1, "Uniform");
		console.log("\nUniform goodness-of-fit test");
		console.log("cannot REJECT Ho at 5% level if P(dof, χ2) > 0.05");
		console.log("χ2 Uniform - total N, mean, sigma, [dof, χ2, P(dof, χ2]): ", myLength,round(mean(testdataU1)[0],rounding),round(stdev(testdataU1)[0],rounding),round(result,rounding));
		

		var myLengthL = 1000;
		var testdataL = Array(2*myLengthL);
		var middle = 10;
		var width = 5;
		for (var i = 0; i < 2*myLengthL; i += 2)
			{
			testdataL[i] = randomLaplace(middle,width);
			testdataL[i+1] = 0;
			};

		// these data (N = 50) are from a Laplace random generator with
		// middle = 10 and b = 7

		testdata = [21.9478, 0, 9.43885, 0, -5.01436, 0, 7.13719, 0, 27.019, 0, 
					29.2715, 0, 11.7373, 0, 10.1227, 0, 3.5013, 0, 4.47865, 0,
					12.252, 0, 11.2333, 0, 9.4548, 0, 8.30638, 0, -0.231563, 0,
					21.8995, 0, 14.4165, 0, 3.97976, 0, 11.08, 0, 10.8615, 0,
					23.0937, 0, 5.68964, 0, 15.8092, 0, 7.31233, 0, 19.0831, 0,
					23.0219, 0, 9.6473, 0, 8.37414, 0, 3.72041, 0, 8.32278, 0,
					5.34293, 0, 8.55167, 0, 13.3784, 0, 16.2084, 0, 9.00131, 0,
					7.8496, 0, 0.504579, 0, 13.2438, 0, -4.18604, 0, 17.1178, 0,
					20.7244, 0, 15.6459, 0, 30.1795, 0, 14.1204, 0, 3.60343, 0,
					6.91392, 0, 8.26113, 0, 9.67093, 0, 13.9476, 0, 12.7806, 0];
		myLength = testdata.length/2;
		
		var result = chiSquaredGoodnessOfFit(testdataL,"Laplace");
		console.log("\nLaplace goodness-of-fit test");
		console.log("cannot REJECT Ho at 5% level if P(dof, χ2) > 0.05");
		console.log("Laplace χ2 - total N, [dof, χ2, P(dof, χ2)]: ", myLengthL, round(result,rounding));

		var result = chiSquaredGoodnessOfFit(testdataL,"Gauss");
		console.log("\nLaplace goodness-of-fit test");
		console.log("cannot REJECT Ho at 5% level if P(dof, χ2) > 0.05");
		console.log("Laplace data & Gauss model χ2 - total N, [dof, χ2, P(dof, χ2)]: ", myLengthL, round(result,rounding));

// return;  // for debugging
		console.log("\nStop testChiSquared");
		console.log("");
		};
	};


// ****************************************************************
// Based upon Mathematica's TTest function and "Numerical Recipes in C",
// Second Edition, Section 14.2. There is a  hypothesis that data collected
// are from a certain probability distribution with mean mu. This is the
// "one-sample test concerning means". This routine uses the Student-t test
// to see if this hypothesis, Ho, should be rejected.

// The t statistic uses the incomplete beta function Numerical Recipes 6.4
// which points to the procedure "betai" and "betacf".
//
// Checked Tuesday, 28 February 2017

function betai(a, b, x) 
	{
	if ( (x < 0.0) || (x > 1) ) throw("routine betai: bad x");
	var bt = 0.0;
	if ( (x === 0.0) || (x === 1.0) )
		bt = 0.0;
	else
		bt = Math.exp(gammln(a+b) - gammln(a) - gammln(b) +
		 a*Math.log(x) + b*Math.log(1-x));
	if ( x < (a+1)/(a+b+2) )
		return bt*betacf(a,b,x)/a;  // direct use of continued fraction exp.
	else
		return 1-bt*betacf(b,a,1-x)/b;  // symmetry xform then continued fraction exp.
	};
	
	
function betacf(a, b, x)
	{  // Note: no checking of variables
	var qab = a+b;
	var qap = a+1;
	var qam = a-1;
	var c = 1;
	var d = 1 - (qab*x/qap);
	if (Math.abs(d) < fPMIN) d = fPMIN;
	d = 1/d;
	var h = d;
	var m2 = 0;
	var aa = 0;
	var del = 0;

	for (var m = 1; m <= iTMAX; m++)
		{
		m2 = 2*m;
		aa = m*(b-m)*x/((qam+m2)*(a+m2));
		d = 1+(aa*d);  // even step of the recurrence
		if (Math.abs(d) < fPMIN) d = fPMIN;
		c = 1+(aa/c);
		if (Math.abs(c) < fPMIN) c = fPMIN;
		d = 1/d;
		h = h*(d*c);
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d = 1+(aa*d);  // odd step of the recurrence
		if (Math.abs(d) < fPMIN) d = fPMIN;
		c = 1+(aa/c);
		if (Math.abs(c) < fPMIN) c = fPMIN;
		d = 1/d;
		del = d*c;
		h *= del;
		if (Math.abs(del-1) < ePS) break;  // are we done?
		};

	if (m > iTMAX) throw("betacf: m or b too large or iTMAX too small");
	return h;
	};


function tTest(data, mu)
	{
	'use strict';
	// Estimate mean & standard deviation from the sample data
	var avg = mean(data)[0];  // ML estimate of average
	var std = stdev(data)[0];  // ML estimate of standard deviation
	var dof = (data.length/2) - 2;  // degrees-of-freedom
	
	var t = (avg - mu)*Math.sqrt(dof)/std;
	
	// We want probability > alpha if the null hypothesis
	// (avg ≈ mu) is not to be rejected at the alpha percent level.
	var probT = betai(dof/2, 0.5, dof/(dof + t*t));

	return [dof, t, probT]
	};

function testTTest(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testTTest");
		
		var alpha = 0.05;  // significance level
		var string = [];
		
		// test betai and betacf,
		var dof = 20;
		var a = dof/2;
		var b = 0.5;
		var t = 0.2;
		var x = dof/(dof + t*t);
		var y = betai(a,b,x);
		console.log("\nbetai - a,b,t,x,betai(a,b,x),1-betai(): ", a,b,t,round(x,rounding),round(y,rounding),round(1-y,rounding));


		// these data (N = 50) are from a Gaussian random generator with
		// mu = 10 and sigma = 7
		var mu = 10;
		var testdata = [27.1976, 0, 16.2564, 0, 23.5566, 0, 9.60064, 0, 12.2811, 0,
		 16.3955, 0, 13.6815, 0, 0.408964, 0, 4.46461, 0, 12.0228, 0, 7.18286, 0,
		 5.88196, 0, 10.0249, 0, 3.83276, 0, 10.0125, 0, 1.3949, 0, 11.4326, 0,
		 2.11597, 0, 4.50059, 0, 8.2486, 0, 3.79871, 0, 15.123, 0, 9.02716, 0,
		 10.5584, 0, 17.3638, 0, 5.19686, 0, 11.0076, 0, 8.92813, 0, 13.9632, 0, 
		 11.466, 0, 14.1743, 0, 1.06765, 0, 9.55211, 0, 19.4067, 0, 16.6116, 0,
		 20.7735, 0, 2.08838, 0, 18.9789, 0, 14.9563, 0, 4.00338, 0, 21.7476, 0,
		 6.14579, 0, 14.4046, 0, 8.75238, 0, -3.07091, 0, -2.07081, 0, 11.1903, 0,
		 21.5998, 0, 12.6512, 0, 4.90017, 0];
		 
		var myLength = testdata.length/2;
		
		var hypothesis = 9;
		var result = tTest(testdata, hypothesis);
		console.log("\nStudent-t test");
		string = "Cannot REJECT Ho at ";
		string += 100*alpha;
		string += "% level if P(t | dof) > ";
		string += (alpha + ".");
		console.log(string);
		console.log("tTest 1 - total N, mean, sigma, [dof, t, P(dof, t]): ", myLength,round(mean(testdata)[0],rounding),round(stdev(testdata)[0],rounding),round(result,rounding));
		if(result[2] >= alpha)
			{
			string = "The null hypothesis that the mean of the population\nis equal to ";
			string += hypothesis;
			string += " is NOT rejected at the ";
			string += 100*alpha;
			string += "% level based on\nthe Student-t test where dof = ";
			string += round(result[0],rounding);
			string += " and t = ";
			string += round(result[1],rounding);
			string += ". Note\nthe calculated probability = P(t | dof) = ";
			string += (round(result[2],rounding) + ".");
			console.log(string);
			}
		else
			{
			string = "The null hypothesis that the mean of the population\nis equal to ";
			string += hypothesis;
			string += " is REJECTED at the ";
			string += 100*alpha;
			string += "% level based on\nthe Student-t test where dof = ";
			string += round(result[0],rounding);
			string += " and t = ";
			string += round(result[1],rounding);
			string += ". Note\nthe calculated probability = P(t | dof) = ";
			string += round(result[2],rounding);
			console.log(string);
			};


		myLength = 1000;
		var testdataG1 = Array(2*myLength);
		mu = 0;
		var sigma = 1;
		for (var i = 0; i < 2*myLength; i += 2)
			{
			testdataG1[i] = randomGaussian(mu,sigma);
			testdataG1[i+1] = 0;
			}

		alpha = 0.01;
		hypothesis = 0.1;
		result = tTest(testdataG1, hypothesis);
		console.log("\nStudent-t test");
		string = "Cannot REJECT Ho at ";
		string += 100*alpha;
		string += "% level if P(t | dof) > ";
		string += (alpha + ".");
		console.log(string);
		console.log("tTest 1 - total N, mean, sigma, [dof, t, P(dof, t]): ", myLength,round(mean(testdata)[0],rounding),round(stdev(testdata)[0],rounding),round(result,rounding));
		if(result[2] >= alpha)
			{
			string = "The null hypothesis that the mean of the population\nis equal to ";
			string += hypothesis;
			string += " is NOT rejected at the ";
			string += 100*alpha;
			string += "% level based on\nthe Student-t test where dof = ";
			string += round(result[0],rounding);
			string += " and t = ";
			string += round(result[1],rounding);
			string += ". Note\nthe calculated probability = P(t | dof) = ";
			string += (round(result[2],rounding) + ".");
			console.log(string);
			}
		else
			{
			string = "The null hypothesis that the mean of the population\nis equal to ";
			string += hypothesis;
			string += " is REJECTED at the ";
			string += 100*alpha;
			string += "% level based on\nthe Student-t test where dof = ";
			string += round(result[0],rounding);
			string += " and t = ";
			string += round(result[1],rounding);
			string += ". Note\nthe calculated probability = P(t | dof) = ";
			string += round(result[2],rounding);
			console.log(string);
			};


		console.log("\nStop testTTest");
		console.log("");
		};
	};


// ****************************************************************
// Based upon Mathematica's SignTest function. There is a  hypothesis that
// data collected are from a certain probability distribution with median mu.
// This version is a "one-sample test concerning means". The test statistic
// is assumed to follow a BinomialDistribution[n,1/2] where n is the number
// of elements in the data not equal to 0.
//
// The Sign statistic uses the incomplete beta function as described in
// Numerical Recipes 6.4, "Cumulative Binomial Probability Distribution"
// which points to the procedures "betai" and "betacf".
//
// Checked Tuesday, 28 February 2017


function signTest(data, med)
	{
	'use strict';
	// Estimate mean & standard deviation from the sample data
	var middle = median(data)[0];  // ML estimate of average
	var n = data.length/2;  // number of samples
	
	var nPos = 0
	var nNeg = 0;
	for (var i = 0; i < data.length; i += 2)
		{
		if (data[i] > med)
			++nPos;  // how many values above the median
		else if (data[i] < med)
			++nNeg;  // how many values below the median
		else break;
		};
	
	var n = nPos + nNeg;
	var nMax = Math.max(nPos,nNeg);  // See https://en.wikipedia.org/wiki/Sign_test
	
	// We want probability > alpha if the null hypothesis
	// (median ≈ med) is not to be rejected at the alpha percent level.
	// See https://en.wikipedia.org/wiki/Sign_test (again)
	var probS = 2*betai(nMax, n-nMax+1, 1/2);

	return [n, nMax, probS]
	};

function testSignTest(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testSignTest");
		
		var alpha = 0.05;  // significance level
		var string = [];
		
		// test betai and betacf,
		var dof = 20;
		var a = dof/2;
		var b = dof - a + 1;
		var x = 1/2;
		var y = betai(a,b,x);
		console.log("\nbetai - dof,a,b,x,betai(a,b,x): ", dof,a,b,round(x,rounding),y, round(y,rounding));

		// These data (N = 20) are from a Gaussian random generator with
		// mu = 0 and sigma = 1. The result from Mathematica's
		// SignTest = 0.503445. Again, the higher this value the more likely
		// that the hypothesis that mu = 0 cannot be rejected.
		var mu = 0;
		var testdata = 
			[
			  1.10771, 0, 0.630657, 0, -1.08356, 0, -0.479773, 0, -0.691762, 0,
			 -0.727665, 0, -1.90394, 0, -0.541428, 0, -0.193613, 0, 1.11772, 0,
			 -0.528942, 0, -1.8149, 0, 1.06851, 0, 0.512378, 0, 0.686986, 0,
			 -0.321068, 0, -0.658884, 0, 1.02287, 0, -0.594719, 0, 0.126128, 0
			];
		 var myLength = testdata.length/2;
		
		var hypothesis = 0;
		var result = signTest(testdata, hypothesis);
		console.log("\nSign test");
		string = "Cannot REJECT Ho at ";
		string += 100*alpha;
		string += "% level if P(t | dof) > ";
		string += (alpha + ".");
		console.log(string);
		console.log("signTest 1 - total N, mean, sigma, [N, nMax, P(Sign)], P(data), true value): ", myLength,round(mean(testdata)[0],rounding),round(stdev(testdata)[0],rounding),round(result,rounding),result[2],0.503445);
		if(result[2] >= alpha)
			{
			string = "The null hypothesis that the median of the population\nis equal to ";
			string += hypothesis;
			string += " is NOT rejected at the ";
			string += 100*alpha;
			string += "% level based on\nthe Sign test where N = ";
			string += round(result[0],rounding);
			string += " and nMax  = ";
			string += round(result[1],rounding);
			string += ". Note\nthe calculated probability = P(Sign) = ";
			string += (round(result[2],rounding) + ".");
			console.log(string);
			}
		else
			{
			string = "The null hypothesis that the median of the population\nis equal to ";
			string += hypothesis;
			string += " is REJECTED at the ";
			string += 100*alpha;
			string += "% level based on\nthe Sign test where N = ";
			string += round(result[0],rounding);
			string += " and nMax  = ";
			string += round(result[1],rounding);
			string += ". Note\nthe calculated probability = P(Sign) = ";
			string += round(result[2],rounding);
			console.log(string);
			};


		myLength = 1000;
		var testdataG1 = Array(2*myLength);
		mu = 0;
		var sigma = 1;
		for (var i = 0; i < 2*myLength; i += 2)
			{
			testdataG1[i] = randomGaussian(mu,sigma);
			testdataG1[i+1] = 0;
			}
		
		alpha = 0.01;
		hypothesis = 0.07;
		result = signTest(testdataG1, hypothesis);
		console.log("\nSign test");
		string = "Cannot REJECT Ho at ";
		string += 100*alpha;
		string += "% level if P(t | dof) > ";
		string += (alpha + ".");
		console.log(string);
		console.log("signTest 1 - total N, mean, sigma, [N, nMax, P(Sign)]): ", myLength,round(mean(testdata)[0],rounding),round(stdev(testdata)[0],rounding),round(result,rounding));
		if(result[2] >= alpha)
			{
			string = "The null hypothesis that the median of the population\nis equal to ";
			string += hypothesis;
			string += " is NOT rejected at the ";
			string += 100*alpha;
			string += "% level based on\nthe Sign test where N = ";
			string += round(result[0],rounding);
			string += " and nMax  = ";
			string += round(result[1],rounding);
			string += ". Note\nthe calculated probability = P(Sign) = ";
			string += (round(result[2],rounding) + ".");
			console.log(string);
			}
		else
			{
			string = "The null hypothesis that the median of the population\nis equal to ";
			string += hypothesis;
			string += " is REJECTED at the ";
			string += 100*alpha;
			string += "% level based on\nthe Sign test where N = ";
			string += round(result[0],rounding);
			string += " and nMax  = ";
			string += round(result[1],rounding);
			string += ". Note\nthe calculated probability = P(Sign) = ";
			string += round(result[2],rounding);
			console.log(string);
			};


		console.log("\nStop testSignTest");
		console.log("");
		};
	};


// ****************************************************************
// Based upon Mathematica's KolmogorovSmirnovTest function. There is a 
// hypothesis that data collected is from a certain continuous (!) probability
// density function. We examine the absolute distance between the actual data
// cumulative and the cumulative from the hypothesized distribution
// Based upon procedure found in Knuth II Section 3.3.1B and Numerical Recipes
// Second in C, Section 14.3 
//
// Checked Wednesday, 8 February 2017

function cdfGauss(x, mean, std)
	{
	return 0.5 * (1 + erf((x - mean) / (std*Math.SQRT2)));
	};

function erf(x)
	{
	// save the sign of x
	var sign = (x >= 0) ? 1 : -1;
	x = Math.abs(x);

	// constants
	var a1 = 0.254829592;
	var a2 = -0.284496736;
	var a3 = 1.421413741;
	var a4 = -1.453152027;
	var a5 = 1.061405429;
	var p	= 0.3275911;

	// Abramowitz & Stegun formula 7.1.26
	var t = 1.0/(1.0 + p*x);
	var y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);
	return sign * y; // erf(-x) = -erf(x);
	};
	
function cdfLaplace(x, m, b)
	{
	if (x === m)
		y = 1/2;
	else if (x < m)
		y = Math.exp((x-m)/b)/2;
	else
		y = 1 - Math.exp(-(x-m)/b)/2; 
	return y;
	};

function cdfUniform(x, lower, upper)
	{
	return (x-lower)/(upper - lower);
	};

function probks(alam)
	{
	var ePS1 = 0.001;
	var ePS2 = 1e-8;
	var iter = 100;
	
	var fac = 2.0;
	var sum = 0.0;
	
	var term = 0.0;
	var termbf = 0.0;
	
	var a2 = -2.0*alam*alam;
	for (var j = 1; j <= iter; j++)
		{
		term = fac * Math.exp(a2*j*j);
		sum += term;
		if (Math.abs(term) <= ePS1*termbf || Math.abs(term) <= ePS2*sum)
			return sum;
		fac = -fac;  // alternating signs in sum
		termbf = Math.abs(term);
		};
	return 1.0;  // only if failing to converge
	};

// <https://en.wikipedia.org/wiki/Kolmogorov–Smirnov_test#Test_with_estimated_parameters>
// If either the form or the parameters of F(x) are determined from the data, then
// the critical values determined in this way are invalid. In such cases,
// Monte Carlo or other methods may be required, but tables have been prepared.... 
// this means kstwo (ks2GoodnessOfFit) instead of ksone
//
//
// Checked Wednesday, 8 February 2017

function ksone(data, distributionType)
	{
	var fo = 0.0;
	var ff = 0.0;
	var fn = 0.0;
	var dt = 0;
	
	var sorted = sort(data);
	var n = sorted.length/2;
	var d = 0.0;
	
	// Use the hypothesized distribution.
	switch(distributionType)
		{
		case "Gauss":
			var avgData = mean(data)[0];  // ML estimate of mean (Gauss)
			var stdData = stdev(data)[0];  // ML estimate of sigma (Gauss)
			
			for (var j = 1; j <= n; j++)
				{
				fn = j / n;
				ff = cdfGauss(sorted[2*(j-1)], avgData, stdData);
				dt = Math.max(Math.abs(fo-ff),Math.abs(fn-ff));
				if (dt > d) d = dt;
				fo = fn;
				}
			break;
	// See <https://en.wikipedia.org/wiki/Laplace_distribution#Parameter_estimation>
		case "Laplace":
			var middleEst = median(data)[0];  // ML estimate of middle (Laplace)
			var ttotal = 0;
			for (var i = 0; i < n; i++)
				ttotal += Math.abs(data[2*i] - middleEst);
			var bEst = ttotal / n;  // ML estimate of width (Laplace)
			
			for (var j = 1; j <= n; j++)
				{
				fn = j / n;
				ff = cdfLaplace(sorted[2*(j-1)], middleEst, bEst);
				dt = Math.max(Math.abs(fo-ff),Math.abs(fn-ff));
				if (dt > d) d = dt;
				fo = fn;
				}
			break;
		default: throw("Non-existent probability choice");
		};
		var corrKS = Math.sqrt(n) + 0.12 + (0.11/n);
		var prob = probks(corrKS*d);
		return [d, prob];
	};

function ks2GoodnessOfFit(data1, distributionType)
	{
	var fo = 0.0;
	var ff = 0.0;
	var fn = 0.0;
	var dt = 0;
	
	var sorted1 = sort(data1);
	var n1 = sorted1.length/2;
	var dataMC = Array(2*n1);  // MC = Monte Carlo method

	var d = 0.0;
	var prob = 0.0;
	
	// Use the hypothesized distribution.
	switch(distributionType)
		{
		case "Gauss":
			// generate a Gaussian data set with same length, mean, sigma
			var avgData = mean(data1)[0];  // ML estimate of mean (Gauss)
			var stdData = stdev(data1)[0];  // ML estimate of sigma (Gauss)
			
			for (var i = 0; i < 2*n1; i += 2)  // generate MC data
				{
				dataMC[i] = randomGaussian(avgData,stdData);
				dataMC[i+1] = 0;
				};

				
			var sorted2 = sort(dataMC);
			var n2 = n1;
			var j1 = 1;
			var j2 = 1;
			var fn1 = 0.0;
			var fn2 = 0.0;
			var d1 = 0;
			var d2 = 0;
			
			while (j1 <= n1 && j2 <= n2)
				{
				d1 = sorted1[2*(j1-1)];
				d2 = sorted2[2*(j2-1)];
				if (d1 <= d2)
					{
					fn1 = j1 / n1;
					j1 += 1;
					};
				if (d2 <= d1)
					{
					fn2 = j2 / n2;
					j2 += 1;
					};
				dt = Math.abs(fn2 - fn1);
				if (dt > d)
					d = dt;
				};
				
			var n = Math.sqrt((n1*n2)/(n1+n2));
			var corrKS = Math.sqrt(n) + 0.12 + (0.11/n);
			prob = probks(corrKS*d);
			break;
			
	// See <https://en.wikipedia.org/wiki/Laplace_distribution#Parameter_estimation>
		case "Laplace":
			// generate a Laplace data set with same length, middle, width
			var middleEst = median(data1)[0];  // ML estimate of middle (Laplace)
			var ttotal = 0;
			for (var i = 0; i < 2*n1; i += 2)
				ttotal += Math.abs(data1[i] - middleEst);
			var bEst = ttotal / n1;  // ML estimate of width (Laplace)
			
			for (var i = 0; i < 2*n1; i += 2)  // generate MC data
				{
				dataMC[i] = randomLaplace(middleEst,bEst);
				dataMC[i+1] = 0;
				};
				
			var sorted2 = sort(dataMC);
			var n2 = n1;
			var j1 = 1;
			var j2 = 1;
			var fn1 = 0.0;
			var fn2 = 0.0;
			var d1 = 0;
			var d2 = 0;
			var d = 0.0;
			
			while (j1 <= n1 && j2 <= n2)
				{
				d1 = sorted1[j1-1];
				d2 = sorted2[j2-1];
				if (d1 <= d2)
					{
					fn1 = j1 / n1;
					j1 += 1;
					};
				if (d2 <= d1)
					{
					fn2 = j2 / n2;
					j2 += 1;
					};
				dt = Math.abs(fn2 - fn1);
				if (dt > d)
					d = dt;
				};
				
			var n = Math.sqrt((n1*n2)/(n2+n1));
			var corrKS = Math.sqrt(n) + 0.12 + (0.11/n);
			prob = probks(corrKS*d);
			break;
		default: throw("Non-existent probability choice");
		};
	return [d, prob];
	};


function testKS(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testKS");
		
		// Gaussian process data
		var myLengthG = 1000;
		var testdataG = Array(2*myLengthG);
		var mu = 10;
		var sigma = 7;
		
		for (var i = 0; i < 2*myLengthG; i += 2)
			{
			testdataG[i] = randomGaussian(mu,sigma);
			testdataG[i+1] = 0;
			};
		// Use ks2GoodnessOfFit test because we do not know the process parameters
		var result = ks2GoodnessOfFit(testdataG, "Gauss");
		console.log("\nREJECT Ho at 5% level if prob < 0.95");
		console.log("Gauss data & model KS 2 - N, mean, stdev, [d, prob]: ", myLengthG, round(mean(testdataG)[0],rounding),round(stdev(testdataG)[0],rounding), round(result,rounding));


		// Laplace process data
		var myLengthL = 1000;
		var testdataL = Array(2*myLengthL);
		var middle = 73;
		var width = 5;
		for (var i = 0; i < 2*myLengthL; i += 2)
			{
			testdataL[i] = randomLaplace(middle,width);
			testdataL[i+1] = 0;
			};

		// Use ks2GoodnessOfFit test because we do not know the process parameters
		var result = ks2GoodnessOfFit(testdataL, "Laplace");
		console.log("\nREJECT Ho at 5% level if prob < 0.95");
		console.log("Laplace data & model KS 2 - N, median, width, [d, prob]: ", myLengthL, round(median(testdataL)[0],rounding), round(stdev(testdataL)[0]*Math.SQRT1_2,rounding), round(result,rounding));


		// Exponential process data
		var myLengthE = 1000;
		var testdataE = Array(2*myLengthE);
		var lambda = 5;
		for (var i = 0; i < 2*myLengthE; i += 2)
			{
			testdataE[i] = randomExponential(lambda);
			testdataE[i+1] = 0;
			};

		// Use ks2GoodnessOfFit test because we do not know the process parameters
		var result = ks2GoodnessOfFit(testdataE, "Gauss");
		console.log("\nREJECT Ho at 5% level if prob < 0.95");
		console.log("Exponential data & Gauss model KS 2 - N, mean, stdev, [d, prob]: ", myLengthE,round(mean(testdataE)[0],rounding),round(stdev(testdataE)[0],rounding),round(result,rounding));


		console.log("Stop");
		console.log("");
		};
	};

