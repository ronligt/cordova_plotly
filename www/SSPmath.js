// ****************************************************************
// Math, array, and signal processing procedures to be used in the SSP package
// i.t. young
// Tuesday, 5 July 2016
// Wednesday, 21 September 2016: separated fourier code into another file


// ****************************************************************
// checked Sunday, 18 December 2016

function oddQ(number)
	{
	return (number % 2 === 1);
	};

function evenQ(number)
	{
	return (number % 2 === 0);
	};

function testOddQ(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testOddQ");
		
		console.log("oddQ: ",35,oddQ(35));
		console.log("oddQ: ",42,oddQ(42));
		
		console.log("evenQ: ",35,evenQ(35));
		console.log("evenQ: ",42,evenQ(42));
		
		console.log("Stop");
		console.log("");
		};
	};


// ****************************************************************
// The following utility rotates a Fourier spectrum 1D or 2D, so that
// frequency f = 0 is in the middle of the display. If the data are 1D
// and of odd length N, then the rotation is floor(N/2) which is
// exactly in the middle. If the data are of even length then the
// rotation is N/2 which equals floor (N/2).
// If the data are 2D then both dimensions must be a power of two so
// the rotation is N/2 = floor(N/2).
//
// checked Sunday, 14 May 2017

function spectralRotate(data,nn)
	{
	var ndim = nn.length;  // number of dimensions, 1 or 2
	var result = [];
	var tempC = 0;
	if (ndim === 1)
		{
		var lrotate = Math.floor(nn[0]/2);  // for 1D signal use this
		result = data.slice();
		
		for (var k = 0; k < lrotate; k++)
			{
			tempC = result.shift();
			result.push(tempC);
			};
		}
	else if (ndim === 2)
		{
		var rows = nn[0];  // for 2D signal use this
		var cols = nn[1];
		var tempR = [];
	
		result = data.slice();
		// first, rotate rows down 
		for (var i = 0; i < (rows/2); i++)
			{
			tempR = result.shift();
			result.push(tempR);
			};
		 
		// then, rotate columns around 
		for (var i = 0; i < rows; i++)
			{
			tempR = result[i];
			for (var j = 0; j < (columns/2); j++)
				{
				tempC = tempR.shift();
				tempR.push(tempC);
				};
			result[i] = tempR;
			};
		}
	else
		throw "Wrong number of dimensions!";
		
	return result;
	};


// ****************************************************************
// Following are utility routines to convert between formats for 
// complex signals that are either 1D or 2D. This routine was not written
// to handle a general number of dimensions, N.
//
// the "toStandardDisplayFormat" routine converts 1D and 2D signals to the
// special format required by the Plotly.js graphic display procedures
//
// The "split" routine takes data in the form of [r1,i1,r2,i2, ...,rn,in] 
// and converts it to [[r1,r2, ...,rn],[i1,i2, ...,in]] where r1 is 
// the real part of a complex sample and i1 is the imaginary part.Note 
// that this is a 1D example.
//
// The "merge" routine takes data in the form of 
// [[r1,r2, ...,rn],[i1,i2, ...,in]] and converts it to 
// [r1,i1,r2,i2, ...,rn,in] where r1 is the real part of a complex sample 
// and i1 is the imaginary part. Note that this is a 1D example.
//
// checked Wednesday, 17 May 2017

function toStandardDisplayFormat(data,nn)
	{
	// data = vector of even length with alternating real and imaginary values
	// if ndim = 2, nn[rows, columns] is an integer array containing the number of rows
	// and the number of columns 
	
	var standardForm = uRs(data);
	var ndim = nn.length; 	// number of dimensions, 1 or 2
	
	var rows = 1;  			// default is 1D signal
	var cols = 0;
	
	if (ndim === 1)
		{
		cols = nn[0];
		var result = standardForm;
		}
	else if (ndim === 2)
		{
		rows = nn[0];		// but for 2D signal use this
		cols = nn[1];
		var result = createArray(rows,columns);
		for (var i = 0; i < rows; i++)
			for (var j = 0; j < columns; j++)
				result[i][j] = standardForm[i*columns + j];
		}
	else
		throw "Wrong number of dimensions!";
	
	return result;
	};
	
	
function split(data,nn)
	{
	// data = vector of even length with alternating real and imaginary values
	// if ndim = 2, nn[rows, columns] is an integer array containing the number of rows
	// and the number of columns 
	
	var ndim = nn.length;  // number of dimensions, 1 or 2
	var rows = 1;  // default is 1D signal
	var cols = 0;
	if (ndim === 1)
		cols = nn[0];
	else if (ndim === 2)
		{
		rows = nn[0];  // but for 2D signal use this
		cols = nn[1];
		}
	else
		throw "Wrong number of dimensions!";

	var result = createArray(2*rows,cols);
	var m = 0;
	for (var i = 0; i < rows; i++)
		{
		m = 2*i;  // for speed
		var k = 0;
		for (var j = 0; j < cols; j++)
			{
			// It's just a matter of getting the right indices
			k = m*cols + 2*j;  // for speed
			result[m][j] = data[k];
			result[m+1][j] = data[k+1];
			};
		};
	return result;
	};
	
	
function merge(data,nn)
	{
	var ndim = nn.length;
	var rows = 1;  // default is 1D signal
	var cols = 0;
	if (ndim === 1)
		cols = nn[0];
	else if (ndim === 2)
		{
		rows = nn[0];  // but for 2D signal use this
		cols = nn[1];
		}
	else
		throw "Wrong number of dimensions!";
	var result = createArray(2*rows*cols);
		var m = 0;
	for (var i = 0; i < rows; i++)
		{
		m = 2*i;  // for speed
		var k = 0;
		for (var j = 0; j < cols; j++)
			{
			// It's just a matter of getting the right indices
			k = m*cols + 2*j;  // for speed
			result[k] = data[m][j];
			result[k+1] = data[m+1][j];
			};
		};
	return result;
	};


function testMergeSplitSignal(testChoice)
	{
		if(testChoice)
		{
		console.log("\nStart testMergeSplitSignal");
	
		var myLength = 8; // number of complex samples
		var real = createArray(myLength);
		var imag = createArray(myLength);
		var data1D = createArray(2*myLength);
	
		for(var i = 0; i < myLength; i++)
			{
			real[i] = i; imag[i] = i + myLength;
			};
	
		var rcopy = real.slice();
		var icopy = imag.slice();
		var data = [rcopy,icopy];

		var dims = [myLength];
	
		var result = merge(data,dims);
		console.log("merge 1D: ",data, dims, result);
	
		var nextResult = split(result,dims);
		console.log("split 1D: ",result, dims, nextResult);
	
		var p1 = [1,2,3,4,5,6,7,8];  // row 1 reals
		var q1 = [9,10,11,12,13,14,15,16];  // row 1 imags
		var p2 = [17,18,19,20,21,22,23,24];  // row 2 reals
		var q2 = [25,26,27,28,29,30,31,32];  // row 2 imags
		var p3 = [33,34,35,36,37,38,39,40];  // row 3 reals
		var q3 = [41,42,43,44,45,46,47,48];  // row 3 imags
		var p4 = [49,50,51,52,53,54,55,56];  // row 4 reals
		var q4 = [57,58,59,60,61,62,63,64];  // row 4 imags

		var imagex =  // 4 rows x 8 columns
			[
				p1, q1,
				p2, q2,
				p3, q3,
				p4, q4
			];
	
	/*
		p1 = [1,2,3,4];  // row 1 reals
		q1 = [9,10,11,12];  // row 1 imags
		p2 = [17,18,19,20];  // row 2 reals
		q2 = [25,26,27,28];  // row 2 imags

		imagex =  // 2 rows x 4 columns
			[
				p1, q1,
				p2, q2,
			];
	*/
		var tdims = arrayDimension(imagex)
		var dims = [tdims[0]/2,tdims[1]];
	
		var result = merge(imagex,dims);
		console.log("merge 2D: ",imagex,dims, result);
	
		var nextResult = split(result,dims);
		console.log("split 2D: ",result, dims, nextResult);
	
		console.log("Stop Merge & Split");
		console.log("");
		};
	};


// ****************************************************************
// This routine is intended for 1D signals where we want to be able to do
// either convolution or correlation of complex data structures
// checked Sunday, 18 December 2016

function reverse(list)
	{
	var nn = [(list.length)/2];
	var splitData = split(list,nn);
	return merge([splitData[0].reverse(), splitData[1].reverse()],nn)
	};

function testReverse(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testReverse");
		
		var list = ['r1','i1','r2','i2','r3','i3','r4','i4'];
		console.log("reverse start & finish: ",list,reverse(list));
		console.log("Stop testReverse");
		console.log("");
		};
	};
	
	

// ****************************************************************
// Don't forgot that one might want to use getMaxOfArray(abs(cmplxArray)),
// etc.
// checked Tuesday, 10 January 2017

function getMaxOfArray(numArray)
	{
	return Math.max.apply(null, flatten(numArray));
	};


function getMinOfArray(numArray)
	{
	return Math.min.apply(null, flatten(numArray));
	};


function testMinMax(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testMinMax");
		var myLength = 10
		var tlist = Array(myLength); /* temporary list */
		for (var i=0; i < myLength; i++)
			{
			tlist[i] = randomInteger(0,100);
			};
			
		console.log("getMinOfArray List: ",tlist);
		
		var minval = getMinOfArray(tlist);
		console.log("getMinOfArray: ",minval);
	
		var maxval = getMaxOfArray(tlist);
		console.log("getMaxOfArray: ",maxval);
		
		console.log("Stop");
		console.log("");
		};
	};


// ****************************************************************
// checked Sunday, 18 December 2016

function total(data)
	{
	var rTotal = iTotal = 0;
	var n = data.length;
	for(var i = 0; i < n; i += 2)
		{
		rTotal += data[i];  // real part
		iTotal += data[i+1];  // imaginary part
		};
		return [rTotal, iTotal];
	};


function testTotal(testChoice)
	{
	if(testChoice)
		{
		var myLength = 5
		var tlist = createArray(2*myLength); /* temporary list */
		for (var i = 0; i < 2*myLength; i++)
			{
			tlist[i] = randomInteger(1,10);
			};
		
		console.log("\nStart testTotal");
		
		console.log("total input & output: ",tlist,total(tlist));
		
		console.log("Stop testTotal");
		console.log("");
		};
	};


/// ****************************************************************
// Inner product of two arrays (vectors)
// Both 1D and 2D data have the same complex format.
// result for z1 and z2 is z1•conjugate(z2)
// checked Sunday, 18 December 2016

function dot(z1,z2)
	{
	var d1 = dimensions(z1);
	var d2 = dimensions(z2);
	
	if(d1 !== d2)
		throw "Two arrays must have the same length!";

	var result = 0;

	if(d1 === 0)
		result = [(z1*z2),0]; // these are real scalars
	else
		{
		var nn = [(z1.length)/2];
		splitz1 = split(z1,nn);
		var real1 = splitz1[0];
		var imag1 = splitz1[1];
		splitz2 = split(z2,nn);
		var real2 = splitz2[0];
		var imag2 = splitz2[1];
		
		var rsum = 0;
		var isum = 0;
		
		for (var i = 0; i < real1.length; i++)
			{
			var a = real1[i];
			var b = imag1[i];
			var c = real2[i];
			var d = imag2[i];

			rsum += (a*c + b*d);  // conjugate occurs here
			isum += (b*c - a*d);  // conjugate occurs here
			};
		result = [rsum,isum];
		};
	return result; 
	};


function testDot(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testDot");
		console.log("scalar: ",
			dimensions(tdata0),arrayDimension(tdata0),tdata0,dot(tdata0,tdata0));
			
		// test signal
		var tdata1 = [1,2,1,2,1,2,1,2,1,2];
		var tdata2 = [1,0,1,0,1,0,1,0,1,0];
		console.log("complex: ",
			tdata1,dot(tdata1,tdata1),
			tdata1,dot(tdata1,tdata2),
			tdata2,dot(tdata2,tdata2));
			
		// two orthogonal vectors
		var tdata1 = [1,0,3,0,2,0];
		var tdata2 = [1,0,3,0,-5,0];
		console.log("complex: ",
			tdata1,dot(tdata1,tdata1),
			tdata1,dot(tdata1,tdata2),
			tdata2,dot(tdata2,tdata2));

		var myLength = 10000;
		var xlist = createArray(2*myLength); /* temporary list */
		var ylist = createArray(2*myLength); /* temporary list */

		for (var i = 0; i < 2*myLength; i++)
			{
			xlist[i] = randomInteger(0,1);
			ylist[i] = randomInteger(0,1);
			};

		var iterations  = 100;
		var start = performance.now();
		var value = 0;
		for (var i = 0; i < iterations; ++i) value = dot(xlist,ylist);
		var finish = performance.now(); // in ms
		var timed = finish - start;
		console.log("dot - Value, & Duration [ns per vector element]: ",
			value,1e6*timed/(myLength*iterations));

		console.log("Stop testDot");
		console.log("");
		};
	};


// ****************************************************************
// Real numeric data only
// Note you cannot sort complex data. Application to complex data will
// be meaningless.You can sort abs(complex data) or the real part of a complex array
// Checked Monday, 9 January 2017

function sort(data)
	{
	var nn = [(data.length)/2];
	var splitData = split(data,nn);
	var realData = splitData[0].slice();
	var sortData = realData.sort(function(a, b){return a-b});
	var result = [sortData, splitData[1]];
	return merge(result,nn);
	};


function testSort(testChoice)
	{
		if(testChoice)
		{
		console.log("\nStart testSort");

		var real = real2.slice();
		var nn = [real.length];
		var imag = createArray(real.length);
		for (var i = 0; i < real.length; i++)
			imag[i] = 0;
		var tdata = merge([real, imag],nn);
		console.log("sort: ",tdata,sort(tdata));
	

		var myLength = 10;
		var nn = [myLength];
		var treal = createArray(myLength); // temporary list
		var timag = createArray(myLength); // temporary list
		for (var i = 0; i < myLength; i++)
			{
			treal[i] = randomInteger(1,20);
			timag[i] = 0;
			};
		var tdata = merge([treal, timag],nn);
		console.log("sort random numbers: ",tdata,sort(tdata));

	
		console.log("Stop testSort");
		console.log("");
		};
	};


// ****************************************************************
// Checked Monday, 9 January 2017

function flatten(array)
	{
	var result = [];
	for(var i = 0; i < array.length; i++)
		{
		if(Array.isArray(array[i]))
			{
			result = result.concat(flatten(array[i]));
			}
		else
			{
			result.push(array[i]);
			};
		};
	return result;
	};

function testFlatten(testChoice)
	{
		if(testChoice)
		{
		console.log("\nStart testFlatten");
		
		console.log("flatten array (start): ",nlist3D);
		console.log("flatten result: ",nlist3D,flatten(nlist3D));
		console.log("flatten array (finish): ",nlist3D);
		
		console.log("Stop testFlatten");
		console.log("");
		};
	};


// ****************************************************************
// WARNING: Do not use var tarray2 = Array(tdimsi).fill(Array(tdimsj).fill(0));
// as pointers in the "tdimsi" positions point to the last place in the
// "tdimsj" section.

// checked Sunday, 18 December 2016

function createArray(length)
// Taken from <http://stackoverflow.com/questions/966225/how-can-i-create-a-two-dimensional-array-in-javascript/966938#966938>
// Use a createArray(rows,cols,...)
	{
	var arr = Array(length || 0),
	i = length;
	if (arguments.length > 1)
		{
		var args = Array.prototype.slice.call(arguments, 1);
		while(i--) arr[length-1 - i] = createArray.apply(this, args);
		};
	return arr;
	};

// Following Mathematica's protocol replace every element smallter than
// 10^(-10) with 0

// This chop() routine is overkill. All arrays that will be used are 1D
// so the multi-dimensional stuff is unnecessary. But it is working so
// I won't fix it.

function chop(array)
	{
	var dims = arrayDimension(array);
	var i = 0;
	var j = 0;
	var k = 0;
	
	if(dims.length === 1) // 1D real signal
		{
		var tdimsi = dims[0];
		var tarray1 = createArray(tdimsi);
		};
	if(dims.length === 2)  // 1D complex signal or 2D real signal
		{
		var tdimsi = dims[0];
		var tdimsj = dims[1];
		var tarray2 = createArray(tdimsi,tdimsj);
		};
	if(dims.length === 3)  // 2D complex signal
		{
		var tdimsi = dims[0];
		var tdimsj = dims[1];
		var tdimsk = dims[2];
		var tarray3 = createArray(tdimsi,tdimsj,tdimsk);
		};
	
	if(dims.length === 0)  // scalar
		{
		if(array < tooSmall)
			{
			array = 0;
			};
		return array;
		}
	else if(dims.length === 1)  // 1D real signal
		{
		for(var i = 0; i < tdimsi; i++)
			{
			if(Math.abs(array[i]) < tooSmall)
				{
				tarray1[i] = 0;
				}
			else
				{
				tarray1[i] = array[i];
				};
			};
		return tarray1;
		}
	else if(dims.length === 2)  // 1D complex signal or 2D real signal
		{
		for(var i = 0; i < tdimsi; i++)
			{
			for(var j = 0; j < tdimsj; j++)
				{
				if(Math.abs(array[i][j]) < tooSmall)
					{
					tarray2[i][j] = 0;
					}
				else
					{
					tarray2[i][j] = array[i][j];
					};
				};
			};
		return tarray2;
		}
	else  // 2D complex signal
		{
		for(var i = 0; i < tdimsi; i++)
			{
			for(var j = 0; j < tdimsj; j++)
				{
				for(var k = 0; k < tdimsk; k++)
					{
					if(Math.abs(array[i][j][k]) < tooSmall)
						{
						tarray3[i][j][k] = 0;
						}
					else
						{
						tarray3[i][j][k] = array[i][j][k];
						};
					};
				};
			};
		return tarray3;
		};
	};


function testChop(testChoice)
	{
		if(testChoice)
		{
		console.log("\nStart testChop");
		var data = tdata2D;
		console.log("chop: ",arrayDimension(data),data,chop(data));
		console.log("Stop testChop");
		console.log("");
		};
	};


// ****************************************************************
// Again the assumption is that we dealing with the "universal"
// signal format as described in "SSPconstants.js". This means that if
// we ask for "n" complex samples that means 2n number will be taken 
//
// The array represents 1D data then "n" can be positive or negative.
// Positive means that we take the first n value of the input array;
// negative means that we take the last n values of the input array.
//
// checked Sunday, 18 December 2016

function take(array,n)
	{
	if(n === 0)
		{
		// In this (special) case, we just copy the array. 
		var tarray = array.slice();
		}
	else
		{
		// indices matter, especially for speed
		var k = 2*n;
		var tarray = createArray(Math.abs(k));
		var m = array.length + k;
		if (n > 0)
			{
			for (var i = 0; i < k; i += 2)
				{
				tarray[i] = array[i];
				tarray[i+1] = array[i+1];
				};
			}
		else
			{
			for (var i = 0; i < Math.abs(k); i += 2)
				{
				tarray[i] = array[i+m];
				tarray[i+1] = array[i+m+1];
				};
			};
		};
	return tarray;
	};

function testTake(testChoice)
	{
		if(testChoice)
		{
		console.log("\nStart testTake");
		var tlist = nlist4;
		var pos = 0;
		console.log("take: ",tlist,pos,take(tlist,pos));
		pos = 3;
		console.log("take: ",tlist,pos,take(tlist,pos));
		pos = -3;
		console.log("take: ",tlist,pos,take(tlist,pos));
		pos = -1;
		console.log("take: ",tlist,pos,take(tlist,pos));
		console.log("Stop testTake");
		console.log("");
		};
	};


// ****************************************************************
// Again the assumption is that we dealing with the "universal"
// signal format as described in "SSPconstants.js". This means that if
// we ask for "n" complex samples that means 2n number will be taken 
//
// The array represents 2D (image) data with dimensions ddims = [rows, cols]
// the array window contains the portion of the image data that is to be
// taken with coordinates [rtop, rbot, cleft, cright]. Note rbot ≥ rtop and
// cright ≥ cleft.
//
// checked Saturday, 25 February 2017

function imageTake(data,ddims,window)
	{
	
	var rows = ddims[0];
	var cols = ddims[1];
	
	// Note rbot ≥ rtop
	var rtop = window[0];
	var rbot = window[1];
	var trows = rbot - rtop + 1;
	
	// Note cright ≥ cleft
	var cleft = window[2];
	var cright = window[3];
	var tcols = cright - cleft + 1;
	
	if( (trows === rows) && (tcols === cols) )
		{
		// In this (special) case, we just copy the array. 
		var tarray = data.slice();
		}
	else
		{
		var tarray0 = Array(2*trows*tcols);
		tarray0.fill(0);
		
		var tdata = split(data, [2*rows,cols]);
		var tarray1 = split(tarray0, [2*trows,tcols]);
		
		for (var i = 0; i < 2*trows; i++)
			for (var j = 0; j < 2*tcols; j++)
				tarray1[i][j] = tdata[2*rtop+i][j+cleft];
		var tarray = merge(tarray1,[2*trows,tcols]);
		};
	return [[trows,tcols], tarray ] ;
	};

function testImageTake(testChoice)
	{
		if(testChoice)
		{
		console.log("\nStart testImageTake");
		
		var rows = 7;
		var cols = 6;
		var idims = [rows,cols];
		
		var constantImage = createArray(2*rows,cols);
		for (var i = 0; i < cols; i++)
			for (var j = 0; j < 2*rows; j++)
				{
				if (evenQ(j))
					constantImage[j][i] = i;
				else
					constantImage[j][i] = (j-1)/2;
				};
		var inputImage = merge(constantImage, idims);
		
		var window = [1, 5, 2, 5];
		var result = imageTake(inputImage,idims,window);

		console.log("\nimageTake - idims, inputImage, window, result: ",idims, split(inputImage,idims),window,result[0],split(result[1],result[0]));
		
		var window = [0, rows-1, 0, cols-1];
		var result = imageTake(inputImage,idims,window);

		console.log("\nimageTake - the whole image");
		console.log("\nimageTake - idims, inputImage, window, result: ",idims, split(inputImage,idims),window,result[0],split(result[1],result[0]));
		
		console.log("\nimageTake - one row");
		var row1 = rows-3;
		var window = [row1, row1, 0, cols-1];
		var result = imageTake(inputImage,idims,window);

		console.log("\nimageTake - idims, inputImage, window, result: ",idims, split(inputImage,idims),window,result[0],split(result[1],result[0]));
		
		console.log("\nimageTake - one column");
		var col1 = cols-2;
		var window = [0, rows-1, col1, col1];
		var result = imageTake(inputImage,idims,window);

		console.log("\nimageTake - idims, inputImage, window, result: ",idims, split(inputImage,idims),window,result[0],split(result[1],result[0]));
		
		console.log("\nStop testImageTake");
		console.log("");
		};
	};


// ****************************************************************
// Make two arrays have the same length. That length will be the shorter
// of the two.
// checked Sunday, 18 December 2016

function normalizeLength(arrays)
	{
	var array0 = arrays[0].slice();
	var array1 = arrays[1].slice();
	var m0 = arrays[0].length;
	var m1 = arrays[1].length;
	var tarrays = createArray(2);

	if(m0 === m1)
		tarrays =[array0,array1];
	else
		{
		if (m0 > m1)
			tarrays = [take(array0,m1/2),array1.slice()];
		else
			tarrays = [array0.slice(),take(array1,m0/2)];
		};
	return tarrays;
	};

function testNormalizeLength(testChoice)
	{
		if(testChoice)
		{
		console.log("\nStart testNormalizeLength");
		
		var tlists = [nlist1,nlist2];
		console.log("normalizeLength initial: ",nlist1,nlist2);
		result = normalizeLength(tlists);
		console.log("normalizeLength result 1: ",result[0],result[1]);
		console.log("");
		
		tlists = [nlist3,alist1];
		console.log("normalizeLength initial: ",nlist3,alist1);
		result = normalizeLength(tlists);
		console.log("normalizeLength result 2: ",result[0],result[1]);
		console.log("");
		
		tlists = [alist1,nlist3];
		console.log("normalizeLength initial: ",alist1,nlist3);
		result = normalizeLength(tlists);
		console.log("normalizeLength result 3: ",result[0],result[1]);
		
		console.log("Stop testNormalizeLength");
		console.log("");
		};
	};


// ****************************************************************
// starting from a complex array A, return B = A/(max(|A|)
// checked Sunday, 18 December 2016

function ampNormalize(array)
	{
	var n = array.length;
	if (oddQ(n))
		throw("Complex signals must have even length, e.g. 2*17");

	// first find the absolute maximum
	var result = createArray(n);
	var maxval = 0;
	var temp = 0;
	for(var i = 0; i < n; i +=2)
		{
		abs[i] = Math.hypot(array[i],array[i+1]);
		temp = abs[i];
		if (temp > maxval)
			maxval = temp;
		};
		
	// now normalize except if maxval === 0
	if(maxval !== 0)
		{
			for(var i = 0; i < n; i += 2)
				{
					result[i] = array[i]/maxval;
					result[i+1] = array[i+1]/maxval;
				};
		};
	return result;
	};


function testAmpNormalize(testChoice)
	{
		if(testChoice)
		{
		console.log("\nStart testAmpNormalize");

		var z0 = [1,0,0,1,3,4];
		var z2 = merge([real2, imag2],[real2.length]);

		var tlist = ampNormalize(z0);
		console.log("ampNormalize input & output: ",
			round(z0,rounding),round(tlist,rounding));
		
		var tlist = ampNormalize(z2);
		console.log("ampNormalize input & output: ",
			round(z2,rounding),round(tlist,rounding));
	
		console.log("Stop testAmpNormalize");
		console.log("");
		};
	};


// ****************************************************************
// complex 1D array padded to a length that is closest poewer of two and normalized
// such that if A is input then B = A/max(|A|)
// checked Sunday, 18 December 2016

function padAndNormalize(array)
	{
	var n = array.length;
	if (oddQ(n))
		throw("Complex signals must have even length, e.g. 2*17");
	var magic = Math.ceil(Math.log2(n));
	// this number is even
	var samples  = Math.pow(2,magic);
	var iter = samples - n;
		
	var tarray = [];
	for(var i = 0; i < iter; i++)
		{
		tarray[i] = 0;
		};
		
	// This is javascript version of Mathemtica's Join
	var padded = array.concat(tarray);
	
	return ampNormalize(padded);
	};

function testPadAndNormalize(testChoice)
	{
		if(testChoice)
		{
		console.log("\nStart testPadAndNormalize");

		var z0 = [1,0,0,1,3,4];
		var output = padAndNormalize(z0);
		console.log("padAndNormalize input & output lengths & lists:",
			z0.length,z0,output.length,output);
			
		var z1 = z0.concat(z0);
		var output = padAndNormalize(z1);
		console.log("padAndNormalize input & output lengths & lists:",
			z1.length,z1,output.length,output);
			
		console.log("Stop testPadAndNormalize");
		console.log("");
		};
	};


// ****************************************************************
// simple functions for complex arrays
// complex array is z = [real[],imag[]]
// checked Sunday, 18 December 2016

// this functions converts a 'universal' function to a standard format
// example complex 'universal' = [r1, i1, r2, i2, ..., rN, iN]
// example real 'universal' = [r1, 0, r2, 0, ..., rN, 0]
// example imaginary 'universal' = [0, i1, 0, i2, ..., 0, iN]
// example real standard = [r1, r2, ..., rN]
// example imaginary standard = [i1, i2, ..., iN]
//
// the standard format is required by plotting routines so conversion
// is necessary. It's a dirty job but someone ...

// universal real to standard real
function uRs(data)
	{
	var n = data.length;
	var real = createArray(n/2);
	for(var i = 0; i < n; i += 2)  real[i/2] = data[i];
	return real;
	};


// universal imaginary to standard imaginary
function uIs(data)
	{
	var n = data.length;
	var imag = createArray(n/2);
	for(var i = 0; i < n; i += 2)  imag[i/2] = data[i+1];
	return imag;
	};


// standard real to universal real
function suR(data)
	{
	var n = 2*data.length;
	var real = createArray(n);
	for(var i = 0; i < n; i += 2)
		{
		real[i] = data[i/2];
		real[i+1] = 0;
		};
	return real;
	};


// standard imaginary to universal imaginary
function suI(data)
	{
	var n = 2*data.length;
	var imag = createArray(n);
	for(var i = 0; i < n; i += 2)
		{
		imag[i] = 0;
		imag[i+1] = data[i/2];
		};
	return imag;
	};


// these functions work with the 'universal' format
function re(data)
	{
	var n = data.length;
	var real = createArray(n);
	for(var i = 0; i < n; i += 2)
		{
		real[i] = data[i];
		real[i+1] = 0;
		};
	return real;
	};


function im(data)
	{
	var n = data.length;
	var imag = createArray(n);
	for(var i = 0; i < n; i += 2)
		{
		imag[i] = 0;
		imag[i+1] = data[i+1];
		};
	return imag;
	};


function arg(data)
	{
	var n = data.length;
	var phase = createArray(n);
	for(var i = 0; i < n; i += 2)
		{
		phase[i] = Math.atan2(data[i+1],data[i]);
		phase[i+1] = 0;
		};
	return phase;
	};


function abs(data)
	{
	var n = data.length;
	var mag = createArray(n);
	for(var i = 0; i < n; i += 2)
		{
		mag[i] = Math.hypot(data[i],data[i+1]);
		mag[i+1] = 0;
		};
	return mag;
	};


function abssq(data)
	{ // absolute value SQUARED
	var n = data.length;
	var mag2 = createArray(n);
	for(var i = 0; i < n; i += 2)
		{
		mag2[i] = (data[i]**2) + (data[i+1]**2);
		mag2[i+1] = 0;
		};
	return mag2;
	};


function conjugate(z)
	{
	var n = z.length;
	var zstar = createArray(n);
	for(var i = 0; i < n; i += 2)
		{
		zstar[i] = z[i];
		zstar[i+1] = -z[i+1];
		};
	return zstar;
	};


// Converts list of points in the (r, theta) coordinates to (x,y) coordinates
// Note in 2D, r ≥ 0 and -π ≤ theta ≤ +π
// Input data structure polarArray is [r1,theta1,r2,theta2,r3,theta3...,rN,thetaN]
// Output data structure xyArray is [x1,y1,x2,y2,x3,y3...,xN,yN]

function fromPolarCoordinates(polarArray)
	{
	var n = polarArray.length;
	var xyArray = createArray(n);

	for(var i = 0; i < n; i += 2)
		{
		xyArray[i] = polarArray[i]*Math.cos(polarArray[i+1]);
		xyArray[i+1] = polarArray[i]*Math.sin(polarArray[i+1]);
		};
	return xyArray;
	};


// Converts list of points in (x,y) coordinates to (r, theta) coordinates
// Note in 2D, r ≥ 0 and -π ≤ theta ≤ +π
// Input data structure xyArray is [x1,y1,x2,y2,x3,y3...,xN,yN]
// Output data structure polarArray is [r1,theta1,r2,theta2,r3,theta3...,rN,thetaN]

function toPolarCoordinates(xyArray)
	{
	var n = xyArray.length;
	var polarArray = createArray(n);

	for(var i = 0; i < n; i += 2)
		{
		polarArray[i] = Math.hypot(xyArray[i],xyArray[i+1]);
		polarArray[i+1] = Math.atan2(xyArray[i+1],xyArray[i]);
		};
	return polarArray;
	};


function testComplex(testChoice)
	{
		if(testChoice)
		{
		console.log("\nStart testComplex");
		var temp = [];
		var z1 = merge([real1, imag1],[real1.length]);
		var z2 = merge([real2, imag2],[real2.length]);
	
		console.log("Complex # = ",z2, z2.length);
	
		console.log("Real part = ",re(z2));
	
		console.log("Imaginary part = ",im(z2));
	
		temp = arg(z2);
		for(var i = 0; i < temp.length; i++)
			{
			temp[i] *= (180/Math.PI);
			};
		console.log("phase displayed in degrees: ",temp);
	
		temp = abs(z2);
		console.log("magnitude: ",temp);
	
		temp = conjugate(z2);
		console.log("conjugate input & output: ",z2,temp);
	
		var polar = toPolarCoordinates(z2);
		console.log("toPolarCoordinates: ",round(z2,rounding),round(polar,rounding));
		console.log("fromPolarCoordinates output: ",round(polar,rounding),round(fromPolarCoordinates(polar),rounding));
		console.log("compare: ",round(z2,rounding),round(fromPolarCoordinates(polar),rounding));
	
		console.log("Stop testComplex");
		console.log("");
		};
	};

// ****************************************************************
// log function for complex signal arrays
// note that this returns an array with the log in the real positions
// and the angle in the imaginary positions
// 
// checked Monday, 9 January 2017

function log10(array)
	{
	var temp = toPolarCoordinates(array);
	var n = temp.length;
	var x = 0;

	for(var i = 0; i < n; i += 2)
		{
		x = temp[i];
		if(x > 0)
			temp[i] = Math.log10(x);
		else if (x === 0)
			temp[i] = -Infinity;
		else
			temp[i] = NaN;
		};
	return temp;
	};


function testLog(testChoice)
	{
		if(testChoice)
		{
		console.log("\nStart testLog");
		var tdata = [0,-10,1,1,1,0,5,0,10,0,500,0,100000,0,2,0,Math.PI,0,0,0,-6,0];
	
		console.log("log10: ",tdata,log10(tdata));
	
		console.log("Stop testLog");
		console.log("");
		};
	};


// ****************************************************************
// simple functions for real signal arrays
// checked Sunday, 18 December 2016

function clip(array,xmin,xmax,ymin,ymax)
//clip elements of input array to limits of ymin and ymax
	{
	var temp = [];
	var n = array.length;
	var slope = (ymax - ymin) / (xmax - xmin);
	var intercept = (ymin*xmax - ymax*xmin) / (xmax - xmin);

	for(var i = 0; i < n; i += 2)
		{
		if(array[i] < xmin)
			temp[i] = ymin;
		else
			{
			if(array[i] > xmax)
				temp[i] = ymax;
			else
				temp[i] = slope*array[i]+intercept;
			}
		temp[i+1] = 0;
		};
	return temp;
	};


function round(data,rounding)
// Round a floating point array of numbers to n digits after the decimal point
// n is defined as 10^(-n), e.g. 0.001 for n = 3 as in round(temp,0.001)
	{
	if (dimensions(data) === 0)
		return rounding*Math.round(data/rounding);
	
	var temp = [];
	var n = data.length;
	for(var i = 0; i < n; i++)
		{
		temp[i] = rounding*Math.round(data[i]/rounding);
		};
	return temp;
	};


function testReal(testChoice)
	{
		if(testChoice)
		{
		console.log("\nStart testReal");
		var tdata = [0.1,0,1,0,5,0,10,0,500,0,100000,0,2,0,Math.PI,0,0,0,-6,0];
	
		console.log("clip: ",3,6.5,4,12,nlist1,clip(nlist1,3,6.5,4,12));
	
		var temp = log10([0.1,1,5,10,500,100000,2,Math.PI]);
		console.log("round: ",temp,round(temp,rounding));
	
		console.log("Stop testReal");
		console.log("");
		};
	};


// ****************************************************************
// Array dimension checker
// arrayDimension Returns:
//		false when array dimensions are different
//		an Array when is rectangular 0d (i.e. an object) or >=1d
// dimensions returns the number of dimensions in the array
// arrayDimension source: 
// <http://stackoverflow.com/questions/13814621/how-can-i-get-the-dimensions-of-a-multidimensional-javascript-array>
// checked Sunday, 18 December 2016

function arrayDimension(a)
	{
	// Make sure it is an array
	if (a instanceof Array)
		{
		// First element is an array
		var sublength = arrayDimension(a[0]);
		if (sublength === false)
			{
			// Dimension is different
			return false;
			}
		else
			{
			// Compare every element to make sure they are of the same dimensions
			for (var i = 1; i < a.length; i++)
				{
				var _sublength = arrayDimension(a[i]);
				// HACK: compare arrays...
				if (_sublength === false || sublength.join(",") != _sublength.join(","))
					{
					// If the dimension is different (i.e. not rectangular)
					return false;
					}
				};
			// OK now it is "rectangular" (could you call 3d "rectangular"?)
			return [a.length].concat(sublength);
			}
		}
	else
		{
		// Not an array
		return [];
		};
	};


function dimensions(array)
	{
	return arrayDimension(array).length;
	};


function testDimensions(testChoice)
	{
		if(testChoice)
		{
		console.log("\nStart testDimensions");
		
		console.log("scalar: ",tdata0,
			dimensions(tdata0),arrayDimension(tdata0));
		console.log("1D real: ",tdata1D,
			dimensions(tdata1D),arrayDimension(tdata1D));
		console.log("1D complex: ",tdata1Dc,
			dimensions(tdata1Dc),arrayDimension(tdata1Dc));
		console.log("2D real: ",tdata2D,
			dimensions(tdata2D),arrayDimension(tdata2D));
			
		console.log("Stop testDimensions");
		};
	};


// ****************************************************************
// Remember that due to the data format of complex numbers arrays, the amount
// of places to be rotated or shifted should be EVEN.
// checked Sunday, 18 December 2016

function rotateRight(array,rotate)
// assuming arrays are either 1D or 2D signals
	{
	let n = array.length;
	let atemp = [];
	for(var i = 0; i < n; i++)
		atemp[(i + rotate) % n] = array[i];
	return atemp;
	};


function rotateLeft(array,rotate)
	{
	let n = array.length;
	let atemp = [];
	for(var i = 0; i < n; i++)
		atemp[i] = array[(i + rotate) % n];
	return atemp;
	};


function testRotate(testChoice)
	{
		if(testChoice)
		{
		console.log("\nStart testRotate");

		var rotate = 2;
		console.log("rotateRight input, rotate & output &: ",
			alist1,rotate,rotateRight(alist1,rotate));
	
		var z2 = merge([real2, imag2],[real2.length]);

		var rotate = 2;
		console.log("rotateLeft input, rotate & output &: ",
			z2,rotate,rotateRight(z2,rotate));
	
		console.log("Stop testRotate");
		console.log("");
		};
	};


// ****************************************************************
// Filter and signal definitions:
//	  1D - the complex coefficients
//	  2D - the filter size and its complex coefficients, that is,
//		   filter[0] = [rows,cols]
//		   filter[1] = [complex coefficients]
//
// checked Sunday, 18 December 2016

const identityFilter = 
	[ // Notice slightly MODIFIED complex data format
		[3,3],
		[0, 0, 0, 0, 0, 0,
		 0, 0, 1, 0, 0, 0, 
		 0, 0, 0, 0, 0, 0]
	];

const crossFilter = 
	[ // Notice slightly MODIFIED complex data format
		[3,3],
		[1, 0, 0, 0, -1, 0,
		 0, 0, 0, 0, 0, 0,
		-1, 0, 0, 0, 1, 0]
	];

const sobelFilter = 
	[ // Notice slightly MODIFIED complex data format
		[3,3],
		[-1, 0, 0, 0, 1, 0,
		 -2, 0, 0, 0, 2, 0,
		-1, 0, 0, 0, 1, 0]
	];

const youngLPFFilter = 
	[
		0.00115333,0,0.0335179,0,0.0051661,0,-0.00569677,0,-0.0137668,0,-0.00892154,0,
		0.00641936,0,0.0182781,0,0.0133235,0,-0.00705987,0,-0.0248268,0,-0.0203065,0,
		0.00759625,0,0.0354256,0,0.0328779,0,-0.00799214,0,-0.0573142,0,-0.0633473,0,
		0.00825774,0,0.141242,0,0.271113,0,0.324989,0,0.271113,0,0.141242,0,
		0.00825774,0,-0.0633473,0,-0.0573142,0,-0.00799214,0,0.0328779,0,0.0354256,0,
		0.00759625,0,-0.0203065,0,-0.0248268,0,-0.00705987,0,0.0133235,0,0.0182781,0,
		0.00641936,0,-0.00892154,0,-0.0137668,0,-0.00569677,0,0.0051661,0,0.0335179,0,
		0.00115333,0
	];

const youngBPFFilter = 
	[
		0.0152455,0,0.0372391,0,0.0167843,0,0.0188983,0,0.00176545,0,-0.0079223,0,
		-0.0135453,0,-0.00775344,0,0.00419412,0,0.0156094,0,0.0165927,0,0.00195226,0,
		-0.0274134,0,-0.0623665,0,-0.0899674,0,-0.0974055,0,-0.0781374,0,-0.0336456,0,
		0.0257054,0,0.0850691,0,0.128413,0,0.144367,0,0.128413,0,0.0850691,0,
		0.0257054,0,-0.0336456,0,-0.0781374,0,-0.0974055,0,-0.0899674,0,-0.0623665,0,
		-0.0274134,0,0.00195226,0,0.0165927,0,0.0156094,0,0.00419412,0,-0.00775344,0,
		-0.0135453,0,-0.0079223,0,0.00176545,0,0.0188983,0,0.0167843,0,0.0372391,0,
		0.0152455,0
	];

const youngHPFilter = 
	[
		-0.00301318,0,-0.0334262,0,-0.00655957,0,-0.00213755,0,0.00724907,0,0.0144713,0
		,0.0155009,0,0.00841314,0,-0.00480427,0,-0.0184909,0,-0.0254566,0,-0.0202915,0,
		-0.00254728,0,0.0219084,0,0.0420527,0,0.0455796,0,0.023733,0,-0.0241999,0,
		-0.0896421,0,-0.156564,0,-0.206513,0,0.775007,0,-0.206513,0,-0.156564,0,
		-0.0896421,0,-0.0241999,0,0.023733,0,0.0455796,0,0.0420527,0,0.0219084,0,
		-0.00254728,0,-0.0202915,0,-0.0254566,0,-0.0184909,0,-0.00480427,0,0.00841314,0,
		0.0155009,0,0.0144713,0,0.00724907,0,-0.00213755,0,-0.00655957,0,-0.0334262,0,
		-0.00301318,0
	];

const youngNotchFilter = 
	[
		0.0000862854,0,-0.00235865,0,-0.0047307,0,-0.00216514,0,0.0043166,0,0.00576716,0,
		0.000753313,0,0.00230128,0,0.0150072,0,0.0173841,0,-0.00766623,0,-0.0345907,0,
		-0.0258923,0,0.000800227,0,-0.00931434,0,-0.0459417,0,-0.0176896,0,0.105425,0,
		0.191674,0,0.088882,0,-0.146544,0,0.728244,0,-0.146544,0,0.088882,0,
		0.191674,0,0.105425,0,-0.0176896,0,-0.0459417,0,-0.00931434,0,0.000800227,0,
		-0.0258923,0,-0.0345907,0,-0.00766623,0,0.0173841,0,0.0150072,0,0.00230128,0,
		0.000753313,0,0.00576716,0,0.0043166,0,-0.00216514,0,-0.0047307,0,-0.00235865,0,
		0.0000862854,0
	];

const youngDerivFilter = 
	[
		0.0701528,0,-0.0862785,0,0.0251206,0,-0.0159789,0,0.0134707,0,-0.012734,0,
		0.0126871,0,-0.0130001,0,0.0135774,0,-0.0143798,0,0.0154115,0,-0.0167443,0,
		0.018335,0,-0.0204765,0,0.0233466,0,-0.0270312,0,0.0322177,0,-0.0400577,0,
		0.0532393,0,-0.0797112,0,0.159218,0,0,0,-0.159218,0,0.0797112,0,
		-0.0532393,0,0.0400577,0,-0.0322177,0,0.0270312,0,-0.0233466,0,0.0204765,0,
		-0.018335,0,0.0167443,0,-0.0154115,0,0.0143798,0,-0.0135774,0,0.0130001,0,
		-0.0126871,0,0.012734,0,-0.0134707,0,0.0159789,0,-0.0251206,0,0.0862785,0,
		-0.0701528,0
	];


function testFilters(testChoice)
	{
		if(testChoice)
		{
		console.log("\nStart testChoice");
		
		console.log("identityFilter = ",
			identityFilter.length,
			identityFilter[0],
			identityFilter[1],
			total(identityFilter[1]));

		console.log("crossFilter = ",
			crossFilter.length,
			crossFilter[0],
			crossFilter[1],
			total(crossFilter[1]));

		console.log("sobelFilter = ",
			sobelFilter.length,
			sobelFilter[0],
			sobelFilter[1],
			total(sobelFilter[1]));

		console.log("youngLPFFilter = ",
			dimensions(youngLPFFilter),
			arrayDimension(youngLPFFilter),
			total(youngLPFFilter),
			youngLPFFilter);

		console.log("youngBPFFilter = ",
			dimensions(youngBPFFilter),
			arrayDimension(youngBPFFilter),
			total(youngBPFFilter),
			youngBPFFilter);

		console.log("youngHPFilter = ",
			dimensions(youngHPFilter),
			arrayDimension(youngHPFilter),
			total(youngHPFilter),
			youngHPFilter);

		console.log("youngNotchFilter = ",
			dimensions(youngNotchFilter),
			arrayDimension(youngNotchFilter),
			total(youngNotchFilter),
			youngNotchFilter);

		console.log("youngDerivFilter = ",
			dimensions(youngDerivFilter),
			arrayDimension(youngDerivFilter),
			total(youngDerivFilter),
			youngDerivFilter);
			
		console.log("Stop testChoice");
		};
	};

// ****************************************************************
// For 1D signals
// checked Wednesday, 21 December 2016

function signalChoiceList(length,type)
	{
	var myLength = oddQ(length) ? 2*(length + 1) : 2*length;
	var signal = createArray(myLength);
	var signala = createArray(myLength/2);
	var signalb = createArray(myLength/2);
	const choices = 6;

	switch(type)
		{
			case 1:  // constant
				for(var i = 0; i < myLength; i += 2)
					{
					signal[i] = 1;
					signal[i+1] = 0;
					};
				break;
			case 2:  // ramp(n)
				for(var i = 0; i < myLength; i += 2)
					{
					signal[i] = i/2;
					signal[i+1] = 0;
					};
				break;
			case 3:  // sqrt(n)
				for(var i = 0; i < myLength; i += 2)
					{
					signal[i] = Math.sqrt(i/2);
					signal[i+1] = 0;
					};
				break;
			case 4:  // n^2
				for(var i = 0; i < myLength; i += 2)
					{
					signal[i] = i*i/4;
					signal[i+1] = 0;
					};
				break;
			case 5:  // e^(-n)
				for(var i = 0; i < myLength; i += 2)
					{
					signal[i] = Math.exp(-i/2);
					signal[i+1] = 0;
					};
				break;
			case 6:  // ramp+discontinuity(n)
				for(var i = 0; i < myLength/2; i += 2)
					{
					signala[i] = 1;
					signala[i+1] = 0;
					};
				for(var i = 0; i < myLength/2; i += 2)
					{
					signalb[i] = i/2+4;
					signalb[i+1] = 0;
					};
				signal = signala.concat(signalb);
				break;
			default:
				throw("Non-existent signal choice");
		};
	return signal;
	};


function testSignalChoiceList(testChoice)
	{
		if(testChoice)
		{
		console.log("\nStart testSignalChoiceList");
		
		var tlist = [];
		var localLength = 10;
		var signalType = 6;

		tlist = signalChoiceList(localLength,signalType);
		console.log("signalChoiceList: ",tlist.length,tlist);

		console.log("Stop testSignalChoiceList");
		console.log("");
		};
	};


// ****************************************************************
// Compute running sum of absolute value of a signal |x[n]| to identify
// where the real signal begins. Even if the signal is 2D it will be
// treated as 1D.
// checked Monday, 9 January 2017

function runningSum(data)
	{
	var myData = abs(flatten(data));
	var rSum = createArray(myData.length);
	
	rSum[0] = myData[0];
	rSum[1] = 0;
	for(var i = 2; i < rSum.length; i += 2)
		{
		rSum[i] = rSum[i-2] + myData[i];
		rSum[i+1] = 0;
		};
	var newMax = getMaxOfArray(rSum);
	for(var i = 0; i < rSum.length; i += 2)
		rSum[i] /= newMax;

	return rSum;
	};


function testRunningSum(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testRunningSum");
		
		var tdata = [-3,0,1,0,50,0,12,0,0,-2];
		console.log("runningSum - samples, input, output: ",
			tdata.length,tdata,runningSum(tdata));
	
		var tdata = nlist5.slice();
		console.log("runningSum - samples, input, output: ",
			tdata.length,tdata,runningSum(tdata));
	
		var myLength = 100
		var tdata = Array(myLength); // temporary list
		var lambda = 9;
		for (var i=0; i<myLength; i++)
			{
			tdata[i] = randomRayleigh(lambda); // both real and imag values
			};
		console.log("runningSum - samples, input, output: ",
			tdata.length,tdata,runningSum(tdata));
			
		console.log("Stop testRunningSum");
		console.log("");
		};
	};


// ****************************************************************
// Starting at the nth position in the list find the first value that ≥ p.
// then return the position i and the value data[i]. The position
// is based upon a 1D signal. The signal will first be processed with
// runningSum(). Note 1 ≤ n ≤ list.length and 0 ≤ p ≤ 1
// checked Monday, 9 January 2017

function select(data,n,p)
	{
	var temp = runningSum(data);
	var plocal = 0;
	if(p < 0 )
		{ plocal = 0; }
	else if(p > 1 )
		{ plocal = 1; }
	else plocal = p;
	
	for(var i = 2*(n - 1); i < 2*temp.length; i += 2)
		if(temp[i] >= plocal) return [i/2,data[i]];
	return false;
	};


function testSelect(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testSelect");
		
		var tdata = [1,0,1,0,1,0,1,0,1,0,1,0];
		var start = 1;
		var p = 0.29;
		var result = select(tdata,start,p);
		console.log("select - samples, input, running, p, position, value[position]: ",
			tdata.length,tdata,runningSum(tdata),p,result[0],result[1] );
	
		var tdata = [-3,0,1,0,50,0,12,0,0,-2];
		var start = 1;
		var p = 0.69;
		var result = select(tdata,start,p);
		console.log("select - samples, input, running, p, position, value[position]: ",
			tdata.length,tdata,runningSum(tdata),p,result[0],result[1] );
	
		var myLength = 8
		var tlist = Array(myLength); /* temporary list */
		var lambda = 9;
		for (var i=0; i<myLength; i++)
			tlist[i] = randomRayleigh(lambda); // both real and imag values
		var p = 0.50;
		var result = select(tlist,start,p);
		console.log("select - samples, input, running, p, position, value[position]: ",
			tlist.length,tlist,runningSum(tlist),p,result[0],result[1] );
			
		console.log("Stop testSelect");
		console.log("");
		};
	};


// ****************************************************************
// Starting at the nth position in the list find the first value in
// the 1D real signal that matches the "element" provided.
// Note 1 ≤ n ≤ list.length
// checked Monday, 9 January 2017

function position(data,n,element)
	{
	var temp = data.slice();
	
	for(var i = 2*(n - 1); i < 2*temp.length; i += 2)
		if(temp[i] === element) return i/2;
	return false;
	};


function testPosition(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testPosition");
		
		var tdata = [-3,0,1,0,50,0,12,0,0,-2];
		var start = 1;
		var pos = 2;
		var elem = tdata[2*pos];
		
		console.log("position - samples, input, element, output: ",
			tdata.length,tdata,elem,position(tdata,start,elem));
		
		start = 4;
		console.log("position - samples, input, element, output: ",
			tdata.length,tdata,elem,position(tdata,start,elem));
		
		tdata = ["a","b","c","d",Math.PI,-6];
		start = 1;
		pos = 1;
		elem = tdata[2*pos];
		
		console.log("position - samples, input, element, output: ",
			tdata.length,tdata,elem,position(tdata,start,elem));
	
		start = 1;
		var myLength = 8
		var tlist = Array(myLength); // temporary list
		var lambda = 9;
		for (var i=0; i<myLength; i++)
			{
			tlist[i] = randomRayleigh(lambda);
			};
			
		elem = tlist[2];
		console.log("position - samples, input, element, output: ",
			tlist.length,tlist,elem,position(tlist,start,elem));

			
		console.log("Stop testPosition");
		console.log("");
		};
	};


// ****************************************************************
// Starting at the 1st pixel in the 2D REAL image, find all the values
// in the image that match the "element" provided.
// Return the list of positions
// checked Tuesday, 10 January 2017

function pixelValuePositions(data,nn,element)
	{
	var rows = nn[0];  // but for 2D signal use this
	var cols = nn[1];
	var temp = data.slice();  // make a copy
	
	var coords = [];
	var rtemp = 0;
	var ctemp = 0;
	for(var i = 0; i < temp.length; i += 2)
		if(temp[i] === element) 
			{ // how to turn linear index into 2D index
			ctemp = (i/2) % cols;
			rtemp = Math.floor((i/2) / cols);
			coords.push([rtemp,ctemp]);
			};
	return coords;
	};


function testPixelValuePositions(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testPixelValuePositions");
		
		var nn = [3,5]
		var myLength = nn[0]*nn[1]; // number of complex samples
		var data2D = createArray(2*myLength);
		
		for(var i = 0; i< 2*myLength; i += 2)
			{
			data2D[i] = randomInteger(1,5); data2D[i+1] = 0;
//			data2D[i] = i; data2D[i+1] = 0;
			};
		
		var elem = 3;
		var result = pixelValuePositions(data2D,nn,elem);
		console.log("pixelValuePosition: ",nn,elem,data2D,result);
		
		elem = 2;
		result = pixelValuePositions(data2D,nn,elem);
		console.log("pixelValuePosition: ",nn,elem,data2D,result);
		
		console.log("Stop testPixelValuePositions");
		console.log("");
		};
	};


// ****************************************************************
// get a list of prime factors of a number
// based upon <http://www.javascripter.net/math/primes/factorization.htm>
// checked Tuesday, 10 January 2017

function factor(n,result)
	{
	if (isNaN(n) || !isFinite(n) || n%1!=0 || n==0)
		return [0];
	if (n < 0)
		return factor(-n,result);
	var minFactor = leastFactor(n);
	var myLength = result.push(minFactor);
	if (n == minFactor)
		{
		return result;
		}
	else
		{
		return factor(n/minFactor,result);
		};
	};

// find the least factor in n by trial division
function leastFactor(n)
	{
	if (isNaN(n) || !isFinite(n)) return NaN;  
	if (n==0) return 0;  
	if (n%1 || n*n<2) return 1;
	if (n%2==0) return 2;  
	if (n%3==0) return 3;  
	if (n%5==0) return 5;  
	var m = Math.sqrt(n);
	for (var i=7;i<=m;i+=30)
		{
		if (n%i==0)      return i;
		if (n%(i+4)==0)  return i+4;
		if (n%(i+6)==0)  return i+6;
		if (n%(i+10)==0) return i+10;
		if (n%(i+12)==0) return i+12;
		if (n%(i+16)==0) return i+16;
		if (n%(i+22)==0) return i+22;
		if (n%(i+24)==0) return i+24;
		}
	return n;
	}

// Based upon Mathematica's FactorInteger function. Given a number n itemp
// returns the prime factors and the number of times each factor is involved
// For example, if n = 90 then calling factorInteger(factor(n,result)) returns
// the array is primes = [[2,1],[3,2],[5,1]]
function factorInteger(factors)
	{
	var uniqPrimes = [];
	var powerPrimes = [];
	var uniq = factors[0];
	var primes = [];
	var i = 0;
	var count = 1;
	
	for(i = 1; i < factors.length; i++)
		{
		if(factors[i] == uniq)
			{
			count++;
			}
		else
			{
			uniqPrimes.push(uniq);
			powerPrimes.push(count);
			uniq = factors[i];
			count = 1;
			};
		};
		
	// close down gracefully
	i--;
	uniqPrimes.push(factors[i]);
	if(factors[i] == factors[i-1])
		{
		powerPrimes.push(count);
		}
	else
		{
		powerPrimes.push(1);
		};

	// reproduce Mathematica format
	for(var i = 0; i < uniqPrimes.length; i++)
		{
		primes.push([uniqPrimes[i],powerPrimes[i]]);
		};
	return primes;
	};

function testFactorInteger(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testFactorInteger");
		var n = 9*5*2;
		var result=[];
		var endresult = factor(n,result);
		var endprimes = factorInteger(endresult);
		console.log("factorInteger n, primes, endprimes: ",
			n,endresult,endprimes);
		
		n = 2*7*7;
		result=[];
		endresult = factor(n,result);
		var endprimes = factorInteger(endresult);
		console.log("factorInteger n, primes, endprimes: ",
			n,endresult,endprimes);
		
		n = 9*17*25*Math.pow(2,9);
		result=[];
		endresult = factor(n,result);
		var endprimes = factorInteger(endresult);
		console.log("factorInteger n, primes, endprimes: ",
			n,endresult,endprimes);
		
		n = Math.pow(2,47)+1;
		result=[];
		endresult = factor(n,result);
		var endprimes = factorInteger(endresult);
		console.log("factorInteger n, primes, endprimes: ",
			n,endresult,endprimes);
		
		console.log("Stop testFactorInteger");
		console.log("");
		};
	};


// ****************************************************************
// Window functions for spectral estimation
// window functions are even and non-zero over the interval |t| <= 1
// For windows where it is relevant a typical value of sigma is 1/4.
// checked Tuesday, 10 January 2017

// Block (rectangular) window
function block(t,sigma)
	{
	if ( t < -1 || t > 1)
		return 0;
	else
		return 1;
	};
	
// Bartlett (triangular) window
function bartlett(t,sigma)
	{
	if ( t < -1 || t > 1)
		return 0;
	else
		return 1 - Math.abs(t);
	};
	
// Parzen (parabolic) window
function parzen(t,sigma)
	{
	if ( t < -1 || t > 1)
		return 0;
	else
		return 1 - (t*t);
	};
	
// Truncated Gaussian window
function truncGauss(t,sigma)
	{
	if ( t < -1 || t > 1)
		return 0;
	else
		return Math.exp(-(t*t)/(2*sigma*sigma));
	};
	
// Truncated Verbeek window
function truncVerbeek(t,sigma)
	{
	if ( t < -1 || t > 1)
		return 0;
	else
		{
		var z = (t*t)/(2*sigma*sigma);
		return (1 + z + (z*z) / 2) * Math.exp(-z);
		};
	};
	
// Tukey window
function tukey(t,sigma)
	{
	if ( t < -1 || t > 1)
		return 0;
	else
		return (1 + Math.cos(t * Math.PI)) / 2;
	};
	
// Hamming window
function hamming(t,sigma)
	{
	if ( t < -1 || t > 1)
		return 0;
	else
		return (0.54 + 0.46 * Math.cos(t * Math.PI));
	};
	
// Hann window
function hann(t,sigma)
	{
	if ( t < -1 || t > 1)
		return 0;
	else
		return Math.cos(t * Math.PI / 2);
	};
	
// Young window
function young(t,sigma)
	{
	if ( t < -1 || t > 1)
		return 0;
	else
		return Math.E * Math.exp(-1 / (1 - (t * t)));
	};
	
// Ideal Gaussian window
function idealGauss(t,sigma)
	{
	return Math.exp(-(t*t)/(2*sigma*sigma));
	};
	

function testWindows(testChoice)
	{
		if(testChoice)
		{
		console.log("\nStart testWindows");
		
		var t = 1/2  // play with this
		var sigma = 1/4;  // standard value
		
		console.log("Windows: t, value ");
		console.log("	Block: t, value ", t, block(t,sigma));
		console.log("	Bartlett: t, value ", t, bartlett(t,sigma));
		console.log("	Parzen: t, value ", t, parzen(t,sigma));
		console.log("	Truncated Gaussian: t, value ", t, truncGauss(t,sigma));
		console.log("	Truncated Verbeek: t, value ", t, truncVerbeek(t,sigma));
		console.log("	Tukey: t, value ", t, tukey(t,sigma));
		console.log("	Hamming: t, value ", t, hamming(t,sigma));
		console.log("	Young: t, value ", t, young(t,sigma));
		console.log("	Ideal Gauss: t, value ", t, idealGauss(t,sigma));
		
		console.log("Stop testWindows");
		console.log("");
		};
	};


// ****************************************************************
// Following is based upon Mathematica's procedure RecurrenceFilter.
// The data is filtered by a linear, constant coefficient difference equation.
// The left coefficients represent the recursion part on the output "y". Th
// right coefficients represent the FIR part on the input "x". Note that setiing
// leftCoeffs = [1] means that no recursion will take place and the entire
// procedure is just a FIR filter.
//
// The following assumptions are made: 1) this procedure is intended for 1D data;
// 2) x is causal; 3) leftCoeffs[0] ≠ 0; 4) y[n] = 0 for n < 0, 5) the coefficients and the data are real.
// checked Thursday, 22 December 2016


function recurrenceFilter(coefficients, data)
	{
	var leftCoeffs = coefficients[0].slice();
	var rightCoeffs = coefficients[1].slice();
	var x = data.slice();
	var leftLength = leftCoeffs.length;
	var rightLength = rightCoeffs.length;
	var xLength = x.length;
	
	// make coefficient of y[n] = 1
	var normCoeff = leftCoeffs[0];
	for (var i = 0; i < leftLength; i += 2)
		leftCoeffs[i] /= normCoeff;
	for (var i = 0; i < rightLength; i += 2)
		rightCoeffs[i] /= normCoeff;
	
	// this is the output data array
	var y = createArray(xLength);
	for (var i = 0; i < xLength; i++)  y[i] = 0;
		
	var causal = 0;
	var temp = 0;
	if(leftLength == 2) // This is a simple FIR filter, no recursion
		{
		for (var n = 0; n < xLength; n += 2)
			{
			for (var i = 0; i < rightLength; i += 2)
				{
				causal = n-i;
				if (causal < 0)
					temp = 0
				else
					temp = x[causal];
				y[n] += (rightCoeffs[i] * temp);
				};
			};
		}
	else // This is the whole nine yards, i.e. recursion
		{
		for (var n = 0; n < xLength; n += 2)
			{
			for (var i = 0; i < rightLength; i += 2)
				{
				causal = n-i;
				if (causal < 0)
					temp = 0
				else
					temp = x[causal];
				y[n] += (rightCoeffs[i] * temp);
				};
			for (var i = 2; i < leftLength; i += 2)
				{
				causal = n-i;
				if (causal < 0)
					temp = 0
				else
					temp = y[causal];
				y[n] -= (leftCoeffs[i] * temp);
				};
			};
		};
	return y;
	};


function testRecurrenceFilter(testChoice)
	{
		if(testChoice)
		{
		console.log("\nStart testRecurrenceFilter");
		
		var leftCoeffs = [1,0];
		var rightCoeffs = [1,0,-1/2,0];
		var coefficients = [ leftCoeffs, rightCoeffs ];
		
		var x = [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0];
		var y = recurrenceFilter(coefficients, x);
		console.log("recurrenceFilter: x, coefficients, y",
			x, coefficients, y);
		
		leftCoeffs = [1,0,-1/2,0];
		rightCoeffs = [1,0];
		coefficients = [ leftCoeffs, rightCoeffs ];
		
		x = [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
		y = recurrenceFilter(coefficients, x);
		console.log("recurrenceFilter: x, coefficients, y",
			x, coefficients, y);
			
		x = [0,0,0,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0];
		y = recurrenceFilter(coefficients, x);
		console.log("recurrenceFilter: x, coefficients, y",
			x, coefficients, y);
			
		leftCoeffs = [1,0];
		rightCoeffs = youngBPFFilter;
		coefficients = [ leftCoeffs, rightCoeffs ];
		x = signalChoiceList(100,2);
		y = recurrenceFilter(coefficients, x);
		console.log("recurrenceFilter BPF: x, coefficients, y",
			x, coefficients, round(y,rounding));
			
		leftCoeffs = [1,0];
		rightCoeffs = [1,0,-1,0];
		coefficients = [ leftCoeffs, rightCoeffs ];
		x = signalChoiceList(10,4);
		y = recurrenceFilter(coefficients, x);
		console.log("recurrenceFilter n^2: x, coefficients, y",
			x, coefficients, y);
			
		console.log("Stop testRecurrenceFilter");
		console.log("");
		};
	};


// ****************************************************************
// Following is based upon Mathematica's procedure GaussianMatrix.
// a matrix is produce of dimensions (2n+1) x (2n+1) with the center
// at (n+1, n+1). values are the bivariate Gaussian with a standard
// deviation of s. Normally, s = n. The values are chosen so that the
// area under the bivariate Gaussian is 1.

// Normalized Gaussians and derivatives
// checked Thursday, 22 December 2016


function normGauss(t,sigma)
	{
	return Math.exp(-(t*t)/(2*sigma*sigma))/(Math.sqrt(2*Math.PI)*sigma);
	};
	
function bivariateGauss(x,sigmaX, y,sigmaY)
	{
	return normGauss(x,sigmaX)*normGauss(y,sigmaY);
	};
	
function dXGauss(x,sigmaX, y,sigmaY)  // derivative in x direction
	{
	return -x*bivariateGauss(x,sigmaX, y,sigmaY)/(sigmaX*sigmaX);
	};
	
function dYGauss(x,sigmaX, y,sigmaY)  // derivative in y direction
	{
	return -y*bivariateGauss(x,sigmaX, y,sigmaY)/(sigmaY*sigmaY);
	};
	
function dXYGauss(x,sigmaX, y,sigmaY)  // derivative in x and y directions
	{
	return x*y*bivariateGauss(x,sigmaX, y,sigmaY)/(sigmaX*sigmaX*sigmaY*sigmaY);
	};
	
function gaussianMatrix(rows , cols, s, derivatives)
	{
	if (evenQ(rows) || evenQ(cols)) throw("ROWS and COLS must both be odd.");
	var mx = (cols+1)/2;
	var my = (rows+1)/2;
	// Note column-to-column is the x direction,
	// row-to-row is the y direction
	var gaussImage = createArray(cols,rows);
	// the Gauss image is real but we must create an imaginary part
	var zeroes = createArray(cols*rows);
	zeroes.fill(0);
	
	switch(derivatives)
		{
		case "None":  // no derivative
			var ttotal = 0;
			var temp = 0;
			for (var i = 0; i < rows; i++)
				for (var j = 0; j < cols; j++)
					{
					temp = bivariateGauss((i-mx),s,(j-my),s);
					gaussImage[j][i] = temp;
					ttotal += temp;
					};
			// normalize area under Gauss to one (following Mathematica)
			for (var i = 0; i < rows; i++)
				for (var j = 0; j < cols; j++)
					gaussImage[j][i] /= ttotal;
			break;
		case "x":  // derivative in x direction
			for (var i = 0; i < rows; i++)
				for (var j = 0; j < cols; j++)
					gaussImage[j][i] = dXGauss((i-mx),s,(j-my),s);
			break;
		case "y":  // derivative in y direction
			for (var i = 0; i < rows; i++)
				for (var j = 0; j < cols; j++)
					gaussImage[j][i] = dYGauss((i-mx),s,(j-my),s);
			break;
		case "xy":  // derivative in x & y directions
			for (var i = 0; i < rows; i++)
				for (var j = 0; j < cols; j++)
					gaussImage[j][i] = dXYGauss((i-mx),s,(j-my),s);
			break;
		default:
			throw("Non-existent derivative choice");
		};
	var dims = [rows,cols];
	var vFormat = merge([flatten(gaussImage),zeroes],[zeroes.length]);
	var kernFormat = [dims, vFormat];
	return kernFormat;
	};


function testGaussianMatrix(testChoice)
	{
		if(testChoice)
		{
		console.log("\nStart testGaussianMatrix");
		
		var rows = 7;  // must be odd
		var cols = 7;  // must be odd
		var sigma = 3/2
		
		var derivatives = "None";
		var result = gaussianMatrix(rows, cols, sigma, derivatives);
		console.log("gaussianMatrix: rows, cols, sigma, derivatives, result, total(result)", rows,cols,sigma,derivatives,result[0],round(result[1],rounding),round(total(result[1]),rounding));
		
		rows = 9;  // must be odd
		cols = 7;  // must be odd
		derivatives = "x";
		result = gaussianMatrix(rows, cols, sigma, derivatives);
		console.log("gaussianMatrix: rows, cols, sigma, derivatives, result, total(result)", rows,cols,sigma,derivatives,result[0],round(result[1],rounding),round(total(result[1]),rounding),split(round(result[1],rounding),[rows,cols]));
		
		derivatives = "y";
		result = gaussianMatrix(rows, cols, sigma, derivatives);
		console.log("gaussianMatrix: rows, cols, sigma, derivatives, result, total(result)", rows,cols,sigma,derivatives,result[0],round(result[1],rounding),round(total(result[1]),rounding));
		
		rows = 21;
		cols = 21;
		sigma = 10/2
		derivatives = "xy";
		result = gaussianMatrix(rows, cols, sigma, derivatives);
		console.log("gaussianMatrix: rows, cols, sigma, derivatives, result, total(result)", rows,cols,sigma,derivatives,result[0],round(result[1],rounding),round(total(result[1]),rounding));
		
		console.log("Stop testGaussianMatrix");
		console.log("");
		};
	};


/// ****************************************************************
// Transpose of a row vector (1D array of nn = [1,n]) into a column vector
// array of nn = [n,1]. Because the elements of the vector are complex,
// the length of the vector must be EVEN.
// For applicability, see test program below
// checked Tuesday, 10 January 2017

function vTranspose(v)
	{
	var vLength = (v.length)/2;
	return [v,[vLength,1]]; 
	};


function testVTranspose(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testVTranspose");
		
		xlist = [1,0,2,0,3,0];
		ylist = ["a","b","c","d"];
		
		var result = vTranspose(xlist);
		console.log("vTranspose - v, vt: ", xlist, result[0], result[1]);

		result = vTranspose(ylist);
		console.log("vTranspose - v, vt: ", ylist, result[0], result[1]);
		
		console.log("Stop testVTranspose");
		console.log("");
		};
	};


/// ****************************************************************
// Outer product of a column vector vt (N x 1) and a row vector u (1 x M)
// to produce a matrix that is (N x M)
// type can be 1D arrays of real numbers
// checked Tuesday, 10 January 2017

function vOuter(v,u)
	{
	var vLength = v.length;
	var uLength = u.length;
	var matrix = [];
	var mNN = [vLength/2, uLength/2];
	
	for (var i = 0; i < vLength; i += 2)
		{
		for (var j = 0; j < uLength; j += 2)
			{
			var r1 = v[i];
			var i1 = v[i+1];
			var r2 = u[j];
			var i2 = u[j+1];
			rx = (r1*r2 - i1*i2);
			ix = (i1*r2 + r1*i2);
			matrix.push(rx,ix);
			};
		};
	return [matrix, mNN]; 
	};


function testVOuter(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testVOuter");
		
		var xlist = [1,0,2,0,3,0,4,0];
		var ylist = [5,0,7,0,11,0];
		
		var result = vOuter(xlist,ylist);
		console.log("vTranspose - xlist, ylist, matrix: ",
			xlist, ylist, result[0], result[1]);
		
		xlist = [1,1,2,3,3,6,4,8];
		ylist = [5,-1,7,0,11,5];
		
		result = vOuter(xlist,ylist);
		console.log("vTranspose - xlist, ylist, matrix: ",
			xlist, ylist, result[0], result[1]);
		
		console.log("Stop testVOuter");
		console.log("");
		};
	};


// ****************************************************************
// Following is based upon Mathematica's procedures ListConvolve and
// ListCorrelate. The length of the output equals the length of the input signal
// but zero-padding and cyclic indexing are used to get the intuitive result.
// The following assumptions are made: 1) this procedure is intended for 1D data;
// where the input is causal.
//
// checked Tuesday, 10 January 2017

function listCorrelate(kernel, data)
	{
	return listConvolve(reverse(kernel), data)
	};
	
	
function listConvolve(kernel, data)
	{
	var filter = kernel.slice();
	var fLength = filter.length;
	var input = data.slice();
	var iLength = input.length;
	// zero-pad the input array with the length of the filter array
	// the new input length = original input length + filter length
	for (var i = 0; i < fLength; i++)
		iLength = input.push(0);
	// but don't zero-pad the output array
	var output = Array(iLength-fLength);
	var outLength = output.length;

	var tempr = 0;
	var tempi = 0;
	var index = 0;
	var ir = 0;
	var ii = 0;
	var fr = 0;
	var fi = 0;
	// This is the convolution
	for (var j = 0; j < outLength; j += 2)
		{
		output[j] = 0;
		output[j+1] = 0;
		tempr = 0;
		tempi = 0;
		for (var i = 0; i < fLength; i += 2)
			{
			index = (j-i) % iLength;  // cyclic index
			if ( index < 0 ) index += iLength;
			
			ir = input[index];
			ii = input[index+1];
			
			fr = filter[i];
			fi = filter[i+1];
			
			tempr += ( ir*fr - ii*fi );
			tempi += ( ii*fr + ir*fi );
			};
		output[j] = tempr;
		output[j+1] = tempi;
		};
	return output;
	};


function testListConvolve(testChoice)
	{
		if(testChoice)
		{
		console.log("\nStart testListConvolve");
		
		var kernel = [1,0,2,0,3,0];
		var x = [0,0,0,0,1,0,0,0,0,0,0,0,0,0];
		var y = listConvolve(kernel, x);
		console.log("listConvolve: kernel, x, y", kernel, x, y);
			
		kernel = [1,0,2,0,3,0];
		var x = [1,0,0,0,1,-Math.PI,-1,1,0,0,Math.E,0,0,0];
		var y = listConvolve(kernel, x);
		console.log("listConvolve: kernel, x, y", kernel, x, y);
		// OK compared to Mathematica Tuesday, 10 January 2017
			
		var x = [1,0,0,0,0,0,0,0,0,0];
		var y = listConvolve(kernel, x);
		console.log("listConvolve: kernel, x, y", kernel, x, y);
			
		var x = [0,0,0,0,0,0,1,0,2,0,3,0,0,0,0,0,0,0];
		var y = listConvolve(kernel, x);
		console.log("listConvolve: kernel, x, y", kernel, x, y);
			
		var x = [0,0,0,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0];
		var y = listConvolve(kernel, x);
		console.log("listConvolve: kernel, x, y", kernel, x, y);
			
		kernel = [1,0,-1,0];
		var x = signalChoiceList(10,4);
		var y = listConvolve(kernel, x);
		console.log("listConvolve: kernel, x, y", kernel, x, y);
			
		kernel = [1,-1,-1,1];
		var x = signalChoiceList(10,4);
		var y = listConvolve(kernel, x);
		console.log("listConvolve: kernel, x, y", kernel, x, y);
			
		console.log("Stop testListConvolve");
		console.log("");
		};
	};

function testListCorrelate(testChoice)
	{
		if(testChoice)
		{
		console.log("\nStart testListCorrelate");
		
		var kernel = [1,0,2,0,3,0];
		var x = [0,0,0,0,1,0,0,0,0,0,0,0,0,0];
		var y = listCorrelate(kernel, x);
		console.log("listCorrelate: kernel, x, y", kernel, reverse(kernel), x, y);
			
		var x = [1,0,0,0,0,0,0,0,0,0];
		var y = listCorrelate(kernel, x);
		console.log("listCorrelate: kernel, x, y", kernel, reverse(kernel), x, y);
			
		var x = [0,0,0,0,1,0,2,0,3,0,0,0,0,0];
		var y = listCorrelate(kernel, x);
		console.log("listCorrelate: kernel, x, y", kernel, reverse(kernel), x, y);
			
		var x = [0,0,0,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0];
		var y = listCorrelate(kernel, x);
		console.log("listCorrelate: kernel, x, y", kernel, reverse(kernel), x, y);
			
		kernel = [1,0,-1,0];
		var x = signalChoiceList(10,4);
		var y = listCorrelate(kernel, x);
		console.log("listCorrelate: kernel, x, y", kernel, reverse(kernel), x, y);
			
		console.log("Stop testListCorrelate");
		console.log("");
		};
	};


// ****************************************************************
// Following is based upon Mathematica's procedure ArrayPad. Padding
// is added to the END of a 1D (complex) signal. For a 1D signal, the "signal"
// is just a list of complex values in the "universal" format.
// "padding" refers to the number of complex values to be appended to the
// kernel and "value" is the complex value that is to be appended.

// checked Thursday, 23 February 2017

function arrayPad(signal, padding, padvalue)
	{
	var tail = [];
	for (var i = 0; i < padding; i++)
		tail.push(padvalue[0],padvalue[1]);  // build a tail of constant complex values
		
	return signal.concat(tail);  // return padded array in 1D universal format
	};


function testArrayPad(testChoice)
	{
		if(testChoice)
		{
		console.log("\nStart testArrayPad");
		
		var x = [1,2,3,4,5,6,7,8];
		
		var padding = 3;  // number of complex signal values to be added, e.g. 3
		var padvalue = [0, 0];
		var result = arrayPad(x, padding, padvalue);
		console.log("\ntestArrayPad: signal,padding, padvalue, result", x,padding, padvalue, result);
		
		var padding = 3;  // number of complex signal values to be added, e.g. 3
		var padvalue = [-1, -2];
		var result = arrayPad(x, padding, padvalue);
		console.log("\ntestArrayPad: signal,padding, padvalue, result", x,padding, padvalue, result);
		
		var padding = 2;  // number of complex signal values to be added, e.g. 3
		var padvalue = ["a", "b"];
		var result = arrayPad(x, padding, padvalue);
		console.log("\ntestArrayPad: signal,padding, padvalue, result", x,padding, padvalue, result);
		
		console.log("\nStop testArrayPad");
		console.log("");
		};
	};


// ****************************************************************
// Following is based upon Mathematica's procedures ListConvolve and
// ListCorrelate. The length of the output equals the length of the input
// signal but padding and cyclic indexing are used to get the intuitive
// result. The following assumptions are made:
// 1) This procedure is intended for 2D data in the "universal" format;
// 2) The padded value is complex "padvalue", usually [0, 0].
// 3) Padding is to the right and bottom of the 2D array.
// 4) The size of the padding is derived from the kernel size which is
// (2n+1) x (2n+1) with the center at (n+1, n+1).
//
// Remember, "kernel" is a package
//		   kernel[0] = [rows,cols]
//		   kernel[1] = [complex coefficients in a 1D list]
//
// checked Friday, 24 February 2017

function listCorrelate2D(kernel, xdims, x)
	{
	var fdims = kernel[0].slice();
	var frows = fdims[0];
	var fcols = fdims[1];
	var h0 = Array(2*frows*fcols);  // flipped correlation filter
	h0.fill(0);  // initialize
	var h1 = split(h0,fdims);  // reformat to special 2D
	var filter = split(kernel[1].slice(),fdims);  // reformat to special 2D

	// flipping from h[i,j] to h[-i,-j]
	for (var j = 0; j < fcols; j++)
		{
		for (var i = 0; i < 2*frows; i += 2)
			h1[i][j] = filter[2*frows-i-2][fcols-j-1];
		for (var i = 1; i < 2*frows; i += 2)
			h1[i][j] = filter[2*frows-i][fcols-j-1];
		};
	var h = merge(h1,fdims);  // flipped correlation filter
	
	return listConvolve2D([fdims,h], xdims, x);
	};
	
	
function listConvolve2D(kernel, xdims, x)
	{
	var fdims = kernel[0].slice();
	var frows = fdims[0];
	var fcols = fdims[1];
	var filter = split(kernel[1].slice(),fdims);

	
	var xrows = xdims[0];
	var xcols = xdims[1];
	var input = split(x,xdims);  // reformat to special 2D
	
	// zero-pad the input array with 2D size of filter array
	var padvalue = [0, 0];
	var padded = imagePad(kernel, xdims, x, padvalue);
	var pdims = padded[0];
	var pInput = split(padded[1],pdims);  // return to special 2D format
	
	// but don't zero-pad the output array
	var out = Array(2*xrows*xcols);
	out.fill(0);  // initialize
	var odims = [xrows, xcols];
	var output = split(out, odims);  // reformat to special 2D

	var tempr = 0;
	var tempi = 0;
	var rindex = 0;  // rows index
	var cindex = 0;  // column index
	var pr = 0;  // padded real
	var pi = 0;  // padded imaginary
	var fr = 0;  // filter real
	var fi = 0;  // filter imaginary
	
	// This is the convolution
	for (var r = 0; r < 2*xrows; r += 2)  // through complex image rows
		{
		for (var c = 0; c < xcols; c++)  // through image columns
			{
			tempr = 0;
			tempi = 0;
			for (var i = 0; i < 2*frows; i += 2)  // through filter rows
				{
				for (var j = 0; j < fcols; j++)  // through filter columns
					{
					rindex = (r-i) % (2*xrows);  // cyclic index
					if ( rindex < 0 ) rindex += 2*xrows;  // special 2D complex format!
					cindex = (c-j) % xcols;  // cyclic index
					if ( cindex < 0 ) cindex += xcols;
					
					pr = pInput[rindex][cindex];
					pi = pInput[rindex+1][cindex];
					
					fr = filter[i][j];
					fi = filter[i+1][j];
					
					tempr += ( pr*fr - pi*fi );
					tempi += ( pi*fr + pr*fi );
					};
				};
			output[r][c] = tempr;
			output[r+1][c] = tempi;
			};
		};
	return merge(output,odims);
	};


function imagePad(kernel, xdims, x, padvalue)
	{
	var fdims = kernel[0].slice();
	var frows = fdims[0];
	var fcols = fdims[1];
	var filter = split(kernel[1].slice(),fdims);
	
	var xrows = xdims[0];
	var xcols = xdims[1];
	var input = split(x,xdims);  // reformat to special 2D
	
	// zero-pad the input array with 2D size of filter array
	var pLength = (xrows + frows)*(xcols + fcols);  // the new size
	var padded = Array(2*pLength);  // create array
	
	for (var j = 0; j < 2*pLength; j += 2)
		{
		padded[j] = padvalue[0];  // intialize real values
		padded[j+1] = padvalue[1];  // intialize imaginary values
		};
		
	var pdims = [(xrows + frows),(xcols + fcols)];  // 2D dimensions
	var pInput = split(padded, pdims);  // reformat to special 2D
	
	// the reformatting was done to make the following addressing easier
	for (var j = 0; j < 2*xrows; j++)
		for (var i = 0; i < xcols; i++)
			pInput[j][i] = input[j][i];  // copy input array to padded array
	
	return [pdims, merge(pInput,pdims)];  // return universal format
	};


function testImagePad(testChoice)
	{
		if(testChoice)
		{
		console.log("\nStart testImagePad");
		
		var rows = 7;
		var cols = 6;
		var idims = [rows,cols];
		var value = [3, 4];
		
		var constantImage = createArray(2*rows,cols);
		for (var i = 0; i < cols; i++)
			for (var j = 0; j < 2*rows; j++)
				{
				if (evenQ(j))
					constantImage[j][i] = 1;
				else
					constantImage[j][i] = 2;
				};
		var inputImage = merge(constantImage, idims);
		
	const localKernel0 =
		[ // Notice complex data package format
			[3,3],
			[0, 0, 0, 0, 0, 0,
			 -1, 0, 0, 0, 1, 0,
			 0, 0, 0, 0, 0, 0]
		];
		
		var kernel = localKernel0;  // remember, "kernel" is a package

		var result = imagePad(kernel, idims, inputImage, value);
		console.log("\ntestImagePad: kernel,idims, inputImage, result", split(kernel[1],kernel[0]),idims,split(inputImage,idims),result[0],split(result[1],result[0]),result[1]);
		
		var value = ["pr", "pi"];
		var constantImage = createArray(2*rows,cols);
		for (var i = 0; i < cols; i++)
			for (var j = 0; j < 2*rows; j++)
				{
				if (evenQ(j))
					constantImage[j][i] = "xr";
				else
					constantImage[j][i] = "xi";
				};
		var inputImage = merge(constantImage, idims);
		
		var result = imagePad(kernel, idims, inputImage, value);
		console.log("\ntestImagePad: kernel,idims, inputImage, result", split(kernel[1],kernel[0]),idims,split(inputImage,idims),result[0],split(result[1],result[0]),result[1]);
		
		console.log("\nStop testImagePad");
		console.log("");
		};
	};


	function testListConvolve2D(testChoice)
	{
		if(testChoice)
		{
		console.log("\nStart testListConvolve2D");
		var rows = 9;  // must be odd
		var cols = 7;  // must be odd
		var idims = [rows,cols];
		
		var impulseImage = createArray(2*rows,cols);
		for (var i = 0; i < cols; i++)
			for (var j = 0; j < 2*rows; j++)
				impulseImage[j][i] = 0;
		impulseImage[rows-1][(cols-1)/2] = 1;
		impulseImage[rows][(cols-1)/2] = 0;
		
		var stepImage = createArray(2*rows,cols);
		for (var i = 0; i < cols; i++)
			for (var j = 0; j < 2*rows; j++)
				{
				if (oddQ(j))
					stepImage[j][i] = 0;
				else if ( i < Math.floor(cols/2))
					stepImage[j][i] = 1;
				else
					stepImage[j][i] = 0;
				};
				
		var triImage = createArray(2*rows,cols);
		for (var i = 0; i < cols; i++)
			for (var j = 0; j < 2*rows; j++)
				{
				if (oddQ(j))
					triImage[j][i] = 0;
				else if ( i < cols-2)
					triImage[j][i] = i;
				else
					triImage[j][i] = 0;
				};
				
		var inputImage = impulseImage;
		
	const localKernel0 =
		[ // Notice complex data package format
			[3,3],
			[0, 0, 0, 0, 0, 0,
			 -1, 0, 0, 0, 1, 0,
			 0, 0, 0, 0, 0, 0]
		];
		
	const localKernel1 =
		[ // Notice complex data package format
			[3,3],
			[1, 0, 2, 0, 3, 0,
			 0, 0, 0, 0, 0, 0,
			 0, 0, 0, 0, 0, 0]
		];
		
		var kernel = localKernel1;  // remember, "kernel" is a package

		var x = merge(inputImage, idims);
		
		var y = listConvolve2D(kernel, idims, x);
		console.log("\nlistConvolve2D Identity: kernel,split(kernel[1],kernel[0]), split(x,idims), split(y,idims)", kernel,split(kernel[1],kernel[0]),split(x,idims), split(y,idims));
		

		kernel = sobelFilter;
		var y = listConvolve2D(kernel, idims, x);
		console.log("\nlistConvolve2D Sobel: kernel, split(kernel[1],kernel[0]), split(x,idims), split(y,idims)", kernel,split(kernel[1],kernel[0]),split(x,idims), split(y,idims));
		
		kernel = crossFilter;
		var y = listConvolve2D(kernel, idims, x);
		console.log("\nlistConvolve2D Cross: kernel, split(kernel[1],kernel[0]), split(x,idims), split(y,idims)", kernel,split(kernel[1],kernel[0]),split(x,idims), split(y,idims));
		

		console.log("\nStop testListConvolve2D");
		console.log("");
		};
	};


function testListCorrelate2D(testChoice)
	{
		if(testChoice)
		{
		console.log("\nStart testListCorrelate2D");
		var rows = 9;  // must be odd
		var cols = 7;  // must be odd
		var idims = [rows,cols];
		var impulseImage = createArray(2*rows,cols);
		for (var i = 0; i < cols; i++)
			for (var j = 0; j < 2*rows; j++)
				impulseImage[j][i] = 0;
		impulseImage[rows-1][(cols-1)/2] = 1;
		impulseImage[rows][(cols-1)/2] = 0;
		
		var stepImage = createArray(2*rows,cols);
		for (var i = 0; i < cols; i++)
			for (var j = 0; j < 2*rows; j++)
				{
				if (oddQ(j))
					stepImage[j][i] = 0;
				else if ( i < Math.floor(cols/2))
					stepImage[j][i] = 1;
				else
					stepImage[j][i] = 0;
				};
				
		var triImage = createArray(2*rows,cols);
		for (var i = 0; i < cols; i++)
			for (var j = 0; j < 2*rows; j++)
				{
				if (oddQ(j))
					triImage[j][i] = 0;
				else if ( i < cols-2)
					triImage[j][i] = i;
				else
					triImage[j][i] = 0;
				};
				
				
		var inputImage = triImage;
		
		const testFilter = 
			[ // Notice slightly MODIFIED complex data format
				[3,3],
				[1, 0, 0, 0, 0, 0,
				 0, 0, 1, 0, 0, 0, 
				 0, 0, 0, 0, 0, -1]
			];
		
	const localKernel0 =
		[ // Notice complex data package format
			[3,3],
			[0, 0, 0, 0, 0, 0,
			 -1, 0, 0, 0, 1, 0,
			 0, 0, 0, 0, 0, 0]
		];
		
	const localKernel1 =
		[ // Notice complex data package format
			[3,3],
			[1, 0, 2, 0, 3, 0,
			 0, 0, 0, 0, 0, 0,
			 0, 0, 0, 0, 0, 0]
		];
		
		var kernel = localKernel1;  // remember, "kernel" is a package
		var x = merge(inputImage, idims);
		var y = listCorrelate2D(kernel, idims, x);

		console.log("\nlistCorrelate2D Identity: kernel,split(kernel[1],kernel[0]), split(x,idims), split(y,idims)", kernel,split(kernel[1],kernel[0]),split(x,idims), split(y,idims));
		
		kernel = sobelFilter;
		var y = listCorrelate2D(kernel, idims, x);
		console.log("\nlistCorrelate2D Sobel: kernel, split(kernel[1],kernel[0]), split(x,idims), split(y,idims)", kernel,split(kernel[1],kernel[0]),split(x,idims), split(y,idims));
		
		kernel = crossFilter;
		var y = listCorrelate2D(kernel, idims, x);
		console.log("\nlistCorrelate2D Cross: kernel, split(kernel[1],kernel[0]), split(x,idims), split(y,idims)", kernel,split(kernel[1],kernel[0]),split(x,idims), split(y,idims));
		
		console.log("\nStop testListCorrelate2D");
		console.log("");
		};
	};


// ****************************************************************
// We want to multiply two Fourier spectra to achieve frequency domain
// filtering. Or we want to multiply two time (space) signals to achieve
// modulation. This procedure realizes this for complex arrays. These
// arrays must have the same length, niet waar?
//
// Checked Tuesday, 10 January 2017

function filterMult(array1, array2)
	{
	var nn1 = arrayDimension(array1);
	var nn2 = arrayDimension(array2);
	
	var ndim1 = nn1.length;
	var ndim2 = nn2.length;
	if ((ndim1 !== 1) || (ndim2 !== 1))
		throw "Arrays must be in universal format!";

	var length1 = nn1[0];
	var length2 = nn2[0];
	if (length1 !== length2) throw "Arrays must have same length!";
	
	var output = Array(length1)
	
	var a1r = 0;  // array 1 real
	var a1i = 0;  // array 1 imaginary
	var a2r = 0;  // array 2 real
	var a2i = 0;  // array 2 imaginary

	for (var i = 0; i < length1; i += 2)
		{
		a1r = array1[i];
		a1i = array1[i+1];
		
		a2r = array2[i];
		a2i = array2[i+1];
		
		output[i] = a1r*a2r - a1i*a2i;
		output[i+1] = a1i*a2r + a1r*a2i;
		};
	return output;
	};
	
function testFilterMult(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testFilterMult");
		
		// test if computation is done properly
		test1 = [7,0,2,-2,3,-3,0,1];
		test2 = [7,0,1/2,0,3,+3,0,1];
		console.log("\nfilterMult: ",test1,test2,filterMult(test1,test2));
		
		// test if length check is executed properly
		test1 = [1,0,2,0,3,0,4,0];
		test2 = [1,0,2,0,3,0];
//		console.log("\nfilterMult: ",test1,test2,filterMult(test1,test2));
		
		// test if format check is executed properly
		var test1 = nlist2D;
		var test2 = nlist2D;
//		console.log("\nfilterMult: ",test1,test2,filterMult(test1,test2));
		
		console.log("\nStop testFilterMult");
		console.log("");
		};
	};


// ****************************************************************
// To produce of copy a multi-dimensional array we use a recursive
// routine clone found at <http://stackoverflow.com/questions/2294703/multidimensional-array-cloning-using-javascript>
// correct usage is: var output = input.clone()
//
// Checked Tuesday, 10 January 2017

function clone(inputArray)
	{
	var nn = arrayDimension(inputArray);
	var ndim = nn.length;
	var rows = 1;  // default is 1D signal
	var cols = 0;
	if (ndim === 1)
		cols = nn[0];
	else if (ndim === 2)
		{
		rows = nn[0];  // but for 2D signal use this
		cols = nn[1];
		}
	else
		throw "Wrong number of dimensions!";

	var outputArray = createArray(rows,cols)
	for (var i = 0; i < rows; i++)
		{
		for (var j = 0; j < cols; j++)
			{
			outputArray[i][j] = inputArray[i][j];
			};
		};
	return outputArray;
	};
	
function testCloneCopy(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testCloneCopy");
		
		test = tdata2D1.slice();
		var resultClone = clone(tdata2D1);
		console.log("Clone: ",tdata2D1,test,resultClone);
		
		test[test.length/2] = 'a';
		
		console.log("Clone: ",tdata2D1,test,resultClone);
		};
	};


// ****************************************************************
// To fit data (x) to a function with parameters (a,b), (e.g. f(x) = ax + b).
// we use the Levenberg-Marquardt algorithm. This implementation comes from:
// <https://github.com/mljs/levenberg-marquardt>.
// 
// correct usage is: var output = input.clone()
//
// Checked 


// const Matrix = require('ml-matrix');

// Difference of the matrix function over the parameters
//
// @param {{x:Array<number>, y:Array<number>}} data - Array of points to
// fit in the format [x1, x2, ... ], [y1, y2, ... ]
// @param {Array<number>} params - Array of previous parameter values
// @param {number} gradientDifference - Adjustment for decrease the damping parameter
// @param {function} paramFunction - The parameters and returns a function
// with the independent variable as a parameter
// @return {Matrix}
//

function gradientFunction(data, params, gradientDifference, paramFunction)
	{
	'use strict';
	const n = paramFunction.length;
	const m = data.x.length;

	var ans = new Array(n);
	const func = paramFunction(...params);

	for (var param = 0; param < n; param++)
		{
		ans[param] = new Array(m);
		
		var auxParams = params.concat();
		auxParams[param] += gradientDifference;
		var funcParam = paramFunction(...auxParams);
		
		for (var point = 0; point < m; point++)
			{
			ans[param][point] = func(data.x[point]) - funcParam(data.x[point]);
			};
		};
	return new Matrix(ans);
	};


// Matrix function over the samples
// @ignore
// @param {{x:Array<number>, y:Array<number>}} data - Array of points to
// fit in the format [x1, x2, ... ], [y1, y2, ... ]
// @param {Array<number>} params - Array of previous parameter values
// @param {function} paramFunction - The parameters and returns a function
// with the independent variable as a parameter
// @return {Matrix}
//

function matrixFunction(data, params, paramFunction)
	{
	'use strict';
	const m = data.x.length;

	var ans = new Array(m);
	const func = paramFunction(...params);

	for (var point = 0; point < m; point++)
		{
		ans[point] = data.y[point] - func(data.x[point]);
		};
	return new Matrix([ans]);
	};


// Iteration for Levenberg-Marquardt
// @ignore
// @param {{x:Array<number>, y:Array<number>}} data - Array of points to
// fit in the format [x1, x2, ... ], [y1, y2, ... ]
// @param {Array<number>} params - Array of previous parameter values
// @param {number} damping - Levenberg-Marquardt parameter
// @param {number} gradientDifference - Adjustment for decrease the damping parameter
// @param {function} parameterizedFunction - The parameters and returns a function
// with the independent variable as a parameter
// @return {Array<number>}
//

function step(data, params, damping, gradientDifference, parameterizedFunction)
	{
	'use strict';
	var identity = Matrix.eye(parameterizedFunction.length)
		.mul(damping * gradientDifference * gradientDifference);
	var gradientFunc = gradientFunction
		(data, params, gradientDifference, parameterizedFunction);
	var matrixFunc = matrixFunction(data, params, parameterizedFunction).transpose();
	params = new Matrix([params]);

	var inverse = Matrix.inv(identity.add(gradientFunc.mmul(gradientFunc.transposeView())));
	params = params.sub
		(
((inverse.mmul(gradientFunc)).mmul(matrixFunc).mul(gradientDifference)).transposeView()
		);
	return params.to1DArray();
	};


// Calculate current error
// @ignore
// @param {{x:Array<number>, y:Array<number>}} data - Array of points to fit
// in the format [x1, x2, ... ], [y1, y2, ... ]
// @param {Array<number>} parameters - Array of current parameter values
// @param {function} parameterizedFunction - The parameters and returns a function
// with the independent variable as a parameter
// @return {number}
//

function errorCalculation(data, parameters, parameterizedFunction)
	{
	'use strict';
	var error = 0;
	const func = parameterizedFunction(...parameters);

	for (var i = 0; i < data.x.length; i++)
		{
		error += Math.abs(data.y[i] - func(data.x[i]));
		};
	return error;
	};


const defaultOptions =
	{
	damping: undefined,
	gradientDifference: 10e-2,
	initialValues: undefined,
	maxIterations: 100,
	errorTolerance: 10e-3
	};


// Curve fitting algorithm
// @param {{x:Array<number>, y:Array<number>}} data - Array of points to fit
// in the format [x1, x2, ... ], [y1, y2, ... ]
// @param {function} parameterizedFunction - The parameters and returns a function
// with the independent variable as a parameter
// @param {object} [options] - Options object
// @param {number} [options.damping = undefined] - Levenberg-Marquardt parameter
// @param {number} [options.gradientDifference = 10e-2] - Adjustment for decrease
// of the damping parameter
// @param {Array<number>} [options.initialValues = undefined] - Array of initial
// parameter values
// @param {number} [options.maxIterations = 100] - Maximum of allowed iterations
// @param {number} [options.errorTolerance = 10e-3] - Minimum uncertainty allowed
// for each point
// @return {{parameterValues: Array<number>, parameterError: number, iterations: number}}
//


function levenbergMarquardt(data, parameterizedFunction, options)
	{
	'use strict';
	// verify that damping is not undefined
	if ((!options) || (!options.damping) || (options.damping <= 0))
		{
		throw new TypeError('The damping option should be a positive number');
		};
	
	// assign default values
	options = Object.assign({}, defaultOptions, options);
	
	// fill with default value for initialValues
	if (!options.initialValues)
		{
		options.initialValues = new Array(parameterizedFunction.length);
		
		for (var i = 0; i < parameterizedFunction.length; i++)
			{
			options.initialValues[i] = 1;
			};
		};

	// check that the data have the correct format
	if (!data.x || !data.y)
		{
		throw new TypeError('The data parameter should have a x and y elements');
		}
		else if ((data.x.constructor !== Array) || (data.x.length < 2) ||
			   (data.y.constructor !== Array) || (data.y.length < 2))
			   	{
				throw new TypeError('The data parameter elements should be an array with more than 2 points');
				}

	const dataLen = data.x.length;
	if (dataLen !== data.y.length)
		{
		throw new RangeError('The data parameter elements should have the same size');
		}
	
	// initial parameters
	var parameters = options.initialValues;

	// check errorCalculation
	var error = errorCalculation(data, parameters, parameterizedFunction);
	var converged = error <= options.errorTolerance;

	for (var iteration = 0; (iteration < options.maxIterations) && !converged; iteration++)
		{
		// step function
		parameters = step(data, parameters, options.damping, 
				options.gradientDifference, parameterizedFunction);
		
		// reevaluate errorCalculation
		error = errorCalculation(data, parameters, parameterizedFunction);
		converged = error <= options.errorTolerance;
		};

	// return example
	return	{
			parameterValues: parameters,
			parameterError: error,
			iterations: iteration
			};
	};


function testLevenbergMarquardt(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testLevenbergMarquardt");
		
		function sinFunction(a, b)
			{
			return (t) => (a * Math.sin(b * t));
			};
		
		var length = 20;
		var data =
			{
			x: Array(length),
			y: Array(length)
			};

		var param0 = [2, 2];
		var sampleFunction = sinFunction(param0[0], param0[1]);
		
		for (var i = 0; i < length; i++)
			{
			data.x[i] = i;
			data.y[i] = sampleFunction(i);
			};
		
		const options =
			{
			damping: 0.001,
			initialValues: [3, 3]
			};

		var ans = LM(data, sinFunction, options);
//		var ans = levenbergMarquardt(data, sinFunction, options);
//		var ans = ["p1","p2"];
		console.log("Levenberg-Marquardt - data, parameters, result ",data, param0,ans);
		
		console.log("Stop testLevenbergMarquardt");
		};
	};


function testErrorCalculation(testChoice)
	{
	if(testChoice)
		{
		console.log("\nStart testErrorCalculation");
		
		function sinFunction(a, b, c)
			{
			return (t) => (a * Math.sin(b*t) + Math.exp(c*t));
			}
		const length = 20;
		let data =
			{
			x: Array(length),
			y: Array(length)
			};
		var param0 = [2,2,-3];
		let sampleFunction = sinFunction(param0[0], param0[1], param0[2]);
		for (let i = 0; i < length; i++)
			{
			data.x[i] = i;
			data.y[i] = sampleFunction(i);
			}
		var param1 = [2.1,2.1,-1];
		var ans0 =
			errorCalculation(data, param0, sinFunction);  // should be ≈ (0, 10e-3)
		var ans1 =
			errorCalculation(data, param1, sinFunction);  // should not be.≈(0, 10e-3)
		console.log("testErrorCalculation - true, estimate",param0, ans0);
		console.log("testErrorCalculation - true, estimate",param1, ans1);
		
		console.log("Stop testErrorCalculation");
		};
	};
/*
describe('errorCalculation test', function () {

    it('Simple case', function () {
		function sinFunction(a, b) {
			return (t) => (a * Math.sin(b * t));
		}


		errorCalculation(data, [2, 2], sinFunction).should.be.approximately(0, 10e-3);
		errorCalculation(data, [4, 4], sinFunction).should.not.be.approximately(0, 10e-3);
	});
});
*/


// ****************************************************************
// Image measurements based upon the command ImageMeasurements in
// Mathematica v 11.
// This version is for 1D "images" and 2D images. Note that although all
// computed values are scalars, they are returned in the "universal"
// format meaning: [real, imaginary].
//
// Further, for 2D signals there is a restriction. The number of rows MUST
// be a power of two, e.g. rows = 2^7 and the number of columns MUST also be
// a power of two, e,g.columns =  2^8. (See the routine "fourND" in
// "SSPfourier.js" for an explanation.). This means the data.length =
// 2*rows*columns = 2*(2^7)*(2^8) = 2^16. Caveat emptor!
//
// checked: Wednesday, 12 April 2017

function imageMeasurements(image, nn, dataType, drempel)

	{
	var imageInfo =
		{
		type: undefined,				// 'real', 'imaginary', or 'complex'
		rows: undefined,				// e.g. 512
		cols: undefined,				// e.g. 128
		maxIntensity: undefined,		// only for 'real' or 'imaginary'
		minIntensity: undefined,		// only for 'real' or 'imaginary'
		meanIntensity: undefined,		// may be 'complex'
		stndDevIntensity: undefined,	// may be 'complex'
		totalIntensity: undefined,		// may be 'complex'
		energyIntensity: undefined,		// may be 'complex'
		medianIntensity: undefined,		// only for 'real' or 'imaginary'
		threshold: undefined,			// only for 'real'
		area: undefined					// only for 'real'
		};
	
	imageInfo.type = dataType;
	imageInfo.threshold = drempel;
	
	var dims = nn.length;
	imageInfo.rows = 1;
	imageInfo.cols = nn[0];
	if (dims === 2)
		{
		imageInfo.rows = nn[0]
		imageInfo.cols = nn[1];
		};
	
	var dlength = 2*imageInfo.rows*imageInfo.cols;
	var max = -Infinity;
	var min = Infinity;
	
	imageInfo.meanIntensity = mean(image);
	imageInfo.stndDevIntensity = stdev(image);
	imageInfo.totalIntensity = total(image);
	
	if (imageInfo.type === 'real')
		{
		for (let i = 0; i < dlength; i += 2)
			{
			if (image[i] > max) max = image[i];
			if (image[i] < min) min = image[i];
			}
		imageInfo.maxIntensity = [max,0];
		imageInfo.minIntensity = [min,0];
		imageInfo.meanIntensity[1] = 0;
		imageInfo.totalIntensity[1] = 0;
		imageInfo.medianIntensity = median(image);
		};
		
	if (imageInfo.type === 'imaginary')
		{
		for (let i = 1; i < dlength; i += 2)
			{
			if (image[i] > max) max = image[i];
			if (image[i] < min) min = image[i];
			}
		imageInfo.maxIntensity = [0,max];
		imageInfo.minIntensity = [0,min];
		imageInfo.meanIntensity[0] = 0;
		imageInfo.totalIntensity[0] = 0;
		
		let tempdata = image.slice();
		let tempswap = 0;
		for (let i = 0; i < tempdata.length; i += 2)
				{
				// swap real and imaginary parts
				tempswap = tempdata[i];
				tempdata[i] = tempdata[i+1];
				tempdata[i+1] = tempswap;
				}
		let swappedmedian = median(tempdata)
		imageInfo.medianIntensity = [swappedmedian[1], swappedmedian[0]];
		};
		
	if (imageInfo.type === 'complex')
		{
		var value = 0;
		for (let i = 0; i < dlength; i += 2)
			{
			value = Math.hypot(image[i],image[i+1]);
			if (value > max) max = value;
			if (value < min) min = value;
			}
		imageInfo.maxIntensity = [max,0];
		imageInfo.minIntensity = [min,0];
		imageInfo.medianIntensity = undefined;
		};
	
	// energy is the sum of the magnitude-squared complex values
	var energy = 0;
	for (let i = 0; i < dlength; i += 2)
		energy += (Math.pow(image[i],2)+Math.pow(image[i+1],2));
	imageInfo.energyIntensity = [energy,0];
	
	// area = number of samples (complex) whose absolute value > threshold
	var oppervlakte = 0;
	for (let i = 0; i < dlength; i += 2)
		{
		oppervlakte += (Math.hypot(image[i],image[i+1]) > drempel) ? 1 : 0;
		};
	imageInfo.area = [oppervlakte, 0];
	
	return imageInfo;
	};


// ****************************************************************
// Root finding algorithm using Brent Method as described in Section 9.3 of
// "Numerical Recipes in C: but implemented at:
// GitHub <https://gist.github.com/borgar/3317728>
//
// Searches the interval from "lowerLimit" to "upperLimit"
// for a root (i.e., zero) of the function "func" with respect to 
// its first argument using Brent's method root-finding algorithm.
//
// @param {function} function for which the root is sought.
// @param {number} the lower point of the interval to be searched.
// @param {number} the upper point of the interval to be searched.
// @param {number} the desired accuracy (convergence tolerance).
// @param {number} the maximum number of iterations.
// @returns an estimate for the root within accuracy.
//
// Translated from zeroin.c in http://www.netlib.org/c/brent.shar.

/*
function uniroot ( func, lowerLimit, upperLimit, errorTol, maxIter )
	{
	var a = lowerLimit
	, b = upperLimit
	, c = a
	, fa = func(a)
	, fb = func(b)
	, fc = fa
	, s = 0
	, fs = 0
	, tol_act		// Actual tolerance
	, new_step		// Step at this iteration
	, prev_step		// Distance from the last but one to the last approximation
	, p				// Interpolation step is calculated in the form p/q;
					// division is delayed until the last moment
	, q
	;

	errorTol = errorTol || 0;
	maxIter  = maxIter  || 1000;

	while ( maxIter-- > 0 )
		{
		prev_step = b - a;
		if ( Math.abs(fc) < Math.abs(fb) )
			{
			// Swap data for b to be the best approximation
			a = b, b = c, c = a;
			fa = fb, fb = fc, fc = fa;
			}
		tol_act = 1e-15 * Math.abs(b) + errorTol / 2;
		new_step = (c-b)/2;
		if ( Math.abs(new_step) <= tol_act || fb === 0 )
			{
			return b; // Acceptable approx. is found
			}
		
		// Decide if the interpolation can be tried
		if ( Math.abs(prev_step) >= tol_act && Math.abs(fa) > Math.abs(fb) )
			{
			// If prev_step was large enough and was in true direction,
			// Interpolatiom may be tried
			var t1, cb, t2;
			cb = c - b;
			
			if ( a === c )
				{ // If we have only two distinct points linear
				// interpolation can only be applied
				t1 = fb / fa;
				p = cb * t1;
				q = 1.0 - t1;
				}
			else
				{ // Quadric inverse interpolation
				q = fa / fc, t1 = fb / fc, t2 = fb / fa;
				p = t2 * (cb * q * (q - t1) - (b - a) * (t1 - 1));
				q = (q - 1) * (t1 - 1) * (t2 - 1);
				}
			
			if ( p > 0 )
				{
				q = -q;  // p was calculated with the opposite sign; make p positive
				}
			else
				{
				p = -p;  // and assign possible minus to q
				}
			
			if ( p < ( 0.75 * cb * q - Math.abs( tol_act * q ) / 2 ) &&
				 p < Math.abs( prev_step * q / 2 ) )
					{ 
					// If (b + p / q) falls in [b,c] and isn't too large it is accepted
					new_step = p / q;
					}
			  // If p/q is too large then the bissection procedure can reduce [b,c] range to more extent
			}
		
		if ( Math.abs( new_step ) < tol_act )
			{ // Adjust the step to be not less than tolerance
			new_step = ( new_step > 0 ) ? tol_act : -tol_act;
			}
		
		a = b, fa = fb; // Save the previous approx.
		b += new_step, fb = func(b);  // Do step to a new approxim.
		
		if ( (fb > 0 && fc > 0) || (fb < 0 && fc < 0) )
			{
			c = a, fc = fa; // Adjust c for it to have a sign opposite to that of b
			}
		}
	};



var test_counter;
function f1 (x) { test_counter++; return (Math.pow(x,2)-1)*x - 5; }
function f2 (x) { test_counter++; return Math.cos(x)-x; }
function f3 (x) { test_counter++; return Math.sin(x)-x; }
function f4 (x) { test_counter++; return (x + 3) * Math.pow(x - 1, 2); }
[
  [f1, 2, 3],
  [f2, 2, 3],
  [f2, -1, 3],
  [f3, -1, 3], 
  [f4, -4, 4/3]
].forEach(function (args)
	{
  test_counter = 0;
  var root = uniroot.apply( pv, args );
  ;;;console.log( 'uniroot:', args.slice(1), root, test_counter );
	})

*/

// ****************************************************************
// ****************************************************************
// Various Javascript programs

// ****************************************************************
// Starting from three, equal-length, 1D lists {i, j, x}, produce a 2D array
// jovin such that jovin[i[k],j[k]] = v[k] for k = 1 to n.
//

function jovinArray(i, j, v, jovin)
	{
	// Note computation is in place
	// previous contents of array at [r,c] are overwritten
	var r = 0;
	var c = 0;
	for (let k = 0; k < i.length; k++)
		{
		r = i[k];
		c = j[k]; 
		jovin[r][c] = v[k];
		};
	return;
	};

function testJovin(testChoice,theLength) {
	if (testChoice) {
		console.log("\nStart testJovin");
		
		// final 2D array size
		var jrows = 2**12;
		var jcols = 2**12;
		// create 2D array of proper size
		var jovin = createArray(jrows, jcols);
		// initialize 2D array
		for (var r = 0; r < jrows; r++)
			for (var c = 0; c < jcols; c++)
				jovin[r][c] = 0;
		
		// number of sample values
		var myLength = theLength;
		var i = Array(myLength);  // row position
		var j = Array(myLength);  // column position
		var v = Array(myLength);  // integer datum value
		for (var k = 0; k < myLength; k++)
			{
			i[k] = randomInteger(0,jcols-1);  // random row coordinate
			j[k] = randomInteger(0,jrows-1);  // random column coordinate
			v[k] = randomInteger(0,(2**16)-1);  // random datum value (integer)
			};
		
		var iterations  = 1000;
		var start = performance.now();  // in ms
		for (var t = 0; t < iterations; ++t)
			jovinArray(i, j, v, jovin);
		var finish = performance.now();  // in ms
		var timed = finish - start;  // in ms
		console.log("testJovin:",
		"\t\t number of data points = "+myLength,
		"\t\t Duration [µs] = "+round(1e3*timed/iterations, rounding),
		"\t\t Duration [µs per vector element] = "+round(1e3*timed/(myLength*iterations), rounding));
		
//		jovinArray(i, j, v, jovin);
		console.log("testJovin - i, j, v, jovin: ", i, j, v, jovin);
		
		console.log("Stop testJovin");
		console.log("");
	};
};
