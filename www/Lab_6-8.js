// ****************************************************************
// Laboratory experiment 6.8
// i.t. young
// Thursday, 28 September 2017
// 
// ****************************************************************
// checked Tuesday, 3 October 2017

// global defined here (for now)
var dispState = 'SH';
var digiState = 'ANALOG';
var experiment = '6.8a';
var deltaN = 0;  // change in number of samples (zoom)
var power = 0;
var autoCorrData = [];
var autoFiltCorrData = [];
var rpsd = [];
var origFilter = [];
var nbins = 1;
var newSamples = 1;

var overkill = 2**3;  // power of two ≥ 4
var binningFactor = overkill*layoutCT.width;
var lengthMax = 2**Math.floor(Math.log2(100000/binningFactor));
lengthMax *= binningFactor;
var startSample = 0;
var stopSample = lengthMax;
var middleSample = (stopSample+startSample)/2;

var newSamples = lengthMax;
var sentence = '';
var sentence1 = '';
var sentence2 = '';

var rData = Array(2*lengthMax);
rData.fill(0);  // initialize
var yData =Array(lengthMax);
var yfData =Array(lengthMax);
var xDataDT = Array(lengthMax);
var tauDT = Array(lengthMax);
var xDataCT = Array(lengthMax);
var tauCT = Array(lengthMax);
var origFilter = Array(lengthMax);
var zeroes = Array(lengthMax);

var nn = [lengthMax]  // number of complex values
var ndim = 1;  // 1-D signal

var nLabels = 5;
var freqTicks = Array(nLabels);
var freqLabels = Array(nLabels);
var scaleLabel = 2;

var leftCoeffs = [1,0];
var rightCoeffs = youngBPFFilter;
var coefficients = [ leftCoeffs, rightCoeffs ];

var thisZoom = 0;
var thisAlg = 't';  // 't'=filter time domain & 'f'=filter freq. domain

// for normal (Gaussian) distribution
var myMean = 0;
var sigma = 0.25;
// for Uniform distribution
var upper = sigma*Math.sqrt(3);
var lower = - upper;
// for Exponential distribution
var lambda = sigma;

// templates in SSPplotting.js
// r=raw, f=filtered, s=signal, h=histogram, c=correlation, p=psd
// t=continuous time, d = discrete time
var rstg = cloneObj(CTPlot);  
var rsdg = cloneObj(DTPlot);
var rctg = cloneObj(CTPlot);
var rcdg = cloneObj(DTPlot);

var fstg = cloneObj(CTPlot);
var fsdg = cloneObj(DTPlot);
var fctg = cloneObj(CTPlot);
var fcdg = cloneObj(DTPlot);

var rstl = cloneObj(layoutCT);
var rsdl = cloneObj(layoutDT);
var rctl = cloneObj(layoutCT);
var rcdl = cloneObj(layoutDT);

var fstl = cloneObj(layoutCT);
var fsdl = cloneObj(layoutDT);
var fctl = cloneObj(layoutCT);
var fcdl = cloneObj(layoutDT);

var rhg = cloneObj(Histogram);
var fhg = cloneObj(Histogram);

var rhl = cloneObj(layoutH);
var fhl = cloneObj(layoutH);

var rpg = cloneObj(SpectrumPlot);
var fpg = cloneObj(SpectrumPlot);

var rpl = cloneObj(layoutSpect);
var fpl = cloneObj(layoutSpect);

// Define 12 display packages
var rst = [[rstg], rstl];
var rsd = [[rsdg], rsdl];

var fst = [[fstg], fstl];
var fsd = [[fsdg], fsdl];

var rct = [[rctg], rctl];
var rcd = [[rcdg], rcdl];

var fct = [[fctg], fctl];
var fcd = [[fcdg], fcdl];

var rp = [[rpg], rpl];
var fp = [[fpg], fpl];

var rh = [[rhg], rhl];
var fh = [[fhg], fhl];


// *****************************************************************
// Functions defined here

function process()
	{
	// Note use of global definitions for thisAlg, rdata, etc.
	// compute normalized autocorrelation GG through Fourier domain
	var rSpect = fft(rData, nn, ndim, FORWARD);  // signal spectrum
	var Spect,fpsdTemp,fpsd;
	var rpsdTemp = abssq(rSpect);  // power spectral density
	rpsd = uRs(rpsdTemp);
	autoCorrData = uRs(fft(rpsdTemp, nn, ndim, BACKWARD));
	
	if (thisAlg === 't')
		{
		//  filter in t-domain
		let fData = recurrenceFilter(coefficients, rData);  // t-domain filtering
		yfData = uRs(fData);  // data is real
		// compute normalized autocorrelation FF through Fourier domain
		Spect = fft(fData, nn, ndim, FORWARD);  // signal spectrum
		}
	else if (thisAlg === 'f')
		{
		//  filter in f-domain
		//  first, zero pad the filter array
		let zpLength = (rData.length - coefficients[1].length)/2;
		let zpFilter = arrayPad(coefficients[1], zpLength, [0,0])
		Spect = fft(zpFilter, nn, ndim, FORWARD);  // filter spectrum
		Spect = filterMult(rSpect, Spect);  // f-domain filtering
		yfData = uRs(fft(Spect, nn, ndim, BACKWARD));
		}
	else
		throw('Houston, we have a problem in process.');

	// compute power spectrum, normalized autocorrelation, |H(Omega)|**2;
	fpsdTemp = abssq(Spect);  // power spectral density
	fpsd = uRs(fpsdTemp);  // power spectral density
	autoFiltCorrData = uRs(fft(fpsdTemp, nn, ndim, BACKWARD));
	for(var i = 0; i < rpsd.length; i++)
		origFilter[i] = fpsd[i]/rpsd[i];

	// normalize autocorrelations & power spectra for display
	autoCorrData = normalizeDisp(autoCorrData);
	autoFiltCorrData = normalizeDisp(autoFiltCorrData);
	rpsd = normalizeDisp(rpsd);
	origFilter = normalizeDisp(origFilter);
	};

function rotateSymmetrics()
	{
	// Note use of global definitions for autoCorrData, etc.
	autoCorrData = rotateRight(autoCorrData,autoCorrData.length/2);
	autoFiltCorrData = rotateRight(autoFiltCorrData,autoFiltCorrData.length/2);
	rpsd = rotateRight(rpsd,rpsd.length/2);
	origFilter = rotateRight(origFilter,origFilter.length/2);
	};

function getData(target)
	{
	// activate displays
	if (target === 'acquire')
		{
		rstg.visible = true;
		rsdg.visible = true;
		rhg.visible = true;
		rctg.visible = true;
		rcdg.visible = true;
		rpg.visible = true;
		}
	else if (target === 'filter')
		{
		fstg.visible = true;
		fsdg.visible = true;
		fhg.visible = true;
		fctg.visible = true;
		fcdg.visible = true;
		fpg.visible = true;
		}
	else
		throw('Houston, we have a problem in getData.');
	myFunctionDisp(dispState,digiState,wA,wB,wC,wD);
	};

// Choose zoom factor with slider
function myFunctionZoom(val)
	{
	thisZoom = val;
	power = Math.pow(2,thisZoom);
	newSamples = Math.floor(lengthMax/power);
	if (oddQ(newSamples)) newSamples--;  // make it even
	sentence1 =  "N = "+newSamples+" samples";
	document.querySelector('#placeN1').value = "N = "+newSamples+" samples";
	sentence2 = "= "+d1round(1000*newSamples/sampFreq)+" [ms]";
	document.querySelector('#placeN2').value = sentence2;
	deltaN = newSamples/2;

	startSample = middleSample - deltaN;
	stopSample = middleSample + deltaN;

	// all of the following are real
	var xCDataNew = createArray(newSamples);  // time index
	var xDDataNew = createArray(newSamples);  // time index
	var yDataNew = createArray(newSamples);  // data
	var zeroesNew = createArray(newSamples);
	var zFdN = createArray(newSamples);  // zoom filt data
	var freqTicksNew = createArray(nLabels);
	// the following involve real, rotated arrays
	var tauCTnew = createArray(newSamples);  // correlation index
	var tauDTnew = createArray(newSamples);  // correlation index
	var zAuCoDN = createArray(newSamples);  // zoom raw phi
	var zAuFiCoDN = createArray(newSamples);  // zoom filt phi
	var zRpsdN = createArray(newSamples);  // zoom raw psd
	var zfilter = createArray(newSamples);  // zoom filter

	// here is where the "zooming" occurs
	for (var i = 0; i < newSamples; i++)
		{
		// first the signals
		xCDataNew[i] = xDataCT[i+startSample];
		xDDataNew[i] = xDataDT[i+startSample];
		yDataNew[i] = yData[i+startSample];  // raw data
		zFdN[i] = yfData[i+startSample];  // filtered data
		zeroesNew[i] = 0;
		// now the correlation functions and the power spectral densities
		zAuCoDN[i] = autoCorrData[i+startSample];  // raw ø
		zAuFiCoDN[i] = autoFiltCorrData[i+startSample];  // filt ø
		zRpsdN[i] = rpsd[i+startSample];  // raw psd
		zfilter[i] = origFilter[i+startSample];  // filt psd
		tauCTnew[i] = tauCT[i+startSample];
		tauDTnew[i] = tauDT[i+startSample];
		};

	rstg.x = binning(xCDataNew,binningFactor);
	fstg.x = rstg.x;
	rstg.y = binning(yDataNew,binningFactor);
	fstg.y = binning(zFdN,binningFactor);

	rsdg.x = xDDataNew;
	fsdg.x = rsdg.x;
	rsdg.y = yDataNew;
	fsdg.y = zFdN;
	
	rsdg.error_y.array = zeroesNew;
	fsdg.error_y.array = rsdg.error_y.array;
	rsdg.error_y.arrayminus = yDataNew;
	fsdg.error_y.arrayminus = zFdN;
	
	rctg.x = binning(tauCTnew,binningFactor);
	fctg.x = rctg.x;
	rctg.y = binning(zAuCoDN,binningFactor);
	fctg.y = binning(zAuFiCoDN,binningFactor);

	rcdg.x = tauDTnew;
	fcdg.x = rcdg.x;
	rcdg.y = zAuCoDN;
	fcdg.y = zAuFiCoDN;
	
	rcdg.error_y.array = zeroesNew;
	fcdg.error_y.array = rcdg.error_y.array;
	rcdg.error_y.arrayminus = zAuCoDN;
	fcdg.error_y.arrayminus = zAuFiCoDN;

	fpl.yaxis4.tickmode = 'log';
	if (newSamples <= maxDTdisplay)
		{
		if (newSamples > changeDTparams)
			{
			rsdg.error_y.thickness = DTlineThick;
			rsdg.marker.size = 3*DTlineThick;
		
			fsdg.error_y.thickness = DTlineThick;
			fsdg.marker.size = 3*DTlineThick;
		
			rcdg.error_y.thickness = DTlineThick;
			rcdg.marker.size = 3*DTlineThick;
		
			fcdg.error_y.thickness = DTlineThick;
			fcdg.marker.size = 3*DTlineThick;
			}
		else
			{
			rsdg.error_y.thickness = 2*DTlineThick;
			rsdg.marker.size = 5*DTlineThick;
		
			fsdg.error_y.thickness = 2*DTlineThick;
			fsdg.marker.size = 5*DTlineThick;
		
			rcdg.error_y.thickness = 2*DTlineThick;
			rcdg.marker.size = 5*DTlineThick;
		
			fcdg.error_y.thickness = 2*DTlineThick;
			fcdg.marker.size = 5*DTlineThick;
			};
//		rhl.yaxis3.tickmode = 'linear';
//		fhl.yaxis3.tickmode = 'linear';
		digiState = 'DIGITAL';
		}
	else
		{
//		rhl.yaxis3.tickmode = 'auto';
//		fhl.yaxis3.tickmode = 'auto';
		digiState = 'ANALOG';
		};

	rhg.x = yDataNew;
	nbins = Math.floor(Math.sqrt(rhg.x.length));
	rhg.nbinsx = nbins;
	rhg.xbins.start = Math.min(rhg.x);
	rhg.xbins.end = Math.max(rhg.x);

	fhg.x = zFdN;
	fhg.nbinsx = nbins;
	fhg.xbins.start = Math.min(fhg.x);
	fhg.xbins.end = Math.max(fhg.x);

	fhl.title = 'Filtered Amplitude Histogram';

	// adjust freqTicks & freqLabels
	scaleLabel = power;
	rpg.x = binning(xDDataNew,binningFactor);
	let rpgLength = rpg.x.length;
	
	for (var i = 0; i < nLabels; i++)
		{
		// choose label positions from binned data
		freqTicksNew[i] = rpg.x[Math.floor(i*rpgLength/(nLabels-1))];
		if (i === (nLabels-1)) freqTicksNew[i] = rpg.x[rpgLength-1];

		if (i < (nLabels-1)/2)
			{
			freqLabels[i] = negPiSymbol + '/'+(scaleLabel)
			if (oddQ(i))  freqLabels[i] = negPiSymbol + '/'+(2*scaleLabel)
			}
		else if (i === (nLabels-1)/2)
			freqLabels[i] = 0
		else
			{
			freqLabels[i] = posPiSymbol + '/'+(scaleLabel)
			if (oddQ(i))  freqLabels[i] = posPiSymbol + '/'+(2*scaleLabel)
			};
		};
	if(power === 1)  // special formatting case to avoid '/' symbol
		{
		freqLabels[0] = negPiSymbol;
		freqLabels[nLabels-1] = posPiSymbol;
		};

	rpg.y = normalizeDisp(binning(zRpsdN,binningFactor));
	rpl.xaxis4.tickvals = freqTicksNew;
	rpl.xaxis4.ticktext = freqLabels;

	fpg.x = rpg.x;
	fpg.y = binning(zfilter,binningFactor);  // no normalization
	fpl.xaxis4.tickvals = freqTicksNew;
	fpl.xaxis4.ticktext = freqLabels;

	myFunctionDisp(dispState,digiState,wA,wB,wC,wD);
	};

// *****************************************************************
// test signal synthesized here
// For the sampling frequency, see: SSPconstants.js
//
function prepareLab_6_8( )
	{
	sentence1 =  "N = "+newSamples+" samples";
	document.querySelector('#placeN1').value = "N = "+newSamples+" samples";
	sentence2 = "= "+d1round(1000*newSamples/sampFreq)+" [ms]";
	document.querySelector('#placeN2').value = sentence2;
	document.querySelector('#placeT').value = d1round(newSamples/sampFreq);
	
	// initialize data in complex format
	for (var i = 0; i < 2*lengthMax; i += 2)
		{
		rData[i] = randomGaussian(myMean,sigma);
		};

	yData = uRs(rData);  // data is real
	
	for (var i = startSample; i < stopSample; i++)
		{
		xDataDT[i] = i;  // discrete time index
		tauDT[i] = i - middleSample;
		xDataCT[i] = 1000*i/sampFreq;  // continuous time in [ms]
		tauCT[i] = 1000*(i - middleSample)/sampFreq;
		zeroes[i] = 0;
		};
	
	for (var i = 0; i < nLabels; i++)
		{
		freqTicks[i] = Math.floor(i*newSamples/(nLabels-1));
		if (i === (nLabels-1)) freqTicks[i] = newSamples-1;
		
		if (i < (nLabels-1)/2)
			freqLabels[i] = evenQ(i) ? negPiSymbol : negPiSymbol+'/'+scaleLabel
		else if (i === (nLabels-1)/2)
			freqLabels[i] = 0
		else
			freqLabels[i] = evenQ(i) ? posPiSymbol : posPiSymbol+'/'+scaleLabel
		};
		
		process();
		

// signal processing initiating ends here
// *****************************************************************
// now start displays
//
		// raw data: signal & amplitude histogram
		
		rstg.x = binning(xDataCT,binningFactor);
		rstg.y = binning(yData,binningFactor);
		rstl.annotations[0].text = 't [ms]';
		

		rhg.x = yData;
		nbins = Math.floor(Math.sqrt(rhg.x.length));
		rhg.nbinsx = nbins;
		rhg.xbins.start = Math.min(rhg.x);
		rhg.xbins.end = Math.max(rhg.x);
		rhl.annotations[0].text = 'signal amplitude';
		

		// filtered data: signal & amplitude histogram
		fstg.x = rstg.x;
		fstg.y = binning(yfData,binningFactor);
		fstl.title = 'Filtered Gaussian Noise';
		fstl.annotations[0].text = 't [ms]';


		fhg.x = yfData;
		fhg.nbinsx = nbins;
		fhg.xbins.start = Math.min(fhg.x);
		fhg.xbins.end = Math.max(fhg.x);
		fhl.title = 'Filtered Amplitude Histogram';
		fhl.annotations[0].text = 'signal amplitude';
		
		// raw & filtered data, corr. function, power spectral densities & filter
		
		var rotSamples = Math.floor(autoCorrData.length/2);

		// these "symmetrics" require binning because they are ANALOG
		rctg.y = rotateRight(binning(autoCorrData,binningFactor),binningFactor/2);
		rctg.x = binning(tauCT,binningFactor);
		rctl.title = 'Normalized '+phi+'<sub>GG</sub>('+tau+' = kT)';
		rctl.annotations[0].text = tau+' [ms]';
		autoCorrData = rotateRight(autoCorrData,rotSamples);

		fctg.y = rotateRight(binning(autoFiltCorrData,binningFactor),binningFactor/2);
		fctg.x = rctg.x;
		fctl.title = 'Normalized '+phi+'<sub>FF</sub>('+tau+' = kT)';
		fctl.annotations[0].text = rctl.annotations[0].text;
		autoFiltCorrData = rotateRight(autoFiltCorrData,rotSamples);

		rpg.y = rotateRight(binning(rpsd,binningFactor),binningFactor/2);
		rpg.x = binning(xDataDT,binningFactor);
		rpl.title = 'Normalized S<sub>GG</sub>('+Omega+')';
		freqTicks[nLabels-1] = rpg.x[rpg.x.length-1];
		rpl.xaxis4.tickvals = freqTicks;
		rpl.xaxis4.ticktext = freqLabels;
		rpl.annotations[0].text = blanks12+Omega;
		rpsd = rotateRight(rpsd,rotSamples);

		fpg.y = rotateRight(binning(origFilter,binningFactor),binningFactor/2);
		fpg.x = rpg.x;
		fpl.title = 'S<sub>FF</sub>('+Omega+')/S<sub>GG</sub>('+Omega+')';
		fpl.xaxis4.tickvals = rpl.xaxis4.tickvals;
		fpl.xaxis4.ticktext = rpl.xaxis4.ticktext;
		fpl.annotations[0].text = rpl.annotations[0].text;
		origFilter = rotateRight(origFilter,rotSamples);

		// these do not require binning because they are DIGITAL
		rsdg.x = xDataDT;
		rsdg.y = yData;
		rsdg.error_y.array = zeroes;
		rsdg.error_y.arrayminus = yData;
		rsdl.annotations[0].text = 'n';

		fsdg.x = rsdg.x;
		fsdg.y = yfData;
		fsdg.error_y.array = rsdg.error_y.array;
		fsdg.error_y.arrayminus = yfData;
		fsdl.title = 'Filtered Gaussian Noise';
		fsdl.annotations[0].text = rsdl.annotations[0].text;

		rcdg.x = tauDT;
		rcdg.y = autoCorrData;
		rcdl.title = 'Normalized '+phi+'<sub>GG</sub>[k]';
		rcdl.annotations[0].text = 'k';

		fcdg.x = rcdg.x;
		fcdg.y = autoFiltCorrData;
		fcdl.title = 'Normalized '+phi+'<sub>FF</sub>[k]';
		fcdl.annotations[0].text = rcdl.annotations[0].text;

		deltaN = 0;  // change in number of samples (zoom)
		power = 0;

// Await interaction

	};
