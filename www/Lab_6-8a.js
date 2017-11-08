// ****************************************************************
// Laboratory experiment 6.8 (a)
// i.t. young
// Thursday, 28 September 2017
// 
// ****************************************************************
// checked Tuesday, 3 October 2017


function executeLab_6_8a( )
	{
	// critical values defined here (for now)
	dispState = 'SH';
	digiState = 'ANALOG';

	thisAlg = 't';
	leftCoeffs = [1,0];
	rightCoeffs = youngBPFFilter;
	coefficients = [ leftCoeffs, rightCoeffs ];
	rotateSymmetrics();
	myFunctionZoom(thisZoom);

// signal processing initiating ends here
// *****************************************************************
// now start displays

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
	
	// this is where the first (empty) display gets placed
	myFunctionDisp(dispState,digiState,wA,wB,wC,wD);
	

	// to be used below when switching dispState
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

	var deltaN = 0;  // change in number of samples (zoom)
	var power = 0;

// Await interaction
	};
