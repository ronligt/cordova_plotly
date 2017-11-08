// ****************************************************************
// SSPplotting.js
//
// Plotting routines to be used with SSP package
// This includes constants and data structures associated with plotly.js
// i.t. young
//
// Sunday, 13 August 2017
// Friday, 18 August 2017
// Tuesday, 19 September 2017
// Tuesday, 26 September 2017

// The "universal" signal format is characterized by an array name (e.g. "data"),
// an array "nn", and a number of dimensions "ndim". The "data" array takes data
// in the form of [r1,i1,r2,i2, ...,rn,in] where "r1" is the real part of a
// complex sample and "i1" is the imaginary part.
// 
// For a 1D signal, "nn" is an array with one element that gives the number of
// complex signal samples (e.g. nn[0] = 17). The number of dimensions
// ndim = nn.length. In this example ndim = 1. Note that because the data are
// complex data.length = 2*17 = 34. Signals that are 1D can have ANY length.

// For a 2D signal (image), "nn" is an array with two elements that gives
// the number of rows and columns where each sample is complex, that is,
// nn = [rows,columns]. In this case, the number of dimensions
// ndim = nn.length = 2. Obviously, rows - n[0] and columns = nn[1].

// For 2D signals there is a restriction. The number of rows MUST be a power
// of two, e.g. rows = 2^7 and the number of columns MUST also be a power of
// two, e,g.columns =  2^8. (See the routine "fourND" in "SSPfourier.js" for
// an explanation.). This means the data.length = 2*rows*columns = 2*(2^7)*(2^8)
// = 2^16. Caveat emptor!

// Added a formatting subroutine to be used with plotly.js. This seemed like
// the best place as it is not signal processing
//
// To make a textual annotation to display on a plotly.js graph, we assume
// the graph is in a normalized picture coordinates 0 ≤ x,y ≤ 1. Then:
//      theText is the annotation string. e.g. 'foo'
//      xCoord where the annotation text should be placed, e.g. x = 0
//      xPosition is the positioning of the text relative to x as in 'center'
//      yCoord where the annotation text should be placed, e.g. y = 0.5 or y = 1.1
//      yPosition is the positioning of the text relative to y as in 'top'
//

// ****************************************************************

const goldenRatio = 2/(1+Math.sqrt(5));

// INDEX to number of images / graphs in a row
const oneImage = 1-1;
const twoImages = 2-1;
const threeImages = 3-1;

const noMode = {displayModeBar: false};

// For Greek letters and other stuff, see: <https://brajeshwar.github.io/entities/>
const tau = '\u03C4';
const phi = '<i>\u03C6</i>';
const Omega = '\u03A9';
const blanks12 = '            ';

// For labeling axes
const negPiSymbol = '–π';
const posPiSymbol = '+π';


// ****************************************************************
// Parameters for various plotly.js routines
//
// get display dimensions
var availableWidth = window.innerWidth ||
	document.documentElement.clientWidth ||
	document.body.clientWidth;

// not used as yet but maybe in the future
var availableHeight = window.innerHeight ||
	document.documentElement.clientHeight ||
	document.body.clientHeight;

// also not used but monitor - (menus + docks)
var windowWidth = window.screen.availWidth;
var windowHeight = window.screen.availHeight;

var baseSize1 = Math.floor(0.80*availableWidth);  // one graph or image per row
var baseSize2 = Math.floor(0.68*availableWidth/2);  // two graphs or images per row
var baseSize3 = Math.floor(0.50*availableWidth/3);  // three graphs or images per row
var graphWidth = [baseSize1, baseSize2, baseSize3];  // graphs or images per row
var graphHeight = [baseSize1, baseSize2, baseSize3];  // graphs or images per row

var d1round = d3.format('.1f');
var d2round = d3.format('.2f');
var d3round = d3.format('.3f');

// font names for iOS can be found at <http://iosfonts.com>
// this includes bold and italic
var baseFontSize = 25;
var myTitleFont = 'Helvetica, Arial, GillSans, san-serif';
var myTitleSize = 0.8*baseFontSize;
var myTitleColor = 'black';
var annotateFontSize = 0.8*baseFontSize;

var axisLabelFont = 'Georgia, Palatino, Bodoni, serif';
var axisLabelSize = 0.6*baseFontSize;  // size of dimension label, i.e. 'time'
var axisLabelColor = 'grey';

var axisTitleFont = 'Georgia-Italic, Palatino-Italic, BodoniSvtyTwoITCTT-BookIta';
var axisTitleSize = 0.7*baseFontSize;  // size of tick numbers, i.e. '0.2'
var axisTitleColor = 'black';

var axisTickFont = 'Georgia, Palatino, Bodoni, serif';
var axisTickSize = 0.6*baseFontSize;  // size of tick numbers, i.e. '0.2'
var axisTickColor = 'grey';
var axisTicksPosition = 'outside';

var titleAboveAxis = 1.15;
var titleRightAxis = 1.05;
var titleBelowAxis = -0.2;

var colorBarTitlePosition = 'right';
var colorBarTitleFont = 'Georgia-Italic, Palatino-Italic, BodoniSvtyTwoITCTT-BookIta';
var colorBarTitleSize = 0.7*baseFontSize;
var colorBarTitleColor = 'black';

var DTmarker = 'circle';
var DTlineThick = 1;
var DTmarkerSize = 3*DTlineThick;
var DTmarkerColor = 'maroon';
var DTlineColor = '#8c8c8c';
var maxDTdisplay = 192;  // max samples in DT format, otherwise screen crowded
var changeDTparams = 48;  // max samples in "fat" format, otherwise screen crowded

var zeroline = 1;
var zeroLineThick = 1;
var scatterMarkerSize = 7;
var signalLineColor = 'steelblue';
var histoColor = '#009999';
var thinHistoBar = 1;
var thickHistoBar = 4;
var spectLineColor = 'DarkGreen';
var scatterMarkerColor = 'maroon';
var scatterLineColor = 'green';

// **************************************************************
// utility program to clone object. See:
// <https://stackoverflow.com/questions/728360/how-do-i-correctly-clone-a-javascript-object>
//
// Note this utility is recursive
//
//  ity Tuesday, 15 August 2017

function cloneObj(obj)
	{
	var copy;
	// Handle the 3 simple types, and null or undefined
	if (null == obj || "object" != typeof obj) return obj;
	
	// Handle Date
	if (obj instanceof Date)
		{
		copy = new Date();
		copy.setTime(obj.getTime());
		return copy;
		}
	// Handle Array
	if (obj instanceof Array)
		{
		copy = [];
		for (var i = 0, len = obj.length; i < len; i++)
			{
			copy[i] = cloneObj(obj[i]);
			}
		return copy;
		}
	// Handle Object
	if (obj instanceof Object)
		{
		copy = {};
		for (var attr in obj)
			{
			if (obj.hasOwnProperty(attr)) copy[attr] = cloneObj(obj[attr]);
			}
		return copy;
		}
	throw new Error("Unable to copy obj! Its type isn't supported.");
	};
	
var annotateTemplate =
	{
		xref: 'paper',
		yref: 'paper',
		x: 0.5,
		xanchor: 'center',
		y: titleBelowAxis,
		yanchor: 'top',
		text: '',
		font:
			{
			family: axisTitleFont,
			size: annotateFontSize,
			color: axisTitleColor
			},
		showarrow: false
	}

// **************************************************************
// Continuous-Time (CT) plotting
// Left vert. axis label: makeAnnotation(0,'center',1.05,'bottom','G(t)')
// Right horiz. axis label: Annotate(0.5,'center',titleBelowAxis,'top',' t [ms]'),
// Center top data: Annotate(0.5,'center',1.15,'top','N = '+xDataCT.length+' samples')

// ity Sunday, 13 August 2017


var CTPlot =  // continuous-time plot
	{
		type: 'scatter',
		visible: false,
		mode: 'lines',
		line: {shape: 'spline'},
		hoverinfo: 'none',
		marker: { size: 3, color: signalLineColor },
		xaxis: 'x1',
		yaxis: 'y1',
		x: [],
		y: [],
	};

var layoutCT =
	{
	showlegend: false,
	
	title :'Gaussian Noise',
	titlefont:
		{
		family: myTitleFont,
		size: myTitleSize,
		color: myTitleColor
		},
	
	xaxis1:
		{
		titlefont:
			{
			family: axisTitleFont,
			size: axisTitleSize,
			color: axisTitleColor,
			},
		zeroline: false,
		zerolinewidth: zeroLineThick,
		showline: true,
		showticklabels: true,
		ticks: axisTicksPosition,
		tickfont:
			{
			family: axisTickFont,
			size: axisTickSize,
			color: axisTickColor
			},
		autotick: true,
		},
	
	yaxis1:
		{
		titlefont:
			{
			family: axisTitleFont,
			size: axisTitleSize,
			color: axisTitleColor,
			},
		zeroline: true,
		zerolinewidth: zeroLineThick,
		showline: true,
		showticklabels: true,
		ticks: axisTicksPosition,
		tickfont:
			{
			family: axisTickFont,
			size: axisTickSize,
			color: axisTickColor
			},
		type: '',
		autorange: true,
		},
		
	width: Math.floor(1.3*graphWidth[twoImages]),
	height: Math.floor(1.3*graphHeight[twoImages]),

	annotations:
		[annotateTemplate, annotateTemplate, annotateTemplate],
	};


// **************************************************************
// Discrete-Time (DT) Layout
// Left vert. axis label: Annotate(0,'center',1.05,'bottom','G[n]')
// Right horiz. axis label: Annotate(0.5,'center',1.5*titleBelowAxis,'bottom','n'),
// Center top data: Annotate(0.5,'center',1.15,'top','N = ' +xDataDT.length + ' samples')
//
// ity Sunday, 13 August 2017

var DTPlot =  // discrete-time plot
	{
	type: 'scatter',
	visible: false,
	mode: 'markers',
	hoverinfo: 'none',
	marker: {symbol: DTmarker, size: DTmarkerSize, color: DTmarkerColor},
	error_y: {
		type: 'data',
		thickness: DTlineThick,
		color: DTlineColor,
		symmetric: false,
		array: [], // use zero-filled array
		arrayminus: [], // re-use y values
		width: 0 // hide cross-bars at end of error bars
		},
	xaxis: 'x2',
	yaxis: 'y2',
	x: [],
	y: [],
	};

var layoutDT =
	{
		showlegend: false,
		
		title :'Gaussian Noise',
		titlefont:
			{
			family: myTitleFont,
			size: myTitleSize,
			color: myTitleColor
			},
		
		xaxis2:
			{
			titlefont:
				{
				family: axisTitleFont,
				size: axisTitleSize,
				color: axisTitleColor,
				},
			zeroline: false,
			zerolinewidth: zeroLineThick,
			showline: true,
			showticklabels: true,
			ticks: axisTicksPosition,
			tickfont:
				{
				family: axisTickFont,
				size: axisTickSize,
				color: axisTickColor
				},
			tickformat: '.0f',
			autotick: true,
			},
		
		yaxis2:
			{
			titlefont:
				{
				family: axisTitleFont,
				size: axisTitleSize,
				color: axisTitleColor,
				},
				
			zeroline: true,
			zerolinewidth: zeroLineThick,
			showline: true,
			showticklabels: true,
			ticks: axisTicksPosition,
			tickfont:
				{
				family: axisTickFont,
				size: axisTickSize,
				color: axisTickColor
				},
			type: '',
			autorange: true,
			},
			
		width: Math.floor(1.3*graphWidth[twoImages]),
		height: Math.floor(1.3*graphHeight[twoImages]),
		
		annotations:
			[annotateTemplate, annotateTemplate, annotateTemplate],
		};


// **************************************************************
// Histogram Layout
//
// Left vert. axis label: Annotate(0,'center',1.05,'bottom','n(G)')
// Right horiz. axis label: Annotate(titleRightAxis,'left',0,'top',' G')
//
// ity Sunday, 13 August 2017

var Histogram =
	{
	type: 'histogram',
	visible: false,
	cumulative: { enabled : false },
	histfunc: 'count',
	histnorm: '',
	hoverinfo: 'none',
	marker: {color: histoColor, line: {color: histoColor, width: 1} },
	xaxis: 'x3',
	yaxis: 'y3',
	autobinx: false,
	nbinsx: 256,
	xbins: { start: 0, end: 256, },
	x: [],
	};


var layoutH =
		{
		showlegend: false,
		
		title :'Amplitude Histogram',
		titlefont:
			{
			family: myTitleFont,
			size: myTitleSize,
			color: myTitleColor
			},
		
		bargap: 0.3,
		
		xaxis3:
			{
			titlefont:
				{
				family: axisTitleFont,
				size: axisTitleSize,
				color: axisTitleColor,
				},
			zeroline: false,
			zerolinewidth: zeroLineThick,
			showline: true,
			showticklabels: true,
			ticks: axisTicksPosition,
			tickfont:
				{
				family: axisTickFont,
				size: axisTickSize,
				color: axisTickColor
				},
			},
		
		yaxis3:
			{
			titlefont:
				{
				family: axisTitleFont,
				size: axisTitleSize,
				color: axisTitleColor,
				},
			zeroline: false,
			zerolinewidth: zeroLineThick,
			showline: true,
			showticklabels: true,
			ticks: axisTicksPosition,
			tickmode: 'auto',
			tickfont:
				{
				family: axisTickFont,
				size: axisTickSize,
				color: axisTickColor
				},
			},
			
		width: Math.floor(1.3*graphWidth[twoImages]),
		height: Math.floor(1.3*graphHeight[twoImages]),
		
		annotations:
			[annotateTemplate, annotateTemplate, annotateTemplate]
	};


// **************************************************************
// Frequency domain plotting
// Left vert. axis label: makeAnnotation(0,'center',1.05,'bottom','|G(f)|')
// Right horiz. axis label: Annotate(0.5,'center',titleBelowAxis,'top',' f [Hz]'),
// Center top data: Annotate(0.5,'center',1.15,'top','N = '+xDataCT.length+' samples')

// ity Sunday, Tuesday, 15 August 2017


var SpectrumPlot =  // continuous-time plot
	{
		type: 'scatter',
		visible: false,
		mode: 'lines',
		line: {shape: 'spline'},
		hoverinfo: 'none',
		marker: { size: 3, color: signalLineColor },
		xaxis: 'x4',
		yaxis: 'y4',
		x: [],
		y: [],
	};

var layoutSpect =
	{
		showlegend: false,
		
		title :'Fourier Spectrum',
		titlefont:
			{
			family: myTitleFont,
			size: myTitleSize,
			color: myTitleColor
			},
		
		xaxis4:
			{
			titlefont:
				{
				family: axisTitleFont,
				size: axisTitleSize,
				color: axisTitleColor,
				},
			zeroline: false,
			zerolinewidth: zeroLineThick,
			showline: true,
			showticklabels: true,
			ticks: axisTicksPosition,
			tickfont:
				{
				family: axisTickFont,
				size: axisTickSize,
				color: axisTickColor
				},
			tickvals: [],
			ticktext: [],
			type: '',
			},
		
		yaxis4:
			{
			titlefont:
				{
				family: axisTitleFont,
				size: axisTitleSize,
				color: axisTitleColor,
				},
			zeroline: true,
			zerolinewidth: zeroLineThick,
			showline: true,
			showticklabels: true,
			ticks: axisTicksPosition,
			tickfont:
				{
				family: axisTickFont,
				size: axisTickSize,
				color: axisTickColor
				},
			type: 'log',
			},
			
		width: Math.floor(1.3*graphWidth[twoImages]),
		height: Math.floor(1.3*graphHeight[twoImages]),
	
		annotations:
			[annotateTemplate, annotateTemplate, annotateTemplate]
	};


// **************************************************************
// Blank plot
// Left vert. axis label: makeAnnotation(0,'center',1.05,'bottom','|G(f)|')
// Right horiz. axis label: Annotate(0.5,'center',titleBelowAxis,'top',' f [Hz]'),
// Center top data: Annotate(0.5,'center',1.15,'top','N = '+xDataCT.length+' samples')

// ity Wednesday, 13 September 2017


var BlankPlot =  // discrete-time plot
	{
		type: 'scatter',
		visible: false,
		hoverinfo: 'none',
		xaxis: 'x5',
		yaxis: 'y5',
		x: [],
		y: [],
	};

var layoutBlank =
	{
		showlegend: false,
	
		title :'',
		titlefont:
			{
			family: myTitleFont,
			size: myTitleSize,
			color: myTitleColor
			},
	
		xaxis5:
			{
			zeroline: false,
			zerolinewidth: 0,
			showline: false,
			showticklabels: false,
			showspikes: false,
			showgrid: false,
			},
	
		yaxis5:
			{
			zeroline: false,
			zerolinewidth: 0,
			showline: false,
			showticklabels: false,
			showgrid: false,
			},
		
		width: Math.floor(1.3*graphWidth[twoImages]),
		height: Math.floor(1.2*graphHeight[twoImages]),
	
		annotations:
			[annotateTemplate, annotateTemplate, annotateTemplate]
	};


// **************************************************************
// Bar Plot Layout
//
// Left vert. axis label: Annotate(0,'center',1.05,'bottom','n(G)')
// Right horiz. axis label: Annotate(titleRightAxis,'left',0,'top',' G')
//
// ity Wednesday, 13 September 2017


var BarPlot =  // bar plot of two counts, n(H) and n(T)
	{
		type: 'bar',
		visible: false,
		text: [],
		textposition: 'auto',
		textfont: {color: 'white'},
		hoverinfo: 'none',
		marker: {
			color: ['rgb(70, 130, 180)', 'rgb(180, 70, 70))'],
			line: { color: 'black', width: 1.5 }, 
			},
		xaxis: 'x6',
		yaxis: 'y6',
		x: ['Heads', 'Tails'],
		y: [],
	};


var layoutBar =
	{
		showlegend: false,
	
		title :'Histogram of Heads & Tails',
		titlefont:
			{
			family: myTitleFont,
			size: myTitleSize,
			color: myTitleColor
			},
	
		xaxis6:
			{
			titlefont:
				{
				family: axisTitleFont,
				size: axisTitleSize,
				color: axisTitleColor,
				},
			tickfont:
				{
				family: axisTickFont,
				size: axisTickSize,
				color: axisTickColor
				},
			},
	
		yaxis6:
			{
			titlefont:
				{
				family: axisTitleFont,
				size: axisTitleSize,
				color: axisTitleColor,
				},
			
			zeroline: true,
			zerolinewidth: zeroLineThick,
			showline: true,
			showticklabels: true,
			ticks: axisTicksPosition,
			tickfont:
				{
				family: axisTickFont,
				size: axisTickSize,
				color: axisTickColor
				},
			autotick: true,
			},
		
		width: Math.floor(1.3*graphWidth[twoImages]),
		height: Math.floor(1.2*graphHeight[twoImages]),
		
		annotations:
			[annotateTemplate, annotateTemplate, annotateTemplate]
	};


// **************************************************************
// Special Image Plot Layout
//
// Left vert. axis label: Annotate(0,'center',1.05,'bottom','n(G)')
// Right horiz. axis label: Annotate(titleRightAxis,'left',0,'top',' G')
//
// ity Wednesday, 13 September 2017


var ImagePlot =  // bar plot of two counts, n(H) and n(T)
	{
		type: 'scatter',
		visible: false,
		mode: 'markers',
		hoverinfo: 'none',
		marker: {
			symbol: 'circle',
			size: 9*DTmarkerSize,
			color: 'rgba(70, 130, 180,0.75)'
			},
		xaxis: 'x7',
		yaxis: 'y7',
		x: [],
		y: [],
	};


var layoutImage =
	{
		showlegend: false,
		
		title :'Show me the money',
		titlefont:
			{
			family: myTitleFont,
			size: myTitleSize,
			color: myTitleColor
			},
		
		xaxis7:
			{
			titlefont:
				{
				family: axisTitleFont,
				size: axisTitleSize,
				color: axisTitleColor,
				},
			zeroline: false,
			zerolinewidth: zeroLineThick,
			showline: false,
			showticklabels: false,
			showspikes: false,
			ticks: axisTicksPosition,
			tickfont:
				{
				family: axisTickFont,
				size: axisTickSize,
				color: axisTickColor
				},
			tickvals: [0,1,2,3,4,5,6,7,8],
			},
		
		yaxis7:
			{
			titlefont:
				{
				family: axisTitleFont,
				size: axisTitleSize,
				color: axisTitleColor,
				},
				
			zeroline: false,
			zerolinewidth: zeroLineThick,
			showline: false,
			showticklabels: false,
			ticks: axisTicksPosition,
			tickfont:
				{
				family: axisTickFont,
				size: axisTickSize,
				color: axisTickColor
				},
			
			autorange: true,
			tickvals: [0,1,2,3,4,5,6,7,8],
		},
		
		width: Math.floor(1.3*graphWidth[twoImages]),
		height: Math.floor(1.2*graphHeight[twoImages]),
		
		annotations:
			[annotateTemplate, annotateTemplate, annotateTemplate]
		};


// **************************************************************
// 3D Surface Plot
//
// Left vert. axis label: Annotate(0,'center',1.05,'bottom','n(G)')
// Right horiz. axis label: Annotate(titleRightAxis,'left',0,'top',' G')
//
// ity Monday, 18 September 2017



var surfacePlot =  // surface plot
	{
		type: 'surface',
		visible: false,
		z: [],
		opacity: 1,
		colorscale: [],
		hoverinfo: 'none',
		showscale: true,
		colorbar: 
			{
			title: '',
			titleside: colorBarTitlePosition,
			titlefont:
				{
				family: colorBarTitleFont,
				size: colorBarTitleSize,
				color: colorBarTitleColor
				},
			tickfont:
				{
				family: axisTickFont,
				size: axisTickSize,
				color: axisTickColor
				},
			},
		};

var surfaceLayout =
	{
		showlegend: false,
		title:'',
		titlefont:
			{
			family: myTitleFont,
			size: myTitleSize,
			color: myTitleColor
			},
		
		scene:
			{
			xaxis:
				{
				title: '',
				titlefont:
					{
					family: axisTitleFont,
					size: axisTitleSize,
					color: axisTitleColor,
					},
				zeroline: true,
				zerolinewidth: zeroline,
				showline: true,
				showticklabels: true,
				ticks: axisTicksPosition,
				tickfont:
					{
					family: axisTickFont,
					size: axisTickSize,
					color: axisTickColor
					},
				tickvals: [],
				ticktext: [],
				},
			yaxis:
				{
				title: '',
				titlefont:
					{
					family: axisTitleFont,
					size: axisTitleSize,
					color: axisTitleColor,
					},
				zeroline: true,
				zerolinewidth: zeroline,
				showline: true,
				showticklabels: true,
				ticks: axisTicksPosition,
				tickfont:
					{
					family: axisTickFont,
					size: axisTickSize,
					color: axisTickColor
					},
				tickvals: [],
				ticktext: [],
				},
			zaxis:
				{
				title: '',
				type: '',
				titlefont:
					{
					family: axisTitleFont,
					size: axisTitleSize,
					color: axisTitleColor,
					},
				zeroline: true,
				zerolinewidth: zeroline,
				showticklabels: true,
				ticks: axisTicksPosition,
				tickfont:
					{
					family: axisTickFont,
					size: axisTickSize,
					color: axisTickColor
					},
				showline: true,
				},
			},
		
		width: Math.floor(graphWidth[oneImage]),
		height: Math.floor(graphHeight[oneImage]),
		
		annotations:
			[annotateTemplate, annotateTemplate, annotateTemplate]
	};


// **************************************************************
// Bin data for display.
//
// ity Tuesday, 29 August 2017
//
// The number of pixels in the horizontal (screen) direction will be limited
// (e.g. 120). The number of data points could be large (e.g. 6000). 
// We do not (necessarily) need to reduce the 6000 points to 120; the
// plotting program can do that for us. To visually emphasize the excess of data
// points over pixel width we use "overkill", (e.g. overkill = 8). The 
// graphing program (e.g. plotly.js) will take care of the rest. We,
// therefore, subsample the original "data" to produce the data samples
// that will be passed to the plotting program. This reduced number of
// data samples is targetLength = "overkill*layout.width",
// (e.g. targetLength = 8*120 = 960). The reason for all of this is that
// binning with the graphing program can be excruciatingly slow.
// So we do the heavy lifting.


			function binning(data,targetLength)
				{
				let step = Math.floor(data.length/targetLength);
				if (step > 1)
					{
					let result = createArray(targetLength);
					for (var i = 0; i < targetLength; i++ )
						result[i] = data[step*i];  // sampling step
					return result;
					}
				else
					{
					return data;
					};
				};


// **************************************************************
// Normalize data for display.
//
// ity Tuesday, 29 August 2017
//
// display data are simple lists. Numerical lists are normalized so
// that the maximum value is 1

function normalizeDisp(array)
	{
	var n = array.length;
	// first find the absolute maximum
	let result = Array(n);
	let maxval = 0;
	let temp = 0;
	for(var i = 0; i < n; i++)
		{
		temp = Math.abs(array[i]) ;
			if (temp > maxval)  maxval = temp;
		};
		
	// now normalize except if maxval === 0
	if(maxval !== 0)
		{
		for(var i = 0; i < n; i++)
			result[i] = array[i]/maxval;
		};
	return result;
	};


// **************************************************************
// Choose display mode with radio buttons
// wA = window A, wB = window B, etc.

function myFunctionDisp(val,target,wA,wB,wC,wD)
	{
	// display CT signal and histogram
	dispState = val;
	if (digiState === 'ANALOG' && dispState === 'SH')
		{
		Plotly.newPlot(wA, rst[0], rst[1], noMode);  // CT signal
		Plotly.newPlot(wB, rh[0], rh[1], noMode);  // histogram
		Plotly.newPlot(wC, fst[0], fst[1], noMode);  // CT signal
		Plotly.newPlot(wD, fh[0], fh[1], noMode);  // histogram
		}
	// display DT signal and histogram
	else if (digiState === 'DIGITAL' && dispState === 'SH')
		{
		Plotly.newPlot(wA, rsd[0], rsd[1], noMode);  // DT signal
		Plotly.newPlot(wB, rh[0], rh[1], noMode);  // histogram
		Plotly.newPlot(wC, fsd[0], fsd[1], noMode);  // DT signal
		Plotly.newPlot(wD, fh[0], fh[1], noMode);  // histogram
		}
	// display CT correlation and power spectrum
	else if (digiState === 'ANALOG' && dispState === 'CP')  
		{
		Plotly.newPlot(wA, rct[0], rct[1], noMode);  // CT correlation
		Plotly.newPlot(wB, rp[0], rp[1], noMode);  // power spectral density
		Plotly.newPlot(wC, fct[0], fct[1], noMode);  // CT correlation
		Plotly.newPlot(wD, fp[0], fp[1], noMode);  // power spectral density
		}
	// display DT correlation and power spectrum
	else if (digiState === 'DIGITAL' && dispState === 'CP')
		{
		Plotly.newPlot(wA, rcd[0], rcd[1], noMode);  // DT correlation
		Plotly.newPlot(wB, rp[0], rp[1], noMode);  // power spectral density
		Plotly.newPlot(wC, fcd[0], fcd[1], noMode);  // DT correlation
		Plotly.newPlot(wD, fp[0], fp[1], noMode);  // power spectral density
		}
	else
		throw('Houston, we have a problem in myFunctionDisp.');
	};


// **************************************************************

// Play sound with button
function myFunctionPlaySound(val,target)
	{
	if (target === 'PlayOrig')  // display signal & histogram
		{
		let valS = 'PlayOrig';
		// call to play sound goes here
		console.log('Playback button pressed = ', valS);
		}
	else if (target === 'PlayFilt')  // display correlation and power spectrum
		{
		let valS = 'PlayFilt';
		// call to play sound goes here
		console.log('Playback button pressed = ', valS);
		}
	else 
		throw('Houston, we have a problem in myFunctionPlaySound.');
	};

