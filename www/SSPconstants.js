// ****************************************************************
// SSPconstants.js
//
// Math and signal processing constants to be used in the SSP package
// i.t. young
// Tuesday, 5 July 2016
// Tuesday, 13 December 2016
// Monday, 29 May 2017
// Tuesday, 27 June 2017
// Wednesday, 28 June 2017
// Tuesday, 4 July 2017
// Sunday, 13 August 2017
// 

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


// ****************************************************************

var testMe = true; /* default */
var rounding = 0.001;

var FORWARD = 1;
var BACKWARD = -1;

// Digital filters
var leftCoeffs = [];
var rightCoeffs = [];
var dummyFilter = [ leftCoeffs, rightCoeffs ];

// Global constants
const tooSmall = 1e-10;
const tooLarge = 1e+10;
const perc1 = 0.01;
const perc5 = 0.05;
const perc10 = 0.10;
const epsilon = 0.0001;

const iTMAX = 100;    // Numerical Recipes gamma functions
const ePS = 3e-07;    // Numerical Recipes gamma functions
const fPMIN = 1e-30;  // Numerical Recipes gamma functions

// For the sampling frequency, see: 
// <https://documentation.apple.com/en/soundtrackpro/usermanual/index.html#chapter=B%26section=2%26tasks=true>
//
const sampFreq = 44100;  // Hz
const timeWindow = 3.0;  // max length of time signal in [s]

var lengthMax = timeWindow*sampFreq;  // 3.0s x 44.1 kHz sampling rate

// ****************************************************************

