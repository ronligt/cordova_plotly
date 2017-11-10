// Default values
var max_time =         1.5; // maximum recording time [seconds]
var default_time =     0.1 * max_time; // default time
var step_time =        0.001 * max_time; // slider steps
var sample_rate =      44100; // samples/second
var max_num_elem =     Math.floor(max_time * sample_rate); // maximum number of samples
var default_num_elem = Math.floor(default_time * sample_rate);
var step_sample =      0.001 * max_num_elem;

// Adjusting values in html
// - slider with time
// $("#slide1").attr('min', step_time);
// $("#slide1").attr('step', step_time);
// $("#slide1").attr('max', max_time);
// $("#slide1").attr('value', default_time);
// $("#slide1value").html(default_time);
// - slider with samples
$("#slide1").attr('min', step_sample);
$("#slide1").attr('step', step_sample);
$("#slide1").attr('max', max_num_elem);
$("#slide1").attr('value', default_num_elem);
$("#slide1value").html(default_num_elem);

// SSP stuff
var startSample = 0;
var stopSample = max_num_elem;
var middleSample = (stopSample+startSample)/2;

var nn = [max_num_elem*2];  // number of complex values
var ndim = 1;  // 1-D signal

var leftCoeffs = [1,0];
var rightCoeffs = youngBPFFilter;
var coefficients = [ leftCoeffs, rightCoeffs ];

var domRect = document.querySelector('#rawdata').getBoundingClientRect();
var div_width = domRect.width;

// Initialise data
var time = Array(max_num_elem);
// var timeBuffer = new ArrayBuffer(max_num_elem*Float32Array.BYTES_PER_ELEMENT);
// var timeView = new Float32Array(timeBuffer);
var tau = Array(max_num_elem);
var raw_data = Array(max_num_elem);
// var rawBuffer = new ArrayBuffer(max_num_elem*Float32Array.BYTES_PER_ELEMENT);
// var rawView = new Float32Array(rawBuffer);
var complex_data = Array(max_num_elem*2);

// var trace_raw = {
//   name: 'raw',
//   x: time,
//   y: raw_data,
//   type: 'scatter',
//   mode: 'lines',
//   line: {
//     smoothing: 0,
//   },
//   opacity: 0.8,
// }

// Calculate filter in time domain
var complex_filtered_data = recurrenceFilter(coefficients, complex_data);  // t-domain filtering
var real_filtered_data = uRs(complex_filtered_data);  // data is real

// var trace_filter = {
//   name: 'filter',
//   x: time,
//   y: real_filtered_data,
//   type: 'scatter',
//   mode: 'lines',
//   line: {
//     smoothing: 0,
//   },
//   opacity: 0.8,
// }

// // Calculate filter in frequency domain
// var rSpect = fft(complex_data, nn, ndim, FORWARD);  // signal spectrum
// // var Spect,fpsdTemp,fpsd;
// var rpsdTemp = abssq(rSpect);  // power spectral density
// // var rpsd = uRs(rpsdTemp);
// var autoCorrData = uRs(fft(rpsdTemp, nn, ndim, BACKWARD));
// console.log(rSpect.length,rpsdTemp.length,autoCorrData.length);

// Initialise plots
var layout_data = {
  title: "Data",
  xaxis: {
    title: "time [ms]",
    // type: 'lineair',
    // range: [1,4],
    // type: "date",
    // rangeslider: {},
  },
  showLegend: false,
  margin: {
    l: 30,
    r: 30,
  },
};
var layout_histogram = {
  title: "Histogram",
  showLegend: false,
  margin: {
    l: 30,
    r: 30,
  },
};

// Draw plots with default number of elements
// - slider with time
graph(default_time);
// - slider with samples
graph(default_num_elem);

// Redraw plot after slider change
$('#slide1').change(function(){
  graph($(this).val());
});

// Draw plot function
// - slider with time
// function graph(end_time) {
// - slider with samples
function graph(num_elem) {
  // var slice_time = time.slice(0,num_elem);
  // var slice_time = timeView.slice(0,num_elem);
  // var slice_time = timeView.subarray(0, num_elem);
  // var slice_time = Array.from(slice_time);
  // var slice_raw_data = raw_data.slice(0,num_elem);
  // var slice_raw = rawView.subarray(0, num_elem);
  // var slice_raw = Array.from(slice_raw);
  // var slice_real_filtered_data = real_filtered_data.slice(0,num_elem);

  var slice_time=[], slice_raw_data=[], slice_real_filtered_data=[];

  var step = 1;
  var max_width = 4*num_elem/max_num_elem*div_width;
  if (num_elem > max_width) {
    step = Math.ceil(num_elem/max_width);
  }
  var j = 0;
  for (var i=0; i<num_elem; i += step) {
    slice_time[j] = time[i];
    slice_raw_data[j] = raw_data[i];
    slice_real_filtered_data[j++] = real_filtered_data[i];
  }

  console.log(max_width)

  var trace_raw_data = {
    name: "raw",
    x: slice_time,
    y: slice_raw_data,
    type: 'scatter',
    opacity: 0.8,
  };
  // var layout_data = {
  //   title: "Data",
  //   xaxis: {
  //     title: "time [ms]",
  //     type: 'lineair',
  //     range: [0,end_time],
  //     // type: "date",
  //     // rangeslider: {},
  //   },
  //   showLegend: false,
  //   margin: {
  //     l: 30,
  //     r: 30,
  //   },
  // };
  var trace_filtered_data = {
    name: "filtered",
    x: slice_time,
    y: slice_real_filtered_data,
    type: 'scatter',
    opacity: 0.8,
  }
  var trace_raw_histogram = {
    name: "raw",
    x: raw_data.slice(0,num_elem),
    type: 'histogram',
    opacity: 0.8,
  }
  var trace_filtered_histogram = {
    name: "filtered",
    x: real_filtered_data.slice(0,num_elem),
    type: 'histogram',
    opacity: 0.8,
  }

  Plotly.newPlot('rawdata', [trace_raw_data, trace_filtered_data], layout_data, {staticPlot: true});

  Plotly.newPlot('histogram', [trace_raw_histogram, trace_filtered_histogram], layout_histogram, {staticPlot: true});

  console.log("width", layout_data.width)

}
