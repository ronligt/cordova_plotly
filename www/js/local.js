// Default values
var max_num_elem = 72000;
var sample_rate = 44; // samples/ms
var default_num_elem = 512;

// SSP stuff
var leftCoeffs = [1,0];
var rightCoeffs = youngBPFFilter;
var coefficients = [ leftCoeffs, rightCoeffs ];

// Initialise data
var time = Array(max_num_elem);
var raw_data = Array(max_num_elem);
var complex_data = Array(max_num_elem*2);

// max_x = Array.apply(null, Array(max_num_elem)).map(function(_,i){return i;});
// max_y = Array.apply(null, Array(max_num_elem)).map(function(_,i){return Math.floor(Math.random()*255);});
// max_y = Array.apply(null, Array(max_num_elem)).map(function(_,i){return randomGaussian(0.0, 1.0);});

for (var i=0; i<max_num_elem; i++) {
  time[i] = i/sample_rate; // ms
  raw_data[i] = randomGaussian(0.0, 1.0);
  complex_data[i*2] = raw_data[i];
  complex_data[i*2+1] = 0;
}
var layout_data = {
  title: "Data",
  xaxis: {
    title: "time [ms]",
  },
  showLegend: false,
  margin: {
    l: 10,
    r: 10,
  },
};
var layout_histogram = {
  title: "Histogram",
  showLegend: false,
  margin: {
    l: 30,
    r: 10,
  },
};

graph(default_num_elem);

$('#slide1').change(function(){
  graph($(this).val());
});

function graph(num_elem) {
  var slice_time = time.slice(0,num_elem);
  var slice_raw_data = raw_data.slice(0,num_elem);
  var slice_complex_data = complex_data.slice(0,num_elem*2)

  var complex_filtered_data = recurrenceFilter(coefficients, slice_complex_data);  // t-domain filtering
  var real_filtered_data = uRs(complex_filtered_data);  // data is real

  console.log(raw_data.length, " ", complex_filtered_data.length, " ", real_filtered_data.length);

  var trace_raw_data = {
    name: "raw",
    x: slice_time,
    y: slice_raw_data,
    type: 'scatter',
    opacity: 0.8,
  };
  var trace_filtered_data = {
    name: "filtered",
    x: slice_time,
    y: real_filtered_data,
    type: 'scatter',
    opacity: 0.8,
  }
  var trace_raw_histogram = {
    name: "raw",
    x: slice_raw_data,
    type: 'histogram',
    opacity: 0.8,
  }
  var trace_filtered_histogram = {
    name: "filtered",
    x: real_filtered_data,
    type: 'histogram',
    opacity: 0.8,
  }

  Plotly.newPlot('rawdata', [trace_raw_data, trace_filtered_data], layout_data, {staticPlot: true}, {displayModeBar: false});

  Plotly.newPlot('histogram', [trace_raw_histogram, trace_filtered_histogram], layout_histogram, {staticPlot: true});

}
