// Default values
var max_num_elem = 72000;
var sample_rate = 44; // samples/ms
var default_num_elem = 512;

// Initialise data
var max_x = Array(max_num_elem);
var max_y = Array(max_num_elem);

// max_x = Array.apply(null, Array(max_num_elem)).map(function(_,i){return i;});
// max_y = Array.apply(null, Array(max_num_elem)).map(function(_,i){return Math.floor(Math.random()*255);});
// max_y = Array.apply(null, Array(max_num_elem)).map(function(_,i){return randomGaussian(0.0, 1.0);});

for (var i=0; i<max_num_elem; i++) {
  max_x[i] = i/sample_rate; // ms
  max_y[i] = randomGaussian(0.0, 1.0);
}
var layout_rawdata = {
  title: "Raw data",
  xaxis: {
    title: "time [ms]",
  },
  showLegend: false,
};
var layout_histogram = {
  title: "Histogram",
  showLegend: false,
};

graph(default_num_elem);

$('#slide1').change(function(){
  graph($(this).val());
});

function graph(num_elem) {
  var trace_rawdata = {
    x: max_x.slice(0,num_elem),
    y: max_y.slice(0,num_elem),
    type: 'scatter',
  };
  var trace_histogram = {
    x: max_y.slice(0,num_elem),
    type: 'histogram',
  }

  Plotly.newPlot('rawdata', [trace_rawdata], layout_rawdata, {staticPlot: true});

  Plotly.newPlot('histogram', [trace_histogram], layout_histogram, {staticPlot: true});

}
