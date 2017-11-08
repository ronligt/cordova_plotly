// Default values
var max_num_elem = 2048;
var default_num_elem = 512;

// Initialise data
max_x = Array.apply(null, Array(max_num_elem)).map(function(_,i){return i;});
max_y = Array.apply(null, Array(max_num_elem)).map(function(_,i){return Math.floor(Math.random()*255);});
var layout1 = {
  title: "plaatje 1",
  type: 'scatter',
  showLegend: false,
};

graph(default_num_elem);

$('#slide1').change(function(){
  graph($(this).val());
});

function graph(num_elem) {
  var trace1 = {
    x: max_x.slice(0,num_elem),
    y: max_y.slice(0,num_elem),
  };

  Plotly.newPlot('myDiv1', [trace1], layout1, {staticPlot: true});
}
