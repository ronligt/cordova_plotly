// Default values
var max_time =         1.5; // maximum recording time [seconds]
var default_time =     0.1 * max_time; // default time
var step_time =        0.001 * max_time; // slider steps
var sample_rate =      44100; // samples/second
var max_num_elem =     Math.floor(max_time * sample_rate); // maximum number of samples
var default_num_elem = Math.floor(default_time * sample_rate);
var step_sample =      Math.floor(0.001 * max_num_elem);

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

var domRect = document.querySelector('#data').getBoundingClientRect();
var canvas_width = parseInt(domRect.width);
console.log(canvas_width);

// Initialise data
var timeBuffer = new ArrayBuffer(max_num_elem*Float32Array.BYTES_PER_ELEMENT);
var timeView = new Float32Array(timeBuffer);
var tau = Array(max_num_elem);
var tauBuffer = new ArrayBuffer(max_num_elem*Float32Array.BYTES_PER_ELEMENT)
var tauView = new Float32Array(tauBuffer)
var rawBuffer = new ArrayBuffer(max_num_elem*Float32Array.BYTES_PER_ELEMENT);
var rawView = new Float32Array(rawBuffer);
// var complex_data = Array(max_num_elem*2);
var complexBuffer = new ArrayBuffer(2*max_num_elem*Float32Array.BYTES_PER_ELEMENT)
var complexView = new Float32Array(complexBuffer)

// Data generation
for (var i=startSample; i<stopSample; i++) {
  timeView[i] = i/sample_rate; // ms
  tauView[i] = (i-middleSample)/sample_rate; // ms
  rawView[i] = randomGaussian(0.0, 1.0);
  complexView[i*2] = rawView[i]
  complexView[i*2+1] = 0
}

// Calculate filter in time domain
var complexFilter = recurrenceFilter(coefficients, complexView);  // t-domain filtering
var complexFilterView = Float32Array.from(complexFilter);
var realFilterBuffer = new ArrayBuffer(max_num_elem*Float32Array.BYTES_PER_ELEMENT)
var realFilterView = new Float32Array(realFilterBuffer)
realFilterView = complexFilterView.filter(function (elem, index){return (index % 2) == 0})

// // Calculate filter in frequency domain
// var rSpect = fft(complex_data, nn, ndim, FORWARD);  // signal spectrum
// // var Spect,fpsdTemp,fpsd;
// var rpsdTemp = abssq(rSpect);  // power spectral density
// // var rpsd = uRs(rpsdTemp);
// var autoCorrData = uRs(fft(rpsdTemp, nn, ndim, BACKWARD));
// console.log(rSpect.length,rpsdTemp.length,autoCorrData.length);

var canvas_data = document.getElementById("data");
var canvas_histogram = document.getElementById("histogram");

graph(default_num_elem);
activateSlider();

function graph(num_elem) {
  var xdata, ydata;
  var xfilter, yfilter;
  if (num_elem > canvas_width) {
    var modulo = Math.ceil(num_elem / (canvas_width*(1+4*num_elem/max_num_elem)))
    xdata = timeView.subarray(0,num_elem).filter(filterModulo)
    var minY, maxY
    ydata = rawView.subarray(0,num_elem).filter(filterModulo)
    yfilter = realFilterView.subarray(0,num_elem).filter(function (elem, index) {return (index % modulo) == 0 })
  } else {
    xdata = timeView.subarray(0, num_elem)
    ydata = rawView.subarray(0, num_elem)
    yfilter = realFilterView.subarray(0, num_elem)
    var minY=ydata[0], maxY=ydata[0]
    for (var i=0; i<num_elem; i++) {
      minY = (ydata[i] < minY?ydata[i]:minY)
      maxY = (ydata[i] > maxY?ydata[i]:maxY)
    }
  }

  var handle_data = fig.figure(canvas_data, "time [s]", "", 0, xdata[xdata.length-1], minY, maxY);
  fig.plot(handle_data, xdata, ydata, 'line', 'blue', []);
  fig.plot(handle_data, xdata, yfilter, 'line', 'red', [])
  var histogram = fig.hist(canvas_histogram, ydata,100*(1+4*num_elem/max_num_elem));

  function filterModulo(elem, index) {
    if ((index % modulo) == 0) {
      minY = (minY === undefined?elem:(elem < minY?elem:minY))
      maxY = (maxY === undefined?elem:(elem > maxY?elem:maxY))
      return true
    }
    return false
  }
}

var view_data, view_histogram;
var options_data, options_histogram;
var chart_data, chart_histogram;
// google.charts.load('current', {'packages':['corechart']});
// google.charts.load('current', {'packages':['line']});
// google.charts.setOnLoadCallback(prepareChart);

function prepareChart() {
  // setup Data
  // var data = new google.visualization.DataTable();
  // data.addColumn('number', 'time');
  // data.addColumn('number', 'raw');
  // data.addColumn('number', 'filtered');
  // for (var i=startSample; i<stopSample; i++) {
  //   data.addRows([[time[i], raw_data[i], real_filtered_data[i]]]);
  // };

  // setup Views
  // view_data = new google.visualization.DataView(data);
  // view_histogram = new google.visualization.DataView(data);
  // view_histogram.setColumns([1,2]);

  // setup options
  options_data = {
    title: 'Data',
    enableInteractivity: false,
    legend: {
      position: 'in',
    },
    hAxis: {
      title: 'time [s]'
    },
  };
  options_histogram = {
    title: 'Histogram',
    enableInteractivity: false,
    legend: {
      position: 'in',
    },
    bar: {
      gap: 0,
    },
    dataOpacity: 0.5,
    histogram: {
      maxNumBuckets: 200,
    },
  };

  // setup charts
  // chart_data = new google.visualization.LineChart(document.getElementById('rawdata'));
  // chart_histogram = new google.visualization.Histogram(document.getElementById('histogram'));

  // drawChart(default_num_elem);
  // activateSlider();
}
function drawChart(num_elem) {
  // view_data.setRows(0, num_elem-1);

  if (num_elem > div_width) {
    var filter_num_elem = div_width * (1 + Math.floor(num_elem / max_num_elem * 4));
  } else {
    var filter_num_elem = num_elem;
  }
  var filter = Array.from(new Array(filter_num_elem),(val,index)=>Math.floor(index/filter_num_elem*num_elem));
  view_data.setRows(filter);
  // chart_data.draw(view_data, options_data);

  // var chart = new google.charts.Line(document.getElementById('rawdata'));
  // chart.draw(view, google.charts.Line.convertOptions(options));
  options_histogram = {
    title: 'Histogram',
    enableInteractivity: false,
    legend: {
      position: 'in',
    },
    bar: {
      gap: 0,
    },
    histogram: {
      maxNumBuckets: filter_num_elem/4,
    },
    hAxis: {
      ticks: Array.from(new Array(10),(val,index)=>index*(5+5)/10-5),
    }
  };

  // view_histogram.setRows(filter);
  // chart_histogram.draw(view_histogram, options_histogram);
}
function activateSlider() {
  // Redraw plot after slider change
  $('#slide1').change(function(){
    // drawChart(parseInt($(this).val()));
    graph(parseInt($(this).val()));
  });
}
