function hist_aves(data_file, container){
    $.getJSON( data_file, function( data ) {
        var plot_data = data;
        console.log(plot_data[0]);
        Highcharts.chart(container, {
            title: {
                text: 'Histogram over index averages (all winter days, all location, 1951 - 2005)'
            },
            xAxis: [{
                title: { text: 'Data' },
                alignTicks: false
            }, {
                title: { text: 'Histogram' },
                alignTicks: false,
                opposite: true
            }],

            yAxis: [{
                title: { text: 'Data' }
            }, {
                title: { text: 'Histogram' },
                opposite: true
            }],

            series: [{
                name: 'Histogram',
                type: 'histogram',
                binsNumber: 20,
                xAxis: 1,
                yAxis: 1,
                baseSeries: 's1',
                zIndex: 1,
            }, {
                name: 'Data',
                visible: false,
                type: 'scatter',
                data: plot_data,
                id: 's1',
                zIndex: -1,
                marker: {
                    radius: 1.5
                }
            }]
        });
    });
}

function hist_sums(data_file, container){
    $.getJSON(data_file, function( data ) {
        var plot_data = data;
        Highcharts.chart(container, {
            title: {
                text: 'Histogram over index sums over winter (all location, 1951 -2005)'
            },
            xAxis: [{
                title: { text: 'Data' },
                alignTicks: false
            }, {
                title: { text: 'Histogram' },
                alignTicks: false,
                opposite: true
            }],

            yAxis: [{
                title: { text: 'Data' }
            }, {
                title: { text: 'Histogram' },
                opposite: true
            }],

            series: [{
                name: 'Histogram',
                type: 'histogram',
                binsNumber: 20,
                xAxis: 1,
                yAxis: 1,
                baseSeries: 's1',
                zIndex: 1
            }, {
                name: 'Data',
                visible: false,
                type: 'scatter',
                data: plot_data,
                id: 's1',
                zIndex: -1,
                marker: {
                    radius: 1.5
                }
            }]
        });
    });
}

function index_sums_years(data_file, container) {
    $.getJSON(data_file, function( ts_data ) {
        Highcharts.chart(container, {
            title: {
                text: 'Sum of indices over season (DJF) for each year, averaged over all gridpoints'
            },
            subtitle: {
                text: 'Index = number of Deg C below 5th percentile'
            },
            yAxis: {
                title: {
                    text: 'Average sum of Index'
                }
            },
            legend: {
                layout: 'vertical',
                align: 'right',
                verticalAlign: 'middle'
            },
            plotOptions: {
                series: {
                    label: {
                        connectorAllowed: false
                    },
                    pointStart: 1951
                }
            },
            series: [{
                name: 'Average Sum of Index',
                data: ts_data
            }]/*,
            responsive: {
                rules: [{
                    condition: {
                        maxWidth: 500
                    },
                    chartOptions: {
                        legend: {
                            layout: 'horizontal',
                            align: 'center',
                            verticalAlign: 'bottom'
                        }
                    }
                }]
            }
            */
        });
    });
}

function index_aves_doys(data_file, container){
    $.getJSON( data_file, function( ts_data ) {
        var data = ts_data;
        Highcharts.chart(container, {
            title: {
                text: 'Average over all years (1951 - 2011) for each day in season (Dec 1 = 1, Feb 28 = 90), averaged over all gridpoints'
            },
            subtitle: {
                text: 'Index = number of Deg C below 5th percentile'
            },
            yAxis: {
                title: {
                    text: 'Average sum of Index'
                }
            },
            legend: {
                layout: 'vertical',
                align: 'right',
                verticalAlign: 'middle'
            },
            plotOptions: {
                series: {
                    label: {
                        connectorAllowed: false
                    },
                    pointStart: 1
                }
            },
            series: [{
                name: 'Average Index',
                data: data
            }]/*,
            responsive: {
                rules: [{
                    condition: {
                        maxWidth: 500
                    },
                    chartOptions: {
                        legend: {
                            layout: 'horizontal',
                            align: 'center',
                            verticalAlign: 'bottom'
                        }
                    }
                }]
            }
            */
        });
    });
}

function box_plot_sums_for_locations(data_file, container){
    $.getJSON( data_file, function( data ) {
        var locs = [], raw_data = [], plot_data = [];
        var max, min, median, lq, uq; 
        $.each(data, function(key, val) {
            locs.push(String(key));
            raw_data.push(val.sort(function(a, b){return a-b}));
        });
        $.each(raw_data, function(idx, lst) {
            max = lst.max();
            min = lst.min();
            median = lst.median();
            lq = percentile(lst, 25);
            uq = percentile(lst, 75);
            plot_data.push([min, lq, median, uq, max]);
        }); 
        Highcharts.chart(container, {
            chart: {
                type: 'boxplot'
            },
            title: {
                text: 'Box Plot Sums at 5 locations'
            },
            legend: {
                enabled: false
            },
            xAxis: {
                categories: locs,
                title: {
                    text: 'Locations'
                }
            },
            yAxis: {
                title: {
                    text: 'Index'
                }
            },
            series: [{
                name: 'Index over season',
                data: plot_data,
                tooltip: {
                    headerFormat: '<em>Location {point.key}</em><br/>'
                }
            }]
        });
    });
}

function box_plot_aves_for_locations(data_file, container){
    $.getJSON( data_file, function( data ) {
        var locs = [], raw_data = [], plot_data = [];
        var max, min, median, lq, uq;
        $.each(data, function(key, val) {
            locs.push(String(key));
            raw_data.push(val.sort(function(a, b){return a-b}));
        });
        $.each(raw_data, function(idx, lst) {
            max = lst.max();
            min = lst.min();
            median = lst.median();
            lq = percentile(lst, 25);
            uq = percentile(lst, 75);
            plot_data.push([min, lq, median, uq, max]);
        });
        Highcharts.chart(container, {
            chart: {
                type: 'boxplot'
            },
            title: {
                text: 'Box Plot Averages at 5 locations'
            },
            legend: {
                enabled: false
            },
            xAxis: {
                categories: locs,
                title: {
                    text: 'Locations'
                }
            },
            yAxis: {
                title: {
                    text: 'Index'
                }
            },
            series: [{
                name: 'Locations',
                data: plot_data,
                tooltip: {
                    headerFormat: '<em>Location {point.key}</em><br/>'
                }
            }]
        });
    });
}


function box_plot_num_days(data_file, container){
    $.getJSON( data_file, function( data ) {
        var years = [], raw_data = [], plot_data = [];
        var max, min, median, lq, uq;
        $.each(data, function(key, val) {
            years.push(String(key));
            raw_data.push(val.sort(function(a, b){return a-b}));
        });
        $.each(raw_data, function(idx, lst) {
            max = lst.max();
            min = lst.min();
            median = lst.median();
            lq = percentile(lst, 25);
            uq = percentile(lst, 75);
            plot_data.push([min, lq, median, uq, max]);
        });
        Highcharts.chart(container, {
            chart: {
                type: 'boxplot'
            },
            title: {
                text: 'Number of days below 5th percentile'
            },
            legend: {
                enabled: false
            },
            xAxis: {
                categories: years,
                title: {
                    text: 'Years'
                }
            },
            yAxis: {
                title: {
                    text: 'Number of days'
                }
            },
            series: [{
                name: 'Years',
                data: plot_data,
                tooltip: {
                    headerFormat: '<em>Year {point.key}</em><br/>'
                }
            }]
        });
    });
}

function pca_ts(data_file, container, comp_idx){
    $.getJSON( data_file, function( ts_data ) {
        var data = [], i, d_utc, ymd_str, y, m, d;
        var doy = 0;
        for (i = 0; i < ts_data.length; i++){
            doy+=1;
            ymd_str = ts_data[i][0].split('-');
            y = parseInt(ymd_str[0]);
            m = parseInt(ymd_str[1]) - 1;
            d = parseInt(ymd_str[2]);
            //console.log(String(y) + '-' + String(m) + '-' + String(d))
            d_utc = Date.UTC(y,m,d);
            val = ts_data[i][1] / 1000.0;
            //data.push([d_utc, val]);
            data.push([doy, val])
        }
        Highcharts.chart(container, {
            chart: {
                zoomType: 'x'
            },
            title: {
                //text: 'PCA time series ' + String(comp_idx) + ' component'
                text:''
            },
            subtitle: {
                text: 'Click and drag in the plot area to zoom in'
            },
            xAxis: {
                /*
                type: 'datetime',
                dateTimeLabelFormats: { // don't display the dummy year
                    month: '%e. %b',
                    year: '%b'
                },
                */
                title: {
                    text: 'Time Step'
                },
            },
            
            yAxis: {
                title: {
                    text: 'Standard deviations'
                }
            },
            legend: {
                layout: 'vertical',
                align: 'right',
                verticalAlign: 'middle'
            },
            series: [{
                name: 'PCA_' + String(comp_idx),
                data: data
            }]
        });
    });
}
