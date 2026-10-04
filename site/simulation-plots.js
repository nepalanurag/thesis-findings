/* FDR figure for the simulation study page.
   Data quoted from Chapter 4 of the thesis (averaged over
   18 scenarios x 50 runs). Lower FDR is better. */

(function () {
  "use strict";

  var rows = [
    { m: "Early fusion, adaptive", fdr: 0.010, c: "#1B9E77" },
    { m: "Late fusion, adaptive", fdr: 0.387, c: "#377EB8" },
    { m: "Baseline LASSO", fdr: 0.526, c: "#9a9a9a" },
    { m: "Early fusion, standard", fdr: 0.720, c: "#8fd0b8" },
    { m: "Late fusion, standard", fdr: 0.810, c: "#a8c6e2" }
  ];

  var trace = {
    x: rows.map(function (d) { return d.fdr; }),
    y: rows.map(function (d) { return d.m; }),
    type: "bar",
    orientation: "h",
    marker: { color: rows.map(function (d) { return d.c; }),
              line: { color: "#333", width: 1 } },
    text: rows.map(function (d) { return "FDR " + d.fdr.toFixed(2); }),
    hoverinfo: "text+y",
    showlegend: false
  };

  Plotly.newPlot("sim-fdr-plot", [trace], {
    title: { text: "False discovery rate by method (lower is better)",
             font: { size: 16, family: "Georgia, serif" } },
    xaxis: { title: "FDR (mean over 18 scenarios)", range: [0, 0.9] },
    margin: { l: 200, r: 20, t: 60, b: 50 },
    plot_bgcolor: "#faf8f3", paper_bgcolor: "#ffffff"
  }, { responsive: true, displayModeBar: false });
})();
