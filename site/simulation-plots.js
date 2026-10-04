/* Figures for the simulation study page.
   Data quoted from Chapter 4 of the thesis (Figure 4.2 heatmap and text,
   averaged over 18 scenarios x 50 runs). */

(function () {
  "use strict";

  var methods = ["Early fusion, adaptive", "Late fusion, adaptive", "Baseline LASSO",
                 "Early fusion, standard", "Late fusion, standard"];
  var colors = ["#1B9E77", "#377EB8", "#9a9a9a", "#8fd0b8", "#a8c6e2"];

  // Figure 1: False discovery rate (lower is better)
  var fdr = [0.010, 0.387, 0.526, 0.720, 0.810];
  Plotly.newPlot("sim-fdr-plot", [{
    x: fdr,
    y: methods,
    type: "bar",
    orientation: "h",
    marker: { color: colors, line: { color: "#333", width: 1 } },
    text: fdr.map(function (d) { return "FDR " + d.toFixed(2); }),
    hoverinfo: "text+y",
    showlegend: false
  }], {
    title: { text: "False discovery rate by method (lower is better)",
             font: { size: 16, family: "Georgia, serif" } },
    xaxis: { title: "FDR (mean over 18 scenarios)", range: [0, 0.9] },
    margin: { l: 200, r: 20, t: 60, b: 50 },
    plot_bgcolor: "#faf8f3", paper_bgcolor: "#ffffff"
  }, { responsive: true, displayModeBar: false });

  // Figure 2: Selection stability, mean pairwise Jaccard (higher is better)
  var jac = [0.590, 0.389, 0.058, 0.240, 0.330];
  Plotly.newPlot("sim-jac-plot", [{
    x: jac,
    y: methods,
    type: "bar",
    orientation: "h",
    marker: { color: colors, line: { color: "#333", width: 1 } },
    text: jac.map(function (d) { return "Jaccard " + d.toFixed(2); }),
    hoverinfo: "text+y",
    showlegend: false
  }], {
    title: { text: "Selection stability by method (higher is better)",
             font: { size: 16, family: "Georgia, serif" } },
    xaxis: { title: "Mean pairwise Jaccard across runs", range: [0, 0.7] },
    margin: { l: 200, r: 20, t: 60, b: 50 },
    plot_bgcolor: "#faf8f3", paper_bgcolor: "#ffffff"
  }, { responsive: true, displayModeBar: false });
})();
