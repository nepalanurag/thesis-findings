/* F1 figure for the simulation study page.
   Data from the rendered adaptive_permutation_group_lasso notebook
   (mean and SD over 5 runs, n = 100, 30 true causal features). */

(function () {
  "use strict";

  var rows = [
    { m: "Early fusion, standard", f1: 0.960, sd: 0.089, c: "#1B9E77" },
    { m: "Late fusion, standard", f1: 0.774, sd: 0.185, c: "#377EB8" },
    { m: "Early fusion, adaptive", f1: 0.398, sd: 0.213, c: "#8fd0b8" },
    { m: "Late fusion, adaptive", f1: 0.345, sd: 0.227, c: "#a8c6e2" },
    { m: "Baseline LASSO", f1: 0.130, sd: 0.104, c: "#9a9a9a" }
  ];

  var trace = {
    x: rows.map(function (d) { return d.f1; }),
    y: rows.map(function (d) { return d.m; }),
    error_x: { type: "data", array: rows.map(function (d) { return d.sd; }),
               visible: true, color: "#555", thickness: 1.4, width: 8 },
    type: "bar",
    orientation: "h",
    marker: { color: rows.map(function (d) { return d.c; }),
              line: { color: "#333", width: 1 } },
    text: rows.map(function (d) {
      return "F1 " + d.f1.toFixed(2) + " (+/- " + d.sd.toFixed(2) + ")";
    }),
    hoverinfo: "text+y",
    showlegend: false
  };

  Plotly.newPlot("sim-f1-plot", [trace], {
    title: { text: "F1 against known truth by method (mean +/- SD, 5 runs)",
             font: { size: 16, family: "Georgia, serif" } },
    xaxis: { title: "F1 score", range: [0, 1.1] },
    margin: { l: 200, r: 20, t: 60, b: 50 },
    plot_bgcolor: "#faf8f3", paper_bgcolor: "#ffffff"
  }, { responsive: true, displayModeBar: false });
})();
