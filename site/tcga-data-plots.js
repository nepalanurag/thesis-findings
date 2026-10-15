/* Descriptive figures for the TCGA-BRCA data itself: PCA, control-gene
   violins, top-genes heatmap. Data from site/data/tcga-data.json, built by
   ~/workspace/thesis-validation/plot_data_tcga.R (n = 1043 primary tumors). */

(function () {
  "use strict";

  var ERP = "#1B9E77", ERN = "#d95f02", DARK = "#333333";
  var TITLE_FONT = { size: 16, family: "Georgia, serif" };
  var BG = { plot_bgcolor: "#faf8f3", paper_bgcolor: "#ffffff" };

  fetch("data/tcga-data.json").then(function (r) { return r.json(); }).then(function (D) {

    // ---------- PCA scatter colored by ER status ----------
    var pca = D.pca;
    var traces = ["ER+", "ER-"].map(function (grp) {
      var pts = pca.filter(function (d) { return d.ER === grp; });
      return {
        x: pts.map(function (d) { return d.PC1; }),
        y: pts.map(function (d) { return d.PC2; }),
        mode: "markers",
        marker: { color: grp === "ER+" ? ERP : ERN, size: 5, opacity: 0.55,
                  line: { width: 0 } },
        name: grp + " (n=" + pts.length + ")",
        hoverinfo: "name"
      };
    });

    Plotly.newPlot("pca-plot", traces, {
      title: { text: "The ER program is visible before any model: PCA of 1,500 genes",
               font: TITLE_FONT },
      xaxis: { title: "PC1 (" + D.pca_var[0] + "% variance)", zeroline: false },
      yaxis: { title: "PC2 (" + D.pca_var[1] + "% variance)", zeroline: false },
      margin: { l: 60, r: 20, t: 60, b: 55 },
      legend: { x: 1.02, y: 0.9 },
      plot_bgcolor: BG.plot_bgcolor, paper_bgcolor: BG.paper_bgcolor,
      hovermode: "closest"
    }, { responsive: true, displayModeBar: false });

    // ---------- Control-gene violins by ER status ----------
    var genes = ["ESR1", "PGR", "FOXA1", "GATA3"];
    var vtraces = [];
    genes.forEach(function (g) {
      ["ER+", "ER-"].forEach(function (grp) {
        var vals = D.violins
          .filter(function (d) { return d.gene === g && d.ER === grp; })
          .map(function (d) { return d.expr; });
        vtraces.push({
          type: "violin",
          y: vals,
          x0: g,
          name: grp,
          legendgroup: grp,
          showlegend: g === "ESR1",
          scalegroup: g,
          side: grp === "ER+" ? "positive" : "negative",
          line: { color: DARK, width: 1 },
          fillcolor: grp === "ER+" ? ERP : ERN,
          opacity: 0.75,
          box: { visible: true },
          meanline: { visible: true },
          hoverinfo: "y+name"
        });
      });
    });

    Plotly.newPlot("violin-plot", vtraces, {
      title: { text: "The positive controls separate ER+ from ER- on their own",
               font: TITLE_FONT },
      yaxis: { title: "log2(RSEM + 1) expression" },
      xaxis: { title: "" },
      violinmode: "overlay",
      margin: { l: 60, r: 20, t: 60, b: 50 },
      legend: { x: 1.02, y: 0.9 },
      plot_bgcolor: BG.plot_bgcolor, paper_bgcolor: BG.paper_bgcolor,
      hovermode: "closest"
    }, { responsive: true, displayModeBar: false });

    // ---------- Heatmap: top 50 variable genes, samples ordered ER+/ER- ----------
    var H = D.heatmap;
    var erColors = H.er.map(function (e) { return e === "ER+" ? ERP : ERN; });

    Plotly.newPlot("heatmap-plot", [
      {
        z: H.er.map(function (e) { return [e === "ER+" ? 1 : 0]; }),
        x: ["ER"],
        type: "heatmap",
        colorscale: [[0, ERN], [1, ERP]],
        showscale: false,
        hoverinfo: "skip",
        xaxis: "x2", yaxis: "y"
      },
      {
        z: H.z,
        x: H.genes,
        type: "heatmap",
        colorscale: "RdBu",
        reversescale: true,
        zmin: -3, zmax: 3,
        showscale: true,
        colorbar: { title: "z-score", titleside: "right" },
        hoverinfo: "x",
        xaxis: "x", yaxis: "y"
      }
    ], {
      title: { text: "Top 50 variable genes: the ER program as a visible block",
               font: TITLE_FONT },
      xaxis: { title: "", tickangle: -45, tickfont: { size: 9 },
               domain: [0.07, 1] },
      xaxis2: { title: "", domain: [0, 0.05], showticklabels: false },
      yaxis: { title: "Samples (ordered: ER+, then ER-)", showticklabels: false },
      margin: { l: 60, r: 60, t: 60, b: 110 },
      height: 520,
      plot_bgcolor: BG.plot_bgcolor, paper_bgcolor: BG.paper_bgcolor
    }, { responsive: true, displayModeBar: false });

  }).catch(function (err) {
    console.error("tcga-data plots failed:", err);
  });
})();
