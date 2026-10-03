/* Interactive figures for the TCGA-BRCA results page.
   Data reproduced from plot_stability_bubble.R and five_arm_summary_table.csv
   (real computed results, Sep 2026). */

(function () {
  "use strict";

  var EF = "#1B9E77", LF = "#377EB8";

  // ---------- Figure 1: stability bubble plot ----------
  var arms = [
    { short: "EF current", family: "EF-Adap", jaccard: 0.229, auc: 0.930, sd: 0.017, nsel: 236.0,
      sig_pam50: "5/5", sig_onco: "4/5", controls: "ESR1 5/5, PGR 4/5, FOXA1 3/5, GATA3 3/5" },
    { short: "EF S1", family: "EF-Adap", jaccard: 0.482, auc: 0.940, sd: 0.008, nsel: 206.6,
      sig_pam50: "5/5", sig_onco: "5/5", controls: "ESR1, PGR, FOXA1, GATA3 all 5/5" },
    { short: "EF S2", family: "EF-Adap", jaccard: 0.483, auc: 0.940, sd: 0.008, nsel: 206.2,
      sig_pam50: "5/5", sig_onco: "5/5", controls: "ESR1, PGR, FOXA1, GATA3 all 5/5" },
    { short: "EF S3", family: "EF-Adap", jaccard: 0.680, auc: 0.948, sd: 0.006, nsel: 122.8,
      sig_pam50: "5/5", sig_onco: "5/5", controls: "ESR1, PGR, FOXA1, GATA3 all 5/5" },
    { short: "LF current", family: "LF-Adap", jaccard: 0.327, auc: 0.924, sd: 0.015, nsel: 498.0,
      sig_pam50: "5/5", sig_onco: "2/5", controls: "ESR1 5/5, PGR 1/5, FOXA1 5/5, GATA3 5/5" },
    { short: "LF S1", family: "LF-Adap", jaccard: 0.535, auc: 0.946, sd: 0.012, nsel: 437.2,
      sig_pam50: "5/5", sig_onco: "2/5", controls: "ESR1 5/5, PGR 5/5, FOXA1 4/5, GATA3 5/5" },
    { short: "LF S2", family: "LF-Adap", jaccard: 0.557, auc: 0.945, sd: 0.013, nsel: 427.8,
      sig_pam50: "5/5", sig_onco: "2/5", controls: "ESR1 5/5, PGR 5/5, FOXA1 4/5, GATA3 5/5" },
    { short: "LF S3", family: "LF-Adap", jaccard: 0.588, auc: 0.944, sd: 0.016, nsel: 206.8,
      sig_pam50: "5/5", sig_onco: "1/5", controls: "ESR1 2/5, PGR 1/5, FOXA1 2/5, GATA3 2/5" }
  ];

  function hoverText(a) {
    return "<b>" + a.short + "</b><br>" +
      "Jaccard: " + a.jaccard.toFixed(3) + "<br>" +
      "AUC: " + a.auc.toFixed(3) + " (+/- " + a.sd.toFixed(3) + ")<br>" +
      "Features selected: " + a.nsel.toFixed(1) + "<br>" +
      "PAM50 enriched: " + a.sig_pam50 + "<br>" +
      "Oncotype DX enriched: " + a.sig_onco + "<br>" +
      "Controls: " + a.controls;
  }

  var bubbleTraces = [];
  ["EF-Adap", "LF-Adap"].forEach(function (fam) {
    var pts = arms.filter(function (a) { return a.family === fam; });
    var col = fam === "EF-Adap" ? EF : LF;
    // dashed progression path current -> S1 -> S2 -> S3
    bubbleTraces.push({
      x: pts.map(function (a) { return a.jaccard; }),
      y: pts.map(function (a) { return a.auc; }),
      mode: "lines",
      line: { dash: "dash", color: col, width: 1.5 },
      opacity: 0.55,
      hoverinfo: "skip",
      showlegend: false
    });
    // bubbles with error bars
    bubbleTraces.push({
      x: pts.map(function (a) { return a.jaccard; }),
      y: pts.map(function (a) { return a.auc; }),
      error_y: {
        type: "data",
        array: pts.map(function (a) { return a.sd; }),
        visible: true, color: "#8a8a8a", thickness: 1.2, width: 6
      },
      mode: "markers",
      marker: {
        size: pts.map(function (a) { return a.nsel; }),
        sizemode: "area",
        sizeref: 2 * 498 / (64 * 64),
        color: pts.map(function (a, i) { return i === 0 ? col + "55" : col; }),
        line: { color: "#333", width: 1 }
      },
      text: pts.map(hoverText),
      hoverinfo: "text",
      name: fam,
      showlegend: true
    });
  });

  var annotations = arms.map(function (a) {
    return {
      x: a.jaccard, y: a.auc, text: a.short,
      showarrow: false,
      yshift: a.short === "EF S3" ? 26 : (a.short.indexOf("S1") >= 0 ? -26 : 26),
      font: { size: 12, color: "#333", family: "Georgia, serif" }
    };
  });
  annotations.push({
    x: 0.680, y: 0.948, text: "adopted",
    showarrow: false, yshift: 44, font: { size: 12, color: EF, style: "italic" }
  });
  // arrowheads at path ends
  ["EF-Adap", "LF-Adap"].forEach(function (fam) {
    var pts = arms.filter(function (a) { return a.family === fam; });
    var p0 = pts[2], p1 = pts[3];
    annotations.push({
      x: p1.jaccard, y: p1.auc,
      ax: p0.jaccard, ay: p0.auc,
      xref: "x", yref: "y", axref: "x", ayref: "y",
      showarrow: true, arrowhead: 2, arrowsize: 1.2,
      arrowcolor: fam === "EF-Adap" ? EF : LF, opacity: 0.55,
      arrowwidth: 1.5, text: ""
    });
  });

  Plotly.newPlot("bubble-plot", bubbleTraces, {
    title: { text: "Prediction vs reproducibility: stabilizing the adaptive fusion models",
             font: { size: 16, family: "Georgia, serif" } },
    xaxis: { title: "Selection stability (mean pairwise Jaccard across 5 splits)",
             range: [0.15, 0.76] },
    yaxis: { title: "Predictive power (held-out AUC, mean +/- SD)",
             range: [0.908, 0.963] },
    annotations: annotations,
    margin: { l: 70, r: 30, t: 60, b: 60 },
    hovermode: "closest",
    plot_bgcolor: "#faf8f3", paper_bgcolor: "#ffffff",
    legend: { x: 1.02, y: 0.5 }
  }, { responsive: true, displayModeBar: false });

  // ---------- Figure 2: 5-arm head-to-head ----------
  var fiveArm = [
    { arm: "LASSO", nsel: "30 (12)", auc: 0.931, sd: 0.017, jac: 0.201,
      pam50: "2/5", onco: "4/5", ctl: "ESR1 5/5; PGR 0/5; FOXA1 0/5; GATA3 2/5" },
    { arm: "EF-Std", nsel: "23 (51)", auc: 0.899, sd: 0.000, jac: 0.600,
      pam50: "0/5", onco: "0/5", ctl: "all four missed" },
    { arm: "EF-Adap", nsel: "236 (67)", auc: 0.930, sd: 0.017, jac: 0.229,
      pam50: "5/5", onco: "4/5", ctl: "ESR1 5/5; PGR 4/5; FOXA1 3/5; GATA3 3/5" },
    { arm: "LF-Std", nsel: "107 (47)", auc: 0.917, sd: 0.018, jac: 0.363,
      pam50: "0/5", onco: "0/5", ctl: "ESR1 1/5; PGR 0/5; FOXA1 0/5; GATA3 1/5" },
    { arm: "LF-Adap", nsel: "498 (148)", auc: 0.924, sd: 0.015, jac: 0.327,
      pam50: "5/5", onco: "2/5", ctl: "ESR1 5/5; PGR 1/5; FOXA1 5/5; GATA3 5/5" }
  ];
  var barColors = { "LASSO": "#9a9a9a", "EF-Std": "#8fd0b8", "EF-Adap": EF,
                    "LF-Std": "#a8c6e2", "LF-Adap": LF };

  var aucTrace = {
    x: fiveArm.map(function (d) { return d.arm; }),
    y: fiveArm.map(function (d) { return d.auc; }),
    error_y: { type: "data", array: fiveArm.map(function (d) { return d.sd; }),
               visible: true, color: "#555", thickness: 1.4, width: 8 },
    type: "bar",
    marker: { color: fiveArm.map(function (d) { return barColors[d.arm]; }),
              line: { color: "#333", width: 1 } },
    text: fiveArm.map(function (d) {
      return "AUC " + d.auc.toFixed(3) + " (+/- " + d.sd.toFixed(3) + ")<br>" +
        "Jaccard " + d.jac.toFixed(3) + "<br>PAM50 " + d.pam50 +
        "<br>Oncotype " + d.onco + "<br>" + d.ctl;
    }),
    hoverinfo: "text+x",
    showlegend: false
  };

  Plotly.newPlot("fivearm-plot", [aucTrace], {
    title: { text: "Head-to-head: held-out AUC by method (ER-IHC endpoint, 5 splits)",
             font: { size: 16, family: "Georgia, serif" } },
    yaxis: { title: "Held-out AUC (mean +/- SD)", range: [0.86, 0.97] },
    xaxis: { title: "Method arm" },
    margin: { l: 60, r: 20, t: 60, b: 50 },
    plot_bgcolor: "#faf8f3", paper_bgcolor: "#ffffff"
  }, { responsive: true, displayModeBar: false });

  var jacTrace = {
    x: fiveArm.map(function (d) { return d.arm; }),
    y: fiveArm.map(function (d) { return d.jac; }),
    type: "bar",
    marker: { color: fiveArm.map(function (d) { return barColors[d.arm]; }),
              line: { color: "#333", width: 1 } },
    text: fiveArm.map(function (d) { return "Jaccard " + d.jac.toFixed(3); }),
    hoverinfo: "text+x",
    showlegend: false
  };

  Plotly.newPlot("jaccard-plot", [jacTrace], {
    title: { text: "Same arms: selection stability (Jaccard)",
             font: { size: 16, family: "Georgia, serif" } },
    yaxis: { title: "Mean pairwise Jaccard", range: [0, 0.7] },
    xaxis: { title: "Method arm" },
    margin: { l: 60, r: 20, t: 60, b: 50 },
    plot_bgcolor: "#faf8f3", paper_bgcolor: "#ffffff"
  }, { responsive: true, displayModeBar: false });
})();
