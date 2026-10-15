/* Additional TCGA-BRCA figures: selection frequency, enrichment, per-split AUCs.
   Data from site/data/tcga-extra.json, built by
   ~/workspace/thesis-validation/plot_extra_tcga.R (real computed results). */

(function () {
  "use strict";

  var EF = "#1B9E77", GRAY = "#9a9a9a", DARK = "#333333";
  var TITLE_FONT = { size: 16, family: "Georgia, serif" };
  var BASE_LAYOUT = {
    margin: { l: 70, r: 20, t: 60, b: 60 },
    plot_bgcolor: "#faf8f3", paper_bgcolor: "#ffffff",
    hovermode: "closest"
  };

  function layout(extra) {
    var l = {};
    for (var k in BASE_LAYOUT) l[k] = BASE_LAYOUT[k];
    for (var k2 in extra) l[k2] = extra[k2];
    return l;
  }

  fetch("data/tcga-extra.json").then(function (r) { return r.json(); }).then(function (D) {

    // ---------- Figure A: what the adopted model selects ----------
    var sf = D.selection_frequency;
    var genes = sf.map(function (d) { return d.gene; });
    var splits = sf.map(function (d) { return d.splits; });
    var colors = sf.map(function (d) { return d.control ? EF : GRAY; });

    Plotly.newPlot("selfreq-plot", [{
      y: genes,
      x: splits,
      type: "bar",
      orientation: "h",
      marker: { color: colors, line: { color: DARK, width: 0.8 } },
      text: sf.map(function (d) {
        return d.gene + ": selected in " + d.splits + " of 5 splits" +
          (d.control ? " (positive control)" : "");
      }),
      hoverinfo: "text",
      showlegend: false
    }], layout({
      title: { text: "EF-Adap S3: selection frequency across 5 fresh splits (top genes)",
               font: TITLE_FONT },
      xaxis: { title: "Splits selecting the gene (of 5)", range: [0, 5.4], dtick: 1 },
      yaxis: { title: "", automargin: true },
      margin: { l: 90, r: 20, t: 60, b: 50 },
      height: 720,
      annotations: [{
        x: 5, y: genes[genes.length - 1], text: "green = ER positive controls",
        showarrow: false, xanchor: "right", yanchor: "bottom",
        font: { size: 12, color: EF, style: "italic" }
      }]
    }), { responsive: true, displayModeBar: false });

    // ---------- Figure B: signature enrichment across arms ----------
    var en = D.enrichment;
    var armNames = { "Baseline LASSO": "LASSO", "EF-Std": "EF-Std", "EF-Adap": "EF-Adap",
                     "LF-Std": "LF-Std", "LF-Adap": "LF-Adap" };
    var xarms = en.map(function (d) { return armNames[d.arm] || d.arm; });

    function enrichTrace(key, name, color) {
      return {
        x: xarms,
        y: en.map(function (d) { return d[key]; }),
        type: "bar",
        name: name,
        marker: { color: color, line: { color: DARK, width: 1 } },
        text: en.map(function (d, i) {
          var sig = key === "pam50_logp" ? d.pam50_sig : d.onc_sig;
          var ov = key === "pam50_logp" ? d.pam50_overlap : d.onc_overlap;
          return name + ": median overlap " + ov + " genes, significant in " +
                 sig + "/5 splits";
        }),
        hoverinfo: "text+y",
        textposition: "none"
      };
    }

    Plotly.newPlot("enrich-plot",
      [enrichTrace("pam50_logp", "PAM50", "#1B9E77"),
       enrichTrace("onc_logp", "Oncotype DX", "#377EB8")],
      layout({
        title: { text: "Biology check: signature enrichment by arm (-log10 median p)",
                 font: TITLE_FONT },
        yaxis: { title: "-log10(p), median over 5 splits" },
        xaxis: { title: "Method arm" },
        barmode: "group",
        shapes: [{
          type: "line", x0: -0.5, x1: 4.5,
          y0: 1.301, y1: 1.301,
          line: { dash: "dash", color: "#a00", width: 1.5 }
        }],
        annotations: [{
          x: 4.5, y: 1.301, text: "p = 0.05",
          showarrow: false, xanchor: "right", yanchor: "bottom",
          font: { size: 12, color: "#a00" }
        }],
        legend: { x: 1.02, y: 0.9 }
      }), { responsive: true, displayModeBar: false });

    // ---------- Figure C: per-split AUCs (honest spread) ----------
    var strip = D.auc_strip;
    var armOrder = ["LASSO", "EF-Std", "EF-Adap", "LF-Std", "LF-Adap"];
    var armColor = { "LASSO": GRAY, "EF-Std": "#8fd0b8", "EF-Adap": EF,
                     "LF-Std": "#a8c6e2", "LF-Adap": "#377EB8" };
    // jitter x positions so overlapping split dots separate
    var seen = {};
    var traces = armOrder.map(function (arm) {
      var pts = strip.filter(function (d) {
        var a = d.arm.replace("Baseline ", "");
        return a === arm;
      });
      var xj = pts.map(function () {
        var k = arm;
        seen[k] = (seen[k] || 0) + 1;
        return (seen[k] - 3) * 0.07;
      });
      return {
        x: pts.map(function (d, i) { return armOrder.indexOf(arm) + xj[i]; }),
        y: pts.map(function (d) { return d.auc; }),
        mode: "markers",
        marker: { color: armColor[arm], size: 10,
                  line: { color: DARK, width: 1 } },
        text: pts.map(function (d) {
          return arm + " split " + d.split + ": AUC " + d.auc.toFixed(3);
        }),
        hoverinfo: "text",
        name: arm,
        showlegend: false
      };
    });

    Plotly.newPlot("aucstrip-plot", traces, layout({
      title: { text: "Every split, not just the mean: held-out AUC per split",
               font: TITLE_FONT },
      yaxis: { title: "Held-out AUC", range: [0.86, 0.98] },
      xaxis: { title: "Method arm", tickvals: [0, 1, 2, 3, 4], ticktext: armOrder },
      margin: { l: 60, r: 20, t: 60, b: 50 }
    }), { responsive: true, displayModeBar: false });

  }).catch(function (err) {
    console.error("tcga-extra plots failed:", err);
  });
})();
