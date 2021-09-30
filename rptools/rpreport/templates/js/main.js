/**
 * If a float, return it with 3 decimals, or NaN (not a number).
 */
function numberFormatter(params) {
  if (!(isNaN(parseFloat(params.value)))) {
  return parseFloat(params.value).toFixed(3);
  } else {
  return parseFloat(params.value);
  }
}

/**
 * Give the height of the ag-grid column header
 * @return {number} tallestHeaderTextHeight height in pixels.
 */
function headerHeightGetter() {
  const columnHeaderTexts = [
      ...document.querySelectorAll('.ag-header-cell-text'),
  ];
  const clientHeights = columnHeaderTexts.map(
      headerText => headerText.clientHeight
  );
  const tallestHeaderTextHeight = Math.max(...clientHeights);

  return tallestHeaderTextHeight;
}

/**
 * Set the height of the ag-grid column header
 */
function headerHeightSetter() {
  const padding = 20;
  const height = headerHeightGetter() + padding;
  gridOptions.api.setHeaderHeight(height);
}

/* used in displayCards */
let neverSelectedRow = true;

/**
 * Collapse detailed cards according row(s) selected
 * @param {object} event on row selection.
 */
function displayCards(event) {
  /* Getting selected row data */
  let selectedNodes = gridOptions.api.getSelectedNodes();
  let selectedData = selectedNodes.map(node => node.data);

  /* Selecting cards dom elements */
  let selectedRowOnPrintElt = document.getElementById('selectedRowOnPrint');
  let ecCodeCardElt = document.getElementById('ec_code_panel');
  let ecCodeBodyElt = document.getElementById('ec_code_details');
  let thermodynamicCardElt =
      document.getElementById('thermodynamic_chart_panel');
  let ruleScoreCardElt =
      document.getElementById('rule_score_chart_panel');
  let comparativeChartsCardElt =
      document.getElementById('comparative_chart_panel');

  /* Uncomment below for debug purpose
   * alert(`Selected Nodes:\n${JSON.stringify(selectedData)}`); */

  /* Building table of selected rows visible only for printing */
  rowDataTable = "<table class=\"table table-striped\">"
      + "<thead><tr><th scope=\"col\">Pathway</th>"
      + "<th scope=\"col\">Global score</th>"
      + "<th scope=\"col\">FBA</th>"
      + "<th scope=\"col\">ΔGm-1</th>"
      + "<th scope=\"col\">Rule score</th>"
      + "<th scope=\"col\">Steps</th></tr></thead><tbody>";
  for (let i = 0; i < selectedData.length; i++) {
    rowDataTable += "<tr><th scope=\"row\">"
        + selectedData[i]['pathway_name'] + "</th><td>"
        + parseFloat(selectedData[i]['global_score']).toFixed(3)
        + "</td><td>"
        + parseFloat(selectedData[i]['fba_obj_fraction']).toFixed(3)
        + "</td><td>"
        + parseFloat(selectedData[i]['dfG_prime_m']).toFixed(3)
        + "</td><td>"
        + parseFloat(selectedData[i]['mean_rule_score']).toFixed(3)
        + "</td><td>"
        + selectedData[i]['nb_reactions']
        + "</td></tr>";
  }
  rowDataTable += "</tbody></table>";
  selectedRowOnPrintElt.innerHTML = rowDataTable;

  /* If only one row is selected */
  if(selectedData.length==1) {
    /* * Show detail cards... */
    thermodynamicCardElt.classList.remove("d-none");
    thermodynamicCardElt.classList.add("d-block");
    ruleScoreCardElt.classList.remove("d-none");
    ruleScoreCardElt.classList.add("d-block");
    ecCodeCardElt.classList.remove("d-none");
    ecCodeCardElt.classList.add("d-block");
    /* ... and hide comparative card */
    comparativeChartsCardElt.classList.remove("d-block");
    comparativeChartsCardElt.classList.add("d-none");

    // If it's the first time a row is selected, then open all detail cards
    if(neverSelectedRow) {
      neverSelectedRow = false;
      const thermodynamicCardBodyElt =
          document.getElementById('collapseTwo');
      const thermodynamicCardBodyCollapse =
          new bootstrap.Collapse(thermodynamicCardBodyElt, {
            toggle: false
          });
      thermodynamicCardBodyCollapse.show();
      const ruleScoreCardBodyElt =
          document.getElementById('collapseThree');
      const ruleScoreCardBodyCollapse =
          new bootstrap.Collapse(ruleScoreCardBodyElt, {
            toggle: false
          });
      ruleScoreCardBodyCollapse.show();
      const ecCodeCardBodyElt = document.getElementById('collapseOne');
      const ecCodeCardBodyCollapse =
          new bootstrap.Collapse(ecCodeCardBodyElt, {
            toggle: false
          });
      ecCodeCardBodyCollapse.show();
    }

    // Formatting ec_code display
    let RP_indexes = Object.keys(selectedData[0]['reactions']);
    let formated_ec_code = [];
    let thermodynamic_values = [];
    let rule_score_values = [];

    for (let i = 0; i < RP_indexes.length; i++) {
      //building html formated ec-codes
      formated_ec_code[i] = RP_indexes[i] + ": "
      // loop in each reaction, to build ec_code list
      let totalOfEcCodePerReactions =
          selectedData[0]['reactions'][RP_indexes[i]]['ec_code'].length;
      for (let j = 0; j < totalOfEcCodePerReactions; j++) {
        single_ec_code =
            JSON.stringify(selectedData[0]['reactions']
            [RP_indexes[i]]['ec_code'][j]);
        single_ec_code = single_ec_code.replace(/["]+/g, '');
        if(formated_ec_code[i] != RP_indexes[i] + ": ") {
          formated_ec_code[i] = formated_ec_code[i] + ", "
              + "<a href=\"https://www.ebi.ac.uk/intenz/query?cmd=SearchEC&ec="
              + single_ec_code.replace(/.-/g, '') + "\">"
              + single_ec_code + "</a>";
        }
        else
        {
          formated_ec_code[i] = formated_ec_code[i]
              + "<a href=\"https://www.ebi.ac.uk/intenz/query?cmd=SearchEC&ec="
              + single_ec_code.replace(/.-/g, '')
              + "\">"
              + single_ec_code + "</a>";
        }
      }
      formated_ec_code[i] = formated_ec_code[i] + "<br />"

      // building thermodynamic profile data
      thermodynamic_values.push({
          'reaction': RP_indexes[i],
          'dfG_prime_m': JSON.stringify(selectedData[0]['reactions']
              [RP_indexes[i]]['dfG_prime_m'])
      });

      // building rule score data
      rule_score_values.push({
          'reaction': RP_indexes[i],
          'rule_score': JSON.stringify(selectedData[0]['reactions']
              [RP_indexes[i]]['rule_score'])
      });
    }
    ecCodeBodyElt.innerHTML = formated_ec_code.join('');

    // update thermodynamic chart
    chart.series[0].data = thermodynamic_values;

    // update rule score chart
    rule_score_chart.series[0].data = rule_score_values;

  } else if(selectedData.length>1) {
    /* * Hide detail cards... */
    ecCodeCardElt.classList.remove("d-block");
    ecCodeCardElt.classList.add("d-none");
    thermodynamicCardElt.classList.remove("d-block");
    thermodynamicCardElt.classList.add("d-none");
    ruleScoreCardElt.classList.remove("d-block");
    ruleScoreCardElt.classList.add("d-none");
    /* ... and show comparative card */
    comparativeChartsCardElt.classList.remove("d-none");
    comparativeChartsCardElt.classList.add("d-block");

    // Updating data
    globalScoreCompChart.series[0].data = selectedData;
    fbaCompChart.series[0].data = selectedData;
    dfgCompChart.series[0].data = selectedData;
    ruleScoreCompChart.series[0].data = selectedData;
    reactionsNbCompChart.series[0].data = selectedData;
  }
}

const columnDefs = [
  {
    headerName: "Pathway",
    field: "pathway_name",
    filter: 'agTextColumnFilter',
  },
  {
    headerName: "Global score",
    field: "global_score",
    type: 'numericColumn',
    cellClassRules: {
      'text-muted': 'x == null',
    },
    valueFormatter : numberFormatter,
   },
   {
     headerName: "FBA (mmol / gDW / hr)",
     field: "fba_obj_fraction",
     type: 'numericColumn',
     cellClassRules: {
       'text-muted': 'x == null',
     },
     valueFormatter : numberFormatter,
   },
  {
    headerName: "ΔGm-1 (kJ / mol)",
    field: "dfG_prime_m",
    type: 'numericColumn',
    cellClassRules: {
      'text-success': 'x <= 0',
      'text-danger': 'x > 0',
      'text-muted': 'x == null',
    },
    valueFormatter : numberFormatter,
  },
  {
    headerName: "Rule score",
    field: "mean_rule_score",
    type: 'numericColumn',
    cellClassRules: {
      'text-muted': 'x == null',
    },
    valueFormatter : numberFormatter,
   },
  {
    headerName: "Steps",
    field: "nb_reactions",
    type: 'numericColumn',
  }
];

// let the grid know which columns and what data to use
let gridOptions = {
  columnDefs: columnDefs,
  rowData: dataJSON,  /* const dataJSON is defined in data.js file */
  defaultColDef: {
    flex: 1,
    resizable: true,
    sortable: true,
    filter: 'agNumberColumnFilter',
    headerComponentParams: {
      template:
          '<div class="ag-cell-label-container" role="presentation">' +
          '  <span ref="eMenu" class="ag-header-icon ag-header-cell-menu-button"></span>' +
          '  <div ref="eLabel" class="ag-header-cell-label" role="presentation">' +
          '    <span ref="eSortOrder" class="ag-header-icon ag-sort-order"></span>' +
          '    <span ref="eSortAsc" class="ag-header-icon ag-sort-ascending-icon"></span>' +
          '    <span ref="eSortDesc" class="ag-header-icon ag-sort-descending-icon"></span>' +
          '    <span ref="eSortNone" class="ag-header-icon ag-sort-none-icon"></span>' +
          '    <span ref="eText" class="ag-header-cell-text" role="columnheader" style="white-space: normal;"></span>' +
          '    <span ref="eFilter" class="ag-header-icon ag-filter-icon"></span>' +
          '  </div>' +
          '</div>',
    },
  },
  onFirstDataRendered: headerHeightSetter,
  onColumnResized: headerHeightSetter,
  onGridReady: function (params) {
    let defaultSortModel = [
      { colId: 'global_score', sort: 'desc' },
    ];
    params.api.setSortModel(defaultSortModel);
    params.api.sizeColumnsToFit();
  },
  rowSelection: 'multiple',
  onSelectionChanged: displayCards,
};

// setup the grid after the page has finished loading
document.addEventListener('DOMContentLoaded', function() {
  const gridDiv = document.querySelector('#myGrid');
  new agGrid.Grid(gridDiv, gridOptions);

  // create thermodynamic profile chart
  chart = agCharts.AgChart.create({
    container: document.querySelector('#dfG_prime_chart'),
    title: {
      text: 'ΔG\'m of pathway reactions',
      fontWeight: 'normal',
    },
    legend: {
      enabled: false,
    },
    series: [{
      type: 'scatter',
      xKey: 'reaction',
      yKey: 'dfG_prime_m',
      tooltip: {
        renderer: function (params) {
          return {
            content: parseFloat(params.yValue).toFixed(3),
          };
        },
      },
    }],
    axes: [{
      type: 'category',
      position: 'bottom',
      },
      {
      type: 'number',
      position: 'left',
      title: {
        text: 'kJ / mol',
      },
    }],
  });

  // create reactions rule score chart
  rule_score_chart = agCharts.AgChart.create({
    container: document.getElementById('rule_score_chart'),
    title: {
      text: 'Rule scores of pathway reactions',
      fontWeight: 'normal',
    },
    legend: {
      enabled: false,
    },
    series: [{
      type: 'column',
      xKey: 'reaction',
      yKeys: ['rule_score'],
      tooltip: {
        renderer: function (params) {
          return {
            content: parseFloat(params.yValue).toFixed(3),
          };
        },
      },
    }],
    axes: [{
      type: 'category',
      position: 'bottom',
      },
      {
      type: 'number',
      position: 'left',
      }],
  });

  // create global score comparative chart
  globalScoreCompChart = agCharts.AgChart.create({
    container: document.getElementById('global_score_compar_chart'),
    title: {
      text: 'Global score of pathways',
      fontWeight: 'normal',
    },
    legend: {
      enabled: false,
    },
    series: [{
      type: 'column',
      xKey: 'pathway_name',
      yKeys: ['global_score'],
      tooltip: {
        renderer: function (params) {
          return {
            content: parseFloat(params.yValue).toFixed(3),
          };
        },
      },
    }],
    axes: [{
      type: 'category',
      position: 'bottom',
      label: {
        rotation: 335
      }
      },
      {
      type: 'number',
      position: 'left',
      }],
  });

  // create fba_obj_fraction comparative chart
  fbaCompChart = agCharts.AgChart.create({
    container: document.getElementById('fba_compar_chart'),
    title: {
      text: 'FBA of pathways',
      fontWeight: 'normal',
    },
    legend: {
      enabled: false,
    },
    series: [{
      type: 'column',
      xKey: 'pathway_name',
      yKeys: ['fba_obj_fraction'],
      tooltip: {
        renderer: function (params) {
          return {
            content: parseFloat(params.yValue).toFixed(3),
          };
        },
      },
    }],
    axes: [{
      type: 'category',
      position: 'bottom',
      label: {
        rotation: 335
      }
      },
      {
      type: 'number',
      position: 'left',
      }],
  });

  // create dfG_prime_m comparative chart
  dfgCompChart = agCharts.AgChart.create({
    container: document.getElementById('dfg_compar_chart'),
    title: {
      text: 'ΔG\'m of pathways',
      fontWeight: 'normal',
    },
    legend: {
      enabled: false,
    },
    series: [{
      type: 'column',
      xKey: 'pathway_name',
      yKeys: ['dfG_prime_m'],
      tooltip: {
        renderer: function (params) {
          return {
            content: parseFloat(params.yValue).toFixed(3),
          };
        },
      },
    }],
    axes: [{
      type: 'category',
      position: 'bottom',
      label: {
        rotation: 335
      }
      },
      {
      type: 'number',
      position: 'left',
      }],
  });

  // create rule score comparative chart
  ruleScoreCompChart = agCharts.AgChart.create({
    container: document.getElementById('rule_score_compar_chart'),
    title: {
      text: 'Rule score of pathways',
      fontWeight: 'normal',
    },
    legend: {
      enabled: false,
    },
    series: [{
      type: 'column',
      xKey: 'pathway_name',
      yKeys: ['mean_rule_score'],
      tooltip: {
        renderer: function (params) {
          return {
            content: parseFloat(params.yValue).toFixed(3),
          };
        },
      },
    }],
    axes: [{
      type: 'category',
      position: 'bottom',
      label: {
        rotation: 335
      }
      },
      {
      type: 'number',
      position: 'left'
      }],
  });

  // create reactions_nb (steps) comparative chart
  reactionsNbCompChart = agCharts.AgChart.create({
    container: document.getElementById('reactions_nb_compar_chart'),
    title: {
      text: 'Steps of pathways',
      fontWeight: 'normal',
    },
    legend: {
      enabled: false,
    },
    series: [{
      type: 'column',
      xKey: 'pathway_name',
      yKeys: ['nb_reactions'],
      tooltip: {
        renderer: function (params) {
          return {
            content: parseFloat(params.yValue).toFixed(0),
          };
        },
      },
    }],
    axes: [{
      type: 'category',
      position: 'bottom',
      label: {
        rotation: 335
      }
      },
      {
      type: 'number',
      position: 'left',
      }],
  });
});
