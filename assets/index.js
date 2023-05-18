import cytoscape from "cytoscape";
import dagre from "cytoscape-dagre";
(() => {
  const goID = document.body.querySelector("#go-id");
  const goName = document.body.querySelector("#go-name");
  const aspect = document.body.querySelector("#aspect");
  const definition = document.body.querySelector("#definition");
  const pagerank = document.body.querySelector("#pagerank")
  const keywordInputBtn = document.body.querySelector
("#keyword-input-btn");
  const keywordInput = document.body.querySelector("#keyword-input");
  const keywordInputDiv = document.body.querySelector("#keyword-input-div");
  const filterN = document.body.querySelector("#filter-n")
  keywordInputBtn.addEventListener("click", (e) => {
    // GET request to go-viz/keyword.json
    const req = new XMLHttpRequest();
    const loading = document.createElement("p");
    const node = document.createTextNode("Loading...");
    req.onreadystatechange = () => {
      if (req.readyState === XMLHttpRequest.DONE) {
        if (req.status === 200) {
          // console.log(req.responseText);
          // console.log(typeof req.responseText);
          if (req.responseText == '{}') {
            alert("No GO terms matching " + keywordInput.value + ".")
          } else {
            const cytoGraphJson = JSON.parse(req.responseText)
            cytoRender(cytoGraphJson)
          }
        } else {
          alert("Error: Server error.");
        }
        loading.remove();
      }
    };
    req.open("GET", "/go-viz/" + keywordInput.value + ".json" + "/" + filterN.value + "/", true);
    loading.appendChild(node);
    keywordInputDiv.appendChild(loading);
    req.send();
  });
  function cytoRender(cytoJson) {
    cytoscape.use(dagre);
    let cy = cytoscape({
      container: document.getElementById("cy"),
      elements: cytoJson,
      style: [
        // Define node and edge styles here
        {
          selector: "node",
          style: {
            content: "data(id)",
            "text-valign": "center",
            "text-halign": "center",
            "text-wrap": "wrap",
            color: "white",
            "text-outline-width": 2,
            "text-outline-color": "data(color)",
            "background-color": "data(color)",
            label: "data(id)",
            "text-max-width": "100%",
            'width': 'data(size)',
            'height': 'data(size)',
          },
        },
        {
          selector: "edge",
          style: {
            width: 2,
            "line-color": "#b2bec3",
          },
        },
        {
          selector: "node:selected",
          style: {
            "background-color": "#000",
            "text-outline-color": "#000",
          },
        },
      ],
      layout: {
        name: "dagre", // Use the dagre layout for a hierarchical tree layout
        spacingFactor: 1.3,
      },
    });
    // let nodes = cy.nodes().sort((a, b) => {
    //   return a.data('pagerank') - b.data('pagerank')
    // })
    // nodes.slice(0,10).forEach((e)=> {
    //   cy.style.selector
    // })
    // Add a zoom control to the graph
    cy.zoomingEnabled(true);
    // Set up event listeners to highlight nodes and edges when clicked
    cy.on("tap", "node", function (event) {
      var node = event.target;
      // var neighborhood = node.neighborhood().add(node);
      // cy.elements().addClass("faded");
      // neighborhood.removeClass("faded");
      // console.log(node)

      goID.textContent = node.data("id");
      goName.textContent = node.data("name");
      aspect.textContent = node.data("aspect");
      definition.textContent = node.data("definition");
      pagerank.textContent = node.data("pagerank")
    });
  }
  
})();
