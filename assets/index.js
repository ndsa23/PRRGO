import cytoscape from "cytoscape";
import dagre from "cytoscape-dagre";
(() => {
  // elements
  const goID = document.body.querySelector("#go-id");
  const goName = document.body.querySelector("#go-name");
  const aspect = document.body.querySelector("#aspect");
  const definition = document.body.querySelector("#definition");
  const pagerank = document.body.querySelector("#pagerank");
  const keywordInputBtn = document.body.querySelector("#keyword-input-btn");
  const keywordInput = document.body.querySelector("#keyword-input");
  const keywordInputDiv = document.body.querySelector("#keyword-input-div");
  const goForm = document.forms.namedItem("go-form");
  const degTable = document.body.querySelector("#deg-table");
  const degDownloadButton = document.body.querySelector("#download-btn");
  const goImgDownloadButton = document.body.querySelector("#go-btn");
  let degData = {};
  let cy;
  // submit form
  keywordInputBtn.addEventListener("click", (e) => {
    // create "Loading . . ." text
    const loading = document.createElement("p");
    const node = document.createTextNode("Loading...");
    loading.appendChild(node);

    const formData = new FormData(goForm);
    const req = new XMLHttpRequest();

    // AJAX Handler
    req.onreadystatechange = (progress) => {
      if (req.readyState === XMLHttpRequest.DONE) {
        // Success
        if (req.status === 200) {
          // console.log(req.responseText);
          // console.log(typeof req.responseText);
          if (req.responseText == "{}") {
            alert("No GO terms matching " + keywordInput.value + ".");
          } else {
            console.log(req.responseText);
            const parsedData = JSON.parse(req.responseText);
            degData = parsedData["go_data"];
            const cytoGraphData = parsedData["graph"];
            cytoRender(cytoGraphData);
            console.log(degData);
          }

          // Error
        } else {
          alert("Error: Server error.");
        }
        loading.remove();
      }
    };
    // POST request used to transmit form data
    req.open("POST", "/submit/", true);
    keywordInputDiv.appendChild(loading);
    req.send(formData);
    e.preventDefault();
  });

  // renders the GO network graph using the cytoscapeJS JSON
  function cytoRender(cytoJson) {
    cytoscape.use(dagre);
    cy = cytoscape({
      container: document.getElementById("cy"),
      elements: cytoJson,
      style: [
        // Define node and edge styles
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
            width: "data(size)",
            height: "data(size)",
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
    cy.zoomingEnabled(true);

    // Set up event listeners to highlight nodes and edges when clicked
    cy.on("tap", "node", function (event) {
      var node = event.target;
      // add relevant info to GO table
      goID.textContent = node.data("id");
      goName.textContent = node.data("name");
      aspect.textContent = node.data("aspect");
      definition.textContent = node.data("definition");
      pagerank.textContent = node.data("pagerank");
      renderDegTable(node.data("id"));
    });
  }

  // parses DEG data and renders DEG table
  function renderDegTable(goID) {
    // clear DEG table
    while (degTable.firstChild) {
      degTable.removeChild(degTable.firstChild);
    }
    // if there is no goID in degData there were no matching DEGs
    if (!(goID in degData)) {
      let trHead = document.createElement("tr");
      let th = document.createElement("th");
      let thContent = document.createTextNode("No matching DEGs.");
      th.appendChild(thContent);
      trHead.appendChild(th);
      degTable.appendChild(trHead);
      return;
    }

    const goIDDegData = degData[goID];

    // make table header row
    let trHead = document.createElement("tr");
    goIDDegData["columns"].forEach((element) => {
      let th = document.createElement("th");
      let thContent = document.createTextNode(element);
      th.appendChild(thContent);
      trHead.appendChild(th);
    });
    degTable.appendChild(trHead);

    // make other rows of the table
    console.log(goIDDegData["data"]);
    goIDDegData["data"].forEach((row) => {
      let tr = document.createElement("tr");
      row.forEach((element) => {
        let td = document.createElement("td");
        let tdContent = document.createTextNode(element);
        td.appendChild(tdContent);
        tr.appendChild(td);
      });
      degTable.appendChild(tr);
    });
  }
  // DEG Download Button Handler
  // Add a click event listener to the button
  degDownloadButton.addEventListener("click", function () {
    // Create a new XMLHttpRequest object
    const xhr = new XMLHttpRequest();

    // Configure the request
    xhr.open("POST", "/download/", true);
    xhr.responseType = "blob"; // Set the response type to blob
    const formData = new FormData(goForm);
    // Set up the callback function for when the request completes
    xhr.onload = function () {
      if (xhr.status === 200) {
        // Extract the filename from the Content-Disposition header
        const contentDisposition = xhr.getResponseHeader("Content-Disposition");
        console.log(contentDisposition)
        const filenameMatch = contentDisposition.match(/filename="(.+)"/);
        const filename = filenameMatch ? filenameMatch[1] : "downloaded_file";

        // Create a temporary URL for the response blob
        const url = window.URL.createObjectURL(xhr.response);

        // Create a temporary <a> element to initiate the download
        const link = document.createElement("a");
        link.href = url;
        link.download = filename; // Set the file name and extension from the server
        document.body.appendChild(link);
        link.click();

        // Clean up the temporary URL and element
        document.body.removeChild(link);
        window.URL.revokeObjectURL(url);

        console.log("Download request successful");
      } else {
        console.log("Download request failed");
        // add error handling code here
      }
    };

    // Send the request
    xhr.send(formData);
  });

  // GO Network Visualization Download Button
  goImgDownloadButton.addEventListener("click", function() {
      // Export Cytoscape.js graph as jpg
      cy.jpeg({ 
        output: 'blob-promise', // Request the image as a Blob object
        full: true, // Capture the full graph including elements outside the viewport
        scale: 2, // Adjust the scale factor as needed
      }).then(function(blob) {
        // Create a temporary <a> element to initiate the download
        const link = document.createElement('a');
        link.href = URL.createObjectURL(blob);
        link.download = keywordInput.value.replace(' ', '_') + '_graph.png'; // Set the desired file name
    
        // Append the link to the document body and programmatically trigger the download
        document.body.appendChild(link);
        link.click();
    
        // Clean up the temporary URL and element
        document.body.removeChild(link);
        URL.revokeObjectURL(link.href);
      });
  })
})();
