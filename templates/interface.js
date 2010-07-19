function updateVisibility() {
  document.getElementById('file_label').style.display = "none";
  document.getElementById('url_label').style.display = "none";
  document.getElementById('gene_label').style.display = "none";
  document.getElementById('range_label').style.display = "none";

  for (i = 0; i < document.invoer.invoermethode.length; i++) {
    if (document.invoer.invoermethode[i].checked) {
      if (document.invoer.invoermethode[i].value == 'file') {
        file_label.style.display = "";
      }
      else if (document.invoer.invoermethode[i].value == 'url') {
        url_label.style.display = "";
      }
      else if (document.invoer.invoermethode[i].value == 'gene') {
        gene_label.style.display = "";
      }
      else if (document.invoer.invoermethode[i].value == 'chr') {
        range_label.style.display = "";
      }                        
    }//if
  }//for
}//updateVisibility

//Toggle the build option in the batch.html page
function changeBatch(sel) {
    var opt = sel.options[sel.selectedIndex].value;
    if(opt=='ConversionChecker') {
        document.getElementById('build').style.display = "";
    } else {
        document.getElementById('build').style.display = "none";
    }
}

function toggle_visibility(id) {
    var e = document.getElementById(id);
    if (e.style.display == 'block') {
        e.style.display = 'none';
    } else {
        e.style.display = 'block';
    }
}
