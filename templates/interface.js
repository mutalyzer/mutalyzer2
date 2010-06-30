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
