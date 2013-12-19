function updateVisibility() {
  document.getElementById('file_label').style.display = "none";
  document.getElementById('url_label').style.display = "none";
  document.getElementById('gene_label').style.display = "none";
  document.getElementById('range_label').style.display = "none";
  document.getElementById('chrname_label').style.display = "none";

  for (i = 0; i < document.invoer.invoermethode.length; i++) {
    if (document.invoer.invoermethode[i].checked) {
      if (document.invoer.invoermethode[i].value == 'file') {
          document.getElementById('file_label').style.display = "";
      }
      else if (document.invoer.invoermethode[i].value == 'url') {
          document.getElementById('url_label').style.display = "";
      }
      else if (document.invoer.invoermethode[i].value == 'gene') {
          document.getElementById('gene_label').style.display = "";
      }
      else if (document.invoer.invoermethode[i].value == 'chr') {
          document.getElementById('range_label').style.display = "";
      }
      else if (document.invoer.invoermethode[i].value == 'chrname') {
          document.getElementById('chrname_label').style.display = "";
      }
    }//if
  }//for
}//updateVisibility

//Toggle the build option in the batch.html page
function changeBatch(sel) {
    var opt = sel.options[sel.selectedIndex].value;
    if(opt=='PositionConverter') {
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


var doInterval;
function onloadBatch() {
    changeBatch(document.getElementById('batchType'));
    doInterval = setInterval(updatePercentage, 3000);
}

// Get the HTTP Object
function getHTTPObject(){
    if (window.ActiveXObject)
        return new ActiveXObject("Microsoft.XMLHTTP");
    else if (window.XMLHttpRequest)
        return new XMLHttpRequest();
    else {
        alert("Your browser does not support AJAX.");
        return null;
    }
}

function updatePercentage() {
    if (!document.getElementById('resultID')){ return; };
    var id = document.getElementById('resultID').value;
    var total = document.getElementById('totalJobs').value;
    var url = 'progress?resultID='+id+'&totalJobs='+total+'&ajax=1';
    var val = "";

    http = getHTTPObject();
    if (http == null){
        val = "Your browser does not support Ajax, swith to a different ";
        val +="or wait for the email to arrive."
        document.getElementById('percent').innerHTML = val;
        clearInterval(doInterval);
        return null;
    }

    http.open("GET", url, true);
    http.onreadystatechange=function() {
      if(http.readyState == 4) {
          if (http.responseText == "OK"){
              val = "Your job is finished, results can be downloaded from: ";
              val += "<a href='Results_"+id+".txt'>here</a>";
              clearInterval(doInterval);
          } else if (isNaN(http.responseText)){
                  //Something went wrong
                  val = "Updating went wrong please wait for your email";
                  clearInterval(doInterval);
          } else {
                  val = "Your job is in progress and currently at ";
                  val += parseInt(http.responseText);
                  val += "%.";
          }
      }
      document.getElementById('percent').innerHTML = val;
    }
    http.send(null);
}

function clearField(form, fieldName) {
    for (var i = 0; i < form.elements.length; i++) {
        if (form.elements[i].name == fieldName) {
            form.elements[i].value = '';
        }
    }
}
