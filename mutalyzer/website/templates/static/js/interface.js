function updateVisibility() {
  document.getElementById('upload_label').style.display = "none";
  document.getElementById('url_label').style.display = "none";
  document.getElementById('slice_gene_label').style.display = "none";
  document.getElementById('slice_accession_label').style.display = "none";
  document.getElementById('slice_chromosome_label').style.display = "none";

  for (i = 0; i < document.invoer.method.length; i++) {
    if (document.invoer.method[i].checked) {
      if (document.invoer.method[i].value == 'upload') {
          document.getElementById('upload_label').style.display = "";
      }
      else if (document.invoer.method[i].value == 'url') {
          document.getElementById('url_label').style.display = "";
      }
      else if (document.invoer.method[i].value == 'slice_gene') {
          document.getElementById('slice_gene_label').style.display = "";
      }
      else if (document.invoer.method[i].value == 'slice_accession') {
          document.getElementById('slice_accession_label').style.display = "";
      }
      else if (document.invoer.method[i].value == 'slice_chromosome') {
          document.getElementById('slice_chromosome_label').style.display = "";
      }
    }//if
  }//for
}//updateVisibility

//Toggle the build option in the batch.html page
function changeBatch(sel) {
	var opt = $(sel).val();
    // var opt = sel.options[sel.selectedIndex].value;
    if(opt=='position-converter') {
        document.getElementById('assembly_name_or_alias').style.display = "";
    } else {
        document.getElementById('assembly_name_or_alias').style.display = "none";
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

function onloadBatch() {
    changeBatch($('input[name="job_type"]:checked'));
}

function clearField(form, fieldName) {
    for (var i = 0; i < form.elements.length; i++) {
        if (form.elements[i].name == fieldName) {
            form.elements[i].value = '';
        }
    }
}

$(document).ready(function() {
    $('.example-input').on('click', function() {
        $('.example-target').val($(this).text());
    });
    $('.example-input-2').on('click', function() {
        $('.example-target-2').val($(this).text());
    });
})
