// Toggle the build option in the batch.html page.
function changeBatch(sel) {
  var opt = $(sel).val();

  if(opt == 'position-converter') {
    $('#assembly_name_or_alias').show();
  }
  else {
    $('#assembly_name_or_alias').hide();
  }
}

function toggle_visibility(id) {
  $(document.getElementById(id)).toggle();
}

function onloadBatch() {
  changeBatch($('input[name="job_type"]:checked'));
}

function clearField(form, fieldName) {
  var i;

  for (i = 0; i < form.elements.length; i++) {
    if (form.elements[i].name == fieldName) {
      form.elements[i].value = '';
    }
  }
}

$(document).ready(function() {
  $(document).on('click', '.example-input', function(event) {
    var target = document.getElementById($(this).data('for'));

    $(target).val($(this).text());
    event.preventDefault();
    event.stopPropagation();
  });

  $(document).on('change', '.input-select', function(event) {
    var context = document.getElementById(
          $(this).data('context')).getElementsByClassName('subform'),
        target = document.getElementById($(this).data('for')),
        i;

    for (i = 0; i < context.length; i++) {
      context[i].style.display = 'none';
    }
    target.style.display = '';

    event.preventDefault();
    event.stopPropagation();
  });
});
