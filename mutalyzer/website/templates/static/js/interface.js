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
  $('.example-input').on('click', function() {
    var target = document.getElementById($(this).data('for'));

    $(target).val($(this).text());
    return false;
  });

  $('.input-select').on('change', function() {
    var context = document.getElementById(
          $(this).data('context')).getElementsByClassName('subform'),
        target = document.getElementById($(this).data('for')),
        i;

    for (i = 0; i < context.length; i++) {
      context[i].style.display = 'none';
    }
    target.style.display = '';

    return false;
  });
});
