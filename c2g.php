<?php
  $execstr = "python src/c2g.py " . 
             $_GET['LOVD_ver'] . ' ' . $_GET['build'] . ' ' . 
             $_GET['acc'] . ' ' . $_GET['var'];

  @exec($execstr, $aRes, $nRes);
  $sRes = implode("\n", $aRes);

  printf("%s\n", $sRes);
?>
