<?php
  header('Content-type: text/plain; charset=UTF-8');

  $execstr = "python src/c2g.py " . 
             escapeshellarg($_GET['LOVD_ver']) . ' ' . 
             escapeshellarg($_GET['build']) . ' ' . 
             escapeshellarg($_GET['acc']) . ' ' . 
             escapeshellarg($_GET['var']);

  @exec($execstr, $aRes, $nRes);
  $sRes = implode("\n", $aRes);

  printf("%s\n", $sRes);
?>
