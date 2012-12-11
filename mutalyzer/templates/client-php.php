<?php

/*
 Example SOAP client for the Mutalyzer web service in PHP.

 See {path}/webservices

 This code is in the public domain; it can be used for whatever purpose
 with absolutely no restrictions.
*/

$URL = '{path}/services/?wsdl';

?><!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Mutalyzer SOAP client</title>
</head>
<body>

<h1>Mutalyzer SOAP client</h1>

<?php

if (isset($_GET['variant']) && $_GET['variant']) {{

    $variant = $_GET['variant'];

    echo '<h2>Result for '.htmlentities($variant).'</h2>';

    // http://www.dotvoid.com/2008/10/soap-structures-in-php/
    $options = array('features' => SOAP_SINGLE_ELEMENT_ARRAYS);

    $client = new SoapClient($URL, $options);

    $result = $client->checkSyntax(array('variant' => $variant))
                  ->checkSyntaxResult;

    if ($result->valid)
        echo '<p><b>Valid!</b>';
    else
        echo '<p><b>Not</b> valid!';

    if (isset($result->messages->SoapMessage)) {{
        echo '<p>Messages:<ol>';
        foreach ($result->messages->SoapMessage as $message) {{
            echo '<li><code>'.htmlentities($message->errorcode).'</code>: ';
            echo htmlentities($message->message);
        }}
        echo '</ol>';
    }}

}}

?>

<h2>Check the syntax of a variant</h2>

<form action="" method="GET">
    <p>Variant:
    <input name="variant" type="text" value="NM_002001.2:c.1del">
    <input type="submit" value="Check syntax">
</form>

</body>
</html>
