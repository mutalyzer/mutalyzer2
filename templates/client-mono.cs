/*
Compilation instructions:

  wsdl '<tal tal:replace = "path"></tal>/services/?wsdl'
  gmcs /target:library Mutalyzer.cs -r:System.Web.Services
  gmcs /r:Mutalyzer.dll client-mono.cs

Usage:

  ./client-mono.exe 'NM_002001.2:c.1del'
*/

using System;

class MutalyzerClient {

    public static void Main(string [] args) {

        checkSyntax c = new checkSyntax();
        c.variant = args[0];

        Console.WriteLine("Checking " + c.variant);

        Mutalyzer mutalyzer = new Mutalyzer();

        checkSyntaxResponse result = mutalyzer.checkSyntax(c);

        if (result == null)
            Console.WriteLine("[No result]");
        else
            Console.WriteLine(result.checkSyntaxResult);

    }

}
