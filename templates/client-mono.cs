/*
Compilation instructions:

  wsdl '<tal tal:replace = "path"></tal>/service/?wsdl'
  gmcs /target:library MutalyzerService.cs -r:System.Web.Services
  gmcs /r:MutalyzerService.dll client-mono.cs

Usage:

  ./client-mono.exe 'NM_002001.2:c.1del'
*/

using System;

class Mutalyzer {

    public static void Main(string [] args) {

        checkSyntax c = new checkSyntax();
        c.variant = args[0];

        Console.WriteLine("Checking " + c.variant);

        MutalyzerService mutalyzer = new MutalyzerService();

        checkSyntaxResponse result = mutalyzer.checkSyntax(c);

        if (result == null)
            Console.WriteLine("[No result]");
        else
            Console.WriteLine(result.checkSyntaxResult);

    }

}
