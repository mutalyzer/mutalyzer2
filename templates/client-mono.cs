/*
  Example SOAP client for the Mutalyzer webservice in C# for the Mono
  platform.

  See {path}/webservices

  Compilation instructions:
    wsdl '{path}/services/?wsdl'
    gmcs /target:library Mutalyzer.cs -r:System.Web.Services
    gmcs /r:Mutalyzer.dll client-mono.cs

  Usage:
    ./client-mono.exe 'NM_002001.2:c.1del'

  This code is in the public domain; it can be used for whatever purpose
  with absolutely no restrictions.
*/

using System;

class MutalyzerClient {{

    public static void Main(String [] args) {{

        String variant;
        checkSyntax argument;
        Mutalyzer mutalyzer;
        CheckSyntaxOutput result;

        if (args.Length < 1) {{
            Console.WriteLine("Please provide a variant");
            Environment.Exit(1);
        }}

        variant = args[0];

        mutalyzer = new Mutalyzer();

        Console.WriteLine("Checking " + variant + "...");

        argument = new checkSyntax();
        argument.variant = variant;
        result = mutalyzer.checkSyntax(argument).checkSyntaxResult;

        if (result.valid)
            Console.WriteLine("Valid!");
        else
            Console.WriteLine("Not valid!");

        foreach (SoapMessage m in result.messages)
            Console.WriteLine(String.Format("Message ({{0}}): {{1}}",
                                            m.errorcode, m.message));

    }}

}}
