/*
 * Mutalyzer Name Generator v0.1
 * http://www.mutalyzer.nl/
 * © 2010 LUMC
 */

// Regular expressions for Field Content Tests
function isTranscript(x){ 
    var re = new RegExp("^[0-9]+$", "g");
    return re.test(x);
}

function isSequenceType(s){
    var re = new RegExp("^[cgnemrp]$", "g");
    return re.test(s);
}

function isDNASequence(s){
    var re = new RegExp("^[ACTG]+$", "g");
    return re.test(s);
}

function isRNASequence(s){
    var re = new RegExp("^[acug]+$", "g");
    return re.test(s);
}

function isProteinSequence1(s){
    var re = new RegExp("^[GALMFWKQESPVICYHRNDT]+$", "g");
    return re.test(s);
}

function isProteinSequence3(s){
    var re = new RegExp("^(Gly|Ala|Leu|Met|Phe|Trp|Lys|Gln|Glu|Ser|Pro|Val|Ile|Cys|Tyr|His|Arg|Asn|Asp|Thr)+$", "g");
    return re.test(s);
}

function isPosition(p){
    var re = new RegExp("^([-\*]?[0-9]+([-+]{1}[ud]?[0-9]+)?)$", "g");
    return re.test(p);
}

function isReference(s){
    if(isLRG(s))
        reference['tVar'].len = "+";
    else
        reference['tVar'].len = "*";
    var re = new RegExp("^[A-Z_0-9]+(\.[0-9]+)?$", "gi");
    return re.test(s);
}

function isLRG(s){
    var re = new RegExp("^LRG_[0-9]+", "g");
    return re.test(s);
}

function isGeneSymbol(s){
    var re = new RegExp("^[A-Z0-9]+$", "gi");
    return re.test(s);
}

//Reference object 

var reference = {
		'refe':	{'name'		: "Reference",
                 'len'   	: "+",
				 'value' 	: "",
				 'ok'		: false,
				 'check'    : isReference,
				 'errStr'   : "should be of the format \"NM_002001.2\""},
        'seqT': {'name'		: "Sequence Type",
                 'len'   	: "1",
				 'value' 	: "c",
				 'ok'		: false,
				 'check'    : isSequenceType,
			     'errStr'   : "You managed to select an impossible Sequence Type. Muppet!"},
        'gSym': {'name'		: "Gene Symbol",
                 'len'   	: "*",
				 'value' 	: "",
				 'ok'		: false,
				 'check'    : isGeneSymbol,
			     'errStr'   : "should only contain Letters and Numbers"},
        'tVar': {'name'		: "Transcript",
                 'len'   	: "*",
				 'value' 	: "",
				 'ok'		: false,
				 'check'    : isTranscript,
				 'errStr'   : "must be a postive integer"},
		'tlc' : 1,
		'number': "",
		
		getType :  function(){
                    var translate = {'c':'cod', 'g':'gen','n':'non','r':'rna',
                                 'm':'mit', 'e':'est','p1':'pr1','p3':'pr3'};
                    var temp = this.seqT.value;
                    temp += this.seqT.value=="p" ? this.tlc : "";
                    return translate[temp];
		},
	    getFlag  : function(name){return this[name].len;},
        getName  : function(name){return this[name].name;},
        getErr   : function(name){return this[name].errStr;},
        getCheck : function(name){return this[name].check;},
        getHGVS  : function(){
            checkElement(reference);
            var ret = this['refe'].value;
            var gs = ""
            if (isLRG(ret)){
                gs = "t"+this['tVar'].value;
                this['gSym'].value = "";
            }
            if (this['gSym'].value!=""){
                if (this['tVar'].value!=""){
                   gs = "("+this.gSym.value+"_v"+pad(this.tVar.value,3)+")";
                } else {
                   gs = "("+this.gSym.value+")";
                }
            }
            if (this.seqT.value == "e")
                return ret+gs+":";
            else
                return ret+gs+":"+this.seqT.value+".";
        }
};

var variants = new Array(); //global storage of the variants

//Type of regex to use for each sequence Type
var seqTypes = {
		'cod'	: {	'pCheck' 	: isPosition,
					'sCheck'	: isDNASequence,
					'errorStr'	: "must consist of nucleotides [ACTG]"},

		'rna'	: {	'pCheck' 	: isPosition,
					'sCheck'	: isRNASequence,
					'errorStr'	: "must consist of nucleotides [acug]"},

		'gen'	: {	'pCheck' 	: isPosition,
					'sCheck'	: isDNASequence,
					'errorStr'	: "must consist of nucleotides [ACTG]"},

		'mit'	: {	'pCheck' 	: isPosition,
					'sCheck'	: isDNASequence,
					'errorStr'	: "must consist of nucleotides [ACTG]"},

		'non'	: {	'pCheck' 	: isPosition,
					'sCheck'	: isDNASequence,
					'errorStr'	: "must consist of nucleotides [ACTG]"},
				
		'est'	: {	'pCheck' 	: isPosition,
					'sCheck'	: isDNASequence,
					'errorStr'	: "must consist of nucleotides [ACTG]"},

		'pr1'	: {	'pCheck' 	: isPosition,
					'sCheck'	: isProteinSequence1,
					'errorStr'	: "must consist of the Single Letter AminoAcids Code"},

		'pr3'	: {	'pCheck' 	: isPosition,
					'sCheck'	: isProteinSequence3,
					'errorStr'	: "must consist of the Three Letter AminoAcids Code"}
}

//Type of field options for each mutation Type
//char 1 -> Position 1, char 2 -> Position 2
//char 3 -> Sequence 1, char 4 -> Sequence 2
//   + 	 -> Minimal 1	   *   -> optional
//   0   -> Hidden		   1   -> 1 atomic length
var mutTypes = {
    "sub" : {
        "flags" : "+011",
        "S1"    : "Deleted",
        "S2"    : "Inserted"},
    "del" : {
        "flags" : "++*0",
        "S1"    : "Deleted",
        "S2"    : ""},
	"ins" : {
        "flags" : "++0+",
        "S1"    : "",
        "S2"    : "Inserted"},
	"dup" : {
        "flags" : "++*0",
        "S1"    : "Duplicated",
        "S2"    : ""},
	"ind" : {
        "flags" : "++*+",
        "S1"    : "Deleted",
        "S2"    : "Inserted"},
	"inv" : {
        "flags" : "++*0",
        "S1"    : "Inverted",
        "S2"    : ""},

    getType : function(i){
        i = parseInt(i)-1;
        translate = ['sub','del','ins','dup','ind','inv'];
        return translate[i];
    }
};


var VariantField = {
    'P1'  :   { 'value'   : "",
          	    'ok'	  : false,
				'check'   : true,
                'index'   : 0},
    'P2'  :   { 'value'   : "",
          	    'ok'	  : false,
				'check'   : true,
                'index'   : 1},
    'S1'  :   { 'value'   : "",
          	    'ok'	  : false,
				'check'   : true,
                'index'   : 2},
    'S2'  :   { 'value'   : "",
          	    'ok'	  : false,
				'check'   : true,
                'index'   : 3},
    //Set these values on new VariantField
    'mType' : "",
	'number' : 0,

	'ok' :    variantOK,
	'removed' : false,

    getFlag : function(name){
        var i = this[name].index;
        return mutTypes[this["mType"]]["flags"].substring(i,i+1);},
    getName : function(name){
        var S1Name = mutTypes[this["mType"]]["S1"]+" Sequence";
        var S2Name = mutTypes[this["mType"]]["S2"]+" Sequence";
        var defaultnames = ["Start Position", "End Position", S1Name, S2Name];
        return defaultnames[this[name].index];},
    getErr  : function(name){
                if(name.substring(1,0)=="P")
                  return "position notation help";
                else
                  return seqTypes[reference.getType()].errorStr;},
    getCheck : function(name){
                if(name.substring(1,0)=="P")
                  return seqTypes[reference.getType()].pCheck;
                else
                  return seqTypes[reference.getType()].sCheck;},
    getHGVS : function(){
        if(this.removed == true) return "";
        var seqT = reference.getType();
        if (seqT=="pr1" || seqT=="pr3")
            return getProteinHGVS(this);
        else
            return getNormalHGVS(this);
    }
}

function clone(obj){
    if(obj == null || typeof(obj) != 'object')
        return obj;
    var temp = new obj.constructor(); 
    for(var tkey in obj)
        temp[tkey] = clone(obj[tkey]);
    return temp;
}

function newVariantField(mutType){
    var newVarField = clone(VariantField);
    newVarField.mType = mutType;
    newVarField.number = variants.length;
    return newVarField;
}

//getvalue from form
function getValue(id){
	var form = document.forms["mainform"];
	return form[id].value;
}

//Update the HTML?
//Use P1 P2 S1 S2 keys to update the DIV given
function checkElement(elem){
	for (var key in elem){
		var error = "";
		var obj = elem[key];

		//check if obj is checkable
		if (obj.check==undefined) continue;
		
		//Get the value from the form
		var IDt = elem["number"];
        if(IDt === "")
            IDt = key;
        else
            IDt = "V"+IDt+key;
		obj.value = getValue(IDt);

        var flag = elem.getFlag(key);
        var optional = (flag == "*") ? "*" : "";
        var name = elem.getName(key)+optional;
        var check = elem.getCheck(key);
        var errStr = elem.getErr(key);
		
		//check if the value is ok
		obj.ok = check(obj.value);


        if (elem["number"]!==""){
            setName(IDt+"name", name);
            show(IDt+"row");
        }
		switch(flag){
			case "0":
                //if the object is not used, don't check it
                hide(IDt+"row");
				continue;
			case "1":
				obj.ok = (obj.value.length>reference.tlc) ? false : obj.ok;
				if(!obj.ok){
					error += name+" incorrect: substitution must consist of a "+
						"single nucleotide / amino acid<br \>"
				}
				break
			case "*":
				obj.ok = !obj.value ? true : obj.ok;
				break
			case "+":
				if(!obj.value){
					obj.ok = false;
					error += name+" required.<br \>";
				}
		}
		if (!obj.ok){
			error += name + " incorrect: "+ errStr;
		}
		//console.log("checked " + IDt + ": " +error);
		document.getElementById(IDt+"error").innerHTML = error;
	}
}


function update(){
    // check & update reference
	updateReference();
    //Reference Fields
    if (variants.length==0 && reference.refe.ok)
        addVariant();
    for (var i=0; i<variants.length; i++){
        updateVariant(variants[i]);
    }
    var hgvs = generateHGVS();
    var ref = hgvs[0];
    var vari = hgvs[1];
    
    var encVar = encodeURIComponent(vari);
    var url = "checkForward?mutationName="+ref+encVar;
    var link = ref+vari;
    var Output = "<a href=\""+url+"\">"+link+"</a>";


    document.getElementById("output").innerHTML = Output;

}

function updateVariant(variant){
    if (variant.removed == true) return;
    var form = document.forms["mainform"];
    var V = "V"+variant["number"];
    variant.mType = mutTypes.getType(form[V+"mutT"].value);
    variant.P1.value = form[V+"P1"].value;
    variant.P2.value = form[V+"P2"].value;
    variant.S1.value = form[V+"S1"].value;
    variant.S2.value = form[V+"S2"].value;
    //Update Names and visible Fields TODO
    checkElement(variant);
}

function updateReference(){
    form = document.forms["mainform"];
    reference.refe.value = form.refe.value;
    reference.seqT.value = form.seqT.value;
    reference.gSym.value = form.gSym.value;
    if (!reference.gSym.value)
        hide("tVar");
    else
        show("tVar");
    if (isLRG(reference.refe.value)){
        show("tVar");
        hide("gSym");
    }
    else{
        show("gSym");
    }
    if (form.seqT.value == "p")
        show("tlc");
    else
        hide("tlc");
    reference.tVar.value = form.tVar.value;
    var tlc = (form.seqT.value == "p" && form.tlc[1].checked) ? 3 : 1;
    reference.tlc = tlc;
    //Update Visible Fields TODO
    checkElement(reference);
}

function generateHGVS(){
    var hgvs = "";
    var ref = reference.getHGVS();
	//concatonate the variants (and add brackets)
    var hgvsvariants = new Array();
    for (var i=0; i<variants.length; i++){
        if(variants[i].removed == true) continue;
        hgvsvariants[hgvsvariants.length] = variants[i].getHGVS();
    }
    hgvsvariant = hgvsvariants.join(";");
    if(hgvsvariants.length>1)
		hgvsvariant = "["+hgvsvariant+"]"

	return [ref,hgvsvariant];
}



// Object helper functions
function variantOK(){
	var names = "P1,P2,S1,S2".split(',');
	for (var i=0; i < 4; i++) {
		if(this[names[i]].ok == false)
			return false;
	}
	return true;
}

function show(id){
    document.getElementById(id).style.display = "";
}

function hide(id){
    document.getElementById(id).style.display = "none";
}

function pad(i, length) {
    var str = '' + i;
    while (str.length < length) {
        str = '0' + str;
    }
    return str;
}

function setName(id, name){
    document.getElementById(id).innerHTML = name;
}

function getProteinHGVS(variant){
    return "not yet implemented";
}

function getNormalHGVS(variant){
    var S1 = variant.S1.value;
    var S2 = variant.S2.value;
    var P1 = variant.P1.value;
    var P2 = variant.P2.value;
    var mutation = (P1 && P2) ? P1+"_"+P2 : P1;
    switch(variant.mType){
    case "sub":
      mutation += S1+">"+S2;
      break;
    case "del":
      mutation += "del"+S1;
      break;
    case "ins":
      mutation += "ins"+S2;
      break;
    case "dup":
      mutation += "dup"+S1;
      break;
    case "ind":
      if (S1 && S2)
          mutation += "del"+S1+"ins"+S2;
      else
          mutation += "delins"+S2;
      break;
    case "inv":
      mutation += "inv"+S1;
      break;
    default:  //Impossible
      mutation = ""
    }
    return mutation;

}

function addVariant() {
    var inside = document.getElementById('variants');
    var variantnmbr = variants.length;
    var newvar = newVariantField("sub");
	variants[variants.length] = newvar;
    var newvariant = document.createElement('fieldset');
    newvariant.setAttribute("id", "variant"+variantnmbr);

    //populate variant with variant template
    var variantHTML = document.getElementById("varianttemplate").innerHTML;
    var finalHTML = "<legend>Variant "+(variantnmbr+1);
    if(variantnmbr>0){
        finalHTML += " <a onclick=\"removeVariant("+variantnmbr+");\""+
            " class=\"remove\">remove</a> ";
    }
    finalHTML += "</legend>";

    finalHTML += variantHTML.replace(/\{NMBR\}/g, variantnmbr);

    newvariant.innerHTML = finalHTML;

    inside.appendChild(newvariant);
    update();
    window.scrollBy(0,200);
}

function removeVariant(nmbr){
    var inside = document.getElementById('variants');
	var delvar = document.getElementById("variant"+nmbr);
	inside.removeChild(delvar);
    variants[nmbr].removed = true;
    update();
}


//###################

function updateoutput(form){
    updateform(form);
    var errors = formvalidate(form);
    var ref = form.reference.value;
    var gs = "";
    //gene Symbol
    if (form.geneSymbol.value!=""){
        if (form.transcript.value!=""){
            gs = "("+form.geneSymbol.value+"_v"+pad(form.transcript.value,3)+")";
        } else {
            gs = "("+form.geneSymbol.value+")";
        }
    } 

    //Sequence Type
    var seqT = form.sequenceType.value;
    var mutT = form.mutationType.value

    //Start Position
    var startP = form.startP.value;

    //End Position
    var endP = form.endP.value;

    //Create Position
    var loc = (startP && endP) ? startP+"_"+endP : startP;

    //Old Sequence
    var S1 = form.oldSequence.value;
    var S2 = form.newSequence.value;

    var mutation = "";

    //Mutations

    if(seqT=="p"){
        switch(parseInt(mutT)){
            case 1:
              mutation = S1+loc+S2
              break;
            case 2:
              if(endP)
                mutation = S1+startP+"_"+S2+endP+"del";
              else
                mutation = S1+startP+"del";
              break;
            case 3:
              //loc incomplete
              mutation = loc+"ins"+S1;
              break
            case 4:
              mutation = S1+startP+"dup";
            case 5:
            case 6:
            case 7:
            default:
                break;
        }
    } else {
        switch(parseInt(mutT)){
            case 1:   //sub
              mutation = S1+">"+S2;
              break;
            case 2:   //del
              mutation = "del"+S1;
              break;
            case 3:   //ins
              mutation = "ins"+S2;
              break;
            case 4:   //dup
              mutation = "dup"+S1;
              break;
            case 5:   //delins
              if (S1 && S2)
                  mutation = "del"+S1+"ins"+S2;
              else
                  mutation = "delins"+S2;
              break;
            case 6:   //inv
              mutation = "inv"+S1;
              break;
            case 7:   //frameshift
              break;
            case 8:   //microsat
              break;
            default:  //Impossible
              break;
        }
    }


    var refNot = ref+gs+":";
    var vari = ""
    switch(seqT){
       case "e":
           vari = loc+mutation;
           break;
       case "p":
           vari = "p.("+mutation+")";
           break;
       default:
           vari = seqT+"."+loc+mutation;
    }
    var encVar = encodeURIComponent(vari);
    var url = "http://www.mutalyzer.nl/2.0/checkForward?mutationName="+refNot+encVar;
    var link = refNot+vari;
    var Output = "<a href=\""+url+"\">"+link+"</a>";

    document.getElementById("output").innerHTML = Output;
    return form;
}


