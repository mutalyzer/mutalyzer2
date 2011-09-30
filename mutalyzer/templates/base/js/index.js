<!--

// --- GET ELEMENT(LAYERREFERENCE) ---
function getElt ()
{
  if (document.layers)
  { var currentLayer = document.layers[getElt.arguments[0]];
    for (var i = 1; i < getElt.arguments.length && currentLayer; i++)
    { currentLayer = currentLayer.document.layers[getElt.arguments[i]];
    }
    return currentLayer;
  }
  else if (document.all)
  { var elt = eval('document.all.' + getElt.arguments[getElt.arguments.length - 1]);
    return(elt);
  }
  else if (document.getElementById)
  { var elt = document.getElementById(getElt.arguments[getElt.arguments.length - 1]);
    return(elt);
  }
}

// --- IMAGE SWAP ---
function swapImage(imgSrc, imgID, elt)
{
  if (document.getElementsByName)
  { var img = document.getElementsByName(imgID);
    img[0].src = imgSrc;
  }
  else
  { if (swapImage.arguments.length == 3)
    { eval("elt.document." + imgID + ".src = '" + imgSrc + "'");
    }
    else
    { eval("document." + imgID + ".src = '" + imgSrc + "'");
    }
  }
}

// --- IMAGE PRELOAD ---
function preloadImages()
{
  if (document.images)
  { var imgStr = preloadImages.arguments;
    if (!document.preloadArray)
    { document.preloadArray = new Array();
    }

    var n = document.preloadArray.length;
    for (var i = 0; i < preloadImages.arguments.length; i++)
    { document.preloadArray[n] = new Image;
      document.preloadArray[n].src = imgStr[i];
      n++;
    }
  }
}

// -->
