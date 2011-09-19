var a_Act;
var a_navAct;
var timerID;
var timerImg;

function setFont(font, hrefID) {
  document.getElementById(hrefID).className = font;
}//setFont

function setImage(imgSrc, imgID) {
  document.getElementById(imgID).src = imgSrc;
}//setImage

function setSub(state, subID) {
  var subElt = document.getElementById(subID + "_0");
  var i = 1;

  while (subElt) {
    subElt.style.display = state;
    subElt = document.getElementById(subID + "_" + i);
    i++;
  }//while
}//setSub

function setActive(newAct) {
  if (document.getElementById("page_" + newAct)) {
    if (newAct.match("_")) {
      setFont("vertnavsubactive", "page_" + newAct);
      setImage("base/images/bullitlicht2.gif", "b_" + newAct);
    }//if
    else {
      setFont("vertnavactive", "page_" + newAct);
      setImage("base/images/bullitlicht1.gif", "b_" + newAct);
    }//else
  
    setSub("", newAct);
  }//if
}//setActive

function setNavActive(newAct) {
  if (document.getElementById("page_" + newAct)) {
    if (newAct.match("_"))
      setImage("base/images/bullitlicht2.gif", "b_" + newAct);
    else
      setImage("base/images/bullitlicht1.gif", "b_" + newAct);
  
    setSub("", newAct);
  }//if
}//setNavActive

function setDeActive(oldAct) {
  if (document.getElementById("page_" + oldAct)) {
    setFont("vertnavsub", "page_" + oldAct);
  
    if (oldAct.match("_"))
      setImage("base/images/bullitmiddel.gif", "b_" + oldAct);
    else
      setImage("base/images/bullitdonker.gif", "b_" + oldAct);
  }//if
}//setDeActive

function alterActive(a_Alter, method) {
  var Label = a_Alter[0];
  var i;

  for (i = 1; i <= a_Alter.length; i++) {
    if (method)
      setActive(Label);
    else {
      setDeActive(Label);
      setSub("none", Label);
    }//else
    Label = Label + "_" + a_Alter[i];
  }//for
}//alterActive

function alterNavActive(a_Alter, method) {
  var Label = a_Alter[0];
  var i;

  for (i = 1; i <= a_Alter.length; i++) {
    if (method)
      setNavActive(Label);
    else
      setDeActive(Label);
    Label = Label + "_" + a_Alter[i];
  }//for
}//alterNavActive

function navPartDeAct() {
  alterNavActive(a_navAct, 0);
  alterActive(a_Act, 1);
}//navPartDeAct

function swapActive(newAct) {
  var a_newAct = newAct.split("_");

  alterActive(a_Act, 0);
  alterActive(a_newAct, 1);

  a_Act = a_newAct;
  a_navAct = a_newAct;
}//swapActive

function swapNavActive() {
  var a_newAct = timerImg.split("_");

  navPartDeAct();
  alterActive(a_newAct, 1);

  a_navAct = a_newAct;
}//swapNavActive

function navAct(imgSrc, imgID) {
  clearTimeout(timerID);

  timerImg = imgID;
  timerID = setTimeout("swapNavActive();", 500);
}//navAct

function navDeAct(imgSrc, imgID) {
  clearTimeout(timerID);

  timerID = setTimeout("navPartDeAct();", 2000);
}//navDeAct

function initActive() {
  var winLoc;

  winLoc = window.location.href;
  winLoc = winLoc.replace(/https?:\/\/[^\/]*\//, "");

  if (winLoc.match("~"))
    winLoc = winLoc.replace(/[^\/]*\//, "");
  else
    winLoc = winLoc.replace(/[^\/]*\/[^\/]*\//, "");

  winLoc = winLoc.replace(/\.cgi$/, "");
  winLoc = winLoc.replace(/index\.shtml/, "");
  winLoc = winLoc.replace(/\/$/, "");
  winLoc = winLoc.replace("#", "");

  a_Act = winLoc.split("/");
  // This is a Quick Hack (tm)
  if (a_Act[0] == '2.0' || a_Act[0] == 'mutalyzer')
    a_Act.shift();
  a_navAct = a_Act;

  alterActive(a_Act, 1);
}//initActive

function resetAll() {
  var a_Menu = document.getElementsByTagName("tr");
  var i;

  for (i = 0; i < a_Menu.length; i++)
    if (a_Menu[i].id.match("_"))
      a_Menu[i].style.display = "none";

  navPartDeAct();
}//resetAll

function clearForm(form, elementName) {
  for (i = 0; i < form.elements.length; i++) {
    if (form.elements[i].name == elementName) {
      form.elements[i].value = "";
    }//if
  }//for
}//clearForm
