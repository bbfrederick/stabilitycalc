#!/usr/bin/env python

#
# HTML formatting
#


def breaktag(thestring):
    return "\n<br />" + thestring


def paratag(thestring):
    return "<p>" + thestring + "</p>\n"


def imagetag(imagename, imagewidth):
    return '<p><a href="' + imagename + '"><IMG BORDER=0 SRC=' + imagename + ' width="' + str(imagewidth) + '"></a></p>\n'


def tablepropsstarttag(tablewidth, align):
    return '<table rules="none" style="margin-center: 20px; text-align: ' + align + '"border="0" width="' + str(tablewidth) + '" cellpadding="2">\n'


def tablepropsendtag():
    return '</table>\n'


def tablepropstag(thestring, tablewidth, align):
    return tablepropsstarttag(tablewidth, align) + thestring + tablepropsendtag()


def tablerowstarttag():
    return "<tr>"


def tablerowendtag():
    return '</tr>\n'


def tablerowtag(thestring):
    return tablerowstarttag() + thestring + tablerowendtag()


def tableentrytag(thestring):
    return "<td>" + thestring + "</td>"


def valignopt(thealign):
    return 'valign="' + thealign + '" '


def widthopt(thewidth):
    return 'width="' + str(thewidth) + '%" '


def tableentryopttag(thestring, theopts):
    return "<td " + theopts + ">" + thestring + "</td>"


def headertag(thestring):
    return "<h3>" + thestring + "</h3>"


def bigheadertag(thestring):
    return "<h2>" + thestring + "</h3>"


def smallheadertag(thestring):
    return "<h4>" + thestring + "</h4>"


def bigtag(thestring):
    return "<big>" + thestring + "</big>"


def boldtag(thestring):
    return "<b>" + thestring + "</b>"


def greentag(thestring):
    return '<FONT COLOR="00ff00">' + thestring + "</FONT>"


def redtag(thestring):
    return '<FONT COLOR="ff0000">' + thestring + "</FONT>"


def yellowtag(thestring):
    return '<FONT COLOR="ffff00">' + thestring + "</FONT>"


def qualitytag(thestring, thequality):
    if thequality == 0:
        return greentag(thestring)
    if thequality == 1:
        return yellowtag(thestring)
    return redtag(thestring)
