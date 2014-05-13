# -*- coding: utf-8 -*-
# PhylDiag v1.02
# python 2.7
# Copyright © 2013 IBENS/Dyogen : Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France

#display_prog = 'display' # Command to execute to display images.
display_prog = 'firefox' # Command to execute to display images.
from math import sqrt, acos, cos, sin

"""
The following code is a lightweight wrapper around SVG files. The metaphor
is to construct a scene, add objects to it, and then write it to a file
to display it.

Warning, some svgclass may request the drawHomologyMatrix.css file, copied at the end of this file
Rq : this css file is consistent with genomicus notations
"""

def colorstr(rgb): return "#%x%x%x" % (rgb[0]/16,rgb[1]/16,rgb[2]/16) if rgb != 'none' else 'none'

class Scene:
    def __init__(self,name="svg",width=400, height=400):
        self.name = name
        self.items = []
        self.height = height
        self.width = width
        return

    def add(self,item): self.items.append(item)

    def strarray(self):
        #var = ["<?xml version=\"1.0\"?>\n",
        #       "<svg height=\"%d\" width=\"%d\" >\n" % (self.height,self.width),
        #       " <g style=\"fill-opacity:1.0; stroke:black;\n",
        #       "  stroke-width:1;\">\n"]
	var= ['<?xml version="1.0" encoding="utf-8" standalone="no"?>\n',
		'<?xml-stylesheet type="text/css" href="drawHomologyMatrix.css" ?>\n', #Attention, necessite le fichier css avec la charte graphique de genomicus
		'<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n',
		#"<svg height=\"%spt\" version=\"1.1\" viewBox=\"0 0 %s %s\" width=\"%spt\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n" % (self.height, self.width, self.height, self.width),
		"<svg height=\"100%%\" version=\"1.1\" viewBox=\"0 0 %s %s\" width=\"100%%\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n" % (self.width, self.height),
		 '<defs>\n',
		  '<style type="text/css">\n',
			'*{stroke-linecap:square;stroke-linejoin:round;}\n',
		  '</style>\n',
		 '</defs>\n'
		 '<g style="fill-opacity:1.0; stroke:black;\n',
		  'stroke-width:0.25;">\n']

        for item in self.items: var += item.strarray()            
        var += [" </g>\n","</svg>\n"]
        return var

    def write_svg(self,filename=None):
        if filename:
            self.svgname = filename
        else:
            self.svgname = self.name + ".svg"
        file = open(self.svgname,'w')
        file.writelines(self.strarray())
        file.close()
        return

    def display(self,prog=display_prog):
        #os.system("%s %s" % (prog,self.svgname))
        return        
        

class Line:
    def __init__(self,start,end, width=0.03):
        self.start = start #xy tuple
        self.end = end     #xy tuple
        self.width = width    #xy tuple
        return

    def strarray(self):
	    return ["  <line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" style=\"stroke-width:%f\" />\n" %\
                (self.start[0],self.start[1],self.end[0],self.end[1],self.width)]


class Circle:
    def __init__(self,center,radius,color):
        self.center = center #xy tuple
        self.radius = radius #xy tuple
        self.color = color   #rgb tuple in range(0,256)
        return

    def strarray(self):
        return ["  <circle cx=\"%f\" cy=\"%f\" r=\"%f\"\n" %\
                (self.center[0],self.center[1],self.radius),
                "    style=\"fill:%s;\"  />\n" % colorstr(self.color)]

class Rectangle:
    def __init__(self,origin, height, width, fill_opacity=0.65, fill='none', stroke=None, svgClass=None, strokeWidth=0.3):
        self.origin = origin
        self.height = height
        self.width = width
        self.fill_opacity = fill_opacity #value of alpha between 0.0 and 1.0
	self.stroke = stroke
	self.strokeWidth = strokeWidth
	if svgClass == None:
		self.fill = fill
		self.svgClass = None
	else:
		self.fill = None
		self.svgClass = svgClass
        return

    def strarray(self):
	if self.svgClass == None:
		return ["  <rect x=\"%f\" y=\"%f\" height=\"%f\"\n" %\
        	       	(self.origin[0],self.origin[1],self.height),
			"    width=\"%f\" style=\"fill:%s;fill-opacity:%f;stroke:%s;stroke-width:%s\"/>\n" %\
               		(self.width,colorstr(self.fill), self.fill_opacity, self.stroke, self.strokeWidth)]
	else:
		return 	["  <rect x=\"%f\" y=\"%f\" height=\"%f\"\n" %\
        	       	(self.origin[0],self.origin[1],self.height),
			"    width=\"%f\" style=\"fill-opacity:%f;stroke:%s\"\n" %\
	               	(self.width, self.fill_opacity, self.stroke),\
			"    class=\"%s\"/>\n" % self.svgClass ]

class Text:
    def __init__(self,origin, text,size=24, text_anchor="start", fill=(0,0,0), stroke=None, fontWeight="800", fontFamily="Arial", transform=""): # font-family does not work
        self.origin = origin
        self.text = text
        self.size = size
	self.text_anchor = text_anchor #specify where is the anchor : "start" (left), "middle" (center), ...
	self.fill = fill
	self.stroke = stroke
	self.fontWeight = fontWeight
	self.fontFamily = fontFamily
	self.tansform = transform
        return

    def strarray(self):
	    return ["  <text x=\"%f\" y=\"%f\" font-size=\"%f\" text-anchor=\"%s\" fill=\"%s\" stroke=\"%s\" font-weight=\"%s\" font-family=\"%s\" transform=\"%s\" >\n" %\
                (self.origin[0],self.origin[1]+0.5*self.size,self.size, self.text_anchor, colorstr(self.fill), 'none' if self.stroke == None else colorstr(self.stroke), self.fontWeight, self.fontFamily, self.tansform),
                "   %s\n" % self.text,
                "  </text>\n"]


class Ellipse:
    def __init__(self,center,radius_x,radius_y,fill_color,line_color,line_width):
        self.center = center
        self.radiusx = radius_x
        self.radiusy = radius_y
        self.fill_color = fill_color
        self.line_color = line_color
        self.line_width = line_width
    def strarray(self):
        return ["  <ellipse cx=\"%f\" cy=\"%f\" rx=\"%f\" ry=\"%f\"\n" %\
                (self.center[0],self.center[1],self.radius_x,self.radius_y),
                "    style=\"fill:%s;stroke:%s;stroke-width:%f\"/>\n" % (colorstr(self.fill_color),colorstr(self.line_color),self.line_width)]

class Polygon:
    #def __init__(self,points,fill_color=(150,0,0),line_color=(0,0,0),line_width=1):
    #    self.points = points
    #    self.fill_color = fill_color
    #    self.line_color = line_color
    #    self.line_width = line_width
    def __init__(self,points,stroke_width=0.1, SVGclass=None): #fill color peut varier entre 0 et 44 (compris)
        self.points = points
	self.stroke_width = stroke_width
	self.SVGclass = SVGclass
    def strarray(self):
        polygon="<polygon points=\""
        for point in self.points:
            polygon+=" %f,%f" % (point[0],point[1])
        #return [polygon,\
        #       "\" \nstyle=\"fill:%s;stroke:%s;stroke-width:%f\"/>\n" %\
        #       (colorstr(self.fill_color),colorstr(self.line_color),self.line_width)]
        return [polygon,"\"\n",\
			"class=\"%s\" style=\"stroke-width:%s\"" % (self.SVGclass,self.stroke_width) , "/>\n"]

class Gene:
    def __init__(self, start, end, width=0.4, strand=+1, stroke_width=0.1, SVGclass = None, text=None):
	self.start = start
	self.end = end
	self.width = width
	self.strand = strand
	self.stroke_width = stroke_width
	self.SVGclass = SVGclass
	self.text = text
	return

    def strarray(self):
	if self.strand != None:
	    scale = sqrt((self.end[0]-self.start[0])**2 + (self.end[1]-self.start[1])**2)
	    points = [(0.1,-0.5), (0.1,0.5), (0.75,0.5), (0.9, 0), (0.75, -0.5)] if self.strand == +1 else [(0.1,0), (0.25, 0.5), (0.9,0.5), (0.9,-0.5), (0.25, -0.5)] # Defines a gene symbol
	    #rotation
	    theta = acos(float(self.end[0]-self.start[0])/scale) if self.end[1]>self.start[1] else -acos(float(self.end[0]-self.start[0])/scale)
	    points = [(i,float(j*float(self.width))/scale) for (i,j) in points]
	    points = [(i*cos(-theta)+j*sin(-theta), -i*sin(-theta)+j*cos(-theta)) for (i,j) in points]
	    points = [(i*scale+self.start[0],j*scale+self.start[1]) for (i,j) in points]
	    #return Polygon(points,self.color,(0,0,0),1).strarray()
	    if self.text == None:
	 	    return Polygon(points, stroke_width=self.stroke_width, SVGclass=self.SVGclass).strarray() 
	    else:
		    return Polygon(points, stroke_width=self.stroke_width, SVGclass=self.SVGclass).strarray() + Text((float(self.start[0]+self.end[0])/2,float(self.start[1]+self.end[1])/2), self.text,size=self.width*0.8, text_anchor="middle", fill=(0,0,0), stroke=None, fontWeight="800", fontFamily="Arial", transform="").strarray()
        else: 
	    scale = sqrt((self.end[0]-self.start[0])**2 + (self.end[1]-self.start[1])**2)
	    points = [(0.1,-0.5), (0.1,0.5), (0.9,0.5), (0.9, -0.5)] # Defines a rectangle for unoriented genes or TBs
	    #rotation
	    theta = acos(float(self.end[0]-self.start[0])/scale) if self.end[1]>self.start[1] else -acos(float(self.end[0]-self.start[0])/scale)
	    points = [(i,float(j*float(self.width))/scale) for (i,j) in points]
	    points = [(i*cos(-theta)+j*sin(-theta), -i*sin(-theta)+j*cos(-theta)) for (i,j) in points]
	    points = [(i*scale+self.start[0],j*scale+self.start[1]) for (i,j) in points]
	    #return Polygon(points,self.color,(0,0,0),1).strarray()
	    if self.text == None:
		    return Polygon(points, stroke_width=self.stroke_width, SVGclass=self.SVGclass).strarray()
	    else:
		    return Polygon(points, stroke_width=self.stroke_width, SVGclass=self.SVGclass).strarray() + Text((float(self.start[0]+self.end[0])/2,float(self.start[1]+self.end[1])/2), self.text,size=self.width*0.8, text_anchor="middle", fill=(0,0,0), stroke=None, fontWeight="800", fontFamily="Arial", transform="").strarray()

#drawHomologyMatrix.css
"""
# -*- coding: utf-8 -*-
/* css file for drawHomologyMatrix */

svg{
}

.HomologGroup0 {fill:#000}
.HomologGroup1 {fill:#300}
.HomologGroup2 {fill:#500}
.HomologGroup3 {fill:#700}
.HomologGroup4 {fill:#900}
.HomologGroup5 {fill:#B00}
.HomologGroup6 {fill:#D00}
.HomologGroup7 {fill:#F00}
.HomologGroup8 {fill:#F30}
.HomologGroup9 {fill:#F50}
.HomologGroup10 {fill:#F70}
.HomologGroup11 {fill:#F90}
.HomologGroup12 {fill:#FB0}
.HomologGroup13 {fill:#FD0}
.HomologGroup14 {fill:#FF0}
.HomologGroup15 {fill:#DD0}
.HomologGroup16 {fill:#BB0}
.HomologGroup17 {fill:#990}
.HomologGroup18 {fill:#770}
.HomologGroup19 {fill:#690}
.HomologGroup20 {fill:#6B0}
.HomologGroup21 {fill:#6D0}
.HomologGroup22 {fill:#0F0}
.HomologGroup23 {fill:#0B0}
.HomologGroup24 {fill:#080}
.HomologGroup25 {fill:#030}
.HomologGroup26 {fill:#043}
.HomologGroup27 {fill:#065}
.HomologGroup28 {fill:#087}
.HomologGroup29 {fill:#0A9}
.HomologGroup30 {fill:#0BB}
.HomologGroup31 {fill:#0DD}
.HomologGroup32 {fill:#0FF}
.HomologGroup33 {fill:#2DE}
.HomologGroup34 {fill:#4BD}
.HomologGroup35 {fill:#69C}
.HomologGroup36 {fill:#87B}
.HomologGroup37 {fill:#A5A}
.HomologGroup38 {fill:#B38}
.HomologGroup39 {fill:#C07}
.HomologGroup40 {fill:#A07}
.HomologGroup41 {fill:#807}
.HomologGroup42 {fill:#607}
.HomologGroup43 {fill:#407}
.HomologGroup44 {fill:#117}

.NoHomologyInWindow0  {stroke:black; fill-opacity:1.0; fill:#000}
.NoHomologyInWindow1  {stroke:black; fill-opacity:1.0; fill:#111}
.NoHomologyInWindow2  {stroke:black; fill-opacity:1.0; fill:#222}
.NoHomologyInWindow3  {stroke:black; fill-opacity:1.0; fill:#333}
.NoHomologyInWindow4  {stroke:black; fill-opacity:1.0; fill:#444}
.NoHomologyInWindow5  {stroke:black; fill-opacity:1.0; fill:#555}
.NoHomologyInWindow6  {stroke:black; fill-opacity:1.0; fill:#666}
.NoHomologyInWindow7  {stroke:black; fill-opacity:1.0; fill:#777}
.NoHomologyInWindow8  {stroke:black; fill-opacity:1.0; fill:#888}
.NoHomologyInWindow9  {stroke:black; fill-opacity:1.0; fill:#999}
.NoHomologyInWindow10  {stroke:black; fill-opacity:1.0; fill:#AAA}
.NoHomologyInWindow11  {stroke:black; fill-opacity:1.0; fill:#BBB}
.NoHomologyInWindow12  {stroke:black; fill-opacity:1.0; fill:#CCC}
.NoHomologyInWindow13  {stroke:black; fill-opacity:1.0; fill:#DDD}
.NoHomologyInWindow14  {stroke:black; fill-opacity:1.0; fill:#EEE}

.SpeciesSpecificGenes {stroke:black; fill:white; fill-opacity:1.0 }
"""
